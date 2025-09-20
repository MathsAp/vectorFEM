using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net;
using System.Text;
using System.Threading.Tasks;
using FEM;
using static FEM.IProblem;

namespace Core
{
    public class ParabolicProblem : IProblem
    {
        public ParabolicProblem(
            IFiniteElementMesh mesh, 
            ITimeMesh timeMesh, 
            Func<Vector3D, double> initCondition, 
            IDictionary<string, IMaterial> materials, 
            CoordinateSystem system = CoordinateSystem.Cartesian)
        {
            Mesh = mesh;
            TimeMesh = timeMesh;
            InitialCondition = initCondition;
            Materials = materials;
            this.system = system;
        }

        CoordinateSystem system;
        public IDictionary<string, IMaterial> Materials { get; }
        public IFiniteElementMesh Mesh { get; }
        public ITimeMesh TimeMesh { get; }

        Func<Vector3D, double> InitialCondition { get; }
        PardisoSLAE? SLAE { get; set; }

        public void Prepare()
        {
            FemAlgorithms.EnumerateMeshDofs(Mesh);
            SLAE = new PardisoSLAE(new PardisoMatrix(FemAlgorithms.BuildPortraitFirstStep(Mesh), Quasar.Native.PardisoMatrixType.SymmetricIndefinite));
            //SLAE = new(new PardisoNonSymmMatrix(FemAlgorithms.BuildPortraitFirstStep(Mesh), Quasar.Native.PardisoMatrixType.StructurallySymmetric));
            TimeMesh.ChangeCoefs(GetWeightsForInitialCondition());
        }

        double[] GetWeightsForInitialCondition()
        {
            Parallel.ForEach(Mesh.Elements, element =>
            //foreach (var element in Mesh.Elements)
            {
                var material = Materials[element.Material];

                if (material.IsVolume)
                {

                    var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, _ => 1);
                    SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                    var LRP = element.BuildLocalRightPart(Mesh.Vertex, InitialCondition);
                    SLAE?.AddLocalRightPart(element.Dofs, LRP);

                }

            }
            );

            using (PardisoSLAESolver SLAESolver = new PardisoSLAESolver(SLAE!))
            {
                SLAESolver.Prepare();

                var solutionVector = SLAESolver.Solve();

                return solutionVector;
            }
        }

        public void Solve(ISolution result)
        {
            result.AddSolutionVector(TimeMesh[0], TimeMesh.Coefs(1));

            int tN = TimeMesh.Size();
            double[] timeCoefs = new double[3];

            //Func<Vector3D, Vector3D>? velocity = p => new(p.Y / p.X, p.X, 0);
            Func<Vector3D, Vector3D>? velocity = null;
            PardisoSLAESolver? SLAESolver = null;

            for (int iT = 1; iT < tN; ++iT)
            {
                var isChanged = TimeMesh.IsChangedStep(iT);

                SLAE?.Clear();

                var t = TimeMesh[iT];
                double deltaT;
                double deltaT1;
                double deltaT0 = TimeMesh[iT] - TimeMesh[iT - 1];
                int scheme;

                if (iT != 1)
                {
                    scheme = 3;

                    deltaT = TimeMesh[iT] - TimeMesh[iT - 2];
                    deltaT1 = TimeMesh[iT - 1] - TimeMesh[iT - 2];

                    timeCoefs[0] = (deltaT + deltaT0) / (deltaT * deltaT0);
                    timeCoefs[1] = deltaT / (deltaT1 * deltaT0);
                    timeCoefs[2] = -deltaT0 / (deltaT * deltaT1);
                }
                else
                {
                    scheme = 2;
                    timeCoefs[0] = 1 / deltaT0;
                    timeCoefs[1] = 1 / deltaT0;
                }

                Parallel.ForEach(Mesh.Elements, element =>
                //foreach (var element in Mesh.Elements)
                {
                    var material = Materials[element.Material];
                    Func<Vector3D, double> coeff;
                    if (material.IsVolume)
                    {
                        coeff = system is CoordinateSystem.Cylindrical ? p => p.X * material.Lambda!(p) : material.Lambda!;
                        var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Stiffness, coeff);
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                        coeff = system is CoordinateSystem.Cylindrical ? p => p.X * material.Sigma!(p) : material.Sigma!;
                        if (velocity is not null)
                        {
                            LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Convection, coeff, velocity);
                            SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);
                        }

                        //coeff = system is CoordinateSystem.Cylindrical ? p => p.X * material.Sigma!(p) : material.Sigma!;
                        LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, coeff);

                        coeff = system is CoordinateSystem.Cylindrical ? p => p.X * material.F!(p, t) : p => material.F!(p, t);
                        var LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff);
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);


                        double[] LC = new double[element.Dofs.Length];
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM, timeCoefs[0]);
                        for (int i = 1; i < scheme; ++i)
                        {
                            LRP = LinearAlgebraAlgorithms.MultiplyMatrixByVector(LM, GetLocalCoeffs(LC, element.Dofs, TimeMesh.Coefs(i)), timeCoefs[i]);
                            SLAE?.AddLocalRightPart(element.Dofs, LRP);
                        }

                    }
                    else if (material.Is2)
                    {
                        coeff = system is CoordinateSystem.Cylindrical ? p => p.X * material.Theta!(p, t) : p => material.Theta!(p, t);
                        var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, coeff);
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);
                    }
                    else if (material.Is3)
                    {
                        Func<Vector3D, double> betta = system is CoordinateSystem.Cylindrical ? p => p.X * material.Betta!(p) : material.Betta!;
                        var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, betta);
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                        coeff = p => betta(p) * material.UBetta!(p, t);
                        var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, coeff);
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);
                    }
                }
                );

                foreach (var element in Mesh.Elements)
                {
                    var material = Materials[element.Material];

                    if (material.Is1)
                    {
                        var LRP = element.BuildLocalRightPartWithFirstBoundaryConditions(Mesh.Vertex, Coeff => material.Ug!(Coeff, t));
                        SLAE?.AddFirstBoundaryConditions(element.Dofs, LRP);
                    }
                }

                if (isChanged)
                {
                    SLAESolver?.Dispose();
                    SLAESolver = new PardisoSLAESolver(SLAE!);
                    SLAESolver.Prepare();
                }

                var solutionVector = SLAESolver?.Solve();

                result.AddSolutionVector(TimeMesh[iT], solutionVector!);

                TimeMesh.ChangeCoefs(solutionVector!);
            }

            SLAESolver?.Dispose();
        }

        double[] GetLocalCoeffs(double[] lc, int[] dofs, double[] coeffs)
        {
            for (int i = 0; i < dofs.Length; ++i)
                lc[i] = coeffs[dofs[i]];

            return lc;
        }

    }
}
