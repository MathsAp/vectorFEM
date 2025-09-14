using FEM;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public class VectorParabolicProblem : IProblem
    {
        public VectorParabolicProblem(IFiniteElementMesh mesh, ITimeMesh timeMesh, Func<Vector3D, Vector3D> initCondition, IDictionary<string, IMaterial> materials)
        {
            Mesh = mesh;
            TimeMesh = timeMesh;
            Materials = materials;
            InitialCondition = initCondition;
        }


        public IDictionary<string, IMaterial> Materials { get; }

        public IFiniteElementMesh Mesh { get; }

        public ITimeMesh TimeMesh { get; }

        Func<Vector3D, Vector3D> InitialCondition { get; }
        PardisoSLAE? SLAE { get; set; }

        public void Prepare()
        {
            FemAlgorithms.EnumerateMeshDofs(Mesh);
            SLAE = new(new PardisoMatrix(FemAlgorithms.BuildPortraitFirstStep(Mesh), Quasar.Native.PardisoMatrixType.SymmetricIndefinite));
            //SLAE = new(new PardisoNonSymmMatrix(FemAlgorithms.BuildPortraitFirstStep(Mesh), Quasar.Native.PardisoMatrixType.StructurallySymmetric));
            TimeMesh.ChangeCoefs(GetWeightsForInitialCondition());
        }

        double[] GetWeightsForInitialCondition()
        {
           //foreach (var element in Mesh.Elements)
            Parallel.ForEach(Mesh.Elements, element =>
            {
                var material = Materials[element.Material];

                if (material.IsVolume)
                {
                    var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, _ => 1);
                    SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                    var LRP = element.BuildLocalRightPart(Mesh.Vertex, InitialCondition);
                    SLAE?.AddLocalRightPart(element.Dofs, LRP);
                }
            });

            using PardisoSLAESolver solver = new(SLAE!);
            solver.Prepare();

            var solution = solver.Solve();

            return solution;
        }

        public void Solve(ISolution result)
        {
            result.AddSolutionVector(TimeMesh[0], TimeMesh.Coefs(1));

            int tN = TimeMesh.Size();
            double[] timeCoeffs = new double[3];

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

                    timeCoeffs[0] = (deltaT + deltaT0) / (deltaT * deltaT0);
                    timeCoeffs[1] = deltaT / (deltaT1 * deltaT0);
                    timeCoeffs[2] = -deltaT0 / (deltaT * deltaT1);
                }
                else
                {
                    scheme = 2;
                    timeCoeffs[0] = 1 / deltaT0;
                    timeCoeffs[1] = 1 / deltaT0;
                }

                //foreach (var element in Mesh.Elements)
                Parallel.ForEach(Mesh.Elements, element =>
                {
                    var material = Materials[element.Material];

                    if (material.IsVolume)
                    {
                        var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Stiffness, coeff => 1d / material.Mu(coeff));
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                        LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, material.Sigma);

                        var LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff => material.Fv(coeff, t));
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);


                        double[] LC = new double[element.Dofs.Length];
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM, timeCoeffs[0]);
                        for (int i = 1; i < scheme; ++i)
                        {
                            LRP = LinearAlgebraAlgorithms.MultiplyMatrixByVector(LM, GetLocalCoeffs(LC, element.Dofs, TimeMesh.Coefs(i)), timeCoeffs[i]);
                            SLAE?.AddLocalRightPart(element.Dofs, LRP);
                        }
                    }
                    else if (material.Is2)
                    {
                        var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, Mesh.FacePortrait, coeff => material.Htheta(coeff, t));
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);
                    }
                }
                );

                foreach(var element in Mesh.Elements)
                {
                    var material = Materials[element.Material];

                    if (material.Is1)
                    {
                        var LRP = element.BuildLocalRightPartWithFirstBoundaryConditions(Mesh.Vertex, Mesh.FacePortrait, coeff => material.Ag(coeff, t));
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

                result.AddSolutionVector(t, solutionVector!);

                TimeMesh.ChangeCoefs(solutionVector!);
            }

            SLAESolver?.Dispose();
        }

        double[] GetLocalCoeffs(double[] lc, int[] dofs, double[] coeffs)
        {
            for (int i = 0; i < dofs.Length; ++i)
            {
                lc[i] = coeffs[dofs[i]];
            }

            return lc;
        }
    }
}
