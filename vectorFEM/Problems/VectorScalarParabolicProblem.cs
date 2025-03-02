using Core.VectorFiniteElements.FiniteElements3D;
using Core.VectorScalarFiniteElements.FiniteElements2D;
using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static FEM.IFiniteElement;

namespace Core
{
    public class VectorScalarParabolicProblem : IProblem
    {
        public VectorScalarParabolicProblem(IFiniteElementMesh mesh, ITimeMesh timeMesh, IDictionary<string, IMaterial> materials)
        {
            Mesh = mesh;
            TimeMesh = timeMesh;
            //InitialCondition = initCondition;
            Materials = materials;
        }

        public IDictionary<string, IMaterial> Materials { get; }

        public IFiniteElementMesh Mesh { get; }

        public ITimeMesh TimeMesh { get; }

        //Func<Vector3D, Vector3D> InitialCondition { get; }

        PardisoSLAE? SLAE { get; set; }



        public void Prepare()
        {
            FemAlgorithms.EnumerateMeshDofs(Mesh);
            SLAE = new(new PardisoNonSymmMatrix(FemAlgorithms.BuildPortraitFirstStep(Mesh), Quasar.Native.PardisoMatrixType.StructurallySymmetric));

            TimeMesh.ChangeCoefs(GetWeightsForInitialCondition());
        }

        double[] GetWeightsForInitialCondition()
        { 
            return new double[Mesh.NumberOfDofs];
        }

        public void Solve(ISolution result)
        {
            result.AddSolutionVector(TimeMesh[0], TimeMesh.Coefs(1));

            int tN = TimeMesh.Size();
            double[] timeCoeffs = new double[3];

            PardisoSLAESolver? SLAESolver = null;

            for (int iT = 1; iT < tN; ++iT)
            {
                bool isChanged = TimeMesh.IsChangedStep(iT);

                SLAE?.Clear();

                double t = TimeMesh[iT];
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

                Parallel.ForEach(Mesh.Elements, element =>
                {
                    var material = Materials[element.Material];

                    if (material.IsVolume)
                    {
                        if (element.Type == IFiniteElement.ElementType.Vector)
                        {
                            var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Stiffness, coeff => 1 / material.Mu(coeff));
                            SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                            // LinearAlgebraAlgorithms.PrintMatrix(LM);

                            LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, material.Sigma);

                            var LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff => material.Fv(coeff, t));
                            SLAE?.AddLocalRightPart(element.Dofs, LRP);

                            double[] LC = new double[element.Dofs.Length];
                            SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM, timeCoeffs[0]);
                            for (int i = 1; i < scheme; ++i)
                            {
                                LRP = LinearAlgebraAlgorithms.MultiplyMatrixVector(LM, GetLocalCoeffs(LC, element.Dofs, TimeMesh.Coefs(i)), timeCoeffs[i]);
                                SLAE?.AddLocalRightPart(element.Dofs, LRP);
                            }
                        }
                        else if (element.Type == IFiniteElement.ElementType.Scalar)
                        {
                            var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Stiffness, coeff => Constants.Mu0);
                            SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                            //   LinearAlgebraAlgorithms.PrintMatrix(LM);

                            var LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff => material.F(coeff, t));
                            SLAE?.AddLocalRightPart(element.Dofs, LRP);
                        }
                    }
                    else if (material.Is2)
                    {
                        if (element.Type == IFiniteElement.ElementType.Vector)
                        {
                            var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, Mesh.FacePortrait, coeff => material.Htheta(coeff, t));
                            SLAE?.AddLocalRightPart(element.Dofs, LRP);
                        }
                        else if (element.Type == IFiniteElement.ElementType.Scalar)
                        {
                            var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, coeff => material.Theta(coeff, t));
                            SLAE?.AddLocalRightPart(element.Dofs, LRP);
                        }
                    }
                    else if (material.IsInterface)
                    {
                        var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Interface, Mesh.FacePortrait, coeff => 1);
                        SLAE?.Matrix.AddLocal(element.GetDofs(IFiniteElement.DofsType.Scalar), element.GetDofs(IFiniteElement.DofsType.Vector), LM, -1);

                        SLAE?.Matrix.AddLocalTransposed(element.GetDofs(IFiniteElement.DofsType.Vector), element.GetDofs(IFiniteElement.DofsType.Scalar), LM);

                        var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, Mesh.FacePortrait, coeff => material.Hext(coeff, t));
                        SLAE?.AddLocalRightPart(element.GetDofs(IFiniteElement.DofsType.Vector), LRP);

                        //LinearAlgebraAlgorithms.PrintVector(LRP);

                        (int faceS, int faceV) = ((LinearVectorScalarRectangularFiniteElementWithNumInteg)element).GetFaceNumbers(Mesh.FacePortrait);

                        var n = LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(faceV);

                        LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, coeff => -Constants.Mu0 * material.Hext(coeff, t) * n);
                        SLAE?.AddLocalRightPart(element.GetDofs(DofsType.Scalar), LRP);

                        // LinearAlgebraAlgorithms.PrintVector(LRP);
                    }
                });

                foreach (var element in Mesh.Elements)
                {
                    var material = Materials[element.Material];

                    if (material.Is1)
                    {
                        if (element.Type == ElementType.Vector)
                        {
                            var LRP = element.BuildLocalRightPartWithFirstBoundaryConditions(Mesh.Vertex, Mesh.FacePortrait, coeff => material.Ag(coeff, t));
                            SLAE?.AddFirstBoundaryConditions(element.Dofs, LRP);
                        }
                        else if (element.Type == ElementType.Scalar)
                        {
                            var LRP = element.BuildLocalRightPartWithFirstBoundaryConditions(Mesh.Vertex, coeff => material.Ug(coeff, t));
                            SLAE?.AddFirstBoundaryConditions(element.Dofs, LRP);
                        }
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
