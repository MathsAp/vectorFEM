using Core.ScalarFiniteElements.FiniteElements3D;
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
    public class VectorScalarEllipticProblem : IProblem
    {
        public VectorScalarEllipticProblem(IFiniteElementMesh mesh, IDictionary<string, IMaterial> materials)
        {
            Mesh = mesh;
            Materials = materials;
        }

        public IDictionary<string, IMaterial> Materials { get; }

        public IFiniteElementMesh Mesh { get; }

        public ITimeMesh TimeMesh { get; }

        PardisoSLAE? SLAE { get; set; }

        public void Prepare()
        {
            FemAlgorithms.EnumerateMeshDofs(Mesh);
            SLAE = new(new PardisoNonSymmMatrix(FemAlgorithms.BuildPortraitFirstStep(Mesh), Quasar.Native.PardisoMatrixType.StructurallySymmetric));
        }

        public void Solve(ISolution result)
        {
            Parallel.ForEach(Mesh.Elements, element =>
            {
                var material = Materials[element.Material];

                if (material.IsVolume)
                {
                    if (element.Type == IFiniteElement.ElementType.Vector)
                    {
                        var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Stiffness, coeff => 1 / material.Mu(coeff));
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                        LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, material.Sigma);
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                        var LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff => material.Fv(coeff, 1));
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);
                    }
                    else if (element.Type == IFiniteElement.ElementType.Scalar)
                    {
                        var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Stiffness, coeff => Constants.Mu0);
                        SLAE?.Matrix.AddLocal(element.Dofs, element.Dofs, LM);

                        var LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff => material.F(coeff, 1));
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);
                    }
                }
                else if (material.Is2)
                {
                    if (element.Type == IFiniteElement.ElementType.Vector)
                    {
                        var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, Mesh.FacePortrait, coeff => material.Htheta(coeff, 1));
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);
                    }
                    else if (element.Type == IFiniteElement.ElementType.Scalar)
                    {
                        var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, coeff => material.Theta(coeff, 1));
                        SLAE?.AddLocalRightPart(element.Dofs, LRP);
                    }
                }
                else if (material.IsInterface)
                {
                    var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Interface, coeff => 1);
                    SLAE?.Matrix.AddLocal(element.GetDofs(IFiniteElement.DofsType.Scalar), element.GetDofs(IFiniteElement.DofsType.Vector), LM, -1);

                    SLAE?.Matrix.AddLocalTransposed(element.GetDofs(IFiniteElement.DofsType.Vector), element.GetDofs(IFiniteElement.DofsType.Scalar), LM);

                    var LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff => material.Hext(coeff, 1));
                    SLAE?.AddLocalRightPart(element.GetDofs(IFiniteElement.DofsType.Vector), LRP);

                    (int faceS, int faceV) = ((LinearVectorScalarRectangularFiniteElementWithNumInteg)element).GetFaceNumbers(Mesh.FacePortrait);

                    var n = LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(faceV);

                    LRP = element.BuildLocalRightPart(Mesh.Vertex, coeff => -Constants.Mu0 * material.Hext(coeff, 1) * n);
                    SLAE?.AddLocalRightPart(element.GetDofs(DofsType.Scalar), LRP);
                }

            }
            );

            foreach(var element in Mesh.Elements)
            {
                var material = Materials[element.Material];

                if (material.Is1)
                {
                    if (element.Type == ElementType.Vector)
                    {
                        var LRP = element.BuildLocalRightPartWithFirstBoundaryConditions(Mesh.Vertex, Mesh.FacePortrait, coeff => material.Ag(coeff, 1));
                        SLAE?.AddFirstBoundaryConditions(element.Dofs, LRP);
                    }
                    else if (element.Type == ElementType.Scalar)
                    {
                        var LRP = element.BuildLocalRightPartWithFirstBoundaryConditions(Mesh.Vertex, coeff => material.Ug(coeff, 1));
                        SLAE?.AddFirstBoundaryConditions(element.Dofs, LRP);
                    }
                }
            }

            using (PardisoSLAESolver solver = new(SLAE!))
            {
                solver.Prepare();

                var solutionVector = solver.Solve();

                result.AddSolutionVector(1, solutionVector);
            }
        }
    }
}
