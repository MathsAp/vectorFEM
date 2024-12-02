using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public class EllipticProblem : IProblem
    {
        public EllipticProblem(IFiniteElementMesh mesh, IDictionary<string, IMaterial> materials) 
        {
            Materials = materials;
            Mesh = mesh;
        }

        public IDictionary<string, IMaterial> Materials { get; }
        public IFiniteElementMesh Mesh { get; }

        PardisoSLAE? SLAE { get; set; }

        public void Prepare()
        {
            FemAlgorithms.EnumerateMeshDofs(Mesh);
            SLAE = new PardisoSLAE(new PardisoMatrix(FemAlgorithms.BuildPortraitFirstStep(Mesh), Quasar.Native.PardisoMatrixType.SymmetricIndefinite));
        }

        public void Solve(ISolution result)
        {

            Parallel.ForEach(Mesh.Elements, element =>
            {
                var material = Materials[element.Material];

                if (material.IsVolume)
                {
                    var LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Stiffness, material.Lambda);
                    SLAE?.Matrix.AddLocal(element.Dofs, LM);

                    LM = element.BuildLocalMatrix(Mesh.Vertex, IFiniteElement.MatrixType.Mass, material.Sigma);
                    SLAE?.Matrix.AddLocal(element.Dofs, LM);

                    var LRP = element.BuildLocalRightPart(Mesh.Vertex, Coeff => material.F(Coeff, 1));
                    SLAE?.AddLocalRightPart(element.Dofs, LRP);
                }
                else if (material.Is2)
                {
                    var LRP = element.BuildLocalRightPartWithSecondBoundaryConditions(Mesh.Vertex, Coeff => material.Theta(Coeff, 1));
                    SLAE?.AddLocalRightPart(element.Dofs, LRP);
                }
            }
            );

            foreach(var element in Mesh.Elements)
            {
                var material = Materials[element.Material];

                if (material.Is1)
                {
                    var LRP = element.BuildLocalRightPartWithFirstBoundaryConditions(Mesh.Vertex, Coeff => material.Ug(Coeff, 1));
                    SLAE?.AddFirstBoundaryConditions(element.Dofs, LRP);
                }
            }

            using (PardisoSLAESolver solver = new PardisoSLAESolver(SLAE!))
            {
                solver.Prepare();

                var solutionVector = solver.Solve();

                result.AddSolutionVector(1, solutionVector);
            }
        }
    }
}
