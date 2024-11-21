using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    namespace ScalarFiniteElements
    {
        namespace FiniteElements3D
        {
            public class LinearParallelepipedalFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<Vector3D, double>
            {
                public LinearParallelepipedalFiniteElementWithNumInteg(string material, int[] vertexNumber)
                {
                    Material = material;
                    VertexNumber = vertexNumber;
                    MasterElement = CubeMasterElementLinearScalarBasis.GetInstance();
                }

                public IMasterElement<Vector3D, double> MasterElement { get; }

                public string Material { get; }

                public int[] VertexNumber { get; } = new int[8];

                public int NumberOfEdges => 12;

                public int[] Dofs { get; } = new int[12];

                public bool IsVector => false;

                public int NumberOfFaces => 6;

                public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, Func<Vector3D, double> Coeff)
                {
                    int N = Dofs.Length;

                    var LM = new double[N, N];



                    return LM;
                }

                public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, double> F)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Ug)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta)
                {
                    throw new NotImplementedException();
                }

                public static int Mu(int i) => i % 2;
                public static int Nu(int i) => (i / 2) % 2;
                public static int upsilon(int i) => i / 4;

                public int DOFOnEdge(int edge)
                {
                    throw new NotImplementedException();
                }

                public int DOFOnElement()
                {
                    throw new NotImplementedException();
                }

                public (int i, int j) Edge(int edge)
                {
                    throw new NotImplementedException();
                }

                public Vector3D GetGradientAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
                {
                    throw new NotImplementedException();
                }

                public double GetValueAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
                {
                    throw new NotImplementedException();
                }

                public bool IsPointOnElement(Vector3D[] VertexCoords, Vector3D point)
                {
                    throw new NotImplementedException();
                }

                public void SetEdgeDOF(int edge, int n, int dof)
                {
                    throw new NotImplementedException();
                }

                public void SetElementDOF(int n, int dof)
                {
                    throw new NotImplementedException();
                }

                public void SetVertexDOF(int vertex, int dof)
                {
                    throw new NotImplementedException();
                }

                public int[] Face(int face)
                {
                    throw new NotImplementedException();
                }

                public int DOFOnFace(int face)
                {
                    throw new NotImplementedException();
                }

                public void SetFaceDOF(int face, int n, int dof)
                {
                    throw new NotImplementedException();
                }

                public Vector3D GetVectorAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
                {
                    throw new NotImplementedException();
                }

                public Vector3D GetCurlAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, Vector3D> F)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, Vector3D> Ug)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, Vector3D> Theta)
                {
                    throw new NotImplementedException();
                }
            }
        }
    }
}
