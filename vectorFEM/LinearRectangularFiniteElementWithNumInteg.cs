using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static FEM.IFiniteElement;

namespace Core
{
    namespace ScalarFiniteElements
    {
        namespace FiniteElements2D
        {
            public class LinearRectangularFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<Vector2D, double, Vector3D>
            {
                public LinearRectangularFiniteElementWithNumInteg(string material, int[] vertexNumber)
                {
                    Material = material;
                    VertexNumber = vertexNumber;
                    MasterElement = SquareMasterElementLinearScalarBasis.GetInstance();
                }

                public IMasterElement<Vector2D, double, Vector3D> MasterElement { get; }

                public bool IsVector => false;

                public string Material { get; }

                public int[] VertexNumber { get; } = new int[4];

                public int NumberOfEdges => 4;

                public int NumberOfFaces => 1;

                public int[] Dofs { get; } = new int[4];

                public ElementType Type => ElementType.Scalar;

                public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, Func<Vector3D, double> Coeff)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, double> F)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, Vector3D> F)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Ug)
                {
                    int N = Dofs.Length;

                    double[] LRP = new double[N];

                    for (int i = 0; i < N; ++i)
                    {
                        LRP[i] = Ug(VertexCoords[VertexNumber[i]]);
                    }

                    return LRP;
                }

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Ug)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta)
                {
                    int N = Dofs.Length;

                    double[] LRP = new double[N];

                    Vector3D a = VertexCoords[VertexNumber[0]];
                    Vector3D b = VertexCoords[VertexNumber[1]];
                    Vector3D c = VertexCoords[VertexNumber[2]];
                    Vector3D d = VertexCoords[VertexNumber[3]];

                    double S = a.Distance(b) * a.Distance(c);

                    double LocalTheta(Vector2D p) => Theta(a * (1 - p.X) * (1 - p.Y) + b * p.X * (1 - p.Y) +
                                                           c * (1 - p.X) * p.Y       + d * p.X * p.Y);

                    var nodes = MasterElement.QuadratureNodes.Nodes;

                    var localThetaInNodes = CalcLocalFuncInQuadratureNodes(LocalTheta, nodes);

                    for (int i = 0; i < N; ++i)
                    {
                        double sum = 0;
                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            sum += nodes[k].Weight * MasterElement.PsiValues[i, k] * localThetaInNodes[k];
                        }

                        LRP[i] = sum * S;
                    }

                    return LRP;
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Theta)
                {
                    throw new NotImplementedException();
                }

                public int DOFOnEdge(int edge, ElementType type) => 0;

                public int DOFOnElement() => 0;

                public int DOFOnFace(int face, ElementType type) => 0;

                public (int i, int j) Edge(int edge)
                {
                    switch (edge)
                    {
                        case 0:
                            return (0, 1);
                        case 1:
                            return (2, 3);
                        case 2:
                            return (0, 2);
                        case 3:
                            return (1, 3);
                        default:
                            throw new ArgumentException();
                    }
                }

                public int[] Face(int face) => face == 0 ? [0, 1, 2, 3] : throw new ArgumentException();


                public Vector3D GetCurlAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
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

                public Vector3D GetVectorAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
                {
                    throw new NotImplementedException();
                }

                public bool IsPointOnElement(Vector3D[] VertexCoords, Vector3D point)
                {
                    throw new NotImplementedException();
                }

                public void SetEdgeDOF(int edge, int n, int dof, ElementType type)
                {
                    throw new NotImplementedException();
                }

                public void SetElementDOF(int n, int dof)
                {
                    throw new NotImplementedException();
                }

                public void SetFaceDOF(int face, int n, int dof, ElementType type)
                {
                    throw new NotImplementedException();
                }

                public void SetVertexDOF(int vertex, int dof)
                {
                    switch (vertex)
                    {
                        case 0:
                            Dofs[0] = dof;
                            break;
                        case 1: 
                            Dofs[1] = dof;
                            break;
                        case 2:
                            Dofs[2] = dof;
                            break;
                        case 3:
                            Dofs[3] = dof;
                            break;

                        default:
                            throw new ArgumentException();
                    }

                }

                public double[] CalcLocalFuncInQuadratureNodes(Func<Vector2D, double> LocalFunc, Quadratures.QuadratureNode<Vector2D>[] nodes)
                {
                    double[] values = new double[nodes.Length];

                    for (int k = 0; k < nodes.Length; ++k)
                    {
                        values[k] = LocalFunc(nodes[k].Node);
                    }

                    return values;
                }

                public int[] GetDofs(DofsType type)
                {
                    throw new NotImplementedException();
                }

                public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, MatrixType type, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, double> Coeff)
                {
                    throw new NotImplementedException();
                }
            }
        }
    }
}
