using FEM;
using Quadratures;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static FEM.IFiniteElement;

namespace Core
{
    namespace ScalarFiniteElements
    {
        namespace FiniteElements3D
        {
            public class LinearParallelepipedalFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<Vector3D, double, double>
            {
                public LinearParallelepipedalFiniteElementWithNumInteg(string material, int[] vertexNumber)
                {
                    Material = material;
                    VertexNumber = vertexNumber;
                    MasterElement = CubeMasterElementLinearScalarBasis.GetInstance();
                }

                public ElementType Type => ElementType.Scalar;
                public IMasterElement<Vector3D, double, double> MasterElement { get; }

                public string Material { get; }

                public int[] VertexNumber { get; } = new int[8];

                public int NumberOfEdges => 12;

                public int[] Dofs { get; } = new int[8];

                public bool IsVector => false;

                public int NumberOfFaces => 6;

                public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, Func<Vector3D, double> Coeff, Func<Vector3D, Vector3D>? Velocity = null)
                {
                    int N = Dofs.Length;

                    var LM = new double[N, N];

                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double hx = VertexCoords[VertexNumber[1]].X - x0;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double hy = VertexCoords[VertexNumber[2]].Y - y0;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double hz = VertexCoords[VertexNumber[4]].Z - z0;

                    double V = hx * hy * hz;

                    double LocalCoeff(Vector3D p) => Coeff(new(hx * p.X + x0, hy * p.Y + y0, hz * p.Z + z0));

                    var nodes = MasterElement.QuadratureNodes.Nodes;
                    var localCoeffInNodes = CalcLocalFuncInQuadratureNodes(LocalCoeff, nodes);


                    switch(type)
                    {
                        case IFiniteElement.MatrixType.Stiffness:

                            double[,] JTJ = { { 1 / (hx * hx), 0,             0             },
                                              { 0,             1 / (hy * hy), 0             },
                                              { 0,             0,             1 / (hz * hz) } };

                            for (int i = 0; i < N; ++i)
                            {
                                for(int j = 0; j < N; ++j)
                                {
                                    double sum = 0;

                                    for (int k = 0; k < nodes.Length; ++k)
                                    {
                                        sum += nodes[k].Weight * MasterElement.GradValues[i, k] 
                                                            * LinearAlgebraAlgorithms.MultiplyMatrix3By3ByVector(JTJ, MasterElement.GradValues[j, k].AsArray()) 
                                                            * localCoeffInNodes[k];
                                    }

                                    LM[i, j] = sum * V;
                                }
                            }

                            break;

                        case IFiniteElement.MatrixType.Mass:

                            for (int i = 0; i < N; ++i)
                            {
                                for (int j = 0; j < N; ++j)
                                {
                                    double sum = 0;

                                    for (int k = 0; k < nodes.Length; ++k)
                                    {
                                        sum += MasterElement.PsiPsiMatrix[k, i, j] * localCoeffInNodes[k];
                                    }

                                    LM[i, j] = sum * V;
                                }
                            }

                            break;

                        default:
                            throw new ArgumentException();
                    }

                    return LM;
                }

                public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, double> F)
                {
                    int N = Dofs.Length;

                    double[] LRP = new double[N];

                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double hx = VertexCoords[VertexNumber[1]].X - x0;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double hy = VertexCoords[VertexNumber[2]].Y - y0;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double hz = VertexCoords[VertexNumber[4]].Z - z0;

                    double V = hx * hy * hz;

                    double LocalF(Vector3D p) => F(new(hx * p.X + x0, hy * p.Y + y0, hz * p.Z + z0));

                    var nodes = MasterElement.QuadratureNodes.Nodes;
                    var localFInNodes = CalcLocalFuncInQuadratureNodes(LocalF, nodes);

                    for (int i = 0; i < N; ++i)
                    {
                        double sum = 0;
                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            sum += nodes[k].Weight * MasterElement.PsiValues[i, k] * localFInNodes[k];
                        }

                        LRP[i] = sum * V;
                    }

                    return LRP;
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
                public static int Upsilon(int i) => i / 4;

                public int DOFOnEdge(int edge, ElementType type) => 0;

                public int DOFOnElement() => 0;

                public (int i, int j) Edge(int edge)
                {
                    switch (edge)
                    {
                        case 0:
                            return (0, 1);
                        case 1:
                            return (2, 3);
                        case 2:
                            return (4, 5);
                        case 3:
                            return (6, 7);
                        case 4:
                            return (0, 2);
                        case 5:
                            return (1, 3);
                        case 6:
                            return (4, 6);
                        case 7:
                            return (5, 7);
                        case 8:
                            return (0, 4);
                        case 9:
                            return (1, 5);
                        case 10:
                            return (2, 6);
                        case 11:
                            return (3, 7);

                        default:
                            throw new ArgumentException();
                    }
                }

                public Vector3D GetGradientAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point) //
                {
                    int N = Dofs.Length;

                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double hx = VertexCoords[VertexNumber[1]].X - x0;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double hy = VertexCoords[VertexNumber[2]].Y - y0;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double hz = VertexCoords[VertexNumber[4]].Z - z0;

                    double xi = (point.X - x0) / hx;
                    double eta = (point.Y - y0) / hy;
                    double zeta = (point.Z - z0) / hz;

                    double xvalue = 0;
                    double yvalue = 0;
                    double zvalue = 0;

                    double coeff = 0;

                    for (int i = 0; i < N; ++i)
                    {
                        coeff = coeffs[Dofs[i]];

                        xvalue += coeff * LinearBasisDerivatives.Phi[Mu(i)](xi) * (1 / hx) * LinearBasis.Phi[Nu(i)](eta) * LinearBasis.Phi[Upsilon(i)](zeta);
                        yvalue += coeff * LinearBasis.Phi[Mu(i)](xi) * LinearBasisDerivatives.Phi[Nu(i)](eta) * (1 / hy) * LinearBasis.Phi[Upsilon(i)](zeta);
                        zvalue += coeff * LinearBasis.Phi[Mu(i)](xi) * LinearBasis.Phi[Nu(i)](eta) * LinearBasisDerivatives.Phi[Upsilon(i)](zeta) * (1 / hz);
                    }

                    return new Vector3D(xvalue, yvalue, zvalue);
                }

                public double GetValueAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point) //
                {
                    int N = Dofs.Length;

                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double hx = VertexCoords[VertexNumber[1]].X - x0;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double hy = VertexCoords[VertexNumber[2]].Y - y0;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double hz = VertexCoords[VertexNumber[4]].Z - z0;

                    double xi = (point.X - x0) / hx;
                    double eta = (point.Y - y0) / hy;
                    double zeta = (point.Z - z0) / hz;

                    double value = 0;

                    for (int i = 0; i < N; ++i)
                    {
                        value += coeffs[Dofs[i]] * LinearBasis.Phi[Mu(i)](xi) * LinearBasis.Phi[Nu(i)](eta) * LinearBasis.Phi[Upsilon(i)](zeta);
                    }

                    return value;
                }

                public bool IsPointOnElement(Vector3D[] VertexCoords, Vector3D point) //
                {
                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double x1 = VertexCoords[VertexNumber[1]].X;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double y1 = VertexCoords[VertexNumber[2]].Y;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double z1 = VertexCoords[VertexNumber[4]].Z;

                    if (x0 <= point.X && point.X <= x1 &&
                        y0 <= point.Y && point.Y <= y1 &&
                        z0 <= point.Z && point.Z <= z1)
                        return true;

                    return false;
                }

                public void SetEdgeDOF(int edge, int n, int dof, ElementType type)
                {
                    throw new NotImplementedException();
                }

                public void SetElementDOF(int n, int dof)
                {
                    throw new NotImplementedException();
                }

                public void SetVertexDOF(int vertex, int dof)
                {
                    switch(vertex)
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
                        case 4:
                            Dofs[4] = dof;
                            break;
                        case 5:
                            Dofs[5] = dof;
                            break;
                        case 6:
                            Dofs[6] = dof;
                            break;
                        case 7: 
                            Dofs[7] = dof;
                            break;

                        default:
                            throw new ArgumentException();
                    }
                }

                public int[] Face(int face)
                {
                    switch (face)
                    {
                        case 0:
                            return [0, 2, 4, 6];
                        case 1:
                            return [1, 3, 5, 7];
                        case 2:
                            return [0, 1, 4, 5];
                        case 3:
                            return [2, 3, 6, 7];
                        case 4:
                            return [0, 1, 2, 3];
                        case 5:
                            return [4, 5, 6, 7];
                        default:
                            throw new ArgumentException();
                    }
                }

                public int DOFOnFace(int face, ElementType type) => 0;

                public void SetFaceDOF(int face, int n, int dof, ElementType type)
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

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Ug)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Theta)
                {
                    throw new NotImplementedException();
                }

                public double[] CalcLocalFuncInQuadratureNodes(Func<Vector3D, double> LocalFunc, Quadratures.QuadratureNode<Vector3D>[] nodes)
                {
                    double[] values= new double[nodes.Length];

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

                public double CalclIntegralOfSquaredDifference(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Func<Vector3D, Vector3D> u)
                {
                    throw new NotImplementedException();
                }
            }
        }
    }
}
