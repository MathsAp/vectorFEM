using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    namespace VectorFiniteElements
    {
        namespace FiniteElements3D
        {
            public class LinearVectorParallelepipedalFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<Vector3D, Vector3D>
            {
                public LinearVectorParallelepipedalFiniteElementWithNumInteg(string material, int[] vertexNumber)
                {
                    Material = material;
                    VertexNumber = vertexNumber;
                    MasterElement = CubeMasterElementLinearVectorBasis.GetInstance();
                }

                public IMasterElement<Vector3D, Vector3D> MasterElement { get; }

                public string Material { get; }

                public int[] VertexNumber { get; } = new int[8];

                public int NumberOfEdges => 12;

                public int[] Dofs { get; } = new int[12];

                public bool IsVector => true;

                public int NumberOfFaces => 6;

                public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, Func<Vector3D, double> Coeff)
                {
                    int N = Dofs.Length;

                    var LM = new double[N, N];

                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double hx = VertexCoords[VertexNumber[1]].X - x0;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double hy = VertexCoords[VertexNumber[2]].Y - y0;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double hz = VertexCoords[VertexNumber[4]].Z - z0;

                    double LocalCoeff(Vector3D point) => Coeff(new(point.X * hx + x0, point.Y * hy + y0, point.Z * hz + z0));

                    var nodes = MasterElement.QuadratureNodes;

                    switch (type)
                    {
                        case IFiniteElement.MatrixType.Stiffness:
                            for (int i = 0; i < N; ++i)
                            {
                                for (int j = 0; j < N; ++j)
                                {
                                    double sum = 0;

                                    for (int k = 0; k < nodes.Nodes.Length; ++k)
                                    {
                                        var curli = MasterElement.CurlValues[i, k];
                                        var curlj = MasterElement.CurlValues[j, k];

                                        if (curli.X == 0d)
                                            curli = new Vector3D(curli.X, curli.Y / hz, curli.Z / hy);
                                        else if (curli.Y == 0d)
                                            curli = new Vector3D(curli.X / hz, curli.Y, curli.Z / hx);
                                        else if (curli.Z == 0d)
                                            curli = new Vector3D(curli.X / hy, curli.Y / hx, curli.Z);

                                        if (curlj.X == 0d)
                                            curlj = new Vector3D(curlj.X, curlj.Y / hz, curlj.Z / hy);
                                        else if (curlj.Y == 0d)
                                            curlj = new Vector3D(curlj.X / hz, curlj.Y, curlj.Z / hx);
                                        else if (curlj.Z == 0d)
                                            curlj = new Vector3D(curlj.X / hy, curlj.Y / hx, curlj.Z);

                                        sum += curli * curlj * LocalCoeff(nodes.Nodes[k].Node) * nodes.Nodes[k].Weight;
                                    }

                                    LM[i, j] = sum * hx * hy * hz;
                                }
                            }

                            break;

                        case IFiniteElement.MatrixType.Mass:

                            for (int i = 0; i < N; ++i)
                            {
                                for (int j = 0; j < N; ++j)
                                {
                                    double sum = 0;

                                    for (int k = 0; k < nodes.Nodes.Length; ++k)
                                    {
                                        sum += MasterElement.PsiPsiMatrix[k, i, j] * LocalCoeff(nodes.Nodes[k].Node);
                                    }

                                    LM[i, j] = sum * hx * hy * hz; 
                                }
                            }

                            break;
                    }

                    return LM;
                }

                public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, double> F)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, Vector3D> F)
                {
                    int N = Dofs.Length;

                    double[] LRP = new double[N];

                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double hx = VertexCoords[VertexNumber[1]].X - x0;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double hy = VertexCoords[VertexNumber[2]].Y - y0;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double hz = VertexCoords[VertexNumber[4]].Z - z0;

                    Vector3D LocalF(Vector3D point) => F(new(point.X * hx + x0, point.Y * hy + y0, point.Z * hz + z0));

                    var nodes = MasterElement.QuadratureNodes;

                    for (int i = 0; i < N; ++i)
                    {
                        double sum = 0;

                        for (int k = 0; k < nodes.Nodes.Length; ++k)
                        {
                            sum += MasterElement.PsiValues[i, k] * LocalF(nodes.Nodes[k].Node) * nodes.Nodes[k].Weight;
                        }

                        LRP[i] = sum * hx * hy * hz;
                    }

                    return LRP;
                }

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Ug)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Ug)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Theta)
                {
                    throw new NotImplementedException();
                }

                public int DOFOnEdge(int edge) => 1;

                public int DOFOnElement() => 0;

                public (int i, int j) Edge(int edge)
                {
                    switch(edge)
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

                public void SetEdgeDOF(int edge, int n, int dof)
                {
                    switch(edge)
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
                        case 8:
                            Dofs[8] = dof;
                            break;
                        case 9:
                            Dofs[9] = dof;
                            break;
                        case 10:
                            Dofs[10] = dof;
                            break;
                        case 11:
                            Dofs[11] = dof;
                            break;
                        default:
                            throw new ArgumentException();
                    }
                }

                public void SetElementDOF(int n, int dof)
                {
                    throw new NotImplementedException();
                }

                public void SetVertexDOF(int vertex, int dof)
                {
                    return;
                }

                public int[] Face(int face)
                {
                    switch(face)
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

                static public int[] FaceS(int face)
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

                public int DOFOnFace(int face) => 0;

                public void SetFaceDOF(int face, int n, int dof)
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
                    double x0 = VertexCoords[VertexNumber[0]].X;
                    double x1 = VertexCoords[VertexNumber[1]].X;

                    double y0 = VertexCoords[VertexNumber[0]].Y;
                    double y1 = VertexCoords[VertexNumber[2]].Y;

                    double z0 = VertexCoords[VertexNumber[0]].Z;
                    double z1 = VertexCoords[VertexNumber[4]].Z;

                    if (x0 <= point.X && point.X <= x1 
                        && 
                        y0 <= point.Y && point.Y <= y1 
                        && 
                        z0 <= point.Z && point.Z <= z1)
                        return true;

                    return false;
                }
                public Vector3D GetVectorAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
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

                    Vector3D vec = Vector3D.Zero;

                    for (int i = 0; i < N; ++i)
                        vec += coeffs[Dofs[i]] * Linear3DVectorBasis.Phi[i](xi, eta, zeta);

                    return vec;
                }

                public Vector3D GetCurlAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
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

                    Vector3D curl = Vector3D.Zero;

                    for (int i = 0; i < N; ++i)
                    {
                        var curli = Linear3DVectorBasis.curlPhi[i](xi, eta, zeta);

                        if (curli.X == 0d)
                            curli = new Vector3D(curli.X, curli.Y / hz, curli.Z / hy);
                        else if (curli.Y == 0d)
                            curli = new Vector3D(curli.X / hz, curli.Y, curli.Z / hx);
                        else if (curli.Z == 0d)
                            curli = new Vector3D(curli.X / hy, curli.Y / hx, curli.Z);

                        curl += coeffs[Dofs[i]] * curli;
                    }

                    return curl;
                }

                public static Vector3D GetNormal(int face)
                {
                    switch (face)
                    {
                        case 0:
                            return -Vector3D.XAxis;
                        case 1:
                            return Vector3D.XAxis;
                        case 2:
                            return -Vector3D.YAxis;
                        case 3:
                            return Vector3D.YAxis;
                        case 4:
                            return -Vector3D.ZAxis;
                        case 5:
                            return Vector3D.ZAxis;
                        default:
                            throw new ArgumentException();
                    }

                }

                public static int[] LocalEdgesDofsAtFace(int face)
                {
                    switch(face)
                    {
                        case 0:
                            return [4, 6, 8, 10];
                        case 1:
                            return [5, 7, 9, 11];
                        case 2:
                            return [0, 2, 8, 9];
                        case 3:
                            return [1, 3, 10, 11];
                        case 4:
                            return [0, 1, 4, 5];
                        case 5:
                            return [2, 3, 6, 7];

                        default:
                            throw new ArgumentException();
                    }
                }
            }
        }
    }

}
