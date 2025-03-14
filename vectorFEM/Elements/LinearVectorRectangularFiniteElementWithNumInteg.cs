using Core.VectorFiniteElements.FiniteElements3D;
using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.WebSockets;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using static FEM.IFiniteElement;

namespace Core
{
    namespace VectorFiniteElements
    {
        namespace FiniteElements2D
        {
            public class LinearVectorRectangularFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<Vector2D, Vector3D, Vector3D>
            {

                public LinearVectorRectangularFiniteElementWithNumInteg(string material, int[] vertexNumber)
                {
                    Material = material;
                    VertexNumber = vertexNumber;
                    MasterElement = SquareMasterElementLinearVectorBasis.GetInstance();
                }

                public IMasterElement<Vector2D, Vector3D, Vector3D> MasterElement { get; }

                public bool IsVector => true;

                public string Material { get; }

                public int[] VertexNumber { get; } = new int[4];

                public int NumberOfEdges => 4;

                public int NumberOfFaces => 1;

                public int[] Dofs { get; } = new int[4];

                public ElementType Type => ElementType.Vector;

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
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Ug)
                {
                    int N = Dofs.Length;

                    double[] LRP = new double[N];
                    double[,] M = new double[N, N];

                    Vector3D a = VertexCoords[VertexNumber[0]];
                    Vector3D b = VertexCoords[VertexNumber[1]];
                    Vector3D c = VertexCoords[VertexNumber[2]];
                    Vector3D d = VertexCoords[VertexNumber[3]];

                    double S = a.Distance(b) * a.Distance(c);

                    int face = FacePortrait[GetFaceTuple()].Item1.Item2;

                    Vector3D n = LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(face);

                    Vector3D LocalUg(Vector2D p) => Ug(a * (1 - p.X) * (1 - p.Y) + b * p.X * (1 - p.Y) +
                                                       c * (1 - p.X) * p.Y       + d * p.X * p.Y);

                    var nodes = MasterElement.QuadratureNodes.Nodes;

                    for (int i = 0; i < N; ++i)
                    {
                        for (int j = 0; j < N; ++j)
                        {
                            double sum = 0;
                            for (int k = 0; k < nodes.Length; ++k)
                            {
                                sum += MasterElement.FacesPsiNPsiNMatrix[face, k, i, j];
                            }

                            M[i, j] = sum * S;
                        }
                    }

                    var localUgInNodes = CalcLocalFuncInQuadratureNodes(LocalUg, n, nodes);

                    for (int i = 0; i < N; ++i)
                    {
                        double sum = 0;

                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            var node = nodes[k];
                            sum += node.Weight * MasterElement.FacesPsiNValues[face, i, k] * localUgInNodes[k];
                        }

                        LRP[i] = sum * S;
                    }

                    var localDofs = BuildLocalDofsForLocalMatrix(N);

                    var SLAE = new PardisoSLAE(new PardisoMatrix(BuildProfileForLocalMatrix(N), Quasar.Native.PardisoMatrixType.SymmetricIndefinite));
                    SLAE.Matrix.AddLocal(localDofs, localDofs, M);
                    SLAE.AddLocalRightPart(localDofs, LRP);

                    using (PardisoSLAESolver SLAESolver = new PardisoSLAESolver(SLAE))
                    {
                        SLAESolver.Prepare();

                        var solution = SLAESolver.Solve();

                        return solution;
                    }

                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta)
                {
                    throw new NotImplementedException();
                }

                public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Theta)
                {
                    int N = Dofs.Length;

                    double[] LRP = new double[N];

                    Vector3D a = VertexCoords[VertexNumber[0]];
                    Vector3D b = VertexCoords[VertexNumber[1]];
                    Vector3D c = VertexCoords[VertexNumber[2]];
                    Vector3D d = VertexCoords[VertexNumber[3]];

                    double S = a.Distance(b) * a.Distance(c);

                    int face = FacePortrait[GetFaceTuple()].Item1.Item2;

                    Vector3D n = LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(face);

                    Vector3D LocalTheta(Vector2D p) => Theta(a * (1 - p.X) * (1 - p.Y) + b * p.X * (1 - p.Y) +
                                                             c * (1 - p.X) * p.Y       + d * p.X * p.Y);

                    var nodes = MasterElement.QuadratureNodes.Nodes;

                    var localThetaInNodes = CalcLocalFuncInQuadratureNodes(LocalTheta, n, nodes);

                    for (int i = 0; i < N; ++i)
                    {
                        double sum = 0;
                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            var node = nodes[k].Node;
                            sum += nodes[k].Weight * MasterElement.FacesPsiValues[face, i, k] * localThetaInNodes[k]; // Функцию в узлах на конечном элементе можно насчитать один раз
                        }

                        LRP[i] = sum * S;
                    }

                    //for(int i = 0; i < 2; ++i)
                    //{
                    //    double sum = 0;

                    //    if (face == 0 || face == 1)
                    //        for (int k = 0; k < nodes.Length; ++k)
                    //        {
                    //            sum += nodes[k].Weight * new Vector3D(0d, MasterElement.PsiValues[i, k].X, 0d) * LocalTheta(nodes[k].Node);
                    //        }
                    //    else if (face == 2 || face == 3)
                    //        for (int k = 0; k < nodes.Length; ++k)
                    //        {
                    //            sum += nodes[k].Weight * new Vector3D(MasterElement.PsiValues[i, k].X, 0d, 0d) * LocalTheta(nodes[k].Node);
                    //        }
                    //    else
                    //        for (int k = 0; k < nodes.Length; ++k)
                    //        {
                    //            sum += nodes[k].Weight * new Vector3D(MasterElement.PsiValues[i, k].X, 0d, 0d) * LocalTheta(nodes[k].Node);
                    //        }

                    //    LRP[i] = sum * S;
                    //}

                    //for (int i = 2; i < 4; ++i)
                    //{
                    //    double sum = 0;

                    //    if (face == 0 || face == 1)
                    //        for (int k = 0; k < nodes.Length; ++k)
                    //        {
                    //            sum += nodes[k].Weight * new Vector3D(0d, 0d, MasterElement.PsiValues[i, k].Y) * LocalTheta(nodes[k].Node);
                    //        }
                    //    else if (face == 2 || face == 3)
                    //        for (int k = 0; k < nodes.Length; ++k)
                    //        {
                    //            sum += nodes[k].Weight * new Vector3D(0d, 0d, MasterElement.PsiValues[i, k].Y) * LocalTheta(nodes[k].Node);
                    //        }
                    //    else
                    //        for (int k = 0; k < nodes.Length; ++k)
                    //        {
                    //            sum += nodes[k].Weight * new Vector3D(0d, MasterElement.PsiValues[i, k].Y, 0d) * LocalTheta(nodes[k].Node);
                    //        }

                    //    LRP[i] = sum * S;
                    //}

                    return LRP;
                }

                public int DOFOnEdge(int edge, ElementType type) => 1;

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

                public void SetEdgeDOF(int edge, int n, int dof, ElementType type)
                {
                    switch (edge)
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
                    return;
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

                (int, int, int, int) GetFaceTuple()
                {
                    int[] face = [VertexNumber[0], VertexNumber[1], VertexNumber[2], VertexNumber[3]];

                    Array.Sort(face);

                    var faceTuple = (face[0], face[1], face[2], face[3]);

                    return faceTuple;
                }

                SortedSet<int>[] BuildProfileForLocalMatrix(int N)
                {
                    var profile = new SortedSet<int>[N];

                    for (int i = 0; i < N; ++i)
                    {
                        profile[i] = new();

                        for (int j = 0; j < N; ++j)
                            profile[i].Add(j);
                    }

                    return profile;
                }

                int[] BuildLocalDofsForLocalMatrix(int N)
                {
                    var localDofs = new int[N];

                    for (int i = 0; i < N; ++i)
                        localDofs[i] = i;

                    return localDofs;
                }

                public int[] GetDofs(DofsType type)
                {
                    throw new NotImplementedException();
                }

                public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, MatrixType type, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, double> Coeff)
                {
                    throw new NotImplementedException();
                }

                Vector3D[] CalcLocalFuncInQuadratureNodes(Func<Vector2D, Vector3D> LocalFunc, Vector3D n, Quadratures.QuadratureNode<Vector2D>[] nodes)
                {
                    Vector3D[] values = new Vector3D[nodes.Length];

                    for (int k = 0; k < nodes.Length; ++k)
                    {
                        values[k] = Vector3D.Cross(LocalFunc(nodes[k].Node), n);
                    }

                    return values;
                }

                public double CalclIntegralOfSquaredDifference(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Func<Vector3D, Vector3D> u)
                {
                    throw new NotImplementedException();
                }
            }
        }
    }
}
