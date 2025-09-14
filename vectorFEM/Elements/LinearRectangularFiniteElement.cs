using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FEM;
using static FEM.IFiniteElement;

namespace Core
{
    namespace FiniteElements2D
    {
        public class LinearRectangularFiniteElement : IFiniteElement
        {
            public LinearRectangularFiniteElement(string material, int[] vertexNumber)
            {

                Material = material;
                VertexNumber = vertexNumber;
            }

            public ElementType Type => ElementType.Scalar;
            public string Material { get; }
            public int[] VertexNumber { get; } = new int[4];

            public int NumberOfEdges => 4;

            public int[] Dofs { get; } = new int[4];

            public bool IsVector => throw new NotImplementedException();

            public int NumberOfFaces => throw new NotImplementedException();

            public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, Func<Vector3D, double> Coeff, Func<Vector3D, Vector3D>? Velocity = null)
            {
                int N = Dofs.Length;

                var LM = new double[N, N]; // Local Matrix

                double hx = VertexCoords[VertexNumber[1]].X - VertexCoords[VertexNumber[0]].X;
                double hy = VertexCoords[VertexNumber[3]].Y - VertexCoords[VertexNumber[0]].Y;
                double averageCoeff = CalcAverageCoeff(VertexCoords, Coeff);

                var G = LinearLocalMatrices.G;
                var M = LinearLocalMatrices.M;


                switch (type)
                {
                    case IFiniteElement.MatrixType.Stiffness:

                        for (int i = 0; i < N; ++i)
                        {
                            for (int j = 0; j < N; ++j)
                            {
                                LM[i, j] += averageCoeff * ((hy / hx) * G[Mu(i), Mu(j)] * M[Nu(i), Nu(j)] + (hx / hy) * M[Mu(i), Mu(j)] * G[Nu(i), Nu(j)]);
                            }
                        }

                        break;
                    case IFiniteElement.MatrixType.Mass:

                        for (int i = 0; i < N; ++i)
                        {
                            for (int j = 0; j < N; ++j)
                            {
                                LM[i, j] += averageCoeff * hx * M[Mu(i), Mu(j)] * hy * M[Nu(i), Nu(j)];
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

                double[] LRP = new double[N]; // Local Right Part

                double hx = VertexCoords[VertexNumber[1]].X - VertexCoords[VertexNumber[0]].X;
                double hy = VertexCoords[VertexNumber[3]].Y - VertexCoords[VertexNumber[0]].Y;

                var M = LinearLocalMatrices.M;
                double[] LF = CalcLocalVectorF(VertexCoords, F);

                for (int i = 0; i < N; ++i)
                {
                    double el = 0;
                    for (int j = 0; j < N; ++j)
                    {
                        el += hx * M[Mu(i), Mu(j)] * hy * M[Nu(i), Nu(j)] * LF[j];
                    }

                    LRP[i] += el;
                }

                return LRP;
            }
            public int DOFOnEdge(int edge, ElementType type) => 0;
            public int DOFOnElement() => 0;
            public (int i, int j) Edge(int edge)
            {
                switch (edge)
                {
                    case 0:
                        return (0, 1);

                    case 1:
                        return (1, 2);

                    case 2:
                        return (3, 2);

                    case 3:
                        return (0, 3);

                    default:
                        throw new ArgumentException();
                }
            }
            public void SetEdgeDOF(int edge, int n, int dof, ElementType type) => throw new NotSupportedException();

            //public int[] GetDofsOnEdge(int edge)
            //{
            //    var dofs = new int[2];

            //    switch (edge)
            //    {   
            //        case 0:
            //            dofs[0] = Dofs[0];
            //            dofs[1] = Dofs[1];
            //            break;
            //        case 1:
            //            dofs[0] = Dofs[1];
            //            dofs[1] = Dofs[3];
            //            break;
            //        case 2:
            //            dofs[0] = Dofs[2];
            //            dofs[1] = Dofs[3];
            //            break;
            //        case 3:
            //            dofs[0] = Dofs[0];
            //            dofs[1] = Dofs[3];
            //            break;
            //        default:
            //            throw new ArgumentException();
            //    }

            //    return dofs;
            //}

            public void SetElementDOF(int n, int dof) => throw new NotSupportedException();
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
                        Dofs[3] = dof;
                        break;
                    case 3:
                        Dofs[2] = dof;
                        break;

                    default:
                        throw new ArgumentException();
                }
            }

            public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Ug) => throw new NotSupportedException();

            public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta) => throw new NotSupportedException();

            public static int Mu(int i) => i % 2;
            public static int Nu(int i) => i / 2;

            double CalcAverageCoeff(Vector3D[] VertexCoords, Func<Vector3D, double> Coeff)
            {
                double average = 0;

                for (int i = 0; i < VertexNumber.Length; ++i)
                {
                    average += Coeff(VertexCoords[VertexNumber[i]]);
                }

                return average / VertexNumber.Length;
            }

            double[] CalcLocalVectorF(Vector3D[] VertexCoords, Func<Vector3D, double> F)
            {
                double[] LF = new double[Dofs.Length];

                LF[0] = F(VertexCoords[VertexNumber[0]]);
                LF[1] = F(VertexCoords[VertexNumber[1]]);
                LF[2] = F(VertexCoords[VertexNumber[3]]);
                LF[3] = F(VertexCoords[VertexNumber[2]]);

                return LF;
            }

            public bool IsPointOnElement(Vector3D[] VertexCoords, Vector3D point)
            {
                double x0 = VertexCoords[VertexNumber[0]].X;
                double x1 = VertexCoords[VertexNumber[1]].X;
                double y0 = VertexCoords[VertexNumber[0]].Y;
                double y1 = VertexCoords[VertexNumber[3]].Y;


                if (x0 <= point.X && point.X <= x1 && y0 <= point.Y && point.Y <= y1)
                    return true;

                return false;
            }
            public double GetValueAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
            {
                int N = Dofs.Length;
                double value = 0;

                double x0 = VertexCoords[VertexNumber[0]].X;
                double y0 = VertexCoords[VertexNumber[0]].Y;

                double hx = VertexCoords[VertexNumber[1]].X - x0;
                double hy = VertexCoords[VertexNumber[3]].Y - y0;

                for (int i = 0; i < N; ++i)
                    value += coeffs[Dofs[i]] * LinearBasis.Phi[Mu(i)]((point.X - x0) / hx) * LinearBasis.Phi[Nu(i)]((point.Y - y0) / hy);

                return value;
            }
            public Vector3D GetGradientAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
            {
                int N = Dofs.Length;
                double xValue = 0;
                double yValue = 0;
                double zValue = 0;

                double x0 = VertexCoords[VertexNumber[0]].X;
                double y0 = VertexCoords[VertexNumber[0]].Y;

                double hx = VertexCoords[VertexNumber[1]].X - x0;
                double hy = VertexCoords[VertexNumber[3]].Y - y0;

                double coeff = 0;
                for (int i = 0; i < N; ++i)
                {
                    coeff = coeffs[Dofs[i]];
                    xValue += coeff * (1 / hx) * LinearBasisDerivatives.Phi[Mu(i)]((point.X - x0) / hx) * LinearBasis.Phi[Nu(i)]((point.Y - y0) / hy);
                    yValue += coeff * LinearBasis.Phi[Mu(i)]((point.X - x0) / hx) * (1 / hy) * LinearBasisDerivatives.Phi[Nu(i)]((point.Y - y0) / hy);
                }

                return new Vector3D(xValue, yValue, zValue);
            }

            public int[] Face(int face)
            {
                throw new NotImplementedException();
            }

            public int DOFOnFace(int face, ElementType type)
            {
                throw new NotImplementedException();
            }

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
