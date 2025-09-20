using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static FEM.IFiniteElement;

namespace Core.ScalarFiniteElements.FiniteElements2D;

public class QuadraticTriangularFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<Vector2D, double, double>
{
    public QuadraticTriangularFiniteElementWithNumInteg(string material, int[] vertexNumbers)
    {
        Material = material;
        VertexNumber = vertexNumbers;
        MasterElement = TriangularMasterElementQuadraticScalarBasis.GetInstance();
    }

    public IMasterElement<Vector2D, double, double> MasterElement { get; }

    public ElementType Type => ElementType.Scalar;

    public string Material { get; }

    public int[] VertexNumber { get; }

    public int NumberOfEdges => 3;

    public int NumberOfFaces => 1;

    public int[] Dofs { get; } = new int[6];

    public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, Func<Vector3D, double> Coeff, Func<Vector3D, Vector3D>? Velocity = null)
    {
        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];

        double detJ = (c.X - a.X) * (b.Y - a.Y) - (b.X - a.X) * (c.Y - a.Y);
        double detJAbs = Math.Abs(detJ);
        double[,] Jr = { { (b.Y - a.Y) / detJ, (a.X - b.X) / detJ },
                         { (a.Y - c.Y) / detJ, (c.X - a.X) / detJ } };

        double LocalCoeff(Vector2D p) => Coeff(a * (1 - p.X - p.Y) + b * p.Y + c * p.X);

        var nodes = MasterElement.QuadratureNodes.Nodes;
        var localCoeffInNodes = CalcLocalFuncInQuadratureNodes(LocalCoeff, nodes);

        int N = Dofs.Length;
        double[,] LM = new double[N, N];

        switch (type)
        {
            case MatrixType.Stiffness:

                double[,] JTJ = LinearAlgebraAlgorithms.GetMatrixMultipliedByTransposedMatrix(Jr);

                for (int i = 0; i < N; ++i)
                {
                    for (int j = 0; j < N; ++j)
                    {
                        double sum = 0;
                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            sum += nodes[k].Weight * MasterElement.GradValues[i, k]
                                * LinearAlgebraAlgorithms.MultiplyMatrix2By2ByVector(JTJ, MasterElement.GradValues[j, k])
                                * localCoeffInNodes[k];
                        }

                        LM[i, j] = sum * detJAbs;
                    }
                }

                break;

            case MatrixType.Mass:

                for (int i = 0; i < N; ++i)
                {
                    for (int j = 0; j < N; ++j)
                    {
                        double sum = 0;
                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            sum += localCoeffInNodes[k] * MasterElement.PsiPsiMatrix[k, i, j];
                        }

                        LM[i, j] = sum * detJAbs;
                    }
                }

                break;

            case MatrixType.Convection:

                Vector2D LocalVelocity(Vector2D p)
                {
                    if (Velocity is null) return Vector2D.Zero;

                    var res = Velocity(a * (1 - p.X - p.Y) + b * p.Y + c * p.X);

                    return new Vector2D(res.X, res.Y);
                }

                var localVelocityInNodes = CalcLocalFuncInQuadratureNodes(LocalVelocity, nodes);

                for (int i = 0; i < N; ++i)
                {
                    for (int j = 0; j < N; ++j)
                    {
                        double sum = 0;
                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            sum += nodes[k].Weight * localVelocityInNodes[k]
                                * LinearAlgebraAlgorithms.MultiplyMatrix2By2ByVector(Jr, MasterElement.GradValues[j, k], transpose: true)
                                * MasterElement.PsiValues[i, k]
                                * localCoeffInNodes[k];
                        }

                        LM[i, j] = sum * detJAbs;
                    }
                }

                break;
        }

        return LM;
    }

    public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, double> F)
    {
        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];

        double detJ = (c.X - a.X) * (b.Y - a.Y) - (b.X - a.X) * (c.Y - a.Y);
        double detJAbs = Math.Abs(detJ);

        double LocalF(Vector2D p) => F(a * (1 - p.X - p.Y) + b * p.Y + c * p.X);

        var nodes = MasterElement.QuadratureNodes.Nodes;
        var localFInNodes = CalcLocalFuncInQuadratureNodes(LocalF, nodes);

        int N = Dofs.Length;
        double[] LRP = new double[N];
        for (int i = 0; i < N; ++i)
        {
            double sum = 0;
            for (int k = 0; k < nodes.Length; ++k)
            {
                sum += nodes[k].Weight * MasterElement.PsiValues[i, k] * localFInNodes[k];
            }

            LRP[i] = sum * detJAbs;
        }

        return LRP;
    }

    public int DOFOnEdge(int edge, IFiniteElement.ElementType type) => 1;

    public int DOFOnElement() => 0;

    public int DOFOnFace(int face, IFiniteElement.ElementType type) => 0;

    public (int i, int j) Edge(int edge)
    {
        switch(edge)
        {
            case 0: return (0, 1);
            case 1: return (1, 2);
            case 2: return (0, 2);

            default: throw new ArgumentException();
        }
    }

    public int[] Face(int face) => face == 0 ? [0, 1, 2] : throw new ArgumentException();

    public Vector3D GetGradientAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
    {
        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];

        double detJ = (c.X - a.X) * (b.Y - a.Y) - (b.X - a.X) * (c.Y - a.Y);
        double[,] Jr = { { (b.Y - a.Y) / detJ, (a.X - b.X) / detJ },
                         { (a.Y - c.Y) / detJ, (c.X - a.X) / detJ } };

        double xi = ((b.Y - a.Y) * point.X + (a.X - b.X) * point.Y + b.X * a.Y - a.X * b.Y) / detJ;
        double eta = ((a.Y - c.Y) * point.X + (c.X - a.X) * point.Y + a.X * c.Y - c.X * a.Y) / detJ;

        Vector2D lp = new(xi, eta);

        int N = Dofs.Length;
        Vector2D value = Vector2D.Zero;
        for (int i = 0; i < N; ++i)
        {
            value += coeffs[Dofs[i]] * LinearAlgebraAlgorithms.MultiplyMatrix2By2ByVector(Jr, TrianglesQuadraticBasisGradients.Phi[i](lp), transpose: true);
        }

        return value.As3D();
    }

    public double GetValueAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
    {
        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];

        double detJ = (c.X - a.X) * (b.Y - a.Y) - (b.X - a.X) * (c.Y - a.Y);

        double xi = ((b.Y - a.Y) * point.X + (a.X - b.X) * point.Y + b.X * a.Y - a.X * b.Y) / detJ;
        double eta = ((a.Y - c.Y) * point.X + (c.X - a.X) * point.Y + a.X * c.Y - c.X * a.Y) / detJ;

        Vector2D lp = new(xi, eta);

        int N = Dofs.Length;
        double value = 0;
        for (int i = 0; i < N; ++i)
        {
            value += coeffs[Dofs[i]] * TrianglesQuadraticBasis.Phi[i](lp);
        }

        return value;
    }

    public bool IsPointOnElement(Vector3D[] VertexCoords, Vector3D point)
    {
        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];

        Vector3D ba = b - a;
        Vector3D ca = c - a;

        Vector3D n = Vector3D.Cross(ca, ba).Normalize();


        if (Math.Abs(n * (point - a)) > Constants.GeometryEps) return false;

        double c1 = Vector3D.Cross(ca, point - a) * n;
        double c2 = Vector3D.Cross(b - c, point - c) * n;
        double c3 = Vector3D.Cross(a - b, point - b) * n;

        if (c1 >= 0 && c2 >= 0 && c3 >= 0)
            return true;

        return false;
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

            default: throw new ArgumentException();
        }
    }

    public void SetEdgeDOF(int edge, int n, int dof, IFiniteElement.ElementType type)
    {
        switch (edge)
        {
            case 0:
                Dofs[3] = dof;
                break;
            case 1:
                Dofs[4] = dof;
                break;
            case 2:
                Dofs[5] = dof;
                break;

            default: throw new ArgumentException();
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

    public Vector2D[] CalcLocalFuncInQuadratureNodes(Func<Vector2D, Vector2D> LocalFunc, Quadratures.QuadratureNode<Vector2D>[] nodes)
    {
        Vector2D[] values = new Vector2D[nodes.Length];

        for (int k = 0; k < nodes.Length; ++k)
        {
            values[k] = LocalFunc(nodes[k].Node);
        }

        return values;
    }

    public void SetElementDOF(int n, int dof)
    {
        throw new NotImplementedException();
    }

    public void SetFaceDOF(int face, int n, int dof, IFiniteElement.ElementType type)
    {
        throw new NotImplementedException();
    }

    public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, double> Coeff)
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

    public double CalclIntegralOfSquaredDifference(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Func<Vector3D, Vector3D> u) // 
    {
        throw new NotImplementedException();
    }

    public Vector3D GetCurlAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
    {
        throw new NotImplementedException();
    }

    public int[] GetDofs(IFiniteElement.DofsType type)
    {
        throw new NotImplementedException();
    }

    public Vector3D GetVectorAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
    {
        throw new NotImplementedException();
    }

    public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Ug) // /
    {
        throw new NotImplementedException();
    }

    public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta) // /
    {
        throw new NotImplementedException();
    }
}

