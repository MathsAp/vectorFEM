using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static FEM.IFiniteElement;

namespace Core.ScalarFiniteElements.FiniteElements1D;

public class LinearStraightFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<double, double, double>
{
    public LinearStraightFiniteElementWithNumInteg(string material, int[] vertexNumber)
    {
        Material = material;
        VertexNumber = vertexNumber;
        MasterElement = LineSegmentMasterElementLinearScalarBasis.GetInstance();
    }

    public IMasterElement<double, double, double> MasterElement { get; }

    public ElementType Type => ElementType.Scalar;

    public string Material { get; }

    public int[] VertexNumber { get; }

    public int NumberOfEdges => 0;

    public int NumberOfFaces => 0;

    public int[] Dofs { get; } = new int[2];

    public int DOFOnEdge(int edge, IFiniteElement.ElementType type) => 0;

    public int DOFOnFace(int face, IFiniteElement.ElementType type) => 0;

    public int DOFOnElement() => 0;

    public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, Func<Vector3D, double> Coeff, Func<Vector3D, Vector3D>? Velocity = null)
    {
        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];

        double h = a.Distance(b);

        int N = Dofs.Length;
        double[,] LM = new double[N, N];

        double LocalCoeff(double p) => Coeff(a * (1 - p) + b * p);
        var nodes = MasterElement.QuadratureNodes.Nodes;
        var localCoeffInNodes = CalcLocalFuncInQuadratureNodes(LocalCoeff, nodes);

        switch (type)
        {
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

                        LM[i, j] = sum * h;
                    }
                }

                break;

            default:
                throw new ArgumentException();
        }

        return LM;
    }

    public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, double> F) //
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

    public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta)
    {
        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];

        double h = a.Distance(b);

        double LocalTheta(double p) => Theta(a * (1 - p) + b * p);

        var nodes = MasterElement.QuadratureNodes.Nodes;
        var localThetaInNodes = CalcLocalFuncInQuadratureNodes(LocalTheta, nodes);

        int N = Dofs.Length;
        double[] LRP = new double[N];
        for (int i = 0; i < N; ++i)
        {
            double sum = 0;
            for (int k = 0; k < nodes.Length; ++k)
            {
                sum += nodes[k].Weight * MasterElement.PsiValues[i, k] * localThetaInNodes[k];
            }

            LRP[i] = sum * h;
        }

        return LRP;
    }    

    public Vector3D GetGradientAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point) //
    {
        throw new NotImplementedException();
    }

    public double GetValueAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point) //
    {
        throw new NotImplementedException();
    }

    public bool IsPointOnElement(Vector3D[] VertexCoords, Vector3D point) //
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
            default:
                throw new ArgumentException();
        }
    }

    public double[] CalcLocalFuncInQuadratureNodes(Func<double, double> LocalFunc, Quadratures.QuadratureNode<double>[] nodes)
    {
        double[] values = new double[nodes.Length];

        for (int k = 0; k < nodes.Length; ++k)
        {
            values[k] = LocalFunc(nodes[k].Node);
        }

        return values;
    }

    public (int i, int j) Edge(int edge) //
    {
        throw new NotImplementedException();
    }

    public int[] Face(int face) //
    {
        throw new NotImplementedException();
    }    

    public Vector3D GetVectorAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point) //
    {
        throw new NotImplementedException();
    }

    public void SetEdgeDOF(int edge, int n, int dof, IFiniteElement.ElementType type) //
    {
        throw new NotImplementedException();
    }

    public void SetElementDOF(int n, int dof) //
    {
        throw new NotImplementedException();
    }

    public void SetFaceDOF(int face, int n, int dof, IFiniteElement.ElementType type) //
    {
        throw new NotImplementedException();
    }

    public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, double> Coeff) //
    {
        throw new NotImplementedException();
    }

    public double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, Vector3D> F) //
    {
        throw new NotImplementedException();
    }

    public double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Ug) //
    {
        throw new NotImplementedException();
    }

    public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Theta) //
    {
        throw new NotImplementedException();
    }

    public double CalclIntegralOfSquaredDifference(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Func<Vector3D, Vector3D> u) // вернуться
    {
        throw new NotImplementedException();
    }

    public Vector3D GetCurlAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point) //
    {
        throw new NotImplementedException();
    }

    public int[] GetDofs(IFiniteElement.DofsType type) // ?
    {
        throw new NotImplementedException();
    }
}

