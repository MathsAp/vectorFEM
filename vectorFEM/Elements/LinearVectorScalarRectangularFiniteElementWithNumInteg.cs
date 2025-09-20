using Core.VectorFiniteElements.FiniteElements3D;
using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using static FEM.IFiniteElement;

namespace Core.VectorScalarFiniteElements.FiniteElements2D;

public class LinearVectorScalarRectangularFiniteElementWithNumInteg : IFiniteElementWithNumericIntegration<Vector2D, double, Vector3D>//, IFiniteElementWithNumericIntegration<Vector2D, Vector3D, Vector3D>
{
    public LinearVectorScalarRectangularFiniteElementWithNumInteg(string material, int[] vertexNumber)
    {
        Material = material;
        VertexNumber = vertexNumber;
        MasterElement = SquareMasterElementLinearScalarBasis.GetInstance();
        VectorMasterElement = SquareMasterElementLinearVectorBasis.GetInstance();
    }
    public IMasterElement<Vector2D, double, Vector3D> MasterElement { get; }

    public bool IsVector => true;

    public string Material { get; }

    public int[] VertexNumber { get; } = new int[4];

    public int NumberOfEdges => 4;

    public int NumberOfFaces => 1;

    public int[] Dofs { get; } = new int[8];

    public ElementType Type => ElementType.VectorScalar;

    SquareMasterElementLinearVectorBasis VectorMasterElement { get; }

    // IMasterElement<Vector2D, Vector3D, Vector3D> IFiniteElementWithNumericIntegration<Vector2D, Vector3D, Vector3D>.MasterElement => throw new NotImplementedException();

    public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, IFiniteElement.MatrixType type, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, double> Coeff)
    {
        int Ns = 4; // Количество скалярных неизвестных
        int Nv = 4; // Количество векторных неизвестных

        var LM = new double[Ns, Nv];

        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];
        Vector3D d = VertexCoords[VertexNumber[3]];

        double S = a.Distance(b) * a.Distance(c); // S - площадь

        (int faceS, int faceV) = GetFaceNumbers(FacePortrait); // faceS = номер грани у скалярного элемента, faceV - номер грани у векторного элемента.

        double LocalCoeff(Vector2D p) => Coeff(a * (1 - p.X) * (1 - p.Y) + b * p.X * (1 - p.Y) +
                                                c * (1 - p.X) * p.Y + d * p.X * p.Y);

        Vector3D n = LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(faceV);


        var nodes = MasterElement.QuadratureNodes.Nodes;

        // var vnodes = ((IFiniteElementWithNumericIntegration<Vector2D, Vector3D, Vector3D>)this).MasterElement.

        var localCoeffInNodes = CalcLocalFuncInQuadratureNodes(LocalCoeff, nodes);

        switch (type)
        {
            case MatrixType.Interface:

                double[,] J = GetJForFace(VertexCoords, faceV);

                for (int i = 0; i < Ns; ++i)
                {
                    for (int j = 0; j < Nv; ++j)
                    {
                        double sum = 0;
                        for (int k = 0; k < nodes.Length; ++k)
                        {
                            sum += nodes[k].Weight * Vector3D.Cross(LinearAlgebraAlgorithms.MultiplyMatrix3By3ByVector(J, MasterElement.FacesGradValues[faceS, i, k].AsArray()), n) * VectorMasterElement.FacesPsiValues[faceV, j, k] * localCoeffInNodes[k];
                        }

                        LM[i, j] = sum * S;
                    }
                }
                break;

            default:
                throw new ArgumentException();

        }


        return LM;
    }

    public double[,] BuildLocalMatrix(Vector3D[] VertexCoords, MatrixType type, Func<Vector3D, double> Coeff, Func<Vector3D, Vector3D>? Velocity = null)
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
        throw new NotImplementedException();
    }

    public double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta)
    {
        int Ns = 4;

        double[] LRP = new double[Ns];

        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];
        Vector3D d = VertexCoords[VertexNumber[3]];

        double S = a.Distance(b) * a.Distance(c);

        double LocalTheta(Vector2D p) => Theta(a * (1 - p.X) * (1 - p.Y) + b * p.X * (1 - p.Y) +
                                                c * (1 - p.X) * p.Y       + d * p.X * p.Y);

        var nodes = MasterElement.QuadratureNodes.Nodes;

        var localThetaInNodes = CalcLocalFuncInQuadratureNodes(LocalTheta, nodes);

        for (int i = 0; i < Ns; ++i)
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
        int Nv = 4;

        double[] LRP = new double[Nv];

        Vector3D a = VertexCoords[VertexNumber[0]];
        Vector3D b = VertexCoords[VertexNumber[1]];
        Vector3D c = VertexCoords[VertexNumber[2]];
        Vector3D d = VertexCoords[VertexNumber[3]];

        double S = a.Distance(b) * a.Distance(c);

        (int faceS, int faceV) = GetFaceNumbers(FacePortrait);

        Vector3D n = LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(faceV);

        Vector3D LocalTheta(Vector2D p) => Theta(a * (1 - p.X) * (1 - p.Y) + b * p.X * (1 - p.Y) +
                                                    c * (1 - p.X) * p.Y + d * p.X * p.Y);

        var nodes = MasterElement.QuadratureNodes.Nodes;

        var localThetaInNodes = CalcLocalFuncInQuadratureNodes(LocalTheta, n, nodes);

        for (int i = 0; i < Nv; ++i)
        {
            double sum = 0;
            for (int k = 0; k < nodes.Length; ++k)
            {
                var node = nodes[k].Node;
                sum += nodes[k].Weight * VectorMasterElement.FacesPsiValues[faceV, i, k] * localThetaInNodes[k];
            }

            LRP[i] = sum * S;
        }


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

    public int[] Face(int face) => face == 0 ? [0, 1, 2, 3] : throw new ArgumentException();

    public Vector3D GetCurlAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point)
    {
        throw new NotImplementedException();
    }

    public int[] GetDofs(DofsType type)
    {
        switch (type)
        {
            case DofsType.Scalar:
                return Dofs[..4];

            case DofsType.Vector:
                return Dofs[4..];

            default:
                throw new ArgumentException();
        }
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
        switch (edge)
        {
            case 0:
                Dofs[4] = dof;
                break;
            case 1:
                Dofs[5] = dof;
                break;
            case 2:
                Dofs[6] = dof;
                break;
            case 3:
                Dofs[7] = dof;
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

    (int, int, int, int) GetFaceTuple()
    {
        int[] face = [VertexNumber[0], VertexNumber[1], VertexNumber[2], VertexNumber[3]];

        Array.Sort(face);

        var faceTuple = (face[0], face[1], face[2], face[3]);

        return faceTuple;
    }

    public (int, int) GetFaceNumbers(IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait)
    {
        (IFiniteElement? FE1, int scalarNumber) = FacePortrait[GetFaceTuple()].Item1;
        (IFiniteElement? FE2, int vectorNumber) = FacePortrait[GetFaceTuple()].Item2;

        if (FE1!.Type == ElementType.Vector)
            (scalarNumber, vectorNumber) = (vectorNumber, scalarNumber);

        return (scalarNumber, vectorNumber);
    }

    double[,] GetJForFace(Vector3D[] VertexCoords, int face)
    {
        double j11 = 0;
        double j22 = 0;
        double j33 = 0;

        if (face == 0 || face == 1)
        {
            j22 = 1 / (VertexCoords[VertexNumber[1]].Y - VertexCoords[VertexNumber[0]].Y);
            j33 = 1 / (VertexCoords[VertexNumber[2]].Z - VertexCoords[VertexNumber[0]].Z);
        }
        else if (face == 2 || face == 3)
        {
            j11 = 1 / (VertexCoords[VertexNumber[1]].X - VertexCoords[VertexNumber[0]].X);
            j33 = 1 / (VertexCoords[VertexNumber[2]].Z - VertexCoords[VertexNumber[0]].Z);
        }
        else if (face == 4 || face == 5)
        {
            j11 = 1 / (VertexCoords[VertexNumber[1]].X - VertexCoords[VertexNumber[0]].X);
            j22 = 1 / (VertexCoords[VertexNumber[2]].Y - VertexCoords[VertexNumber[0]].Y);
        }

        double[,] J = { { j11,   0,   0 },
                    {   0, j22,   0 },
                    {   0,   0, j33 } };

        return J;
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

    public Vector3D[] CalcLocalFuncInQuadratureNodes(Func<Vector2D, Vector3D> LocalFunc, Vector3D n, Quadratures.QuadratureNode<Vector2D>[] nodes)
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
