using Core;
using FEM;
using Quadratures;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices.Marshalling;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    public class CubeMasterElementLinearVectorBasis : IMasterElement<Vector3D, Vector3D, Vector3D>
    {
        private CubeMasterElementLinearVectorBasis()
        {
            int order = 5;
            int n = 12;
            QuadratureNodes = NumericalIntegration.FactoryQuadratures3D(order, ElemType.Cube);
            PsiValues = MasterElementAlgorithms.CalcVectorPsiValues(n, QuadratureNodes);
            CurlValues = MasterElementAlgorithms.CalcCurlValues(n, QuadratureNodes);
            PsiPsiMatrix = MasterElementAlgorithms.CalcVectorPsiPsiMatrix(n, QuadratureNodes, PsiValues);
        }

        private static Lazy<CubeMasterElementLinearVectorBasis> Instance { get; } = new(() => new CubeMasterElementLinearVectorBasis());

        public Vector3D[,] PsiValues { get; }

        public double[,,] PsiPsiMatrix { get; }

        public Vector3D[,] CurlValues { get; }

        public Vector3D[,] GradValues => throw new NotImplementedException();

        public QuadratureNodes<Vector3D> QuadratureNodes {  get; }

        public Vector3D[,,] FacesPsiValues => throw new NotImplementedException();

        public Vector3D[,,] FacesPsiNValues => throw new NotImplementedException();

        public double[,,,] FacesPsiNPsiNMatrix => throw new NotImplementedException();

        public Vector3D[,,] FacesGradValues => throw new NotImplementedException();

        public static CubeMasterElementLinearVectorBasis GetInstance() => Instance.Value;
    }

    public class CubeMasterElementLinearScalarBasis : IMasterElement<Vector3D, double, double>
    {
        private CubeMasterElementLinearScalarBasis()
        {
            int n = 8;
            int order = 5;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures3D(order, ElemType.Cube);
            PsiValues = MasterElementAlgorithms.CalcScalarPsiValues(n, QuadratureNodes);
            PsiPsiMatrix = MasterElementAlgorithms.CalcScalarPsiPsiMatrix(n, QuadratureNodes, PsiValues);
            GradValues = MasterElementAlgorithms.CalcGradValues(n, QuadratureNodes);
        }

        private static Lazy<CubeMasterElementLinearScalarBasis> Instance { get; } = new(() => new CubeMasterElementLinearScalarBasis());
        public double[,] PsiValues { get; }

        public double[,,] PsiPsiMatrix { get; }

        public double[,] CurlValues => throw new NotImplementedException();

        public Vector3D[,] GradValues { get; }

        public QuadratureNodes<Vector3D> QuadratureNodes { get; }

        public double[,,] FacesPsiValues => throw new NotImplementedException();

        public double[,,] FacesPsiNValues => throw new NotImplementedException();

        public double[,,,] FacesPsiNPsiNMatrix => throw new NotImplementedException();

        public double[,,] FacesGradValues => throw new NotImplementedException();

        public static CubeMasterElementLinearScalarBasis GetInstance() => Instance.Value;
    }

    public class SquareMasterElementLinearVectorBasis : IMasterElement<Vector2D, Vector3D, Vector3D>
    {
        private SquareMasterElementLinearVectorBasis()
        {
            int n = 4;
            int order = 5;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures2D(order, ElemType.Rectangle);
            FacesPsiValues = MasterElementAlgorithms.CalcVectorFacesPsiValues(n, QuadratureNodes);
            FacesPsiNValues = MasterElementAlgorithms.CalcVectorFacesPsiNValues(n, QuadratureNodes);
            FacesPsiNPsiNMatrix = MasterElementAlgorithms.CalcVectorFacesPsiNPsiNMatrix(n, QuadratureNodes, FacesPsiNValues);
        }
        private static Lazy<SquareMasterElementLinearVectorBasis> Instance { get; } = new (() => new SquareMasterElementLinearVectorBasis());
      
        public Vector3D[,] PsiValues => throw new NotImplementedException();

        public Vector3D[,,] FacesPsiValues { get; }
        
        public Vector3D[,,] FacesPsiNValues {  get; }

        public double[,,] PsiPsiMatrix => throw new NotImplementedException();

        public double[,,,] FacesPsiNPsiNMatrix {  get; }

        public Vector3D[,] CurlValues => throw new NotImplementedException();

        public Vector2D[,] GradValues => throw new NotImplementedException();

        public QuadratureNodes<Vector2D> QuadratureNodes { get; }

        public Vector3D[,,] FacesGradValues => throw new NotImplementedException();

        public static SquareMasterElementLinearVectorBasis GetInstance() => Instance.Value;
    }

    public class SquareMasterElementLinearScalarBasis : IMasterElement<Vector2D, double, Vector3D>
    {
        private SquareMasterElementLinearScalarBasis()
        {
            int n = 4;
            int order = 5;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures2D(order, ElemType.Rectangle);
            PsiValues = MasterElementAlgorithms.CalcScalarPsiValues(n, QuadratureNodes);
            PsiPsiMatrix = MasterElementAlgorithms.CalcScalarPsiPsiMatrix(n, QuadratureNodes, PsiValues);
            GradValues = MasterElementAlgorithms.CalcGradValues(n, QuadratureNodes);
            FacesGradValues = MasterElementAlgorithms.CalcFacesGradValues(n, QuadratureNodes);
        }
        private static Lazy<SquareMasterElementLinearScalarBasis> Instance { get; } = new (() => new SquareMasterElementLinearScalarBasis());

        public double[,] PsiValues { get; }

        public double[,,] PsiPsiMatrix { get; }

        public double[,] CurlValues => throw new NotImplementedException();

        public Vector2D[,] GradValues { get; }

        public QuadratureNodes<Vector2D> QuadratureNodes { get; }

        public double[,,] FacesPsiValues => throw new NotImplementedException();

        public double[,,] FacesPsiNValues => throw new NotImplementedException();

        public double[,,,] FacesPsiNPsiNMatrix => throw new NotImplementedException();

        public Vector3D[,,] FacesGradValues { get; }

        public static SquareMasterElementLinearScalarBasis GetInstance() => Instance.Value;
    }
    public class LineSegmentMasterElementLinearScalarBasis : IMasterElement<double, double, double>
    {
        private LineSegmentMasterElementLinearScalarBasis()
        {
            int n = 2;
            int order = 5;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures1D(order, ElemType.StraightLine);
            PsiValues = MasterElementAlgorithms.CalcScalarPsiValues(n, QuadratureNodes);
            PsiPsiMatrix = MasterElementAlgorithms.CalcScalarPsiPsiMatrix(n, QuadratureNodes, PsiValues);
            GradValues = MasterElementAlgorithms.CalcGradValues(n, QuadratureNodes);
        }
        private static Lazy<LineSegmentMasterElementLinearScalarBasis> Instance { get; } = new(() => new LineSegmentMasterElementLinearScalarBasis());

        public double[,] PsiValues { get; }

        public double[,,] FacesPsiValues => throw new NotImplementedException();

        public double[,,] FacesPsiNValues => throw new NotImplementedException();

        public double[,,] FacesGradValues => throw new NotImplementedException();

        public double[,,] PsiPsiMatrix { get; }

        public double[,,,] FacesPsiNPsiNMatrix => throw new NotImplementedException();

        public double[,] CurlValues => throw new NotImplementedException();

        public double[,] GradValues { get; }

        public QuadratureNodes<double> QuadratureNodes { get; }

        public static LineSegmentMasterElementLinearScalarBasis GetInstance() => Instance.Value;
    }

    public class LineSegmentMasterElementQuadraticScalarBasis : IMasterElement<double, double, double>
    {
        private LineSegmentMasterElementQuadraticScalarBasis()
        {
            int n = 3;
            int order = 7;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures1D(order, ElemType.StraightLine);
            PsiValues = MasterElementAlgorithms.CalcScalarPsiValues(n, QuadratureNodes, true);
            PsiPsiMatrix = MasterElementAlgorithms.CalcScalarPsiPsiMatrix(n, QuadratureNodes, PsiValues);
            GradValues = MasterElementAlgorithms.CalcGradValues(n, QuadratureNodes, true);
        }

        private static Lazy<LineSegmentMasterElementQuadraticScalarBasis> Instance { get; } = new(() => new LineSegmentMasterElementQuadraticScalarBasis());

        public double[,] PsiValues { get; }

        public double[,,] FacesPsiValues => throw new NotImplementedException();

        public double[,,] FacesPsiNValues => throw new NotImplementedException();

        public double[,,] FacesGradValues => throw new NotImplementedException();

        public double[,,] PsiPsiMatrix { get; }

        public double[,,,] FacesPsiNPsiNMatrix => throw new NotImplementedException();

        public double[,] CurlValues => throw new NotImplementedException();

        public double[,] GradValues { get; }

        public QuadratureNodes<double> QuadratureNodes { get; }

        public static LineSegmentMasterElementQuadraticScalarBasis GetInstance() => Instance.Value;
    }

    public class TriangularMasterElementQuadraticScalarBasis : IMasterElement<Vector2D, double, double>
    {
        private TriangularMasterElementQuadraticScalarBasis()
        {
            int n = 6;
            int order = 6;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures2D(order, ElemType.Triangle);
            PsiValues = MasterElementAlgorithms.CalcScalarPsiValues(n, QuadratureNodes, true);
            PsiPsiMatrix = MasterElementAlgorithms.CalcScalarPsiPsiMatrix(n, QuadratureNodes, PsiValues);
            GradValues = MasterElementAlgorithms.CalcGradValues(n, QuadratureNodes, true);
        }

        private static Lazy<TriangularMasterElementQuadraticScalarBasis> Instance { get; } = new(() => new TriangularMasterElementQuadraticScalarBasis());

        public double[,] PsiValues { get; }

        public double[,,] FacesPsiValues => throw new NotImplementedException();

        public double[,,] FacesPsiNValues => throw new NotImplementedException();

        public double[,,] FacesGradValues => throw new NotImplementedException();

        public double[,,] PsiPsiMatrix { get; }

        public double[,,,] FacesPsiNPsiNMatrix => throw new NotImplementedException();

        public double[,] CurlValues => throw new NotImplementedException();

        public Vector2D[,] GradValues { get; }

        public QuadratureNodes<Vector2D> QuadratureNodes { get; }

        public static TriangularMasterElementQuadraticScalarBasis GetInstance() => Instance.Value;
    }
}