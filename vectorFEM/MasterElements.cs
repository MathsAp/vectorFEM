using Core;
using Quadratures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices.Marshalling;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    public class CubeMasterElementLinearVectorBasis : IMasterElement<Vector3D, Vector3D>
    {
        private static CubeMasterElementLinearVectorBasis? Instance;

        public Vector3D[,] PsiValues { get; }

        public double[,,] PsiPsiMatrix { get; }

        public Vector3D[,] CurlValues { get; }

        public Vector3D[,] GradValues => throw new NotImplementedException();

        public QuadratureNodes<Vector3D> QuadratureNodes {  get; }

        private CubeMasterElementLinearVectorBasis()
        {
            int order = 5;
            int n = 12;
            QuadratureNodes = NumericalIntegration.FactoryQuadratures3D(order, ElemType.Cube);
            PsiValues = MasterElementAlgorithms.CalcVectorPsiValues(n, QuadratureNodes);
            CurlValues = MasterElementAlgorithms.CalcCurlValues(n, QuadratureNodes);
            PsiPsiMatrix = MasterElementAlgorithms.CalcVectorPsiPsiMatrix(n, QuadratureNodes, PsiValues);
        }


        public static CubeMasterElementLinearVectorBasis GetInstance()
        {
            if (Instance == null)
                Instance = new CubeMasterElementLinearVectorBasis();

            return Instance;
        }
    }

    public class CubeMasterElementLinearScalarBasis : IMasterElement<Vector3D, double>
    {
        private static CubeMasterElementLinearScalarBasis? Instance;
        public double[,] PsiValues { get; }

        public double[,,] PsiPsiMatrix { get; }

        public double[,] CurlValues => throw new NotImplementedException();

        public Vector3D[,] GradValues { get; }

        public QuadratureNodes<Vector3D> QuadratureNodes { get; }

        private CubeMasterElementLinearScalarBasis()
        {
            int n = 8;
            int order = 5;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures3D(order, ElemType.Cube);
            PsiValues = MasterElementAlgorithms.CalcScalarPsiValues(n, QuadratureNodes);
            PsiPsiMatrix = MasterElementAlgorithms.CalcScalarPsiPsiMatrix(n, QuadratureNodes, PsiValues);
            GradValues = MasterElementAlgorithms.CalcGradValues(n, QuadratureNodes);
        }

        public static CubeMasterElementLinearScalarBasis GetInstance()
        {
            if (Instance == null)
                Instance = new CubeMasterElementLinearScalarBasis();

            return Instance;
        }
    }

    public class SquareMasterElementLinearScalarBasis : IMasterElement<Vector2D, double>
    {
        private static SquareMasterElementLinearScalarBasis? Instance;

        public double[,] PsiValues { get; }

        public double[,,] PsiPsiMatrix { get; }

        public double[,] CurlValues => throw new NotImplementedException();

        public Vector2D[,] GradValues { get; }

        public QuadratureNodes<Vector2D> QuadratureNodes { get; }

        private SquareMasterElementLinearScalarBasis()
        {
            int n = 4;
            int order = 5;

            QuadratureNodes = NumericalIntegration.FactoryQuadratures2D(order, ElemType.Rectangle);
            PsiValues = MasterElementAlgorithms.CalcScalarPsiValues(n, QuadratureNodes);
            PsiPsiMatrix = MasterElementAlgorithms.CalcScalarPsiPsiMatrix(n, QuadratureNodes, PsiValues);
            GradValues = MasterElementAlgorithms.CalcGradValues(n, QuadratureNodes);
        }

        public static SquareMasterElementLinearScalarBasis GetInstance()
        {
            if (Instance == null)
                Instance = new SquareMasterElementLinearScalarBasis();

            return Instance;
        }
    }
}
