using Core;
using Quadratures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace FEM
{
    public class CubeMasterElementLinearVectorBasis : IMasterElement<Vector3D, Vector3D, Vector3D>
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

    public class CubeMasterElementLinearScalarBasis : IMasterElement<Vector3D, double, Vector3D>
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

    public static class MasterElementAlgorithms
    {
        public static Vector3D[,] CalcVectorPsiValues(int n, QuadratureNodes<Vector3D> quadratureNodes)
        {
            Vector3D[,] psiValues = new Vector3D[n, quadratureNodes.Nodes.Length];

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    psiValues[i, j] = LinearVectorBasis.Phi[i](node.X, node.Y, node.Z);
                }

            return psiValues;
        }

        public static double[,] CalcScalarPsiValues(int n, QuadratureNodes<Vector3D> quadratureNodes)
        {
            double[,] psiValues = new double[n, quadratureNodes.Nodes.Length];

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    psiValues[i, j] = LinearBasis.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y) * LinearBasis.Phi[upsilon(i)](node.Z);
                }

            return psiValues;
        }

        public static Vector3D[,] CalcCurlValues(int n, QuadratureNodes<Vector3D> quadratureNodes)
        {
            Vector3D[,] curlValues = new Vector3D[n, quadratureNodes.Nodes.Length];

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    curlValues[i, j] = LinearVectorBasis.curlPhi[i](node.X, node.Y, node.Z);
                }

            return curlValues;
        }

        public static double[,,] CalcVectorPsiPsiMatrix(int n, QuadratureNodes<Vector3D> quadratureNodes, Vector3D[,] psiValues)
        {
            double[,,] psiPsiMatrix = new double[quadratureNodes.Nodes.Length, n, n];

            for (int k = 0; k < quadratureNodes.Nodes.Length; ++k)
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        psiPsiMatrix[k, i, j] = psiValues[i, k] * psiValues[j, k];
                    }

            return psiPsiMatrix;
        }

        public static double[,,] CalcScalarPsiPsiMatrix(int n, QuadratureNodes<Vector3D> quadratureNodes, double[,] psiValues)
        {
            double[,,] psiPsiMatrix = new double[quadratureNodes.Nodes.Length, n, n];

            for (int k = 0; k < quadratureNodes.Nodes.Length; ++k)
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        psiPsiMatrix[k, i, j] = psiValues[i, k] * psiValues[j, k];
                    }

            return psiPsiMatrix;
        }

        public static Vector3D[,] CalcGradValues(int n, QuadratureNodes<Vector3D> quadratureNodes)
        {
            Vector3D[,] gradValues = new Vector3D[n, quadratureNodes.Nodes.Length];

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    gradValues[i, j] = new Vector3D(LinearBasisDerivatives.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y) * LinearBasis.Phi[upsilon(i)](node.Z),
                                                    LinearBasis.Phi[Mu(i)](node.X) * LinearBasisDerivatives.Phi[Nu(i)](node.Y) * LinearBasis.Phi[upsilon(i)](node.Z),
                                                    LinearBasis.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y) * LinearBasisDerivatives.Phi[upsilon(i)](node.Z));
                }

            return gradValues;
        }

        public static int Mu(int i) => i % 2;
        public static int Nu(int i) => (i / 2) % 2;
        public static int upsilon(int i) => i / 4;
    }
}
