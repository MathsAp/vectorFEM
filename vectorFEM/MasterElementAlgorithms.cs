using Core.FiniteElements2D;
using Core.ScalarFiniteElements.FiniteElements3D;
using Core;
using Quadratures;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.InteropServices;

namespace FEM
{
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

            int Mu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Mu(i);
            int Nu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Nu(i);
            int Upsilon(int i) => LinearParallelepipedalFiniteElementWithNumInteg.upsilon(i);

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    psiValues[i, j] = LinearBasis.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y) * LinearBasis.Phi[Upsilon(i)](node.Z);
                }

            return psiValues;
        }

        public static double[,] CalcScalarPsiValues(int n, QuadratureNodes<Vector2D> quadratureNodes)
        {
            double[,] psiValues = new double[n, quadratureNodes.Nodes.Length];


            int Mu(int i) => LinearRectangularFiniteElement.Mu(i);
            int Nu(int i) => LinearRectangularFiniteElement.Nu(i);

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    psiValues[i, j] = LinearBasis.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y);
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
            {
                var w = quadratureNodes.Nodes[k].Weight;

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        psiPsiMatrix[k, i, j] = w * psiValues[i, k] * psiValues[j, k];
                    }
            }

            return psiPsiMatrix;
        }

        public static double[,,] CalcScalarPsiPsiMatrix(int n, QuadratureNodes<Vector3D> quadratureNodes, double[,] psiValues)
        {
            double[,,] psiPsiMatrix = new double[quadratureNodes.Nodes.Length, n, n];

            for (int k = 0; k < quadratureNodes.Nodes.Length; ++k)
            {
                var w = quadratureNodes.Nodes[k].Weight;

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        psiPsiMatrix[k, i, j] = w * psiValues[i, k] * psiValues[j, k];
                    }
            }

            return psiPsiMatrix;
        }

        public static double[,,] CalcScalarPsiPsiMatrix(int n, QuadratureNodes<Vector2D> quadratureNodes, double[,] psiValues)
        {
            double[,,] psiPsiMatrix = new double[quadratureNodes.Nodes.Length, n, n];

            for (int k = 0; k < quadratureNodes.Nodes.Length; ++k)
            {
                var w = quadratureNodes.Nodes[k].Weight;

                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                    {
                        psiPsiMatrix[k, i, j] = w * psiValues[i, k] * psiValues[j, k];
                    }
            }

            return psiPsiMatrix;
        }

        public static Vector3D[,] CalcGradValues(int n, QuadratureNodes<Vector3D> quadratureNodes)
        {
            Vector3D[,] gradValues = new Vector3D[n, quadratureNodes.Nodes.Length];

            int Mu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Mu(i);
            int Nu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Nu(i);
            int Upsilon(int i) => LinearParallelepipedalFiniteElementWithNumInteg.upsilon(i);

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    gradValues[i, j] = new Vector3D(LinearBasisDerivatives.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y) * LinearBasis.Phi[Upsilon(i)](node.Z),
                                                    LinearBasis.Phi[Mu(i)](node.X) * LinearBasisDerivatives.Phi[Nu(i)](node.Y) * LinearBasis.Phi[Upsilon(i)](node.Z),
                                                    LinearBasis.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y) * LinearBasisDerivatives.Phi[Upsilon(i)](node.Z));
                }

            return gradValues;
        }

        public static Vector2D[,] CalcGradValues(int n, QuadratureNodes<Vector2D> quadratureNodes)
        {
            Vector2D[,] gradValues = new Vector2D[n, quadratureNodes.Nodes.Length];

            int Mu(int i) => LinearRectangularFiniteElement.Mu(i);
            int Nu(int i) => LinearRectangularFiniteElement.Nu(i);

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    gradValues[i, j] = new Vector2D(LinearBasisDerivatives.Phi[Mu(i)](node.X) * LinearBasis.Phi[Nu(i)](node.Y),
                                                    LinearBasis.Phi[Mu(i)](node.X) * LinearBasisDerivatives.Phi[Nu(i)](node.Y));
                }

            return gradValues;
        }

    }
}
