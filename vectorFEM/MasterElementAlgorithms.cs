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
using Core.VectorFiniteElements.FiniteElements3D;

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
                    psiValues[i, j] = Linear3DVectorBasis.Phi[i](node.X, node.Y, node.Z);
                }

            return psiValues;
        }

        public static Vector2D[,] CalcVectorPsiValues(int n, QuadratureNodes<Vector2D> quadratureNodes)
        {
            Vector2D[,] psiValues = new Vector2D[n, quadratureNodes.Nodes.Length];


            for(int i = 0; i < n; ++i)
                for(int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                {
                    var node = quadratureNodes.Nodes[j].Node;
                    psiValues[i, j] = Linear2DVectorBasis.Phi[i](node.X, node.Y);
                }

            return psiValues;
        }

        public static Vector3D[,,] CalcVectorFacesPsiValues(int n, QuadratureNodes<Vector2D> quadratureNodes)
        {
            Vector3D[,,] psiValues = new Vector3D[6, n, quadratureNodes.Nodes.Length];

            Vector3D[] CubeVertex = [new Vector3D(0, 0, 0), new Vector3D(1, 0, 0), new Vector3D(0, 1, 0), new Vector3D(1, 1, 0),
                                     new Vector3D(0, 0, 1), new Vector3D(1, 0, 1), new Vector3D(0, 1, 1), new Vector3D(1, 1, 1)];

            for (int f = 0; f < 6; ++f)
            {
                var face = LinearVectorParallelepipedalFiniteElementWithNumInteg.FaceS(f);

                Vector3D a = CubeVertex[face[0]];
                Vector3D b = CubeVertex[face[1]];
                Vector3D c = CubeVertex[face[2]];
                Vector3D d = CubeVertex[face[3]];

                var localDofs = LinearVectorParallelepipedalFiniteElementWithNumInteg.LocalEdgesDofsAtFace(f);

                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                    {
                        var node = quadratureNodes.Nodes[j].Node;

                        var newNode = a * (1 - node.X) * (1 - node.Y) + b * node.X * (1 - node.Y) +
                                      c * (1 - node.X) * node.Y       + d * node.X * node.Y;

                        psiValues[f, i, j] = Linear3DVectorBasis.Phi[localDofs[i]](newNode.X, newNode.Y, newNode.Z);
                    }
                }
            }

            return psiValues;
        }

        public static Vector3D[,,] CalcVectorFacesPsiNValues(int n, QuadratureNodes<Vector2D> quadratureNodes)
        {
            Vector3D[,,] psiNValues = new Vector3D[6, n, quadratureNodes.Nodes.Length];

            Vector3D[] CubeVertex = [new Vector3D(0, 0, 0), new Vector3D(1, 0, 0), new Vector3D(0, 1, 0), new Vector3D(1, 1, 0),
                                     new Vector3D(0, 0, 1), new Vector3D(1, 0, 1), new Vector3D(0, 1, 1), new Vector3D(1, 1, 1)];

            for (int f = 0; f < 6; ++f)
            {
                var face = LinearVectorParallelepipedalFiniteElementWithNumInteg.FaceS(f);

                Vector3D a = CubeVertex[face[0]];
                Vector3D b = CubeVertex[face[1]];
                Vector3D c = CubeVertex[face[2]];
                Vector3D d = CubeVertex[face[3]];

                var normal = LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(f);

                var localDofs = LinearVectorParallelepipedalFiniteElementWithNumInteg.LocalEdgesDofsAtFace(f);


                for (int i = 0; i < n; ++i)
                {
                    for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
                    {
                        var node = quadratureNodes.Nodes[j].Node;

                        var newNode = a * (1 - node.X) * (1 - node.Y) + b * node.X * (1 - node.Y) +
                                      c * (1 - node.X) * node.Y + d * node.X * node.Y;

                        psiNValues[f, i, j] = Vector3D.Cross(Linear3DVectorBasis.Phi[localDofs[i]](newNode.X, newNode.Y, newNode.Z), normal);
                    }
                }
            }

            return psiNValues;
        }

        public static double[,] CalcScalarPsiValues(int n, QuadratureNodes<Vector3D> quadratureNodes)
        {
            double[,] psiValues = new double[n, quadratureNodes.Nodes.Length];

            int Mu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Mu(i);
            int Nu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Nu(i);
            int Upsilon(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Upsilon(i);

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
                    curlValues[i, j] = Linear3DVectorBasis.curlPhi[i](node.X, node.Y, node.Z);
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

        public static double[,,] CalcVectorPsiPsiMatrix(int n, QuadratureNodes<Vector2D> quadratureNodes, Vector3D[,] psiValues)
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

        public static double[,,] CalcVectorPsiPsiMatrix(int n, QuadratureNodes<Vector2D> quadratureNodes, Vector2D[,] psiValues)
        {
            double[,,] psiPsiMatrix = new double[quadratureNodes.Nodes.Length, n, n];

            for(int k = 0; k < quadratureNodes.Nodes.Length; ++k)
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

        public static double[,,,] CalcVectorFacesPsiNPsiNMatrix(int n, QuadratureNodes<Vector2D> quadratureNodes, Vector3D[,,] psiNValues)
        {
            double[,,,] psiNPsiNMatrix = new double[6, quadratureNodes.Nodes.Length, n, n];

            for(int f = 0; f < 6; ++f)
            {
                for(int k = 0; k < quadratureNodes.Nodes.Length; ++k)
                {
                    var w = quadratureNodes.Nodes[k].Weight;

                    for (int i = 0; i < n; ++i)
                    {
                        for (int j = 0; j < n; ++j)
                        {
                            psiNPsiNMatrix[f, k, i, j] = w * psiNValues[f, i, k] * psiNValues[f, j, k];
                        }    
                    }
                }
            }

            return psiNPsiNMatrix;
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
            int Upsilon(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Upsilon(i);

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

        public static Vector3D[,,] CalcFacesGradValues(int n, QuadratureNodes<Vector2D> quadratureNodes)
        {
            Vector3D[,,] gradValues = new Vector3D[6, n, quadratureNodes.Nodes.Length];

            Vector3D[] CubeVertex = [new Vector3D(0, 0, 0), new Vector3D(1, 0, 0), new Vector3D(0, 1, 0), new Vector3D(1, 1, 0),
                                     new Vector3D(0, 0, 1), new Vector3D(1, 0, 1), new Vector3D(0, 1, 1), new Vector3D(1, 1, 1)];

            int Mu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Mu(i);
            int Nu(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Nu(i);
            int Upsilon(int i) => LinearParallelepipedalFiniteElementWithNumInteg.Upsilon(i);

            for (int f = 0; f < 6; ++f)
            {
                var face = LinearVectorParallelepipedalFiniteElementWithNumInteg.FaceS(f);

                Vector3D a = CubeVertex[face[0]];
                Vector3D b = CubeVertex[face[1]];
                Vector3D c = CubeVertex[face[2]];
                Vector3D d = CubeVertex[face[3]];

                var nodes = quadratureNodes.Nodes;

                for (int i = 0; i < n; ++i)
                {
                    int ind = face[i];

                    for (int j = 0; j < nodes.Length; ++j)
                    {
                        var node = nodes[j].Node;

                        var newNode = a * (1 - node.X) * (1 - node.Y) + b * node.X * (1 - node.Y) +
                                      c * (1 - node.X) * node.Y       + d * node.X * node.Y;


                        gradValues[f, i, j] = new Vector3D(LinearBasisDerivatives.Phi[Mu(ind)](newNode.X) * LinearBasis.Phi[Nu(ind)](newNode.Y) * LinearBasis.Phi[Upsilon(ind)](newNode.Z),
                                                           LinearBasis.Phi[Mu(ind)](newNode.X) * LinearBasisDerivatives.Phi[Nu(ind)](newNode.Y) * LinearBasis.Phi[Upsilon(ind)](newNode.Z),
                                                           LinearBasis.Phi[Mu(ind)](newNode.X) * LinearBasis.Phi[Nu(ind)](newNode.Y) * LinearBasisDerivatives.Phi[Upsilon(ind)](newNode.Z));
                    }
                }
            }

            return gradValues;
        }
    }

   
}
