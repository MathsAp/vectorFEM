using Core;
using Core.ScalarFiniteElements.FiniteElements1D;
using Core.ScalarFiniteElements.FiniteElements2D;
using Core.ScalarFiniteElements.FiniteElements3D;
using Core.VectorFiniteElements.FiniteElements2D;
using Core.VectorFiniteElements.FiniteElements3D;
using Core.VectorScalarFiniteElements.FiniteElements2D;
using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core2
{
    public class RegularRectangularFiniteElementMesh : IFiniteElementMesh
    {
        public RegularRectangularFiniteElementMesh(string path, Dimension dim)
        {
            this.path = path;
            this.dim = dim;

            if (!Path.Exists(path)) throw new ArgumentException("The path to the mesh folder is incorrect");

            materialNumbers = GetMaterialNumbers();

            InputCalculationArea();
            InputBoundaryConditionsArea();
            InputSplitCalcArea();

            (IXW, X) = CreateOneDimensionalMesh(XW, refineParams[0]);
            (IYW, Y) = CreateOneDimensionalMesh(YW, refineParams[1]);
            if (dim is Dimension.D3) (IZW, Z) = CreateOneDimensionalMesh(ZW, refineParams[2]);

            vertex = CreateVertexArray();
            elements = CreateFiniteElements();
            facePortrait = FemAlgorithms.BuildFacePortrait(this);
        }

        IEnumerable<IFiniteElement> elements;
        public IEnumerable<IFiniteElement> Elements => elements;

        IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> facePortrait;
        public IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait => facePortrait;

        Vector3D[] vertex = [];
        public Vector3D[] Vertex => vertex;

        public int NumberOfDofs { get; set; }


        public enum Dimension { D2, D3 };

        Dimension dim;

        string path;

        IDictionary<int, string> materialNumbers;

        double[] XW = [];
        double[] YW = [];
        double[] ZW = [];

        int[] IXW = [];
        int[] IYW = [];
        int[] IZW = [];

        double[] X = [];
        double[] Y = [];
        double[] Z = [];

        int[,] MW = new int[0, 0];
        int[,] BC = new int[0, 0];

        RefineParams[] refineParams = new RefineParams[3];

        enum Coordinate { x, y, z };

        IDictionary<int, string> GetMaterialNumbers()
        {
            Dictionary<int, string> materialNumbers = new();

            using (StreamReader reader = new(Path.Combine(path, "materials.txt")))
            {
                int n = int.Parse(reader.ReadLine());
                string[] values;
                for (int i = 0; i < n; ++i)
                {
                    values = reader.ReadLine().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                    materialNumbers[int.Parse(values[1])] = values[0];
                }
            }

            return materialNumbers;
        }

        void InputCalculationArea()
        {
            int n = 0;

            using (StreamReader reader = new(Path.Combine(path, "calculationArea.txt")))
            {
                n = int.Parse(reader.ReadLine().Split(' ')[0]);

                XW = new double[n];
                var numbers = reader.ReadLine().Split(' ');

                for (int i = 0; i < n; ++i)
                    XW[i] = double.Parse(numbers[i]);


                n = int.Parse(reader.ReadLine().Split(' ')[0]);

                YW = new double[n];
                numbers = reader.ReadLine().Split(' ');

                for (int i = 0; i < n; ++i)
                    YW[i] = double.Parse(numbers[i]);

                if (dim is Dimension.D3)
                {
                    n = int.Parse(reader.ReadLine().Split(' ')[0]);

                    ZW = new double[n];
                    numbers = reader.ReadLine().Split(' ');

                    for (int i = 0; i < n; ++i)
                        ZW[i] = double.Parse(numbers[i]);
                }

                n = int.Parse(reader.ReadLine().Split(' ')[0]);
                int k = dim is Dimension.D3 ? 8 : 6;
                MW = new int[n, k];

                for (int i = 0; i < n; ++i)
                {
                    numbers = reader.ReadLine().Split(' ');
                    for (int j = 0; j < k; ++j)
                    {
                        MW[i, j] = int.Parse(numbers[j]) - 1;
                    }
                }
            }
        }

        void InputBoundaryConditionsArea()
        {
            int n = 0;

            using (StreamReader reader = new(Path.Combine(path, "boundaryConditionsArea.txt")))
            {
                n = int.Parse(reader.ReadLine().Split(' ')[0]);
                int k = dim is Dimension.D3 ? 8 : 6;
                BC = new int[n, k];

                for (int i = 0; i < n; ++i)
                {
                    var numbers = reader.ReadLine().Split(' ');

                    for (int j = 0; j < k; ++j)
                    {
                        BC[i, j] = int.Parse(numbers[j]) - 1;
                    }
                }
            }
        }

        void InputSplitCalcArea()
        {
            using (StreamReader reader = new(Path.Combine(path, "splitCalcArea.txt")))
            {
                int d = dim is Dimension.D3 ? 3 : 2;
                for (int i = 0; i < d; ++i)
                {
                    refineParams[i] = new();

                    var numbers = reader.ReadLine().Split(' ');

                    int n = numbers.Length / 2;

                    for (int j = 0; j < n; ++j)
                    {
                        refineParams[i].splitCount.Add(int.Parse(numbers[2 * j]));
                        refineParams[i].stretchRatio.Add(double.Parse(numbers[2 * j + 1]));
                    }
                }
            }
        }

        (int[] ImeshW, double[] mesh) CreateOneDimensionalMesh(ReadOnlySpan<double> baseMesh, RefineParams refineParams)
        {
            int n = baseMesh.Length;
            int newN = refineParams.splitCount.Sum() + 1;

            double[] mesh = new double[newN];
            int[] ImeshW = new int[n];

            mesh[0] = baseMesh[0];
            ImeshW[0] = 0;

            int sum = 0;

            for (int i = 0; i < n - 1; ++i)
            {
                int numIntervals = refineParams.splitCount[i];
                double coeff = refineParams.stretchRatio[i];

                double step = 0;
                if (coeff == 1d)
                {
                    step = (baseMesh[i + 1] - baseMesh[i]) / numIntervals;

                    for (int j = 1; j < numIntervals; ++j)
                    {
                        mesh[j + sum] = baseMesh[i] + step * j;
                    }
                }
                else
                {
                    step = (baseMesh[i + 1] - baseMesh[i]) * (1 - coeff) / (1 - Math.Pow(coeff, numIntervals));

                    for (int j = 1; j < numIntervals; ++j)
                    {
                        mesh[j + sum] = baseMesh[i] + step * (1 - Math.Pow(coeff, j)) / (1 - coeff);
                    }
                }

                sum += numIntervals;
                mesh[sum] = baseMesh[i + 1];
                ImeshW[i + 1] = sum;
            }

            return (ImeshW, mesh);
        }

        void ChangeRefineParams(RefineParams refineParams)
        {
            for (int i = 0; i < refineParams.splitCount.Count; ++i)
            {
                refineParams.splitCount[i] *= 2;

                refineParams.stretchRatio[i] = Math.Sqrt(refineParams.stretchRatio[i]);
            }
        }

        Vector3D[] CreateVertexArray()
        {
            int nX = X.Length;
            int nY = Y.Length;
            int nZ = dim is Dimension.D3 ? Z.Length : 1;

            Vector3D[] vertex = new Vector3D[nX * nY * nZ];

            int ind = 0;
            for (int k = 0; k < nZ; ++k)
            {
                for (int j = 0; j < nY; ++j)
                {
                    for (int i = 0; i < nX; ++i)
                    {
                        vertex[ind++] = dim is Dimension.D3 ? new Vector3D(X[i], Y[j], Z[k]) : new Vector3D(X[i], Y[j], 0);
                    }
                }
            }

            return vertex;
        }

        List<IFiniteElement> CreateFiniteElements()
        {
            List<IFiniteElement> elements = [];

            int nAreas = MW.GetLength(0);

            for (int area = 0; area < nAreas; ++area)
            {
                string material = materialNumbers[MW[area, 0]];
                FEM.IFiniteElement.ElementType type = (IFiniteElement.ElementType)MW[area, 1];

                int ix0 = IXW[MW[area, 2]];
                int ix1 = IXW[MW[area, 3]];

                int iy0 = IYW[MW[area, 4]];
                int iy1 = IYW[MW[area, 5]];

                int iz0 = dim is Dimension.D3 ? IZW[MW[area, 6]] : 0;
                int iz1 = dim is Dimension.D3 ? IZW[MW[area, 7]] : 1;

                switch (type)
                {
                    case IFiniteElement.ElementType.Scalar:

                        for (int r = iz0; r < iz1; ++r)
                        {
                            for (int s = iy0; s < iy1; ++s)
                            {
                                for (int p = ix0; p < ix1; ++p)
                                {
                                    IFiniteElement elem = dim is Dimension.D3 
                                        ? new LinearParallelepipedalFiniteElementWithNumInteg(material, Get3DElementVertexNumber(p, s, r)) 
                                        : new LinearRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.z));
                                    elements.Add(elem);
                                }
                            }
                        }

                        break;

                    case IFiniteElement.ElementType.Vector:

                        for (int r = iz0; r < iz1; ++r)
                        {
                            for (int s = iy0; s < iy1; ++s)
                            {
                                for (int p = ix0; p < ix1; ++p)
                                {
                                    elements.Add(new LinearVectorParallelepipedalFiniteElementWithNumInteg(material, Get3DElementVertexNumber(p, s, r)));
                                }
                            }
                        }

                        break;

                    default:
                        throw new ArgumentException("The element type is set incorrectly");
                }
            }

            int nBC = BC.GetLength(0);

            for (int bound = 0; bound < nBC; ++bound)
            {
                string material = materialNumbers[BC[bound, 0]];
                FEM.IFiniteElement.ElementType type = (IFiniteElement.ElementType)BC[bound, 1];

                int ix0 = IXW[BC[bound, 2]];
                int ix1 = IXW[BC[bound, 3]];

                int iy0 = IYW[BC[bound, 4]];
                int iy1 = IYW[BC[bound, 5]];

                int iz0 = dim is Dimension.D3 ? IZW[BC[bound, 6]] : 0;
                int iz1 = dim is Dimension.D3 ? IZW[BC[bound, 7]] : 1;

                if (ix0 == ix1)
                {
                    int p = ix0;

                    for (int r = iz0; r < iz1; ++r)
                    {
                        for (int s = iy0; s < iy1; ++s)
                        {
                            switch (type)
                            {
                                case IFiniteElement.ElementType.Scalar:

                                    IFiniteElement elem = dim is Dimension.D3
                                        ? new LinearRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.x))
                                        : new LinearStraightFiniteElementWithNumInteg(material, Get1DElementVertexNumber(p, s, Coordinate.x));

                                    elements.Add(elem);

                                    break;

                                case IFiniteElement.ElementType.Vector:

                                    elements.Add(new LinearVectorRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.x)));

                                    break;

                                case IFiniteElement.ElementType.VectorScalar:

                                    elements.Add(new LinearVectorScalarRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.x)));

                                    break;

                                default:
                                    throw new ArgumentException("The element type is set incorrectly");
                            }
                        }
                    }
                }
                else if (iy0 == iy1)
                {
                    int s = iy0;

                    for (int r = iz0; r < iz1; ++r)
                    {
                        for (int p = ix0; p < ix1; ++p)
                        {
                            switch (type)
                            {
                                case IFiniteElement.ElementType.Scalar:

                                    IFiniteElement elem = dim is Dimension.D3 
                                        ? new LinearRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.y))
                                        : new LinearStraightFiniteElementWithNumInteg(material, Get1DElementVertexNumber(p, s, Coordinate.y));

                                    elements.Add(elem);

                                    break;

                                case IFiniteElement.ElementType.Vector:

                                    elements.Add(new LinearVectorRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.y)));

                                    break;

                                case IFiniteElement.ElementType.VectorScalar:

                                    elements.Add(new LinearVectorScalarRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.y)));

                                    break;

                                default:
                                    throw new ArgumentException("The element type is set incorrectly");
                            }
                        }
                    }
                }
                else if (iz0 == iz1)
                {
                    int r = iz0;

                    for (int s = iy0; s < iy1; ++s)
                    {
                        for (int p = ix0; p < ix1; ++p)
                        {
                            switch (type)
                            {
                                case IFiniteElement.ElementType.Scalar:

                                    elements.Add(new LinearRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.z)));

                                    break;

                                case IFiniteElement.ElementType.Vector:

                                    elements.Add(new LinearVectorRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.z)));

                                    break;

                                case IFiniteElement.ElementType.VectorScalar:

                                    elements.Add(new LinearVectorScalarRectangularFiniteElementWithNumInteg(material, Get2DElementVertexNumber(p, s, r, Coordinate.z)));

                                    break;

                                default:
                                    throw new ArgumentException("The element type is set incorrectly");
                            }
                        }
                    }
                }
            }

            return elements;
        }

        int[] Get3DElementVertexNumber(int p, int s, int r)
        {
            int[] vertexNumber = new int[8];

            int nX = X.Length;
            int nY = Y.Length;
            int nZ = Z.Length;

            vertexNumber[0] = r * nY * nX + s * nX + p;
            vertexNumber[1] = r * nY * nX + s * nX + p + 1;
            vertexNumber[2] = r * nY * nX + (s + 1) * nX + p;
            vertexNumber[3] = r * nY * nX + (s + 1) * nX + p + 1;
            vertexNumber[4] = (r + 1) * nY * nX + s * nX + p;
            vertexNumber[5] = (r + 1) * nY * nX + s * nX + p + 1;
            vertexNumber[6] = (r + 1) * nY * nX + (s + 1) * nX + p;
            vertexNumber[7] = (r + 1) * nY * nX + (s + 1) * nX + p + 1;

            return vertexNumber;
        }

        int[] Get2DElementVertexNumber(int p, int s, int r, Coordinate coord)
        {
            int[] vertexNumber = new int[4];

            int nX = X.Length;
            int nY = Y.Length;
            int nZ = Z.Length;

            switch (coord)
            {
                case Coordinate.x:

                    vertexNumber[0] = r * nY * nX + s * nX + p;
                    vertexNumber[1] = r * nY * nX + (s + 1) * nX + p;
                    vertexNumber[2] = (r + 1) * nY * nX + s * nX + p;
                    vertexNumber[3] = (r + 1) * nY * nX + (s + 1) * nX + p;

                    break;

                case Coordinate.y:

                    vertexNumber[0] = r * nY * nX + s * nX + p;
                    vertexNumber[1] = r * nY * nX + s * nX + p + 1;
                    vertexNumber[2] = (r + 1) * nY * nX + s * nX + p;
                    vertexNumber[3] = (r + 1) * nY * nX + s * nX + p + 1;

                    break;

                case Coordinate.z:

                    vertexNumber[0] = r * nY * nX + s * nX + p;
                    vertexNumber[1] = r * nY * nX + s * nX + p + 1;
                    vertexNumber[2] = r * nY * nX + (s + 1) * nX + p;
                    vertexNumber[3] = r * nY * nX + (s + 1) * nX + p + 1;

                    break;

            }

            return vertexNumber;
        }

        int[] Get1DElementVertexNumber(int p, int s, Coordinate coord)
        {
            int[] vertexNumber = new int[2];

            int nX = X.Length;

            switch (coord)
            {
                case Coordinate.x:
                    vertexNumber[0] = nX * s + p;
                    vertexNumber[1] = nX * (s + 1) + p;
                    break;
                case Coordinate.y:
                    vertexNumber[0] = nX * s + p;
                    vertexNumber[1] = nX * s + p + 1;
                    break;
            }

            return vertexNumber;
        }
    }
}
