using Core.ScalarFiniteElements.FiniteElements1D;
using Core.ScalarFiniteElements.FiniteElements2D;
using FEM;
using Microsoft.CodeAnalysis;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using Telma.ViewModels;
using static Core.RegularRectangularFiniteElementMesh;

namespace Core;

public class IrregularTriangularFiniteElementMesh : IFiniteElementMesh
{
    public IrregularTriangularFiniteElementMesh(string path)
    {
        if (!Directory.Exists(path)) throw new ArgumentException("The path to the mesh folder is incorrect", nameof(path));

        this.path = path;

        materialNumbers = GetMaterialNumbers();
        InputCalculationArea();
        InputBoundaryConditionsArea();

        Elements = CreateFiniteElements();
        FacePortrait = FemAlgorithms.BuildFacePortrait(this);
    }

    public IEnumerable<IFiniteElement> Elements { get; }

    public IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait { get; }

    Vector3D[] vertex = [];
    public Vector3D[] Vertex => vertex;

    public int NumberOfDofs { get; set; }

    string path;
    IDictionary<int, string> materialNumbers;
    Dictionary<(int i, int j), int[]> edges = [];
    List<((int i, int j) edge, (int intervals, double coeff) splitParams)> edgesWithSplitParams = [];

    int[,] MW = new int[0, 0];
    int[,] BC = new int[0, 0];

    Dictionary<int, string> GetMaterialNumbers()
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

        using StreamReader reader = new(Path.Combine(path, "calculationArea.txt"));
        n = int.Parse(reader.ReadLine().Trim());

        Vector3D[] baseVertex = new Vector3D[n];
        double[] point;
        for (int i = 0; i < n; ++i)
        {
            point = [.. reader.ReadLine().Split(' ').Select(double.Parse)];
            baseVertex[i] = new(point[0], point[1], 0);
        }
        vertex = baseVertex;

        n = int.Parse(reader.ReadLine().Trim());
        for (int i = 0; i < n; ++i)
        {
            int k = int.Parse(reader.ReadLine().Trim());

            string[] splitParams = reader.ReadLine().Split(' ');
            int intervals = int.Parse(splitParams[0]);
            double coeff = double.Parse(splitParams[1]);

            for (int j = 0; j < k; ++j)
            {
                string[] input = reader.ReadLine().Split(' ');

                (int i, int j) edge = (int.Parse(input[0]) - 1, int.Parse(input[1]) - 1);
                edgesWithSplitParams.Add((edge, (intervals, coeff)));
            }
        }

        n = int.Parse(reader.ReadLine().Trim());
        MW = new int[n, 6];
        for (int i = 0; i < n; ++i)
        {
            int[] numbers = [.. reader.ReadLine().Split(' ').Select(int.Parse)];

            for (int j = 0; j < 6; ++j)
                MW[i, j] = numbers[j] - 1;
        }
    }

    void InputBoundaryConditionsArea()
    {

        using StreamReader reader = new(Path.Combine(path, "boundaryConditionsArea.txt"));
        int n = int.Parse(reader.ReadLine().Trim());
        BC = new int[n, 4];

        for (int i = 0; i < n; ++i)
        {
            int[] numbers = [.. reader.ReadLine().Split(' ').Select(int.Parse)];

            for (int j = 0; j < 4; ++j)
                BC[i, j] = numbers[j] - 1;
        }
    }

    void SplitEdges()
    {
        List<Vector3D> newVertex = new(vertex);
        int currVertex = newVertex.Count;

        foreach (var el in edgesWithSplitParams)
        {
            (int i, int j) edge = el.edge;
            int intervals = el.splitParams.intervals;
            double coeff = el.splitParams.coeff;

            Vector3D a = vertex[edge.i];
            Vector3D b = vertex[edge.j];
            int[] vertexNumbers = new int[intervals + 1];
            vertexNumbers[0] = edge.i;
            vertexNumbers[^1] = edge.j;

            double step = 0;
            if (coeff == 1d)
            {
                step = 1d / intervals;

                for (int j = 1; j < intervals; ++j)
                {
                    double p = step * j;
                    newVertex.Add(a * (1 - p) + p * b);
                    vertexNumbers[j] = currVertex++;
                }
            }
            else
            {
                step = (1 - coeff) / (1 - Math.Pow(coeff, intervals));

                for (int j = 1; j < intervals; ++j)
                {
                    double p = step * (1 - Math.Pow(coeff, j)) / (1 - coeff);
                    newVertex.Add(a * (1 - p) + p * b);
                    vertexNumbers[j] = currVertex++;
                }
            }

            if (edge.i > edge.j) edge = (edge.j, edge.i);
            edges[edge] = vertexNumbers;
        }

        vertex = [.. newVertex];
    }

    List<IFiniteElement> CreateFiniteElements()
    {
        SplitEdges();
        List<Vector3D> newVertex = new(vertex);
        List<IFiniteElement> elements = [];

        int nAreas = MW.GetLength(0);

        for (int area = 0; area < nAreas; ++area)
        {
            string material = materialNumbers[MW[area, 0]];
            FEM.IFiniteElement.ElementType type = (IFiniteElement.ElementType)MW[area, 1];

            switch (type)
            {
                case IFiniteElement.ElementType.Scalar:

                    (int[] vetexNumbers, int pN, int sN) = GetVertexNumbersForArea(area, newVertex, newVertex.Count);

                    for (int s = 0; s < sN - 1; ++s)
                    {
                        for (int p = 0; p < pN - 1; ++p)
                        {
                            (int[] tr1, int[] tr2) = GetVertexNumbersForTriangles(newVertex, vetexNumbers, p, s, pN);
                            elements.Add(new QuadraticTriangularFiniteElementWithNumInteg(material, tr1));
                            elements.Add(new QuadraticTriangularFiniteElementWithNumInteg(material, tr2));
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

            (int i, int j) edge = (BC[bound, 2], BC[bound, 3]);
            if (edge.i > edge.j) edge = (edge.j, edge.i);

            int[] vertexNumbers = edges[edge];

            for (int i = 0; i < vertexNumbers.Length - 1; ++i)
            {
                switch (type)
                {
                    case IFiniteElement.ElementType.Scalar:

                        elements.Add(new QuadraticStraightFiniteElement(material, vertexNumbers[i..(i + 2)]));

                        break;

                    default:
                        throw new ArgumentException("The element type is set incorrectly");
                }
            }
        }

        vertex = [.. newVertex];

        return elements;
    }

    (int[], int[]) GetVertexNumbersForTriangles(List<Vector3D> newVertex, int[] vertexNumbers, int p, int s, int pN)
    {
        int[] quadVertexNumber =
        [
            vertexNumbers[s * pN + p],
            vertexNumbers[s * pN + p + 1],
            vertexNumbers[(s + 1) * pN + p],
            vertexNumbers[(s + 1) * pN + p + 1],
        ];

        double dist1 = newVertex[quadVertexNumber[0]].Distance(newVertex[quadVertexNumber[3]]);
        double dist2 = newVertex[quadVertexNumber[1]].Distance(newVertex[quadVertexNumber[2]]);

        int[] tr1 =
        [
            quadVertexNumber[0],
            quadVertexNumber[1],
            quadVertexNumber[dist1 < dist2 ? 3 : 2]
        ];
        int[] tr2 =
        [
            quadVertexNumber[dist1 < dist2 ? 0 : 1],
            quadVertexNumber[3],
            quadVertexNumber[2]
        ];

        //tr1[0] = vertexNumbers[s * pN + p];
        //tr1[1] = vertexNumbers[s * pN + p + 1];
        //tr1[2] = vertexNumbers[(s + 1) * pN + p];

        //tr2[0] = vertexNumbers[s * pN + p + 1];
        //tr2[1] = vertexNumbers[(s + 1) * pN + p + 1];
        //tr2[2] = vertexNumbers[(s + 1) * pN + p];

        return (tr1, tr2);
    }

    (int[], int, int) GetVertexNumbersForArea(int area, List<Vector3D> newVertex, int currVertex)
    {
        (int i, int j) h1 = (MW[area, 2], MW[area, 3]);
        if (h1.i > h1.j) h1 = (h1.j, h1.i);
        (int i, int j) h2 = (MW[area, 4], MW[area, 5]);
        if (h2.i > h2.j) h2 = (h2.j, h2.i);
        (int i, int j) v1 = (MW[area, 2], MW[area, 4]);
        if (v1.i > v1.j) v1 = (v1.j, v1.i);
        (int i, int j) v2 = (MW[area, 3], MW[area, 5]);
        if (v2.i > v2.j) v2 = (v2.j, v2.i);

        int pN = edges[h1].Length;
        int sN = edges[v2].Length;

        int[] vertexNumbers = new int[pN * sN];

        for (int i = 0; i < pN; ++i)
        {
            vertexNumbers[i] = edges[h1][i];
            vertexNumbers[i + pN * (sN - 1)] = edges[h2][i];
        }

        for (int i = 0; i < sN; ++i)
        {
            vertexNumbers[i * pN] = edges[v1][i];
            vertexNumbers[i * pN + (pN - 1)] = edges[v2][i];
        }

        for (int s = 1; s < sN - 1; ++s)
        {
            Vector3D p1 = vertex[edges[v1][s]];
            Vector3D p2 = vertex[edges[v2][s]];
            for (int p = 1; p < pN - 1; ++p)
            {
                Vector3D q1 = vertex[edges[h1][p]];
                Vector3D q2 = vertex[edges[h2][p]];

                Vector3D pp = LinearAlgebraAlgorithms.IntersectTwoSegments(p1, p2, q1, q2);

                newVertex.Add(LinearAlgebraAlgorithms.IntersectTwoSegments(p1, p2, q1, q2));
                vertexNumbers[s * pN + p] = currVertex++;
            }
        }

        return (vertexNumbers, pN, sN);
    }

    public void CreateFilesWithMesh(string path = "")
    {
        using (StreamWriter writer = new StreamWriter(Path.Combine(path, "triangles.txt")))
        {
            foreach (var elem in Elements)
            {
                if (elem.VertexNumber.Length > 2)
                    writer.WriteLine(string.Join(' ', elem.VertexNumber));
            }
        }

        using (StreamWriter writer = new StreamWriter(Path.Combine(path, "vertices.txt")))
        {
            foreach (var v in Vertex)
            {
                writer.WriteLine($"{v.X} {v.Y}");
            }
        }
    }
}
