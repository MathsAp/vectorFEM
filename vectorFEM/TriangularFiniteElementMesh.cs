using Core.ScalarFiniteElements.FiniteElements1D;
using Core.ScalarFiniteElements.FiniteElements2D;
using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TelmaQuasarCommon.Core.Geometry;
using TelmaQuasarCommon.Problems;

namespace Core;

file class GeometryObject(string name, int materialNumber, int[] vertexNumber)
{
    public string Name { get; } = name;

    public int MaterialNumber { get; } = materialNumber;

    public int[] VertexNumber = vertexNumber;
}

public class TriangularFiniteElementMesh : IFiniteElementMesh
{
    public TriangularFiniteElementMesh(string path, bool isTelma = false)
    {
        if (!File.Exists(path)) throw new ArgumentException("The path to the mesh file is incorrect", nameof(path));

        this.path = path;

        (elements, vertex) = ReadMeshFromFile();
        FacePortrait = FemAlgorithms.BuildFacePortrait(this);
    }


    string path;

    List<IFiniteElement> elements;
    public IEnumerable<IFiniteElement> Elements => elements;

    public IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait { get; }

    Vector3D[] vertex;
    public Vector3D[] Vertex => vertex;

    public int NumberOfDofs { get; set; }

    (List<IFiniteElement>, Vector3D[]) ReadMeshFromFile()
    {
        using StreamReader reader = new(path);

        int n = int.Parse(reader.ReadLine());

        Vector3D[] vertex = new Vector3D[n];
        for (int i = 0; i < n; ++i)
        {
            vertex[i] = new Vector3D([.. reader.ReadLine().Split('\t', StringSplitOptions.RemoveEmptyEntries).Select(double.Parse)]);
        }

        n = int.Parse(reader.ReadLine());
        GeometryObject[] geomObjs = new GeometryObject[n];
        for (int i = 0; i < n; ++i)
        {
            string[] input = reader.ReadLine().Split(' ', StringSplitOptions.RemoveEmptyEntries);

            geomObjs[i] = new(input[0], int.Parse(input[3]), [.. input[5..].Select(int.Parse)]);
        }

        n = int.Parse(reader.ReadLine());
        Dictionary<int, string> volumeMaterials = new(n);
        for (int i = 0; i < n; ++i)
        {
            string[] input = reader.ReadLine().Split(' ', StringSplitOptions.RemoveEmptyEntries);
            volumeMaterials[int.Parse(input[0])] = input[1];
        }

        n = int.Parse(reader.ReadLine());
        Dictionary<int, string> boundaryMaterials = new(n);
        for (int i = 0; i < n; ++i)
        {
            string[] input = reader.ReadLine().Split(' ', StringSplitOptions.RemoveEmptyEntries);
            boundaryMaterials[int.Parse(input[0])] = input[1];
        }


        List<IFiniteElement> elems = new(geomObjs.Length);
        foreach(var geomObj in geomObjs)
        {
            IFiniteElement elem = geomObj.Name is "Triangle" 
                ? new QuadraticTriangularFiniteElementWithNumInteg(volumeMaterials[geomObj.MaterialNumber], geomObj.VertexNumber)
                : new QuadraticStraightFiniteElement(boundaryMaterials[geomObj.MaterialNumber], geomObj.VertexNumber);

            elems.Add(elem);
        }

        return (elems, vertex);
    }
}
