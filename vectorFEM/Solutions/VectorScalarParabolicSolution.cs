using Core;
using FEM;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core;

public class VectorScalarParabolicSolution : ISolution
{
    public VectorScalarParabolicSolution(IFiniteElementMesh mesh, ITimeMesh timeMesh, IDictionary<string, IMaterial> materials, string _path = "")
    {
        Mesh = mesh;
        TimeMesh = timeMesh;
        Materials = materials;
        solutionVector = new double[mesh.NumberOfDofs];

        if (_path.Length == 0)
            path = "VectorScalarParabolicProblemWeights";

        if (Directory.Exists(path))
            Directory.Delete(path, true);

        Directory.CreateDirectory(path!);
    }

    double time = -1;
    public double Time
    {
        get { return time; }

        set
        {
            if (TimeMesh[0] <= value && value <= TimeMesh[TimeMesh.Size() - 1])
            {
                if (value != time)
                {
                    int ind = BinarySearch(TimeMesh, value, 0, TimeMesh.Size() - 1);
                    time = TimeMesh[ind];

                    using (StreamReader reader = new StreamReader(Path.Combine(path, time.ToString() + ".txt")))
                    {
                        string? coeff = null;

                        for (int i = 0; (coeff = reader.ReadLine()) != null; ++i)
                        {
                            solutionVector[i] = double.Parse(coeff);
                        }
                    }
                }
            }
        }
    }

    string path = "";
    public ITimeMesh TimeMesh { get; }

    public IFiniteElementMesh Mesh { get; }

    IDictionary<string, IMaterial> Materials { get; }

    double[] solutionVector { get; set; }
    public ReadOnlySpan<double> SolutionVector => solutionVector;

    public void AddSolutionVector(double t, double[] solution)
    {
        using (StreamWriter writer = new StreamWriter(Path.Combine(path, t.ToString() + ".txt"), false))
        {
            foreach (var coeff in solution)
                writer.WriteLine(coeff);
        }
    }

    public Vector3D Gradient(Vector3D point)
    {
        foreach (var elem in Mesh.Elements)
        {
            if (elem.VertexNumber.Length > 4)
            {
                if (elem.IsPointOnElement(Mesh.Vertex, point))
                    return elem.GetGradientAtPoint(Mesh.Vertex, SolutionVector, point);
            }
        }

        return new Vector3D(0, 0, 0);
    }

    public Vector3D Curl(Vector3D point)
    {
        foreach (var element in Mesh.Elements)
        {
            if (element.VertexNumber.Length > 4)
            {
                if (element.IsPointOnElement(Mesh.Vertex, point))
                {
                    return element.GetCurlAtPoint(Mesh.Vertex, solutionVector, point);
                }
            }
        }

        return Vector3D.Zero;
    }

    public double Value(Vector3D point)
    {
        foreach (var elem in Mesh.Elements)
        {
            if (elem.VertexNumber.Length > 4)
            {
                if (elem.IsPointOnElement(Mesh.Vertex, point))
                    return elem.GetValueAtPoint(Mesh.Vertex, SolutionVector, point);
            }
        }

        return 0.0;
    }

    public Vector3D Vector(Vector3D point)
    {
        foreach (var element in Mesh.Elements)
        {
            if (element.VertexNumber.Length > 4)
            {
                if (element.IsPointOnElement(Mesh.Vertex, point))
                {
                    return element.GetVectorAtPoint(Mesh.Vertex, solutionVector, point);
                }
            }
        }

        return Vector3D.Zero;
    }

    public Vector3D H(Vector3D point)
    {
        foreach (var element in Mesh.Elements)
        {
            if (element.VertexNumber.Length > 4)
            {
                if (element.IsPointOnElement(Mesh.Vertex, point))
                {
                    var material = Materials[element.Material];

                    if (element.Type == IFiniteElement.ElementType.Vector)
                    {
                        return (1 / material.Mu(point)) * element.GetCurlAtPoint(Mesh.Vertex, solutionVector, point);
                    }
                    else if (element.Type == IFiniteElement.ElementType.Scalar)
                    {
                        return material.Hext(point, Time) - element.GetGradientAtPoint(Mesh.Vertex, solutionVector, point);
                    }
                }
            }
        }

        return Vector3D.Zero;
    }

    public Vector3D B(Vector3D point)
    {
        foreach (var element in Mesh.Elements)
        {
            if (element.VertexNumber.Length > 4)
            {
                if (element.IsPointOnElement(Mesh.Vertex, point))
                {
                    var material = Materials[element.Material];

                    if (element.Type == IFiniteElement.ElementType.Vector)
                    {
                        return element.GetCurlAtPoint(Mesh.Vertex, solutionVector, point);
                    }
                    else if (element.Type == IFiniteElement.ElementType.Scalar)
                    {
                        return Constants.Mu0 * (material.Hext(point, Time) - element.GetGradientAtPoint(Mesh.Vertex, solutionVector, point));
                    }
                }
            }
        }

        return Vector3D.Zero;
    }

    public double CalcNormL2(Func<Vector3D, Vector3D> u)
    {
        throw new NotImplementedException();
    }

    public static int BinarySearch(ITimeMesh timeMesh, double target, int low, int high)
    {
        int mid = 0;
        double midValue = 0;

        while (low <= high)
        {
            mid = (low + high) / 2;
            midValue = timeMesh[mid];

            if (midValue == target) return mid;
            else if (midValue < target) low = mid + 1;
            else high = mid - 1;
        }

        return mid;
    }
}
