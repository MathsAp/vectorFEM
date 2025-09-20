using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core;

public class VectorScalarEllipticSolution : ISolution
{
    public VectorScalarEllipticSolution(IFiniteElementMesh mesh, IDictionary<string, IMaterial> materials)
    {
        Mesh = mesh;
        Materials = materials;
    }
    public double Time { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
    public ITimeMesh TimeMesh => throw new NotImplementedException();

    public IFiniteElementMesh Mesh { get; }

    IDictionary<string, IMaterial> Materials { get; }

    double[] solutionVector { get; set; } = [];
    public ReadOnlySpan<double> SolutionVector => solutionVector;

    public void AddSolutionVector(double t, double[] solution)
    {
        solutionVector = solution;
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
        foreach(var element in Mesh.Elements)
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
                        return material.Hext(point, 1) - element.GetGradientAtPoint(Mesh.Vertex, solutionVector, point);
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
                        return Constants.Mu0 * (material.Hext(point, 1) - element.GetGradientAtPoint(Mesh.Vertex, solutionVector, point));
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
}
