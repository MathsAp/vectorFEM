using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public class EllipticSolution : ISolution
    {
        public EllipticSolution(IFiniteElementMesh mesh) 
        { 
            Mesh = mesh;
            solutionVector = new double[mesh.NumberOfDofs];
        }
        public double Time { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
        public ITimeMesh TimeMesh => throw new NotImplementedException();

        public IFiniteElementMesh Mesh { get; }

        double[] solutionVector { get; }
        public ReadOnlySpan<double> SolutionVector => solutionVector;

        public void AddSolutionVector(double t, double[] solution)
        {
            int N = Mesh.NumberOfDofs;

            for (int i = 0; i < N; ++i)
            {
                solutionVector[i] = solution[i];
            }
        }

        public Vector3D Gradient(Vector3D point)
        {
            foreach(var elem in Mesh.Elements)
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
                if (element.VertexNumber.Length != 2)
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
            foreach(var elem in Mesh.Elements)
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
                if (element.VertexNumber.Length != 2)
                {
                    if (element.IsPointOnElement(Mesh.Vertex, point))
                    {
                        return element.GetVectorAtPoint(Mesh.Vertex, solutionVector, point);
                    }
                }
            }

            return Vector3D.Zero;
        }
    }
}
