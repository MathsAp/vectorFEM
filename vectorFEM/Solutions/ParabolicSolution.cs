using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FEM;

namespace Core
{
    public class ParabolicSolution : ISolution
    {
        public ParabolicSolution(IFiniteElementMesh mesh, ITimeMesh timeMesh, string _path = "")
        {
            Mesh = mesh;
            TimeMesh = timeMesh;
            solutionVector = new double[mesh.NumberOfDofs];

            if (_path.Length == 0)
                path = "ParabolicProblemWeights";

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
        public IFiniteElementMesh Mesh { get; }
        public ITimeMesh TimeMesh { get; }

        double[] solutionVector { get; }
        public ReadOnlySpan<double> SolutionVector => solutionVector;

        public double Value(Vector3D point)
        {
            foreach (var element in Mesh.Elements)
            {
                if (element.VertexNumber.Length > 4)
                {
                    if (element.IsPointOnElement(Mesh.Vertex, point))
                    {
                        return element.GetValueAtPoint(Mesh.Vertex, solutionVector, point);
                    }
                }
            }

            return 0; // что возвращать, если не попала ни в один элемент?
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

        public Vector3D Gradient(Vector3D point)
        {
            foreach (var element in Mesh.Elements)
            {
                if (element.VertexNumber.Length > 4)
                {
                    if (element.IsPointOnElement(Mesh.Vertex, point))
                    {
                        return element.GetGradientAtPoint(Mesh.Vertex, solutionVector, point);
                    }
                }
            }

            return Vector3D.Zero;
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

        public void AddSolutionVector(double t, double[] solution)
        {
            using (StreamWriter writer = new StreamWriter(Path.Combine(path, t.ToString() + ".txt"), false))
            {
                foreach (var coeff in solution)
                    writer.WriteLine(coeff);
            }
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

        public double CalcNormL2(Func<Vector3D, Vector3D> u)
        {
            throw new NotImplementedException();
        }
    }

    //public class Solution2 : ISolution
    //{
    //    public Solution2(IFiniteElementMesh mesh, ITimeMesh timeMesh)
    //    {
    //        Mesh = mesh;
    //        TimeMesh = timeMesh;
    //    }

    //    public double Time { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

    //    IFiniteElementMesh Mesh { get; }
    //    ITimeMesh TimeMesh { get; }
    //    double[][] SolutionsVectors { get; }

    //    double[] solutionVector {  get; }
    //    public ReadOnlySpan<double> SolutionVector => solutionVector;

    //    public void AddSolutionVector(double t, double[] solution) => throw new NotImplementedException();
    //    public double Value(Vector2D point) => throw new NotImplementedException();
    //}
}
