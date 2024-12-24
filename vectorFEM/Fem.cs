using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata.Ecma335;
using System.Text;
using System.Threading.Tasks;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace FEM
{
    public static class FemAlgorithms
    {
        public static Dictionary<(int, int), int> BuildEdgePortrait(IFiniteElementMesh mesh)
        {
            var dict = new Dictionary<(int, int), int>();
            foreach (var element in mesh.Elements)
            {
                for (int i = 0; i < element.NumberOfEdges; i++)
                {
                    var edge = element.Edge(i);
                    edge = (element.VertexNumber[edge.i], element.VertexNumber[edge.j]);
                    if (edge.i < edge.j) edge = (edge.j, edge.i);
                    var n = element.DOFOnEdge(i, element.Type);
                    if (!dict.TryGetValue(edge, out int c) || c < n) dict[edge] = n;
                }
            }
            return dict;
        }

        public static Dictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> BuildFacePortrait(IFiniteElementMesh mesh)
        {
            var dict = new Dictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))>();


            foreach(var element in mesh.Elements)
            {
                if (element.VertexNumber.Length > 4)
                    for (int i = 0; i < element.NumberOfFaces; ++i)
                    {
                        var face = element.Face(i);

                        face[0] = element.VertexNumber[face[0]];
                        face[1] = element.VertexNumber[face[1]];
                        face[2] = element.VertexNumber[face[2]];
                        face[3] = element.VertexNumber[face[3]];

                        Array.Sort(face);

                        var faceTuple = (face[0], face[1], face[2], face[3]);

                        if (!dict.TryGetValue(faceTuple, out ((IFiniteElement?, int), (IFiniteElement?, int)) value)) dict[faceTuple] = ((element, i), (null, -1));
                        else dict[faceTuple] = (dict[faceTuple].Item1, (element, i));
                    
                    }
            }    

            return dict;
        }

        public static void EnumerateMeshDofs(IFiniteElementMesh mesh)
        {
            int dof = 0;
            int[] VertexDof = new int[mesh.Vertex.Length];

            var nodeNumbers = GetNodeNumbersOfScalarElements(mesh);

            foreach(var num in nodeNumbers)
            {
                VertexDof[num] = dof++;
            }

            foreach (var element in mesh.Elements)
            {
                for (int i = 0; i < element.VertexNumber.Length; i++)
                    element.SetVertexDOF(i, VertexDof[element.VertexNumber[i]]); // для векторных работает вхолостую
            }

            var edges = BuildEdgePortrait(mesh);
            edges = edges.ToDictionary(edgeinfo => edgeinfo.Key, edgeinfo => dof += edgeinfo.Value);

            foreach (var element in mesh.Elements)
            {
                for (int i = 0; i < element.NumberOfEdges; i++)
                {
                    var edge = element.Edge(i);
                    edge = (element.VertexNumber[edge.i], element.VertexNumber[edge.j]);
                    if (edge.i < edge.j) edge = (edge.j, edge.i);
                    var n = element.DOFOnEdge(i, element.Type);
                    var start = edges[edge] - n;
                    for (int j = 0; j < n; j++)
                        element.SetEdgeDOF(i, j, start + j, element.Type);
                }
            }

            foreach (var element in mesh.Elements)
                for (int i = 0; i < element.DOFOnElement(); i++)
                    element.SetElementDOF(i, dof++);

            mesh.NumberOfDofs = dof;
        }
        public static SortedSet<int>[] BuildPortraitFirstStep(IFiniteElementMesh mesh)
        {
            var a = new SortedSet<int>[mesh.NumberOfDofs];
            for (int i = 0; i < mesh.NumberOfDofs; i++) a[i] = new();
            foreach (var element in mesh.Elements)
            {
                for (int i = 0; i < element.Dofs.Length; i++)
                    for (int j = 0; j < element.Dofs.Length; j++)
                        a[element.Dofs[i]].Add(element.Dofs[j]);
            }
            return a;
        }

        static SortedSet<int> GetNodeNumbersOfScalarElements(IFiniteElementMesh mesh)
        {
            var set = new SortedSet<int>();

            foreach(var element in mesh.Elements)
            {
                if (element.Type != IFiniteElement.ElementType.Vector)
                {
                    set.UnionWith(element.VertexNumber);
                }
            }

            return set;
        }
    }
}

