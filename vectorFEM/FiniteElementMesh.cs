using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FEM;

namespace Core;

public class FiniteElementMesh : IFiniteElementMesh
{

    public FiniteElementMesh(IEnumerable<IFiniteElement> elements, Vector3D[] vertex)
    {
        Elements = elements;
        Vertex = vertex;
        FacePortrait = FemAlgorithms.BuildFacePortrait(this);
    }

    public IEnumerable<IFiniteElement> Elements { get; }

    public Vector3D[] Vertex { get; }

    public int NumberOfDofs { get; set; }

    public IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait { get; }

}
