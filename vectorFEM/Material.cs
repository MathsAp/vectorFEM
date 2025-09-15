using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public class Material : IMaterial
    {
        public Material(MaterialType type)
        {
            Type = type;
        }

        public MaterialType Type { get; }

        public bool IsVolume => Type is MaterialType.Volume;

        public bool Is1 => Type is MaterialType.FirstBoundary;

        public bool Is2 => Type is MaterialType.SecondBoundary;

        public bool Is3 => Type is MaterialType.ThirdBoundary;

        public bool IsInterface => Type is MaterialType.Interface;

        public Func<Vector3D, double>? Lambda { get; init; }

        public Func<Vector3D, double>? Sigma { get; init; }

        public Func<Vector3D, double>? Epsilon { get; init; }

        public Func<Vector3D, double>? Mu { get; init; }

        public Func<Vector3D, double>? Betta { get; init; }

        public Func<Vector3D, double, double>? UBetta { get; init; }

        public Func<Vector3D, double, double>? Theta { get; init; }

        public Func<Vector3D, double, Vector3D>? Htheta { get; init; }

        public Func<Vector3D, double, Vector3D>? Hext { get; init; }

        public Func<Vector3D, double, double>? Ug { get; init; }

        public Func<Vector3D, double, Vector3D>? Ag { get; init; }

        public Func<Vector3D, double, double>? F { get; init; }

        public Func<Vector3D, double, Vector3D>? Fv { get; init; }
    }
}
