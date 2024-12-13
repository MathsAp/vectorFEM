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
        public Material(bool isVolume, bool is1, bool is2, bool isInterface, Func<Vector3D, double> lambda, Func<Vector3D, double> sigma, Func<Vector3D, double> epsilon, Func<Vector3D, double> mu, Func<Vector3D, double, double> theta, Func<Vector3D, double, Vector3D> htheta, Func<Vector3D, double, Vector3D> hext, Func<Vector3D, double, double> ug, Func<Vector3D, double, Vector3D> ag, Func<Vector3D, double, double> f, Func<Vector3D, double, Vector3D> fv)
        {
            IsVolume = isVolume;
            Is1 = is1;
            Is2 = is2;
            IsInterface = isInterface;
            Lambda = lambda;
            Sigma = sigma;
            Epsilon = epsilon;
            Mu = mu;
            Theta = theta;
            Htheta = htheta;
            Hext = hext;
            Ug = ug;
            Ag = ag;
            F = f;
            Fv = fv;
        }

        public bool IsVolume { get; }

        public bool Is1 { get; }

        public bool Is2 { get; }

        public bool IsInterface { get; }

        public Func<Vector3D, double> Lambda { get; }

        public Func<Vector3D, double> Sigma { get; }

        public Func<Vector3D, double> Epsilon { get; }

        public Func<Vector3D, double> Mu { get; }

        public Func<Vector3D, double, double> Theta { get; }

        public Func<Vector3D, double, Vector3D> Htheta { get; }

        public Func<Vector3D, double, Vector3D> Hext { get; }

        public Func<Vector3D, double, double> Ug { get; }

        public Func<Vector3D, double, Vector3D> Ag { get; }

        public Func<Vector3D, double, double> F { get; }

        public Func<Vector3D, double, Vector3D> Fv { get; }
    }
}
