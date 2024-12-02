using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public static class Linear2DVectorBasis
    {
        public static readonly Func<double, double, Vector2D>[] Phi =
        {
            (x, y) => new Vector2D(LinearBasis.Phi[0](y), 0d),
            (x, y) => new Vector2D(LinearBasis.Phi[1](y), 0d),

            (x, y) => new Vector2D(0d, LinearBasis.Phi[0](x)),
            (x, y) => new Vector2D(0d, LinearBasis.Phi[1](x))
        };
    }

    public static class Linear3DVectorBasis
    {
        public static readonly Func<double, double, double, Vector3D>[] Phi =
        {
            (x, y, z) => new Vector3D(LinearBasis.Phi[0](y) * LinearBasis.Phi[0](z), 0d, 0d),
            (x, y, z) => new Vector3D(LinearBasis.Phi[1](y) * LinearBasis.Phi[0](z), 0d, 0d),
            (x, y, z) => new Vector3D(LinearBasis.Phi[0](y) * LinearBasis.Phi[1](z), 0d, 0d),
            (x, y, z) => new Vector3D(LinearBasis.Phi[1](y) * LinearBasis.Phi[1](z), 0d, 0d),

            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[0](x) * LinearBasis.Phi[0](z), 0d), //
            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[1](x) * LinearBasis.Phi[0](z), 0d),
            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[0](x) * LinearBasis.Phi[1](z), 0d), //
            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[1](x) * LinearBasis.Phi[1](z), 0d),

            (x, y, z) => new Vector3D(0d, 0d, LinearBasis.Phi[0](x) * LinearBasis.Phi[0](y)), //
            (x, y, z) => new Vector3D(0d, 0d, LinearBasis.Phi[1](x) * LinearBasis.Phi[0](y)),
            (x, y, z) => new Vector3D(0d, 0d, LinearBasis.Phi[0](x) * LinearBasis.Phi[1](y)), //
            (x, y, z) => new Vector3D(0d, 0d, LinearBasis.Phi[1](x) * LinearBasis.Phi[1](y))
        };

        public static readonly Func<double, double, double, Vector3D>[] curlPhi =
        {
            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[0](y) * LinearBasisDerivatives.Phi[0](z), -LinearBasisDerivatives.Phi[0](y) * LinearBasis.Phi[0](z)),
            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[1](y) * LinearBasisDerivatives.Phi[0](z), -LinearBasisDerivatives.Phi[1](y) * LinearBasis.Phi[0](z)),
            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[0](y) * LinearBasisDerivatives.Phi[1](z), -LinearBasisDerivatives.Phi[0](y) * LinearBasis.Phi[1](z)),
            (x, y, z) => new Vector3D(0d, LinearBasis.Phi[1](y) * LinearBasisDerivatives.Phi[1](z), -LinearBasisDerivatives.Phi[1](y) * LinearBasis.Phi[1](z)),

            (x, y, z) => new Vector3D(-LinearBasis.Phi[0](x) * LinearBasisDerivatives.Phi[0](z), 0d, LinearBasisDerivatives.Phi[0](x) * LinearBasis.Phi[0](z)),
            (x, y, z) => new Vector3D(-LinearBasis.Phi[1](x) * LinearBasisDerivatives.Phi[0](z), 0d, LinearBasisDerivatives.Phi[1](x) * LinearBasis.Phi[0](z)),
            (x, y, z) => new Vector3D(-LinearBasis.Phi[0](x) * LinearBasisDerivatives.Phi[1](z), 0d, LinearBasisDerivatives.Phi[0](x) * LinearBasis.Phi[1](z)),
            (x, y, z) => new Vector3D(-LinearBasis.Phi[1](x) * LinearBasisDerivatives.Phi[1](z), 0d, LinearBasisDerivatives.Phi[1](x) * LinearBasis.Phi[1](z)),

            (x, y, z) => new Vector3D(LinearBasis.Phi[0](x) * LinearBasisDerivatives.Phi[0](y), -LinearBasisDerivatives.Phi[0](x) * LinearBasis.Phi[0](y), 0d),
            (x, y, z) => new Vector3D(LinearBasis.Phi[1](x) * LinearBasisDerivatives.Phi[0](y), -LinearBasisDerivatives.Phi[1](x) * LinearBasis.Phi[0](y), 0d),
            (x, y, z) => new Vector3D(LinearBasis.Phi[0](x) * LinearBasisDerivatives.Phi[1](y), -LinearBasisDerivatives.Phi[0](x) * LinearBasis.Phi[1](y), 0d),
            (x, y, z) => new Vector3D(LinearBasis.Phi[1](x) * LinearBasisDerivatives.Phi[1](y), -LinearBasisDerivatives.Phi[1](x) * LinearBasis.Phi[1](y), 0d)

        };
    }
}
