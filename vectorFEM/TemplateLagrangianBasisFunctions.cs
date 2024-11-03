using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public static class LinearBasis
    {
        public static readonly Func<double, double>[] Phi =
        {
            x => 1 - x,
            x => x
        };
    }

    public static class LinearBasisDerivatives
    {
        public static readonly Func<double, double>[] Phi =
        {
            x => -1,
            x => 1
        };
    }

    public static class QuadraticBasis
    {
        public static readonly Func<double, double>[] Phi =
        {
            x => 2 * (x - 1 / 2d) * (x - 1),
            x => -4 * x * (x - 1),
            x => 2 * x * (x - 1 / 2d)
        };
    }

    public static class QuadraticBasisDerivatives
    {
        public static readonly Func<double, double>[] Phi =
        {
            x => 4 * x - 3,
            x => 4 - 8 * x,
            x => 4 * x - 1
        };
    }
    public static class CubicBasis
    {
        public static readonly Func<double, double>[] Phi =
        {
            x => -9 / 2d * (x - 1 / 3d) * (x - 2 / 3d) * (x - 1),
            x => 27 / 2d * x * (x - 2 / 3d) * (x - 1),
            x => -27 / 2d * x * (x - 1 / 3d) * (x - 1),
            x => 9 / 2d * x * (x - 1 / 3d) * (x - 2 / 3d)
        };
    }

    public static class CubicBasisDerivatives
    {
        public static readonly Func<double, double>[] Phi =
        {
            x => 1 / 2d * (36 * x - 27 * x * x - 11),
            x => 1 / 2d * (81 * x * x - 90 * x + 18),
            x => 1 / 2d * (72 * x - 81 * x * x - 9),
            x => 1 / 2d * (27 * x * x - 18 * x + 2)
        };
    }
}
