using ReactiveUI;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core;

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

public static class TrianglesLinearBasis
{
    public static readonly Func<Vector2D, double>[] Phi =
    {
        p => 1 - p.X - p.Y,
        p => p.Y,
        p => p.X
    };
}

public static class TrianglesLinearBasisGradients
{
    public static readonly Func<Vector2D, Vector2D>[] Phi =
    {
        p => new(-1, -1),
        p => new(0, 1),
        p => new(1, 0)
    };
}

public static class TrianglesQuadraticBasis
{
    static readonly Func<Vector2D, double>[] phi = TrianglesLinearBasis.Phi;

    public static readonly Func<Vector2D, double>[] Phi =
    {
        p => phi[0](p) * (2 * phi[0](p) - 1),
        p => phi[1](p) * (2 * phi[1](p) - 1),
        p => phi[2](p) * (2 * phi[2](p) - 1),
        p => 4 * phi[0](p) * phi[1](p),
        p => 4 * phi[1](p) * phi[2](p),
        p => 4 * phi[0](p) * phi[2](p)
    };
}

public static class TrianglesQuadraticBasisGradients
{
    static readonly Func<Vector2D, double>[] phi = TrianglesLinearBasis.Phi;

    public static readonly Func<Vector2D, Vector2D>[] Phi =
    {
        p => new(1 - 4 * phi[0](p), 1 - 4 * phi[0](p)),
        p => new(0, 4 * phi[1](p) - 1),
        p => new(4 * phi[2](p) - 1, 0),
        p => new(-4 * phi[1](p), 4 * (phi[0](p) - phi[1](p))),
        p => new(4 * phi[1](p), 4 * phi[2](p)),
        p => new(4 * (phi[0](p) - phi[2](p)), -4 * phi[2](p))
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
