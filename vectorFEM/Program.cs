using Core;
using Quadratures;
using System.Runtime.CompilerServices;

double func(double x, double y, double z)
{
    return 10 * Math.Pow(x, 9) * 10 * Math.Pow(y, 9) * 10 * Math.Pow(z, 9);
}

QuadratureNodes<Vector3D> QuadratureNodes;

QuadratureNodes = NumericalIntegration.FactoryQuadratures3D(9, ElemType.Cube);

double x0 = -1;
double x1 = 2;

double y0 = 1;
double y1 = 2;

double z0 = -3;
double z1 = 1;

double hx = x1 - x0;
double hy = y1 - y0;
double hz = z1 - z0;

int n = QuadratureNodes.Nodes.Length;

double result = 0; 
for (int i = 0; i < n; ++i)
{
    var node = QuadratureNodes.Nodes[i].Node;
    result += QuadratureNodes.Nodes[i].Weight * func(node.X * hx + x0, node.Y * hy + y0, node.Z * hz + z0);
}

result *= hx * hy * hz;

Console.WriteLine($"Значение интеграла = {result} \n");


