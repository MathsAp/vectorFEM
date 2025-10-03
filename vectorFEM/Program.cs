using Core;
using Core.ScalarFiniteElements.FiniteElements1D;
using Core.ScalarFiniteElements.FiniteElements2D;
using Core.ScalarFiniteElements.FiniteElements3D;
using Core.VectorFiniteElements.FiniteElements2D;
using Core.VectorFiniteElements.FiniteElements3D;
using DynamicData.Aggregation;
using FEM;
using Core.Quadratures;
using System.Drawing;
using System.Globalization;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;
using TelmaQuasarCommon;
using TelmaQuasarCommon.Core.BasisFunction;
using Constants = Core.Constants;
using Vector3D = Core.Vector3D;

CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;

/* Тестирование расчета плотности тока от обмоток */
//var calculator = new TelmaCoilCalculator();
//await calculator.Load("D:\\Telma\\Telma\\CubeTest\\CubeTest.TelmaProject", "my1");

//Console.WriteLine(calculator.Hext(Vector3D.Zero, 10));


// Тестирование интегрирования для нормы
/*
double x00 = -0.2;
double hx = 0.2 - x00;

double y00 = -0.2;
double hy = 0.2 - y00;

double z0 = -0.2;
double hz = 0.2 - z0;

Func<Vector3D, Vector3D> u = x => new Vector3D(x.Y * x.Z + x.Y + x.Z + 10, x.X * x.Z + x.X + x.Z + 2, x.X * x.Y + x.X + x.Y - 6);
Vector3D LocalU(Vector3D point) => u(new(point.X * hx + x00, point.Y * hy + y00, point.Z * hz + z0));

var masterElement = CubeMasterElementLinearVectorBasis.GetInstance();
int N = 12;
double[,] M = new double[N, N];
double[] b = new double[N];
var nodes = masterElement.QuadratureNodes.Nodes;

for (int i = 0; i < N; ++i)
{
    for (int j = 0; j < N; ++j)
    {
        double sum = 0;
        for (int k = 0; k < nodes.Length; ++k)
        {
            sum += masterElement.PsiPsiMatrix[k, i, j];
        }

        M[i, j] = sum * hx * hy * hz;
    }
}

for (int i = 0; i < N; ++i)
{
    double sum = 0;
    for (int k = 0; k < nodes.Length; ++k)
    {
        sum += nodes[k].Weight * masterElement.PsiValues[i, k] * LocalU(nodes[k].Node);
    }

    b[i] = sum * hx * hy * hz;
}


int[] BuildLocalDofsForLocalMatrix(int N)
{
    var localDofs = new int[N];

    for (int i = 0; i < N; ++i)
        localDofs[i] = i;

    return localDofs;
}

SortedSet<int>[] BuildProfileForLocalMatrix(int N)
{
    var profile = new SortedSet<int>[N];

    for (int i = 0; i < N; ++i)
    {
        profile[i] = new();

        for (int j = 0; j < N; ++j)
            profile[i].Add(j);
    }

    return profile;
}

var localDofs = BuildLocalDofsForLocalMatrix(N);

var SLAE = new PardisoSLAE(new PardisoMatrix(BuildProfileForLocalMatrix(N), Quasar.Native.PardisoMatrixType.SymmetricIndefinite));
SLAE.Matrix.AddLocal(localDofs, localDofs, M);
SLAE.AddLocalRightPart(localDofs, b);

double[] sol;

using (PardisoSLAESolver SLAESolver = new PardisoSLAESolver(SLAE))
{
    SLAESolver.Prepare();

    sol = SLAESolver.Solve();
}

double value = 0;

for (int k = 0; k < nodes.Length; ++k)
{
    Vector3D sum = Vector3D.Zero;

    for (int i = 0; i < N; ++i)
    {
        sum += sol[i] * masterElement.PsiValues[i, k];
    }

    Vector3D diff = sum - LocalU(nodes[k].Node);

    value += nodes[k].Weight * diff * diff;
}

value *= hx * hy * hz;

Console.WriteLine($"Value = {Math.Sqrt(value)}");

*/

//RegularRectangularFiniteElementMesh m = new("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh", RegularRectangularFiniteElementMesh.Dimension.D2);

IrregularTriangularFiniteElementMesh m = new("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh\\NewMesh'");
m.CreateFilesWithMesh("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh");
double func(double x, double y, double z)
{
    return 10 * Math.Pow(x, 9) * 10 * Math.Pow(y, 9) * 10 * Math.Pow(z, 9);
}

QuadratureNodes<Vector3D> QuadratureNodes;

QuadratureNodes = NumericalIntegration.FactoryQuadratures3D(9, ElemType.Cube);

//IDictionary<string, IMaterial> mats = MaterialsFactory.CreateMaterials("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh");


//Material mat = new(MaterialType.Volume);
//mat.UBetta = (p, t) => p.Norm + t;

//IMaterial imat = mat;

//Console.WriteLine($"{mats["name"].Lambda(new(1, 2, 3))}, {mats["name"].F(new(1, 2, 3), 10)}");

//(ITimeMesh tm, Func<Vector3D, double> initFunc) = TimeMeshFactory<double>.CreateTimeMesh("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh");

//Console.WriteLine(initFunc(new(1, 2, 3)));

double x00 = 2;
double x1 = -1;

double y00 = 2;
double y1 = 1;

double z0 = 1;
double z1 = -3;

double hx = x1 - x00;
double hy = y1 - y00;
double hz = z1 - z0;

int n = QuadratureNodes.Nodes.Length;

double result = 0;
for (int i = 0; i < n; ++i)
{
    var node = QuadratureNodes.Nodes[i].Node;
    result += QuadratureNodes.Nodes[i].Weight * func(node.X * hx + x00, node.Y * hy + y00, node.Z * hz + z0);
}

result *= hx * hy * hz;

Console.WriteLine($"Значение интеграла = {result} \n ");

//Vector3D vec = Vector3D.Zero;

//var vec2 = new Vector3D(1, 2, 3);

//vec += vec2;

//Console.WriteLine($"Vec = {vec}, VecZero = {Vector3D.Zero}");

//Console.WriteLine($"Vec0 = {LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(0)}\n" +
//                  $"Vec1 = {LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(1)}\n" +
//                  $"Vec2 = {LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(2)}\n" +
//                  $"Vec3 = {LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(3)}\n" +
//                  $"Vec4 = {LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(4)}\n" +
//                  $"Vec5 = {LinearVectorParallelepipedalFiniteElementWithNumInteg.GetNormal(5)}");

//int n = 12;

//var quadratureNodes = NumericalIntegration.FactoryQuadratures2D(5, ElemType.Rectangle);

//Vector3D[,] psiValues = new Vector3D[n, quadratureNodes.Nodes.Length];

//Vector2D[] p = new Vector2D[quadratureNodes.Nodes.Length];

//Vector3D[] CubeVertex = [new Vector3D(0, 0, 0), new Vector3D(1, 0, 0), new Vector3D(0, 1, 0), new Vector3D(1, 1, 0),
//                                     new Vector3D(0, 0, 1), new Vector3D(1, 0, 1), new Vector3D(0, 1, 1), new Vector3D(1, 1, 1)];

//int[] faces = [0, 3];

//foreach (var f in faces)
//{
//    var face = LinearVectorParallelepipedalFiniteElementWithNumInteg.FaceS(f);

//    Vector3D a = CubeVertex[face[0]];
//    Vector3D b = CubeVertex[face[1]];
//    Vector3D c = CubeVertex[face[2]];
//    Vector3D d = CubeVertex[face[3]];

//    var localDofs = LinearVectorParallelepipedalFiniteElementWithNumInteg.LocalEdgesDofsAtFace(f);

//    if (f == 0)
//        for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
//        {
//            var node = quadratureNodes.Nodes[j].Node;

//            var newNode = a * (1 - node.X) * (1 - node.Y) + b * node.X * (1 - node.Y) +
//                            c * (1 - node.X) * node.Y + d * node.X * node.Y;

//            psiValues[0, j] = Linear3DVectorBasis.Phi[10](newNode.X, newNode.Y, newNode.Z);
//        }

//    if (f == 3)
//        for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
//        {
//            var node = quadratureNodes.Nodes[j].Node;

//            var newNode = a * (1 - node.X) * (1 - node.Y) + b * node.X * (1 - node.Y) +
//                            c * (1 - node.X) * node.Y + d * node.X * node.Y;

//            psiValues[1, j] = Linear3DVectorBasis.Phi[10](newNode.X, newNode.Y, newNode.Z);
//        }
//}

//for (int j = 0; j < quadratureNodes.Nodes.Length; ++j)
//{
//    var node = quadratureNodes.Nodes[j].Node;

//    p[j] = Linear2DVectorBasis.Phi[1](node.X, node.Y);
//}

//Vector3D[,,] facesPsiValues = SquareMasterElementLinearVectorBasis.GetInstance().FacesPsiValues;

//Console.WriteLine($"Vec = {vec}, VecZero = {Vector3D.Zero}");

//int[] locDofs = new int[5];
//SortedSet<int>[] profile = new SortedSet<int>[5];

//double[,] M = { { 10d, 3d, 2d, 1d, 0d },
//                { 3d, 10d, 3d, 2d, 1d },
//                { 2d, 3d, 10d, 3d, 2d },
//                { 1d, 2d, 3d, 10d, 3d },
//                { 0d, 1d, 2d, 3d, 10d } };

//double[] LRP = { 26, 45, 60, 69, 70 };

//for (int i = 0; i < 5; ++i)
//{
//    locDofs[i] = i;
//    profile[i] = new();

//    for (int j = 0; j < 5; ++j)
//        profile[i].Add(j);
//}

//var SLAE = new PardisoSLAE(new PardisoMatrix(profile, Quasar.Native.PardisoMatrixType.SymmetricIndefinite));
//SLAE.Matrix.AddLocal(locDofs, M);
//SLAE.AddLocalRightPart(locDofs, LRP);

//using (PardisoSLAESolver SLAESolver = new PardisoSLAESolver(SLAE))
//{
//    SLAESolver.Prepare();

//    var solution = SLAESolver.Solve();

//    Console.WriteLine($"Vec = {vec}, VecZero = {Vector3D.Zero}");
//}




Vector3D[] vertex = { new Vector3D(0, 0, 0), new Vector3D(3, 0, 0), new Vector3D(6, 0, 0), // 0 1 2
                      new Vector3D(0, 3, 0), new Vector3D(3, 3, 0), new Vector3D(6, 3, 0), // 3 4 5
                      new Vector3D(0, 6, 0), new Vector3D(3, 6, 0), new Vector3D(6, 6, 0), // 6 7 8

                      new Vector3D(0, 0, 3), new Vector3D(3, 0, 3), new Vector3D(6, 0, 3), // 9 10 11
                      new Vector3D(0, 3, 3), new Vector3D(3, 3, 3), new Vector3D(6, 3, 3), // 12 13 14
                      new Vector3D(0, 6, 3), new Vector3D(3, 6, 3), new Vector3D(6, 6, 3), // 15 16 17

                      new Vector3D(0, 0, 6), new Vector3D(3, 0, 6), new Vector3D(6, 0, 6), // 18 19 20
                      new Vector3D(0, 3, 6), new Vector3D(3, 3, 6), new Vector3D(6, 3, 6), // 21 22 23
                      new Vector3D(0, 6, 6), new Vector3D(3, 6, 6), new Vector3D(6, 6, 6) }; // 24 25 26

//IFiniteElement[] elements = [ new LinearParallelepipedalFiniteElementWithNumInteg("volume", [0, 1, 3, 4, 9, 10, 12, 13]), new LinearParallelepipedalFiniteElementWithNumInteg("volume", [1, 2, 4, 5, 10, 11, 13, 14]),
//                              new LinearParallelepipedalFiniteElementWithNumInteg("volume", [3, 4, 6, 7, 12, 13, 15, 16]), new LinearParallelepipedalFiniteElementWithNumInteg("volume", [4, 5, 7, 8, 13, 14, 16, 17]),
//                              new LinearParallelepipedalFiniteElementWithNumInteg("volume", [9, 10, 12, 13, 18, 19, 21, 22]), new LinearParallelepipedalFiniteElementWithNumInteg("volume", [10, 11, 13, 14, 19, 20, 22, 23]),
//                              new LinearParallelepipedalFiniteElementWithNumInteg("volume", [12, 13, 15, 16, 21, 22, 24, 25]), new LinearParallelepipedalFiniteElementWithNumInteg("volume", [13, 14, 16, 17, 22, 23, 25, 26]),

//                              new LinearRectangularFiniteElementWithNumInteg("1", [0, 3, 9, 12]), new LinearRectangularFiniteElementWithNumInteg("1", [3, 6, 12, 15]),
//                              new LinearRectangularFiniteElementWithNumInteg("1", [9, 12, 18, 21]), new LinearRectangularFiniteElementWithNumInteg("1", [12, 15, 21, 24]),

//                              new LinearRectangularFiniteElementWithNumInteg("2", [6, 7, 15, 16]), new LinearRectangularFiniteElementWithNumInteg("2", [7, 8, 16, 17]),
//                              new LinearRectangularFiniteElementWithNumInteg("2", [15, 16, 24, 25]), new LinearRectangularFiniteElementWithNumInteg("2", [16, 17, 25, 26]),

//                              new LinearRectangularFiniteElementWithNumInteg("3", [18, 19, 21, 22]), new LinearRectangularFiniteElementWithNumInteg("3", [19, 20, 22, 23]),
//                              new LinearRectangularFiniteElementWithNumInteg("3", [21, 22, 24, 25]), new LinearRectangularFiniteElementWithNumInteg("3", [22, 23, 25, 26]),

//                              new LinearRectangularFiniteElementWithNumInteg("4", [0, 1, 3, 4]), new LinearRectangularFiniteElementWithNumInteg("4", [1, 2, 4, 5]),
//                              new LinearRectangularFiniteElementWithNumInteg("4", [3, 4, 6, 7]), new LinearRectangularFiniteElementWithNumInteg("4", [4, 5, 7, 8]),

//                              new LinearRectangularFiniteElementWithNumInteg("5", [0, 1, 9, 10]), new LinearRectangularFiniteElementWithNumInteg("5", [1, 2, 10, 11]),
//                              new LinearRectangularFiniteElementWithNumInteg("5", [9, 10, 18, 19]), new LinearRectangularFiniteElementWithNumInteg("5", [10, 11, 19, 20]),

//                              new LinearRectangularFiniteElementWithNumInteg("6", [2, 5, 11, 14]), new LinearRectangularFiniteElementWithNumInteg("6", [5, 8, 14, 17]),
//                              new LinearRectangularFiniteElementWithNumInteg("6", [11, 14, 20, 23]), new LinearRectangularFiniteElementWithNumInteg("6", [14, 17, 23, 26]) ];

IFiniteElement[] elements = [ new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [0, 1, 3, 4, 9, 10, 12, 13]), new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [1, 2, 4, 5, 10, 11, 13, 14]),
                              new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [3, 4, 6, 7, 12, 13, 15, 16]), new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [4, 5, 7, 8, 13, 14, 16, 17]),
                              new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [9, 10, 12, 13, 18, 19, 21, 22]), new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [10, 11, 13, 14, 19, 20, 22, 23]),
                              new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [12, 13, 15, 16, 21, 22, 24, 25]), new LinearVectorParallelepipedalFiniteElementWithNumInteg("volume", [13, 14, 16, 17, 22, 23, 25, 26]),

                              new LinearVectorRectangularFiniteElementWithNumInteg("1", [0, 3, 9, 12]), new LinearVectorRectangularFiniteElementWithNumInteg("1", [3, 6, 12, 15]),
                              new LinearVectorRectangularFiniteElementWithNumInteg("1", [9, 12, 18, 21]), new LinearVectorRectangularFiniteElementWithNumInteg("1", [12, 15, 21, 24]),

                              new LinearVectorRectangularFiniteElementWithNumInteg("2", [6, 7, 15, 16]), new LinearVectorRectangularFiniteElementWithNumInteg("2", [7, 8, 16, 17]),
                              new LinearVectorRectangularFiniteElementWithNumInteg("2", [15, 16, 24, 25]), new LinearVectorRectangularFiniteElementWithNumInteg("2", [16, 17, 25, 26]),

                              new LinearVectorRectangularFiniteElementWithNumInteg("3", [18, 19, 21, 22]), new LinearVectorRectangularFiniteElementWithNumInteg("3", [19, 20, 22, 23]),
                              new LinearVectorRectangularFiniteElementWithNumInteg("3", [21, 22, 24, 25]), new LinearVectorRectangularFiniteElementWithNumInteg("3", [22, 23, 25, 26]),

                              new LinearVectorRectangularFiniteElementWithNumInteg("4", [0, 1, 3, 4]), new LinearVectorRectangularFiniteElementWithNumInteg("4", [1, 2, 4, 5]),
                              new LinearVectorRectangularFiniteElementWithNumInteg("4", [3, 4, 6, 7]), new LinearVectorRectangularFiniteElementWithNumInteg("4", [4, 5, 7, 8]),

                              new LinearVectorRectangularFiniteElementWithNumInteg("5", [0, 1, 9, 10]), new LinearVectorRectangularFiniteElementWithNumInteg("5", [1, 2, 10, 11]),
                              new LinearVectorRectangularFiniteElementWithNumInteg("5", [9, 10, 18, 19]), new LinearVectorRectangularFiniteElementWithNumInteg("5", [10, 11, 19, 20]),

                              new LinearVectorRectangularFiniteElementWithNumInteg("6", [2, 5, 11, 14]), new LinearVectorRectangularFiniteElementWithNumInteg("6", [5, 8, 14, 17]),
                              new LinearVectorRectangularFiniteElementWithNumInteg("6", [11, 14, 20, 23]), new LinearVectorRectangularFiniteElementWithNumInteg("6", [14, 17, 23, 26]) ];

//FiniteElementMesh mesh = new(elements, vertex);

//Dictionary<int, string> materialNumbers = new()
//{
//    { 0, "volume" },
//    { 1, "1" },
//    { 2, "2" },
//    { 3, "3" },
//    { 4, "4" },
//    { 5, "5" },
//    { 6, "6" }
//};

//Dictionary<int, string> materialNumbers = new()
//{
//    { 0, "vectorVolume" },
//    { 1, "2" },
//    { 2, "3" },
//    { 3, "4" },
//    { 4, "5" },
//    { 5, "6"},
//    { 6, "7" }
//};

Dictionary<int, string> materialNumbers = new()
{
    { 0, "vectorVolume" },
    { 1, "scalarVolume" },
    { 2, "1" },
    { 3, "2" },
    { 4, "interface" },
};

//RegularParallelepipedalFiniteElementMesh mesh = new("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh", materialNumbers);

//Dictionary<string, IMaterial> Materials = new()
//{
//    { "volume", new Material(true, false, false, x => 4, x => 2, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 2 * (x.X * x.Y * x.Z + x.X * x.Y + x.X * x.Z + x.Y * x.Z + x.X + x.Y + x.Z + 5), (x, t) => Vector3D.Zero) },
//    { "1", new Material(false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => -4 * (x.Y * x.Z + x.Y + x.Z + 1), (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "2", new Material(false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => 4 * (x.X * x.Z + x.X + x.Z + 1), (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "3", new Material(false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => 4 * (x.X * x.Y + x.X + x.Y + 1), (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => x.X * x.Y + x.X + x.Y + 5, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "5", new Material(false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => x.X * x.Z + x.X + x.Z + 5, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "6", new Material(false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 7 * (x.Y * x.Z + x.Y + x.Z) + 11, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) }
//};

//Dictionary<string, IMaterial> Materials = new()
//{
//    { "volume", new Material(true, false, false, x => 0, x => 2, x => 0, x => 10, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D(2 * x.Y, 2 * x.Z, 2 * x.X)) },
//    { "1", new Material(false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => new Vector3D(-t / 10d, -t / 10d, -t / 10d), (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "2", new Material(false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => new Vector3D(-t / 10d, -t / 10d, -t / 10d), (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "3", new Material(false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => new Vector3D(-t / 10d, -t / 10d, -t / 10d), (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D(x.Y * t, 0, x.X * t), (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "5", new Material(false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D(0, x.Z * t, x.X * t), (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "6", new Material(false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D(x.Y * t, x.Z * t, 6 * t), (x, t) => 0, (x, t) => Vector3D.Zero) }
//};

//Dictionary<string, IMaterial> Materials = new()
//{
//    { "volume", new Material(true, false, false, false, x => 0, x => 4, x => 0, x => 5, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D(4 * (x.Y * x.Z + x.Y + x.Z + 5), 4 * (x.X * x.Z + x.X + x.Z + 8), 4 * (x.X * x.Y + x.X + x.Y + 11))) },
//    { "1", new Material(false, false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => new Vector3D(0, 0, 0), (x, t) => Vector3D.Zero,(x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "2", new Material(false, false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => new Vector3D(0, 0, 0), (x, t) => Vector3D.Zero,(x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "3", new Material(false, false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => new Vector3D(0, 0, 0), (x, t) => Vector3D.Zero,(x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D((x.Y + 5) * t, (x.X + 8) * t, (x.X * x.Y + x.X + x.Y + 11) * t), (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "5", new Material(false, true, false, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D((x.Z + 5) * t, (x.X * x.Z + x.X + x.Z + 8) * t, (x.X + 11) * t), (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "6", new Material(false, true, false, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => new Vector3D((x.Y * x.Z + x.Y + x.Z + 5) * t, (7 * x.Z + 14) * t, (7 * x.Y + 17) * t), (x, t) => 0, (x, t) => Vector3D.Zero) }
//};

double I = 2;
double x0 = 0;
double y0 = 0;

Func<Vector3D, double, Vector3D> RealHext = (x, t) => (I / (2 * Math.PI)) * new Vector3D(-(x.Y - y0) / ((x.X - x0) * (x.X - x0) + (x.Y - y0) * (x.Y - y0)), (x.X - x0) / ((x.X - x0) * (x.X - x0) + (x.Y - y0) * (x.Y - y0)), 0);
//Func<Vector3D, double, Vector3D> Hext = (x, t) => calculator.Hext(x, t);

// VKR
//Dictionary<string, IMaterial> Materials = new()
//{
//    { "vectorVolume", new Material(true, false, false, false, x => 0, x => 1e6, x => 0, x => Constants.Mu0 * 1000, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "scalarVolume", new Material(true, false, false, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, Hext, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "1", new Material(false, true, false, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "2", new Material(false, false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
//    { "interface", new Material(false, false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, Hext, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) }
//};

Dictionary<string, IMaterial> Materials = new()
{
    { "vectorVolume", new Material(MaterialType.Volume) },
    { "scalarVolume", new Material(MaterialType.Volume) },
    { "1", new Material(MaterialType.FirstBoundary) },
    { "2", new Material(MaterialType.SecondBoundary) },
    { "interface", new Material(MaterialType.Interface) }
};

// Первый тест
//Dictionary<string, IMaterial> Materials = new()
//{
//    { "vectorVolume", new Material(true, false, false, false, lambda => 0, sigma => 4, epsilon => 0, mu => 5, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => new Vector3D(4 * (fv.Y * fv.Z + fv.Y + fv.Z + 5), 4 * (fv.X * fv.Z + fv.X + fv.Z + 8), 4 * (fv.X * fv.Y + fv.X + fv.Y + 11))) },
//    { "2", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Z + 5, ag.X * ag.Z + ag.X + ag.Z + 8, ag.X + 11), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "3", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Y + 5, ag.X + 8, ag.X * ag.Y + ag.X + ag.Y + 11), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Y * ag.Z + ag.Y + ag.Z + 5, 5 * ag.Z + 12, 5 * (ag.Y + 3)), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "5", new Material(false, false, true, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => Vector3D.Zero) }
//};

// Второй тест
//Dictionary<string, IMaterial> Materials = new()
//{
//    { "vectorVolume", new Material(true, false, false, false, lambda => 0, sigma => 1, epsilon => 0, mu => 2, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (x, t) => new Vector3D(x.Y * x.Y * x.Z * x.Z + x.Y * x.Z + x.Y * x.Y * x.Z + x.Y * x.Z * x.Z + 3, x.X * x.X * x.Z * x.Z + x.X * x.X * x.Z + x.X * x.Z * x.Z + x.X * x.Z + 6, x.X * x.X * x.Y * x.Y + x.X * x.X * x.Y + x.X * x.Y * x.Y + x.X * x.Y + 9)) },
//    { "2", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (x, t) => new Vector3D(x.Z + 5 + x.Z * x.Z, x.X * x.X * x.Z * x.Z + x.X * x.X * x.Z + x.X * x.Z * x.Z + x.X * x.Z + x.X + x.Z + 8 + x.X * x.X + x.Z * x.Z, x.X + 11 + x.X * x.X), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "3", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (x, t) => new Vector3D(x.Y + 5 + x.Y * x.Y, x.X + 8 + x.X * x.X, x.X * x.X * x.Y * x.Y + x.X * x.X * x.Y + x.X * x.Y * x.Y + x.X * x.Y + x.X + x.Y + 11 + x.X * x.X + x.Y * x.Y), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (x, t) => new Vector3D(x.Y * x.Y * x.Z * x.Z + x.Y * x.Y * x.Z + x.Y * x.Z * x.Z + x.Y * x.Z + x.Z + x.Y + 5 + x.Y * x.Y + x.Z * x.Z, 21 * x.Z * x.Z + 21 * x.Z + 28, 21 * x.Y * x.Y + 21 * x.Y + 31), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "5", new Material(false, false, true, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (x, t) => new Vector3D(x.Y - x.Z, x.Z * (x.Y * x.Y + x.Y + 1), - x.Y * (x.Z * x.Z + x.Z + 1)), (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "6", new Material(false, false, true, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (x, t) => new Vector3D((x.Y - 4) * (x.X * x.X + x.X + 1), (4 - x.X) * (x.Y * x.Y + x.Y + 1), 21 * (x.X - x.Y)), (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "7", new Material(false, false, true, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (x, t) => new Vector3D((4 - x.Z) * (x.X * x.X + x.X + 1), 21 * (x.Z - x.X), (x.X - 4) * (x.Z * x.Z + x.Z + 1)), (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => Vector3D.Zero) }
//};


// Третий тест
//Dictionary<string, IMaterial> Materials = new()
//{
//    { "vectorVolume", new Material(true, false, false, false, lambda => 0, sigma => 5, epsilon => 0, mu => 8, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => new Vector3D(5 * (fv.Y * fv.Z + fv.Y + fv.Z + fv.X + 5), 5 * (fv.X * fv.Z + fv.X + fv.Z + fv.Y + 8), 5 * (fv.X * fv.Y + fv.X + fv.Y + fv.Z + 11))) },
//    { "2", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Z + ag.X + 5, ag.X * ag.Z + ag.X + ag.Z + ag.Y + 8, ag.X + ag.Z + 11), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "3", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Y + ag.X + 5, ag.X + ag.Y + 8, ag.X * ag.Y + ag.X + ag.Y + ag.Z + 11), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Y * ag.Z + ag.Y + ag.Z + 9, 5 * ag.Z + ag.Y + 12, 5 * (ag.Y + 3) + ag.Z), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "5", new Material(false, false, true, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => Vector3D.Zero) }
//};

// Четвертый тест
//Dictionary<string, IMaterial> Materials = new()
//{
//    { "vectorVolume", new Material(true, false, false, false, lambda => 0, sigma => sigma.X + sigma.Y + sigma.Z, epsilon => 0, mu => 5, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => new Vector3D((fv.X + fv.Y + fv.Z) * (fv.Y * fv.Z + fv.Y + fv.Z + 5), (fv.X + fv.Y + fv.Z) * (fv.X * fv.Z + fv.X + fv.Z + 8), (fv.X + fv.Y + fv.Z) * (fv.X * fv.Y + fv.X + fv.Y + 11))) },
//    { "2", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Z + 5, ag.X * ag.Z + ag.X + ag.Z + 8, ag.X + 11), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "3", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Y + 5, ag.X + 8, ag.X * ag.Y + ag.X + ag.Y + 11), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Y * ag.Z + ag.Y + ag.Z + 5, 5 * ag.Z + 12, 5 * (ag.Y + 3)), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "5", new Material(false, false, true, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => Vector3D.Zero) }
//};

//// Пятый тест 
//Dictionary<string, IMaterial> Materials = new()
//{
//    { "vectorVolume", new Material(true, false, false, false, lambda => 0, sigma => 4, epsilon => 0, mu => 1 / (mu.X + mu.Y + mu.Z), (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => new Vector3D(4 * fv.Z, 4 * fv.X, 4 * fv.Y)) },
//    { "2", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Z, ag.X, 0), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "3", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(0, ag.X, ag.Y), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "4", new Material(false, true, false, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => Vector3D.Zero, (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => new Vector3D(ag.Z, 4, ag.Y), (f, t) => 0, (fv, t) => Vector3D.Zero) },
//    { "5", new Material(false, false, true, false, lambda => 0, sigma => 0, epsilon => 0, mu => 0, (theta, t) => 0, (htheta, t) => new Vector3D(htheta.X + htheta.Y + htheta.Z, htheta.X + htheta.Y + htheta.Z, htheta.X + htheta.Y + htheta.Z), (hext, t) => Vector3D.Zero, (ug, t) => 0, (ag, t) => Vector3D.Zero, (f, t) => 0, (fv, t) => Vector3D.Zero) }
//};

//double[] t = [0, 0.19];

//RefineParams refineParams = new() { splitCount = { 19 }, stretchRatio = { 1 } };
//TimeMesh timeMesh = new(t, refineParams);

////EllipticProblem problem = new(mesh, Materials);

////VectorParabolicProblem problem = new(mesh, timeMesh, x => Vector3D.Zero, Materials);

////VectorScalarEllipticProblem problem = new(mesh, Materials);

//VectorScalarParabolicProblem problem = new(mesh, timeMesh, Materials);

////VectorEllipticProblem problem = new(mesh, Materials);

//problem.Prepare();
//VectorScalarParabolicSolution solution = new(mesh, timeMesh, Materials);
////VectorScalarEllipticSolution solution = new(mesh, Materials);
//problem.Solve(solution);

//Func<Vector3D, double> RealFunc = x => x.X * x.Y * x.Z + x.X * x.Y + x.X * x.Z + x.Y * x.Z + x.X + x.Y + x.Z + 5;
//Func<Vector3D, double, Vector3D> RealFunc = (x, t) => new Vector3D(x.Y * t, x.Z * t, x.X * t);
//Func<Vector3D, double, Vector3D> RealFunc = (x, t) => new Vector3D(t * (x.Y * x.Z + x.Y + x.Z + 5), t * (x.X * x.Z + x.X + x.Z + 8), t * (x.X * x.Y + x.X + x.Y + 11));

// Первый тест
//Func<Vector3D, Vector3D> RealFunc = x => new Vector3D((x.Y * x.Z + x.Y + x.Z + 5), (x.X * x.Z + x.X + x.Z + 8), (x.X * x.Y + x.X + x.Y + 11));

// Второй тест
//Func<Vector3D, Vector3D> RealFunc = x => new Vector3D((x.Y * x.Z + x.Y + x.Z + 5 + x.Y * x.Y * x.Z * x.Z + x.Y * x.Y * x.Z + x.Y * x.Z * x.Z + x.Y * x.Y + x.Z * x.Z), (x.X * x.Z + x.X + x.Z + 8 + x.X * x.X * x.Z * x.Z + x.X * x.X * x.Z + x.X * x.Z * x.Z + x.X * x.X + x.Z * x.Z), (x.X * x.Y + x.X + x.Y + 11 + x.X * x.X * x.Y * x.Y + x.X * x.X * x.Y + x.X * x.Y * x.Y + x.X * x.X + x.Y * x.Y));

// Третий тест 
//Func<Vector3D, Vector3D> RealFunc = x => new Vector3D((x.Y * x.Z + x.Y + x.Z + x.X + 5), (x.X * x.Z + x.X + x.Z + x.Y + 8), (x.X * x.Y + x.X + x.Y + x.Z + 11));

// Четвертый тест
//Func<Vector3D, Vector3D> RealFunc = x => new Vector3D((x.Y * x.Z + x.Y + x.Z + 5), (x.X * x.Z + x.X + x.Z + 8), (x.X * x.Y + x.X + x.Y + 11));

//// Пятый тест
//Func<Vector3D, Vector3D> RealFunc = x => new Vector3D(x.Z, x.X, x.Y);

//Func<Vector3D, Vector3D> RealGradientFunc = x => new Vector3D(x.Y * x.Z + x.Y + x.Z + 1, x.X * x.Z + x.X + x.Z + 1, x.X * x.Y + x.X + x.Y + 1);

// Первый тест
//Func<Vector3D, double, Vector3D> RealCurlFunc = (x, t) => new Vector3D(0, 0, 0);

// Второй тест
//Func<Vector3D, double, Vector3D> RealCurlFunc = (x, t) => new Vector3D(2 * (x.Y - x.Z) * (x.X * x.X + x.X + 1), 2 * (x.Z - x.X) * (x.Y * x.Y + x.Y + 1), 2 * (x.X - x.Y) * (x.Z * x.Z + x.Z + 1));

// Третий тест
//Func<Vector3D, double, Vector3D> RealCurlFunc = (x, t) => new Vector3D(0, 0, 0);

// Четвертый тест
//Func<Vector3D, double, Vector3D> RealCurlFunc = (x, t) => new Vector3D(0, 0, 0);

//// Пятый тест
//Func<Vector3D, double, Vector3D> RealCurlFunc = (x, t) => new Vector3D(1, 1, 1);



bool flag = false;

//while (flag)
//{
//    Console.WriteLine("Введите x: ");
//    double x = double.Parse(Console.ReadLine()!);

//    Console.WriteLine("Введите y: ");
//    double y = double.Parse(Console.ReadLine()!);

//    Console.WriteLine("Введите z: ");
//    double z = double.Parse(Console.ReadLine()!);

//    Vector3D point = new Vector3D(x, y, z);

//    Console.WriteLine($"Значение в точке численного решения ({x}; {y}; {z}) = " + solution.Value(point));
//    Console.WriteLine($"Значение в точке реального решения ({x}; {y}; {z}) = " + RealFunc(point));


//    Vector3D gradNumerical = solution.Gradient(point);
//    Vector3D gradReal = RealGradientFunc(point);
//    Console.WriteLine($"Градиент численного решения в точке ({x}; {y}; {z}) = " + gradNumerical);
//    Console.WriteLine($"Градиент реального решения в точке ({x}; {y}; {z}) = " + gradReal);


//    Console.WriteLine("Хотите продолжить?");

//    flag = bool.Parse(Console.ReadLine()!);
//}

//var sol = elements[24].BuildLocalRightPartWithFirstBoundaryConditions(vertex, mesh.FacePortrait, x => new Vector3D(0, x.Z, x.X));

//Console.WriteLine($"Norm L2: {solution.CalcNormL2(RealFunc)}");

//while (flag)
//{
//    Console.WriteLine("Введите время: ");
//    double time = double.Parse(Console.ReadLine()!);

//    solution.Time = time;

//    Console.WriteLine("Введите x: ");
//    double x = double.Parse(Console.ReadLine()!);

//    Console.WriteLine("Введите y: ");
//    double y = double.Parse(Console.ReadLine()!);

//    Console.WriteLine("Введите z: ");
//    double z = double.Parse(Console.ReadLine()!);

//    Vector3D point = new Vector3D(x, y, z);

//    //Console.WriteLine($"time = {solution.Time}");
//    Console.WriteLine($"Значение в точке численного решения ({x}; {y}; {z}) = " + solution.Vector(point));
//    //Console.WriteLine($"Значение в точке реального решения ({x}; {y}; {z}) = " + RealFunc(point));

//    Vector3D curlNumerical = solution.Curl(point);
//    //Vector3D curlReal = RealCurlFunc(point, 1);
//    Console.WriteLine($"Ротор численного решения в точке ({x}; {y}; {z}) = " + curlNumerical);
//    //Console.WriteLine($"Ротор реального решения в точке ({x}; {y}; {z}) = " + curlReal);

//    Console.WriteLine("Хотите продолжить?");

//    flag = bool.Parse(Console.ReadLine()!);
//}

//Vector3D[] points =  [ new Vector3D(0.09, 0, 0),       new Vector3D(0.11, 0, 0),       new Vector3D(0.3, 0, 0),
//                       new Vector3D(0, 0, 0.09),       new Vector3D(0, 0, 0.11),       new Vector3D(0, 0, 0.3),
//                       new Vector3D(0.09, 0.09, 0.09), new Vector3D(0.11, 0.11, 0.11), new Vector3D(0.3, 0.3, 0.3) ];

//WriteResultsToFile();

Vector3D[] points = [ new Vector3D(1, 1, 0), new Vector3D(5, 1, 0), new Vector3D(1, 6, 0), new Vector3D(5, 6, 0), new Vector3D(3, 3.5, 0) ];

IFiniteElement[] elems = [new QuadraticTriangularFiniteElementWithNumInteg("scalarVolume", [0, 1, 4]),
            new QuadraticTriangularFiniteElementWithNumInteg("scalarVolume", [1, 3, 4]),
new QuadraticTriangularFiniteElementWithNumInteg("scalarVolume", [3, 2, 4]),
new QuadraticTriangularFiniteElementWithNumInteg("scalarVolume", [2, 0, 4]),
new QuadraticStraightFiniteElement("firstBoundary", [0, 2]),
new QuadraticStraightFiniteElement("secondBoundary", [0, 1]),
new QuadraticStraightFiniteElement("thirdBoundary1", [1, 3]),
new QuadraticStraightFiniteElement("thirdBoundary2", [2, 3])];

FiniteElementMesh mesh2 = new(elems, points);

string path = "C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh";
string path2 = Path.Combine(path, "mesh.txt");

Func<Vector3D, double, double> uReal = (p, t) => p.X + p.Y + t * p.X * p.Y;
Func<Vector3D, double, Vector3D> graduReal = (p, t) => new(1 + t * p.Y, 1 + t * p.X, 0);

IDictionary<string, IMaterial> materials = MaterialsFactory.CreateMaterials(path);
(ITimeMesh tMesh, Func<Vector3D, double> initFunc) = TimeMeshFactory<double>.CreateTimeMesh(path);

ParabolicProblem problem = new(new RegularRectangularFiniteElementMesh(path, RegularRectangularFiniteElementMesh.Dimension.D2), tMesh, initFunc, materials, IProblem.CoordinateSystem.Cylindrical);
//ParabolicProblem problem = new(new TriangularFiniteElementMesh(path2), tMesh, initFunc, materials, IProblem.CoordinateSystem.Cylindrical);
//ParabolicProblem problem = new(mesh2, tMesh, initFunc, materials);
//EllipticProblem problem = new(new RegularRectangularFiniteElementMesh(path, RegularRectangularFiniteElementMesh.Dimension.D2), materials);
problem.Prepare();

ParabolicSolution solution = new(problem.Mesh, tMesh, path);

problem.Solve(solution);

Console.WriteLine("Введите количество точек для проверки: ");
int k = int.Parse(Console.ReadLine());

Random random = new();
for (int i = 0; i < k; ++i)
{
    int time = random.Next(1, 11);
    solution.Time = time;
    time = (int)solution.Time;
    Console.WriteLine($"Установленное время t = {solution.Time}");

    double min = 1;
    double maxX = 5;
    double maxY = 6;

    (double x, double y) = (random.NextDouble() * (maxX - min) + min, random.NextDouble() * (maxY - min) + min);

    Vector3D p = new(x, y, 0);

    Console.WriteLine($"Значение в точке {p} численного u = {solution.Value(p)}");
    Console.WriteLine($"Значение в точке {p} искомого u = {uReal(p, time)}");

    Console.WriteLine($"Значение в точке {p} численного gradu = {solution.Gradient(p)}");
    Console.WriteLine($"Значение в точке {p} искомого gradu = {graduReal(p, time)}");

}

while (!flag) // curr vkr
{
    Console.WriteLine("Введите время: ");
    double time = double.Parse(Console.ReadLine()!);

    solution.Time = time;

    Console.WriteLine($"Установленное время t = {solution.Time}");

    Console.WriteLine("Введите x: ");
    double x = double.Parse(Console.ReadLine()!);

    Console.WriteLine("Введите y: ");
    double y = double.Parse(Console.ReadLine()!);

    Console.WriteLine("Введите z: ");
    double z = double.Parse(Console.ReadLine()!);

    Vector3D point = new Vector3D(x, y, z);

    Console.WriteLine($"Значение в точке {point} численного u = {solution.Value(point)}");
    Console.WriteLine($"Значение в точке {point} искомого u = {uReal(point, time)}");

    Console.WriteLine($"Значение в точке {point} численного gradu = {solution.Gradient(point)}");
    Console.WriteLine($"Значение в точке {point} искомого gradu = {graduReal(point, time)}");

    Console.WriteLine("Хотите продолжить?");

    flag = Console.ReadLine() == "n";
}

//void WriteResultsToFile(string path = "")
//{
//    if (path.Length == 0)
//        path = "ParabolicProblemSolutionAtPoints";

//    if (Directory.Exists(path))
//        Directory.Delete(path, true);

//    Directory.CreateDirectory(path!);

//    foreach(var p in points)
//    {
//        using (StreamWriter writer = new StreamWriter(Path.Combine(path, p.ToString() + ".txt"), false))
//        {
//            for (int i = 0; i < timeMesh.Size(); ++i)
//            {
//                solution.Time = timeMesh[i];
//                writer.WriteLine($"{solution.Time:F2}    {solution.B(p).X:E12}   {solution.B(p).Y:E12}   {solution.B(p).Z:E12}   {solution.B(p).Norm:E12}");
//            }
//        }W
//    }

//    Console.WriteLine("Done");
//}


//IMasterElement<Vector2D, double, Vector3D> ME = SquareMasterElementLinearScalarBasis.GetInstance();
//var MEV = SquareMasterElementLinearVectorBasis.GetInstance();

int[] test = [1, 2, 3, 4, 5, 6, 7, 8];

Console.WriteLine(string.Join(", ", test[..4]));
Console.WriteLine(string.Join(", ", test[4..]));

//Dictionary<int, string> materialNumbers = new()
//{
//    { 0, "volume" },
//    { 1, "1" },
//    { 2, "2" },
//    { 3, "3" },
//    { 4, "4" },
//    { 5, "5" },
//    { 6, "6" }
//};

var test1 = new RegularParallelepipedalFiniteElementMesh("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh", materialNumbers);

Console.WriteLine("Хотите продолжить?");

//double[,] M = { { 1, 2 }, 
//                { 3, 4 } };

//int[] locDofs = new int[4];
//SortedSet<int>[] profile = new SortedSet<int>[4];

//for (int i = 0; i < 4; ++i)
//{
//    locDofs[i] = i;
//    profile[i] = new();

//    for (int j = 0; j < 4; ++j)
//        profile[i].Add(j);
//}

//var SLAE = new PardisoSLAE(new PardisoNonSymmMatrix(profile, Quasar.Native.PardisoMatrixType.StructurallySymmetric));
//SLAE.Matrix.AddLocalTransposed([2, 3], [0, 1], M);
//SLAE.Matrix.AddLocal([0, 1], [2, 3], M, -1);

//Console.WriteLine("Хотите продолжить?");