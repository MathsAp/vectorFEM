using Core;
using Core.ScalarFiniteElements.FiniteElements2D;
using Core.ScalarFiniteElements.FiniteElements3D;
using Core.VectorFiniteElements.FiniteElements2D;
using Core.VectorFiniteElements.FiniteElements3D;
using FEM;
using Quadratures;
using System.Runtime.CompilerServices;

//double func(double x, double y, double z)
//{
//    return 10 * Math.Pow(x, 9) * 10 * Math.Pow(y, 9) * 10 * Math.Pow(z, 9);
//}

//QuadratureNodes<Vector3D> QuadratureNodes;

//QuadratureNodes = NumericalIntegration.FactoryQuadratures3D(9, ElemType.Cube);

//double x0 = -1;
//double x1 = 2;

//double y0 = 1;
//double y1 = 2;

//double z0 = -3;
//double z1 = 1;

//double hx = x1 - x0;
//double hy = y1 - y0;
//double hz = z1 - z0;

//int n = QuadratureNodes.Nodes.Length;

//double result = 0; 
//for (int i = 0; i < n; ++i)
//{
//    var node = QuadratureNodes.Nodes[i].Node;
//    result += QuadratureNodes.Nodes[i].Weight * func(node.X * hx + x0, node.Y * hy + y0, node.Z * hz + z0);
//}

//result *= hx * hy * hz;

//Console.WriteLine($"Значение интеграла = {result} \n ");

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

Dictionary<int, string> materialNumbers = new()
{
    { 0, "vectorVolume" },
    { 1, "scalarVolume" },
    { 2, "1" },
    { 3, "2" },
    { 4, "interface" }
};

RegularParallelepipedalFiniteElementMesh mesh = new("C:\\Users\\bossf\\source\\repos\\vectorFEM\\vectorFEM\\Mesh", materialNumbers);

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
Func<Vector3D, double, Vector3D> Hext = (x, t) => Vector3D.Zero;

Dictionary<string, IMaterial> Materials = new()
{
    { "vectorVolume", new Material(true, false, false, false, x => 0, x => 0, x => 0, x => Constants.Mu0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
    { "scalarVolume", new Material(true, false, false, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, Hext, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
    { "1", new Material(false, true, false, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
    { "2", new Material(false, false, true, false, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) },
    { "interface", new Material(false, false, false, true, x => 0, x => 0, x => 0, x => 0, (x, t) => 0, (x, t) => Vector3D.Zero, Hext, (x, t) => 0, (x, t) => Vector3D.Zero, (x, t) => 0, (x, t) => Vector3D.Zero) }
};

double[] t = [0, 2, 4, 6, 8, 10];

TimeMesh timeMesh = new(t);

//EllipticProblem problem = new(mesh, Materials);

//VectorParabolicProblem problem = new(mesh, timeMesh, x => Vector3D.Zero, Materials);

VectorScalarEllipticProblem problem = new(mesh, Materials);

//VectorEllipticProblem problem = new(mesh, Materials);

problem.Prepare();
//EllipticSolution solution = new(mesh);
VectorScalarEllipticSolution solution = new(mesh, Materials);
problem.Solve(solution);

//Func<Vector3D, double> RealFunc = x => x.X * x.Y * x.Z + x.X * x.Y + x.X * x.Z + x.Y * x.Z + x.X + x.Y + x.Z + 5;
//Func<Vector3D, double, Vector3D> RealFunc = (x, t) => new Vector3D(x.Y * t, x.Z * t, x.X * t);
Func<Vector3D, double, Vector3D> RealFunc = (x, t) => new Vector3D(t * (x.Y * x.Z + x.Y + x.Z + 5), t * (x.X * x.Z + x.X + x.Z + 8), t * (x.X * x.Y + x.X + x.Y + 11));
//Func<Vector3D, Vector3D> RealGradientFunc = x => new Vector3D(x.Y * x.Z + x.Y + x.Z + 1, x.X * x.Z + x.X + x.Z + 1, x.X * x.Y + x.X + x.Y + 1);
Func<Vector3D, double, Vector3D> RealCurlFunc = (x, t) => new Vector3D(0, 0, 0);


bool flag = true;

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

//    Console.WriteLine($"time = {solution.Time}");
//    Console.WriteLine($"Значение в точке численного решения ({x}; {y}; {z}) = " + solution.Vector(point));
//    Console.WriteLine($"Значение в точке реального решения ({x}; {y}; {z}) = " + RealFunc(point, solution.Time));

//    Vector3D curlNumerical = solution.Curl(point);
//    Vector3D curlReal = RealCurlFunc(point, solution.Time);
//    Console.WriteLine($"Ротор численного решения в точке ({x}; {y}; {z}) = " + curlNumerical);
//    Console.WriteLine($"Ротор реального решения в точке ({x}; {y}; {z}) = " + curlReal);

//    Console.WriteLine("Хотите продолжить?");

//    flag = bool.Parse(Console.ReadLine()!);
//}

while (flag) // curr
{
    Console.WriteLine("Введите x: ");
    double x = double.Parse(Console.ReadLine()!);

    Console.WriteLine("Введите y: ");
    double y = double.Parse(Console.ReadLine()!);

    Console.WriteLine("Введите z: ");
    double z = double.Parse(Console.ReadLine()!);

    Vector3D point = new Vector3D(x, y, z);

   // Console.WriteLine($"Значение в точке численного A ({x}; {y}; {z}) = " + solution.Vector(point));
    //Console.WriteLine($"Значение в точке численного V ({x}; {y}; {z}) = " + solution.Value(point));

    Console.WriteLine($"Значение в точке численного H ({x}; {y}; {z}) = " + solution.H(point));
    Console.WriteLine($"Значение в точке реального H ({x}; {y}; {z}) = " + RealHext(point, 1));

    Console.WriteLine($"Значение в точке численного B ({x}; {y}; {z}) = " + solution.B(point));
    Console.WriteLine($"Значение в точке реального B ({x}; {y}; {z}) = " + Constants.Mu0 * RealHext(point, 1));

    Console.WriteLine("Хотите продолжить?");

    flag = bool.Parse(Console.ReadLine()!);
}

IMasterElement<Vector2D, double, Vector3D> ME = SquareMasterElementLinearScalarBasis.GetInstance();
var MEV = SquareMasterElementLinearVectorBasis.GetInstance();

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