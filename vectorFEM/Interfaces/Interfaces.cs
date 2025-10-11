using System.Collections.Generic;
using Core;
using Core.Quadratures;
using Quasar.Native;

namespace FEM;

public interface IFiniteElement
{
    ElementType Type { get; }
    string Material { get; }
    enum MatrixType { Stiffness, Mass, Interface, Convection }
    enum ElementType { Scalar, Vector, VectorScalar };
    enum DofsType { Scalar, Vector };
    int[] VertexNumber { get; } // в порядке локальной нумерации вершин 
    void SetVertexDOF(int vertex, int dof);
    int NumberOfEdges { get; }
    int NumberOfFaces { get; }
    (int i, int j) Edge(int edge);
    int[] Face(int face);
    int DOFOnEdge(int edge, ElementType type);
    int DOFOnFace(int face, ElementType type);
    void SetEdgeDOF(int edge, int n, int dof, ElementType type);
    void SetFaceDOF(int face, int n, int dof, ElementType type);

    int DOFOnElement();
    void SetElementDOF(int n, int dof);
    int[] Dofs { get; }
    int[] GetDofs(DofsType type);
    double[,] BuildLocalMatrix(Vector3D[] VertexCoords, MatrixType type, Func<Vector3D, double> Coeff, Func<Vector3D, Vector3D>? Velocity = null); // у коэффициента первый параметр в локальных координатах элемента - зачем в локальных координатах?
                                                                                                        // Как будем понимать, интегрируем или коэффициент раскладывается, если он не постоянный?
    double[,] BuildLocalMatrix(Vector3D[] VertexCoords, MatrixType type, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, double> Coeff);
    double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, double> F); // у коэффициента первый параметр в локальных координатах элемента  
    double[] BuildLocalRightPart(Vector3D[] VertexCoords, Func<Vector3D, Vector3D> F);
    double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Ug);
    double[] BuildLocalRightPartWithFirstBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Ug);
    double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, Func<Vector3D, double> Theta);
    double[] BuildLocalRightPartWithSecondBoundaryConditions(Vector3D[] VertexCoords, IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait, Func<Vector3D, Vector3D> Theta);

    bool IsPointOnElement(Vector3D[] VertexCoords, Vector3D point); // Проверяет, принадлежит ли точка конечному элементу
    double GetValueAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point); // Получить значение в точке на конечном элементе 
    Vector3D GetVectorAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point); // Получить значение векторного поля в точке на конечном элементе
    Vector3D GetGradientAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point); // Получить градиент в точке на конечном элементе
    Vector3D GetCurlAtPoint(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Vector3D point); // Получить значение ротора векторного поля в точке на конечном элементе

    double CalclIntegralOfSquaredDifference(Vector3D[] VertexCoords, ReadOnlySpan<double> coeffs, Func<Vector3D, Vector3D> u);

}
public interface IFiniteElementWithNumericIntegration<T1, T2, T3> : IFiniteElement
{
    IMasterElement<T1, T2, T3> MasterElement { get; }
}

public interface IFiniteElementMesh
{
    IEnumerable<IFiniteElement> Elements { get; }
    IDictionary<(int, int, int, int), ((IFiniteElement?, int), (IFiniteElement?, int))> FacePortrait { get; } // Как сделать универсальным

    Vector3D[] Vertex { get; } // без повторов
    int NumberOfDofs { get; set; }
}
public interface IMasterElement<T1, T2, T3>
{
    T2[,] PsiValues { get; }
    T2[,,] FacesPsiValues { get; }
    T2[,,] FacesPsiNValues { get; }
    T3[,,] FacesGradValues { get; }
    double[,,] PsiPsiMatrix { get; }
    double[,,,] FacesPsiNPsiNMatrix { get; }
    T2[,] CurlValues { get; }
    T1[,] GradValues {  get; }
    QuadratureNodes<T1> QuadratureNodes { get; }
}
public interface ITimeMesh
{
    double this[int i] { get; } // значение времени по индексу в массиве времени

    int Size(); // количество временных слоев
    double[] Coefs(int i); // массив коэф. для слоев, по идее содержит два слоя(предыдущий (j-1) и пред предыдущий (j-2)). // уточнить по поводу ReadOnlySpan
                           // i будет принимать два значения - 1 и 2, для (j-1) и (j-2) соответственно. 

    void ChangeCoefs(double[] coefs); // добавляем только что насчитанные коэф., заменяя (j-2) на (j-1), а (j-1) на (j)

    //void AddFirstInitialCondition(double[] coefs);
    //void AddSecondInitialCondition(double[] coefs); // Поскольку используем трехслойную схему, то нам требуется два начальных условия,
    //                                                // а у нас только одно начальное условие
    //                                                // Поэтому мы должны получить второе с помощью неявной двухслойной, и добавить к весам, 
    //                                                // Чтобы дальше считать трехслойной
    bool IsChangedStep(int i); // смотрим, поменялся ли шаг по времени

    void DoubleMesh();
}

public interface IProblem
{
    enum CoordinateSystem { Cartesian, Cylindrical }
    IDictionary<string, IMaterial> Materials { get; }
    void Prepare();
    void Solve(ISolution result);
}

public enum MaterialType
{
    Volume = 1,
    FirstBoundary,
    SecondBoundary,
    ThirdBoundary,
    Interface
}

public interface IMaterial
{
    MaterialType Type { get; }

    bool IsVolume { get; }
    bool Is1 { get; }
    bool Is2 { get; }
    bool Is3 { get; }
    bool IsInterface { get; }
    Func<Vector3D, double>? Lambda { get; }
    Func<Vector3D, double>? Sigma { get; }
    Func<Vector3D, double>? Epsilon { get; }
    Func<Vector3D, double>? Mu{ get; }
    Func<Vector3D, double>? Betta { get; }
    Func<Vector3D, double, double>? UBetta { get; }
    Func<Vector3D, double, double>? Theta { get; }
    Func<Vector3D, double, Vector3D>? Htheta { get; }
    Func<Vector3D, double, Vector3D>? Hext { get; }
    Func<Vector3D, double, double>? Ug { get; }
    Func<Vector3D, double, Vector3D>? Ag { get; }
    Func<Vector3D, double, double>? F { get; }
    Func<Vector3D, double, Vector3D>? Fv { get; }
}

//public static class MaterialExtensions
//{
//    public static bool IsVolume(this IMaterial material) =>
//        (material.Type & MaterialType.Volume) != 0;

//    public static bool Is1(this IMaterial material) =>
//        (material.Type & MaterialType.FirstBoundary) != 0;

//    public static bool Is2(this IMaterial material) =>
//        (material.Type & MaterialType.SecondBoundary) != 0;

//    public static bool Is3(this IMaterial material) =>
//        (material.Type & MaterialType.ThirdBoundary) != 0;

//    public static bool IsInterface(this IMaterial material) =>
//        (material.Type & MaterialType.) != 0;
//}

public interface ISolution
{
    double Time { get; set; } // Меняется вектор решения при установке определенного времени
    IFiniteElementMesh Mesh { get; }
    ITimeMesh TimeMesh { get; }
    ReadOnlySpan<double> SolutionVector { get; }
    void AddSolutionVector(double t, double[] solution);
    double Value(Vector3D point);
    Vector3D Gradient(Vector3D point);

    double CalcNormL2(Func<Vector3D, Vector3D> u);
}

public interface IMatrix
{
    int N { get; }
    void SetProfile(SortedSet<int>[] profile);
    void AddLocal(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1d);
    void AddLocalTransposed(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1d);
    void Symmetrize(int dof, double value, double[] RightPart);
    void MultiplyByVector(ReadOnlySpan<double> x, double[] result);
    void Clear();
}

public interface ISLAE
{
    IMatrix Matrix { get; }
    void AddLocalRightPart(int[] dofs, double[] lrp);
    void AddFirstBoundaryConditions(int[] dofs, double[] lrp);
    void Clear();
    void ClearRightPart();
    double[] RightPart { get; }
}

public interface ISLAESolver
{
    ISLAE SLAE { get; }
    double[] Solve();
}
