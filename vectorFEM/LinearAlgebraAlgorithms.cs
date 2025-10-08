using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FEM;
using Quasar.Native;

namespace Core;

public static class LinearAlgebraAlgorithms
{
    public static double[] MultiplyMatrixByVector(double[,] matrix, double[] vector, double coeff = 1d)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        if (cols != vector.Length)
            throw new ArgumentException("Количество столбцов матрицы должно совпадать с длиной вектора.");

        double[] result = new double[rows];

        for (int i = 0; i < rows; i++)
        {
            double sum = 0;
            for (int j = 0; j < cols; j++)
            {
                sum += matrix[i, j] * vector[j];
            }
            result[i] = coeff * sum;
        }


        return result;
    }

    public static double[] MultiplyMatrixByVector(double[,] matrix, ReadOnlySpan<double> vector, double coeff = 1d)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        if (cols != vector.Length)
            throw new ArgumentException("Количество столбцов матрицы должно совпадать с длиной вектора.");

        double[] result = new double[rows];

        for (int i = 0; i < rows; i++)
        {
            double sum = 0;
            for (int j = 0; j < cols; j++)
            {
                sum += matrix[i, j] * vector[j];
            }
            result[i] = coeff * sum;
        }


        return result;
    }

    public static Vector3D MultiplyMatrix3By3ByVector(double[,] matrix, ReadOnlySpan<double> vector, double coeff = 1d)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        if (cols != vector.Length)
            throw new ArgumentException("Количество столбцов матрицы должно совпадать с длиной вектора.");

        if (rows != 3 || cols != 3)
            throw new ArgumentException("Матрица имеет размер, отличный от 3 на 3", nameof(matrix));

        double[] result = new double[rows];

        for (int i = 0; i < rows; ++i)
        {
            double sum = 0;
            for (int j = 0; j < cols; ++j)
            {
                sum += matrix[i, j] * vector[j];
            }

            result[i] = coeff * sum;
        }

        return new Vector3D(result);
    }

    public static Vector3D MultiplyMatrix3By3ByVector(double[,] matrix, Vector3D vector, double coeff = 1d)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        if (rows != 3 || cols != 3)
            throw new ArgumentException("Матрица имеет размер, отличный от 3 на 3", nameof(matrix));

        double[] result = new double[rows];

        for (int i = 0; i < rows; ++i)
        {
            double sum = 0;
            for (int j = 0; j < cols; ++j)
            {
                sum += matrix[i, j] * vector[j];
            }

            result[i] = coeff * sum;
        }

        return new Vector3D(result);
    }

    public static Vector2D MultiplyMatrix2By2ByVector(double[,] matrix, Vector2D vector, double coeff = 1d, bool transpose = false)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        if (rows != 2 || cols != 2)
            throw new ArgumentException("Матрица имеет размер, отличный от 2 на 2", nameof(matrix));

        double[] result = new double[rows];

        for (int i = 0; i < rows; ++i)
        {
            double sum = 0;
            for (int j = 0; j < cols; ++j)
            {
                sum += (transpose ? matrix[j, i] : matrix[i, j]) * vector[j];
            }

            result[i] = coeff * sum;
        }

        return new Vector2D(result);
    }

    public static double[,] GetMatrixMultipliedByTransposedMatrix(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        double[,] result = new double[rows, rows];

        for (int i = 0; i < rows; ++i)
        {
            for (int j = i; j < rows; ++j)
            {
                double sum = 0;
                for (int k = 0; k < cols; ++k)
                {
                    sum += matrix[i, k] * matrix[j, k];
                }

                result[i, j] = sum;
                result[j, i] = sum;
            }
        }

        return result;
    }

    public static Vector3D IntersectTwoSegments(Vector3D p1, Vector3D p2, Vector3D q1, Vector3D q2)
    {
        Vector3D r = p2 - p1;
        Vector3D s = q2 - q1;
        Vector3D w = q1 - p1;

        Vector3D rs = Vector3D.Cross(r, s);
        double t = Vector3D.Cross(w, s) * rs / (rs * rs);

        return p1 + t * r;
    }

    public static void PrintMatrix(double[,] matrix)
    {
        int rows = matrix.GetLength(0);
        int cols = matrix.GetLength(1);

        Console.WriteLine();
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                Console.Write(matrix[i, j] + " ");
            }
            Console.WriteLine();
        }
    }

    public static void PrintVector(ReadOnlySpan<double> vector)
    {
        Console.WriteLine();
        for (int i = 0; i < vector.Length; ++i)
        {
            Console.WriteLine(vector[i]);
        }
        Console.WriteLine();
    }

    public static double[,] SparseMatrixToDense(PardisoNonSymmMatrix matrix)
    {
        int N = matrix.N;

        double[,] newMatrix = new double[N, N];

        for (int i = 0; i < N; ++i)
        {
            int ia0 = matrix.ia[i];
            int ia1 = matrix.ia[i + 1];

            for (int k = ia0; k < ia1; ++k)
            {
                int j = matrix.ja[k];

                newMatrix[i, j] = matrix.values[k];
            }
        }

        return newMatrix;
    }
}
