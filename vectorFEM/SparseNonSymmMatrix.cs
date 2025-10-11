using Core;
using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization.Formatters;
using System.Text;
using System.Threading.Tasks;
using TelmaQuasarCommon;

namespace Core;

public class SparseNonSymmMatrix : IMatrix
{
    public SparseNonSymmMatrix() { }
    public SparseNonSymmMatrix(SortedSet<int>[] profile)
    {
        SetProfile(profile);
    }

    public int N => di.Length;


    int[] ia = [];
    int[] ja = [];
    double[] di = [];
    double[] au = [];
    double[] al = [];

    public void AddLocal(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1)
    {
        for (int i = 0; i < dofsi.Length; ++i)
        {
            for (int j = 0; j < dofsj.Length; ++j)
            {
                int dofi = dofsi[i];
                int dofj = dofsj[j];

                AddElemInMatrix(dofi, dofj, coeff * matrix[i, j]);
            }
        }
    }

    void AddElemInMatrix(int i, int j, double value, bool isSet = false)
    {
        if (i == j)
        {
            if (!isSet)
                LinqExtensions.ThreadSafeAdd(di, i, value);
            else
                LinqExtensions.ThreadSafeSet(di, i, value);
        }
        else if (j < i)
        {
            int ia0 = ia[i];
            int ia1 = ia[i + 1];

            int ind = BinarySearch(ja, j, ia0, ia1 - 1);

            if (ind != -1)
            {
                if (!isSet)
                    LinqExtensions.ThreadSafeAdd(al, ind, value);
                else
                    LinqExtensions.ThreadSafeSet(al, ind, value);
            }

        }
        else
        {
            int ia0 = ia[j];
            int ia1 = ia[j + 1];

            int ind = BinarySearch(ja, i, ia0, ia1 - 1);

            if (ind != -1)
            {
                if (!isSet)
                    LinqExtensions.ThreadSafeAdd(au, ind, value);
                else
                    LinqExtensions.ThreadSafeSet(au, ind, value);
            }
        }
    }

    void SetElemInMatrix(int i, int j, double value) => AddElemInMatrix(i, j, value, true);

    public void AddLocalTransposed(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1)
    {
        throw new NotImplementedException();
    }

    public void Clear()
    {
        Array.Clear(di, 0, di.Length);
        Array.Clear(au, 0, au.Length);
        Array.Clear(al, 0, al.Length);
    }

    public void SetProfile(SortedSet<int>[] profile)
    {
        BuildLowerTriangleRowSparseMatrixPortrait(profile);

        di = new double[profile.Length];
        au = new double[ia.Last()];
        al = new double[ia.Last()];
    }

    void BuildLowerTriangleRowSparseMatrixPortrait(SortedSet<int>[] profile)
    {
        ia = new int[profile.Length + 1];

        for (int i = 0; i < profile.Length; ++i)
            ia[i + 1] = ia[i] + profile[i].Where(j => j < i).Count();

        ja = new int[ia.Last()];

        for (int i = 0; i < profile.Length; ++i)
            profile[i].Where(j => j < i).ToArray().AsSpan().CopyTo(ja.AsSpan(ia[i]));
    }

    public void Symmetrize(int dof, double value, double[] RightPart)
    {
        int ia0 = ia[dof];
        int ia1 = ia[dof + 1];

        for (int j = ia0; j < ia1; ++j)
            LinqExtensions.ThreadSafeSet(al, j, 0);

        for (int j = dof + 1; j < N; ++j)
            SetElemInMatrix(dof, j, 0);

        LinqExtensions.ThreadSafeSet(di, dof, 1);
    }

    static int BinarySearch(int[] arr, int target, int low, int high)
    {
        while (low <= high)
        {
            int mid = (low + high) / 2;
            int midValue = arr[mid];

            if (midValue == target) return mid;
            else if (midValue < target) low = mid + 1;
            else high = mid - 1;
        }

        return -1;
    }

    public void MultiplyByVector(ReadOnlySpan<double> x, double[] result)
    {
        for (int i = 0; i < N; ++i)
        {
            int ia0 = ia[i];
            int ia1 = ia[i + 1];

            int k = ia0;
            double s = 0;
            result[i] = di[i] * x[i];
            for ( ; k < ia1; ++k)
            {
                int j = ja[k];

                s += al[k] * x[j];
                result[j] += au[k] * x[i]; 
            }
            result[i] += s;
        }
    }

    public SparseNonSymmMatrix LUDecomposeMatrix()
    {
        SparseNonSymmMatrix M = new()
        {
            ia = ia,
            ja = ja,
            di = new double[N],
            au = new double[au.Length],
            al = new double[al.Length]
        };

        for (int i = 0; i < N; ++i)
        {
            int i0 = ia[i];
            int i1 = ia[i + 1];

            int k = i0;
            double sumDi = 0;
            for ( ; k < i1; k++)
            {
                int ki = i0;
                int j = ja[k];
                int j0 = ia[j];
                int j1 = ia[j + 1];
                int kj = j0;
                double sumAl = 0;
                double sumAu = 0;
                while (ki < k && kj < j1)
                {
                    int ji = ja[ki];
                    int jj = ja[kj];
                    if (ji > jj) kj++;
                    else if (ji < jj) ki++;
                    else
                    {
                        sumAl += M.al[ki] * M.au[kj];
                        sumAu += M.al[kj] * M.au[ki];
                        ki++;
                        kj++;
                    }
                }

                M.al[k] = (al[k] - sumAl);
                M.au[k] = (au[k] - sumAu) / M.di[j];
                sumDi += M.al[k] * M.au[k];
            }

            M.di[i] = di[i] - sumDi;
        }

        return M;
    }

    public double[] CalculateY(double[] y, ReadOnlySpan<double> b)
    {
        for (int i = 0; i < N; ++i)
        {
            int i0 = ia[i];
            int i1 = ia[i + 1];

            int k = i0;
            double s = 0;
            for ( ; k < i1; ++k)
            {
                int j = ja[k];
                s += al[k] * y[j];
            }
            y[i] = (b[i] - s) / di[i];
        }

        return y;
    }

    public double[] CalculateX(double[] x, ReadOnlySpan<double> y)
    {
        y.CopyTo(x);

        for (int i = N - 1; i >= 0; --i)
        {
            int i0 = ia[i];
            int i1 = ia[i + 1];

            int k = i1 - 1;
            double xi = x[i];
            for ( ; k >= i0; --k)
            {
                int j = ja[k];
                x[j] -= au[k] * xi;
            }
            x[i] = xi;
        }

        return x;
    }
}
