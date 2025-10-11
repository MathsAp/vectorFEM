using System;
using System.Collections.Generic;
using System.Linq;
using System.Linq.Expressions;
using System.Text;
using System.Text.Json.Serialization;
using System.Threading.Tasks;
using FEM;
using Quasar.Native;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace Core;

public class PardisoMatrix : IMatrix, IPardisoMatrix<double>
{
    public PardisoMatrix(SortedSet<int>[] profile, PardisoMatrixType _type)
    {
        SetProfile(profile);

        type = _type;
    }

    public int N => ia.Length - 1;
    int IPardisoMatrix<double>.n => ia.Length - 1;

    int[] ia = [];
    ReadOnlySpan<int> IPardisoMatrix<double>.ia => ia;

    int[] ja = [];
    ReadOnlySpan<int> IPardisoMatrix<double>.ja => ja;

    double[] values = [];
    ReadOnlySpan<double> IPardisoMatrix<double>.a => values;

    PardisoMatrixType type;
    PardisoMatrixType IPardisoMatrix<double>.MatrixType => type;

    public void AddLocal(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1d)
    {
        for (int i = 0; i < dofsi.Length; ++i)
        {
            for (int j = i; j < dofsi.Length; ++j)
            {
                int di = dofsi[i];
                int dj = dofsi[j];

                if (di <= dj)
                {
                    int ia0 = ia[di];
                    int ia1 = ia[di + 1];
                    int ind = BinarySearch(ja, dj, ia0, ia1 - 1);

                    if (ind != -1)
                        LinqExtensions.ThreadSafeAdd(values, ind, coeff * matrix[i, j]);
                }
                else
                {
                    int ia0 = ia[dj];
                    int ia1 = ia[dj + 1];
                    int ind = BinarySearch(ja, di, ia0, ia1 - 1);

                    if (ind != -1)
                        LinqExtensions.ThreadSafeAdd(values, ind, coeff * matrix[i, j]);
                }
            }
        }
    }

    public void AddLocalTransposed(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1d)
    {
        throw new NotImplementedException();
    }

    public void Symmetrize(int dof, double value, double[] RightPart)
    {
        int ia0 = 0;
        int ia1 = 0;

        for (int i = 0; i < dof; ++i)
        {
            ia0 = ia[i];
            ia1 = ia[i + 1];

            int ind = BinarySearch(ja, dof, ia0, ia1 - 1);

            if (ind != -1)
            {
                LinqExtensions.ThreadSafeAdd(RightPart, i, -values[ind] * value);
                LinqExtensions.ThreadSafeSet(values, ind, 0);
            }
        }

        ia0 = ia[dof];
        ia1 = ia[dof + 1];

        LinqExtensions.ThreadSafeSet(values, ia0, 1);

        for (int ind = ia0 + 1; ind < ia1; ++ind)
        {
            int j = ja[ind];
            LinqExtensions.ThreadSafeAdd(RightPart, j, -values[ind] * value);
            LinqExtensions.ThreadSafeSet(values, ind, 0);
        }
    }
    public void Clear()
    {
        Array.Clear(values, 0, values.Length);
    }
    public void SetProfile(SortedSet<int>[] profile)
    {
        BuildUpperTriangleRowSparseMatrixPortrait(profile);

        values = new double[ia.Last()];
    }

    private void BuildUpperTriangleRowSparseMatrixPortrait(SortedSet<int>[] profile)
    {
        ia = new int[profile.Length + 1];

        for (int i = 0; i < profile.Length; ++i)
            ia[i + 1] = ia[i] + profile[i].Where(j => j >= i).Count();

        ja = new int[ia.Last()];

        for (int i = 0; i < profile.Length; ++i)
            profile[i].Where(j => j >= i).ToArray().AsSpan().CopyTo(ja.AsSpan(ia[i]));
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
        throw new NotImplementedException();
    }
}
