using FEM;
using Quasar.Native;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public class PardisoNonSymmMatrix : IMatrix, IPardisoMatrix<double>
    {
        public PardisoNonSymmMatrix(SortedSet<int>[] profile, PardisoMatrixType _type)
        {
            SetProfile(profile);

            type = _type;
        }

        public int N { get => ia.Length - 1; }
        int IPardisoMatrix<double>.n => ia.Length - 1;

        public int[] ia = Array.Empty<int>();
        ReadOnlySpan<int> IPardisoMatrix<double>.ia => ia;

        public int[] ja = Array.Empty<int>();
        ReadOnlySpan<int> IPardisoMatrix<double>.ja => ja;

        public double[] values = Array.Empty<double>();
        ReadOnlySpan<double> IPardisoMatrix<double>.a => values;

        PardisoMatrixType type;
        PardisoMatrixType IPardisoMatrix<double>.MatrixType => type;

        public void AddLocal(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1d)
        {
            for (int i = 0; i < dofsi.Length; ++i)
            {
                for (int j = 0; j < dofsj.Length; ++j)
                {
                    int di = dofsi[i];
                    int dj = dofsj[j];


                    int ia0 = ia[di];
                    int ia1 = ia[di + 1];
                    int ind = BinarySearch(ja, dj, ia0, ia1 - 1);

                    if (ind != -1)
                        LinqExtensions.ThreadSafeAdd(values, ind, coeff * matrix[i, j]);
                  //  else
                      //  throw new Exception();
                }
            }
        }

        public void AddLocalTransposed(int[] dofsi, int[] dofsj, double[,] matrix, double coeff = 1d)
        {
            for (int i = 0; i < dofsi.Length; ++i)
            {
                for (int j = 0; j < dofsj.Length; ++j)
                {
                    int di = dofsi[i];
                    int dj = dofsj[j];


                    int ia0 = ia[di];
                    int ia1 = ia[di + 1];
                    int ind = BinarySearch(ja, dj, ia0, ia1 - 1);

                    if (ind != -1)
                        LinqExtensions.ThreadSafeAdd(values, ind, coeff * matrix[j, i]);
                    //else
                        //throw new Exception();
                }
            }
        }

        public void Symmetrize(int dof, double value, double[] RightPart)
        {
            int ia0 = ia[dof];
            int ia1 = ia[dof + 1];

            int diag = BinarySearch(ja, dof, ia0, ia1 - 1);

            for (int ind = ia0; ind < ia1; ++ind)
            {
                LinqExtensions.ThreadSafeSet(values, ind, 0);
            }

            LinqExtensions.ThreadSafeSet(values, diag, 1);
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
                ia[i + 1] = ia[i] + profile[i].Count();

            ja = new int[ia.Last()];

            for (int i = 0; i < profile.Length; ++i)
                profile[i].ToArray().AsSpan().CopyTo(ja.AsSpan(ia[i]));
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
    }
}
