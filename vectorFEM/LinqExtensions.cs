using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core;

public static class LinqExtensions
{
    public static void ThreadSafeSet(this double[] array, int index, double value)
    {
        if (double.IsNaN(value))
            throw new ArgumentNullException(nameof(value));

        double initialValue, computedValue;
        do
        {
            initialValue = array[index];
            computedValue = value;
        }
        while (initialValue != Interlocked.CompareExchange(ref array[index], computedValue, initialValue));
    }

    public static void ThreadSafeAdd(this double[] array, int index, double value)
    {
        if (double.IsNaN(value))
            throw new ArgumentNullException(nameof(value));

        double initialValue, computedValue;
        do
        {
            initialValue = array[index];
            computedValue = initialValue + value;
        }
        while (initialValue != Interlocked.CompareExchange(ref array[index], computedValue, initialValue));
    }

}
