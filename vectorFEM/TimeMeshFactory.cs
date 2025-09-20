using Core;
using FEM;
using Microsoft.CodeAnalysis.CSharp.Scripting;
using Microsoft.CodeAnalysis.CSharp.Syntax;
using Microsoft.CodeAnalysis.Scripting;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core;

public static class TimeMeshFactory<T>
{
    static ScriptOptions options = ScriptOptions.Default.AddReferences(typeof(Program).Assembly);

    public static (ITimeMesh TimeMesh, Func<Vector3D, T> InitialFunction) CreateTimeMesh(string path)
    {
        if (!Path.Exists(path)) throw new ArgumentException("The path is incorrect", nameof(path));

        using var reader = new StreamReader(Path.Combine(path, "time.txt"));

        int n = int.Parse(reader.ReadLine());
        double[] t = [.. reader.ReadLine().Split(' ', StringSplitOptions.RemoveEmptyEntries).Select(double.Parse)];

        if (t.Length != n) 
            throw new FormatException("The file contains incorrect data (the number of time values does not match the set one)");

        string[] values = reader.ReadLine().Split(' ', StringSplitOptions.RemoveEmptyEntries);

        if (values.Length != 2 * (n - 1))
            throw new FormatException("The file contains incorrect data (the number of pairs to split the intervals should be one less than the number of time values)");

        RefineParams refineParams = new();

        for (int i = 0; i < n - 1; ++i)
        {
            refineParams.splitCount.Add(int.Parse(values[2 * i]));
            refineParams.stretchRatio.Add(double.Parse(values[2 * i + 1]));
        }

        string funcBody = reader.ReadLine().Trim();
        Func<Vector3D, T> initialFunc = Compile<Func<Vector3D, T>>(funcBody, "p => ");

        TimeMesh timeMesh = new(t, refineParams);

        return (timeMesh, initialFunc);
    }

    static TFunc Compile<TFunc>(string body, string prefix) =>
        CSharpScript.EvaluateAsync<TFunc>(prefix + body, options).GetAwaiter().GetResult();
}
