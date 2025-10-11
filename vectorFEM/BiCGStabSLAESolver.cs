using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TelmaQuasarCommon;

namespace Core;

public class BiCGStabSLAESolver : ISLAESolver
{
    public BiCGStabSLAESolver(ISLAE slae, double eps, int maxIterations)
    {
        SLAE = slae;
        if (slae.Matrix is SparseNonSymmMatrix matrix) A = matrix;
        else throw new ArgumentException($"Incorrect matrix type, {nameof(SparseNonSymmMatrix)} was expected");

        this.eps = eps;
        this.maxIterations = maxIterations;

        r = new double[A.N];
        r0 = new double[A.N];
        z = new double[A.N];
        p = new double[A.N];
        y = new double[A.N];
        az = new double[A.N];
        ay = new double[A.N];
        s = new double[A.N];
    }


    public ISLAE SLAE { get; }

    SparseNonSymmMatrix? M;
    SparseNonSymmMatrix A;

    double eps;
    int maxIterations;

    double[] r = [];
    double[] r0 = [];
    double[] z = [];
    double[] p = [];
    double[] y = [];
    double[] az = [];
    double[] ay = [];
    double[] s = [];

    public void Prepare()
    {
        M = A.LUDecomposeMatrix();
    }

    public double[] Solve()
    {
        ClearVectors();
        double[] b = SLAE.RightPart;
        double[] solution = new double[A.N];
        double[] x = solution;

        A.MultiplyByVector(x, az);
        LinearAlgebraAlgorithms.SubtractionVectorsWithCoeffs(b, az, r0);
        double normB = LinearAlgebraAlgorithms.EuclideanNorm(b);
        double relativeResidual = LinearAlgebraAlgorithms.EuclideanNorm(r0) / normB;
        r0.AsSpan().CopyTo(r);
        r0.AsSpan().CopyTo(p);
        double scalarRR0 = LinearAlgebraAlgorithms.ScalarProduct(r, r0);

        int k = 1;
        for ( ; relativeResidual > eps && k < maxIterations; ++k)
        {
            CalculateSLAE(y, p);

            A.MultiplyByVector(y, ay);
            double alpha = scalarRR0 / LinearAlgebraAlgorithms.ScalarProduct(ay, r0);

            LinearAlgebraAlgorithms.AdditionVectorsWithCoeffs(r, ay, s, 1, -alpha);

            CalculateSLAE(z, s);

            A.MultiplyByVector(z, az);
            double omega = LinearAlgebraAlgorithms.ScalarProduct(az, s) / LinearAlgebraAlgorithms.ScalarProduct(az, az);

            LinearAlgebraAlgorithms.AdditionVectorsWithCoeffs(x, y, z, x, 1, alpha, omega);

            LinearAlgebraAlgorithms.AdditionVectorsWithCoeffs(s, az, r, 1, -omega);

            double scalarRR0Prev = scalarRR0;
            scalarRR0 = LinearAlgebraAlgorithms.ScalarProduct(r, r0);
            double beta = (alpha * scalarRR0) / (omega * scalarRR0Prev);

            LinearAlgebraAlgorithms.AdditionVectorsWithCoeffs(r, p, ay, p, 1, beta, -beta * omega);

            relativeResidual = LinearAlgebraAlgorithms.EuclideanNorm(r) / normB;

            Console.WriteLine($"Iteration #{k}, residual: {relativeResidual}");
        }

        A.MultiplyByVector(x, az);
        LinearAlgebraAlgorithms.SubtractionVectorsWithCoeffs(b, az, r);
        relativeResidual = LinearAlgebraAlgorithms.EuclideanNorm(r) / normB;
        Console.WriteLine($"Real residual: {relativeResidual}, Iteration: {k}");

        return solution;
    }

    void CalculateSLAE(double[] x, ReadOnlySpan<double> y)
    {
        M!.CalculateY(x, y);
        M.CalculateX(x, x);
    }

    void ClearVectors()
    {
        Array.Clear(r, 0, A.N);
        Array.Clear(r0, 0, A.N);
        Array.Clear(z, 0, A.N);
        Array.Clear(p, 0, A.N);
        Array.Clear(y, 0, A.N);
        Array.Clear(az, 0, A.N);
        Array.Clear(ay, 0, A.N);
        Array.Clear(s, 0, A.N);
    }
}
