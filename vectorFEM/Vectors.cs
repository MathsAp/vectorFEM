using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.Json.Serialization;
using System.Threading.Tasks;

namespace Core
{
    public enum AngleMeasureUnits { amuRadians = 0, amuDegrees = 1 };

    //[JsonConverter(typeof(Vector2DJsonConverter))]
    public readonly struct Vector2D : IEquatable<Vector2D>
    {
        public static readonly Vector2D Zero = new Vector2D(0, 0);
        public static readonly Vector2D XAxis = new Vector2D(1, 0);
        public static readonly Vector2D YAxis = new Vector2D(0, 1);

        public double X { get; }
        public double Y { get; }
        public double[] AsArray() => new[] { X, Y };

        public Vector2D(double x, double y)
        {
            X = x;
            Y = y;
        }

        public void Deconstruct(out double x, out double y)
            => (x, y) = (X, Y);
        /// <summary>
        ///  из полярных координат в декартовы
        /// </summary>
        public Vector2D(double a, double b, AngleMeasureUnits measure) : this()
        {
            if (measure == AngleMeasureUnits.amuRadians)
            {
                X = a * Math.Cos(b);
                Y = a * Math.Sin(b);
            }
            else
            {
                double c = b * Math.PI / 180;
                X = a * Math.Cos(c);
                Y = a * Math.Sin(c);
            }
        }
        public Vector2D(ReadOnlySpan<double> arr)
        {
#if DEBUG
            if (arr.Length != 2) throw new ArgumentException("Array size error");
#endif
            X = arr[0];
            Y = arr[1];
        }
        public double this[int k]
        {
            get
            {
                return k switch
                {
                    0 => X,
                    1 => Y,
                    _ => throw new Exception("get: Vector2D out of range"),
                };
            }
        }
        public static double Distance(Vector2D a, Vector2D b) => (a - b).Norm;

        public static double SqrDistance(Vector2D a, Vector2D b)
        {
            Vector2D diff = a - b;
            return diff * diff;
        }
        public double Distance(Vector2D b) => Distance(this, b);

        public double SqrDistance(Vector2D b) => SqrDistance(this, b);

        public double Norm => Math.Sqrt(X * X + Y * Y);

        public Vector2D Normalize() => this / Norm;

        public override string ToString() => $"Vec({X}, {Y})";

        public override bool Equals(object? obj) => obj is Vector2D v && Equals(v);

        public override int GetHashCode() => HashCode.Combine(X, Y);

        public bool Equals(Vector2D a) => a.X == X && a.Y == Y;
        public static bool TryParse(string line, out Vector2D res)
        {
            double x, y;
            var words = line.Split(new[] { ' ', '\t', ',', '>', '<', '(', ')' }, StringSplitOptions.RemoveEmptyEntries);
            if (words[0] == "Vec")
            {
                if (words.Length != 3 || !double.TryParse(words[1], out x) || !double.TryParse(words[2], out y))
                {
                    res = Zero;
                    return false;
                }
                else { res = new Vector2D(x, y); return true; }
            }
            if (words.Length != 2 || !double.TryParse(words[0], out x) || !double.TryParse(words[1], out y))
            {
                res = Zero;
                return false;
            }
            else { res = new Vector2D(x, y); return true; }
        }

        public static Vector2D Parse(string line)
        {
            if (!TryParse(line, out Vector2D res))
                throw new FormatException("Can't parse Vector2D!");
            return res;
        }
        public Vector2D Round(int digits) => new Vector2D(Math.Round(X, digits), Math.Round(Y, digits));

        public static Vector2D Vec(double x, double y) => new Vector2D(x, y);
        public Vector3D As3D() => new Vector3D(X, Y, 0);
        #region Static operators

        public static Vector2D operator -(Vector2D a) => new Vector2D(-a.X, -a.Y);

        public static Vector2D operator +(Vector2D a, Vector2D b) => new Vector2D(a.X + b.X, a.Y + b.Y);

        public static Vector2D operator -(Vector2D a, Vector2D b) => new Vector2D(a.X - b.X, a.Y - b.Y);

        public static Vector2D operator /(Vector2D a, double v) => new Vector2D(a.X / v, a.Y / v);

        public static Vector2D operator *(Vector2D a, double v) => new Vector2D(a.X * v, a.Y * v);

        public static Vector2D operator *(double v, Vector2D a) => new Vector2D(v * a.X, v * a.Y);

        public static double operator *(Vector2D a, Vector2D b) => a.X * b.X + a.Y * b.Y;

        public static bool operator ==(Vector2D a, Vector2D b) => a.X == b.X && a.Y == b.Y;

        public static bool operator !=(Vector2D a, Vector2D b) => a.X != b.X || a.Y != b.Y;

        public static Vector2D Cross(Vector2D v1) => new Vector2D(v1.Y, -v1.X);

        public static double Mixed(Vector2D v1, Vector2D v2) => v1.Y * v2.X - v1.X * v2.Y;

        public static Vector2D Sum(Vector2D a, Vector2D b) => new Vector2D(a.X + b.X, a.Y + b.Y);

        #endregion
        #region EqualityComparer

        private class EqualityComparer : IEqualityComparer<Vector2D>
        {
            public int Digits { get; set; }

            public bool Equals(Vector2D v1, Vector2D v2)
            {
                return v1.Round(Digits) == v2.Round(Digits);
            }

            public int GetHashCode(Vector2D obj)
            {
                return obj.Round(Digits).GetHashCode();
            }
        }

        public static IEqualityComparer<Vector2D> CreateComparer(int digits = 7)
        {
            return new EqualityComparer { Digits = digits };
        }


        #endregion
    }

    // [JsonConverter(typeof(Vector3DJsonConverter))]
    public readonly struct Vector3D : INumberBase<Vector3D>
    {
        public static Vector3D Zero { get; } = new Vector3D(0, 0, 0);
        public static Vector3D XAxis { get; } = new Vector3D(1, 0, 0);
        public static Vector3D YAxis { get; } = new Vector3D(0, 1, 0);
        public static Vector3D ZAxis { get; } = new Vector3D(0, 0, 1);
        public static Vector3D[] Axes { get; } = { XAxis, YAxis, ZAxis };

        public double X { get; }
        public double Y { get; }
        public double Z { get; }

        public Vector3D(double x, double y, double z)
        {
            X = x;
            Y = y;
            Z = z;
        }
        public Vector3D(Vector2D vec, double z)
        {
            X = vec.X;
            Y = vec.Y;
            Z = z;
        }
        public Vector3D(ReadOnlySpan<double> arr)
        {
#if DEBUG
            if (arr.Length != 3) throw new ArgumentException();
#endif
            X = arr[0];
            Y = arr[1];
            Z = arr[2];
        }
        public double Distance(Vector3D b) => Distance(this, b);

        public double SqrDistance(Vector3D b) => SqrDistance(this, b);
        public void Deconstruct(out double x, out double y, out double z)
            => (x, y, z) = (X, Y, Z);

        public double this[int k]
        {
            get
            {
                return k switch
                {
                    0 => X,
                    1 => Y,
                    2 => Z,
                    _ => throw new Exception("get: Vector3D out of range"),
                };
            }
        }

        public double[] AsArray() => new[] { X, Y, Z };

        public Vector2D As2D() => new Vector2D(X, Y);

        public double Norm => Math.Sqrt(X * X + Y * Y + Z * Z);

        public double MaxNorm => Math.Max(Math.Abs(X), Math.Max(Math.Abs(Y), Math.Abs(Z)));

        public Vector3D Projection(Vector3D p) => (this * p) * p;

        public Vector3D Normalize() => this / Norm;

        public Vector3D Round(int digits) => new Vector3D(Math.Round(X, digits), Math.Round(Y, digits), Math.Round(Z, digits));

        public override string ToString() => $"Vec({X}, {Y}, {Z})";

        public override bool Equals(object? obj) => obj is Vector3D v && Equals(v);

        public override int GetHashCode() => HashCode.Combine(X, Y, Z);

        public bool Equals(Vector3D a) => a.X == X && a.Y == Y && a.Z == Z;

        public static bool TryParse(string line, out Vector3D res)
        {
            double x, y, z;
            var words = line.Split(new[] { ' ', '\t', ',', '(', ')', '<', '>' }, StringSplitOptions.RemoveEmptyEntries);
            if (words[0] == "Vec")
            {
                if (words.Length != 4 || !double.TryParse(words[1], out x) || !double.TryParse(words[2], out y)
                    || !double.TryParse(words[3], out z))
                {
                    res = Zero;
                    return false;
                }
                res = new Vector3D(x, y, z);
                return true;
            }
            if (words.Length != 3 || !double.TryParse(words[0], out x) || !double.TryParse(words[1], out y)
                || !double.TryParse(words[2], out z))
            {
                res = Zero;
                return false;
            }

            res = new Vector3D(x, y, z);
            return true;
        }
        public static Vector3D Vec(double x, double y, double z) => new Vector3D(x, y, z);

        public static Vector3D Parse(string line)
        {
            Vector3D res;
            if (!TryParse(line, out res))
                throw new FormatException("Can't parse Vector3D!");
            return res;
        }

        #region Static operators

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator -(Vector3D a) => new Vector3D(-a.X, -a.Y, -a.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator +(Vector3D a) => a;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double operator *(Vector3D a, Vector3D b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(double a, Vector3D b) => new Vector3D(a * b.X, a * b.Y, a * b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator *(Vector3D b, double a) => new Vector3D(a * b.X, a * b.Y, a * b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator /(Vector3D a, double v) => new Vector3D(a.X / v, a.Y / v, a.Z / v);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator +(Vector3D a, Vector3D b) => new Vector3D(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D operator -(Vector3D a, Vector3D b) => new Vector3D(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Vector3D a, Vector3D b) => a.X == b.X && a.Y == b.Y && a.Z == b.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Vector3D a, Vector3D b) => a.X != b.X || a.Y != b.Y || a.Z != b.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Cross(Vector3D v1, Vector3D v2) =>
            new Vector3D(v1.Y * v2.Z - v2.Y * v1.Z, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Mixed(Vector3D v1, Vector3D v2, Vector3D v3) =>
            (v1.Y * v2.Z - v2.Y * v1.Z) * v3.X + (v1.Z * v2.X - v1.X * v2.Z) * v3.Y + (v1.X * v2.Y - v1.Y * v2.X) * v3.Z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Sum(Vector3D a, Vector3D b) => a + b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Min(Vector3D a, Vector3D b) =>
            new Vector3D(Math.Min(a.X, b.X), Math.Min(a.Y, b.Y), Math.Min(a.Z, b.Z));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Vector3D Max(Vector3D a, Vector3D b) =>
            new Vector3D(Math.Max(a.X, b.X), Math.Max(a.Y, b.Y), Math.Max(a.Z, b.Z));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Distance(Vector3D a, Vector3D b) => (a - b).Norm;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double SqrDistance(Vector3D a, Vector3D b)
        {
            var diff = a - b;
            return diff * diff;
        }
        #endregion

        #region EqualityComparer

        private class EqualityComparer : IEqualityComparer<Vector3D>
        {
            public int Digits { get; set; }

            public bool Equals(Vector3D v1, Vector3D v2)
            {
                return v1.Round(Digits) == v2.Round(Digits);
            }

            public int GetHashCode(Vector3D obj)
            {
                return obj.Round(Digits).GetHashCode();
            }
        }

        public static IEqualityComparer<Vector3D> CreateComparer(int digits = 7)
        {
            return new EqualityComparer { Digits = digits };
        }
        #endregion

        public static Vector3D Abs(Vector3D value) => throw new NotSupportedException();
        public static bool IsCanonical(Vector3D value) => true;
        public static bool IsComplexNumber(Vector3D value) => false;
        public static bool IsEvenInteger(Vector3D value) => false;
        public static bool IsFinite(Vector3D value) => double.IsFinite(value.X) && double.IsFinite(value.Y) && double.IsFinite(value.Z);
        public static bool IsImaginaryNumber(Vector3D value) => false;
        public static bool IsInfinity(Vector3D value) => double.IsInfinity(value.X) || double.IsInfinity(value.Y) || double.IsInfinity(value.Z);
        public static bool IsInteger(Vector3D value) => false;
        public static bool IsNaN(Vector3D value) => double.IsNaN(value.X) || double.IsNaN(value.Y) || double.IsNaN(value.Z);
        public static bool IsNegative(Vector3D value) => false;
        public static bool IsNegativeInfinity(Vector3D value) => false;
        public static bool IsNormal(Vector3D value) => value.Norm == 1.0;
        public static bool IsOddInteger(Vector3D value) => false;
        public static bool IsPositive(Vector3D value) => false;
        public static bool IsPositiveInfinity(Vector3D value) => false;
        public static bool IsRealNumber(Vector3D value) => false;
        public static bool IsSubnormal(Vector3D value) => false;
        public static bool IsZero(Vector3D value) => value == default;
        public static Vector3D MaxMagnitude(Vector3D x, Vector3D y) => MaxMagnitudeNumber(x, y);
        public static Vector3D MaxMagnitudeNumber(Vector3D x, Vector3D y) => x.X > y.X || (x.X == y.X && x.Y > y.Y) || (x.X == y.X && x.Y == y.Y && x.Z > y.Z) ? x : y;
        public static Vector3D MinMagnitude(Vector3D x, Vector3D y) => MinMagnitudeNumber(x, y);
        public static Vector3D MinMagnitudeNumber(Vector3D x, Vector3D y) => x.X < y.X || (x.X == y.X && x.Y < y.Y) || (x.X == y.X && x.Y == y.Y && x.Z < y.Z) ? x : y;
        public static Vector3D Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider) => Parse(s.ToString(), style, provider);
        public static Vector3D Parse(string s, NumberStyles style, IFormatProvider? provider) => Parse(s);
        public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, [MaybeNullWhen(false)] out Vector3D result) => TryParse(s.ToString(), style, provider, out result);
        public static bool TryParse([NotNullWhen(true)] string? s, NumberStyles style, IFormatProvider? provider, [MaybeNullWhen(false)] out Vector3D result) => TryParse(s!, out result);

        public static Vector3D One => throw new NotSupportedException();
        public static int Radix => throw new NotSupportedException();

        public bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format, IFormatProvider? provider)
        {
            var s = ToString();
            charsWritten = s.Length;
            return s.AsSpan().TryCopyTo(destination);
        }
        public string ToString(string? format, IFormatProvider? formatProvider) => ToString();
        public static Vector3D Parse(ReadOnlySpan<char> s, IFormatProvider? provider) => Parse(s.ToString());
        public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, [MaybeNullWhen(false)] out Vector3D result) => TryParse(s.ToString(), out result);
        public static Vector3D Parse(string s, IFormatProvider? provider) => Parse(s);
        public static bool TryParse([NotNullWhen(true)] string? s, IFormatProvider? provider, [MaybeNullWhen(false)] out Vector3D result) => TryParse(s!, out result);

        public static Vector3D AdditiveIdentity => Zero;

        public static Vector3D operator --(Vector3D value) => throw new NotSupportedException();
        public static Vector3D operator /(Vector3D left, Vector3D right) => throw new NotSupportedException();
        public static Vector3D operator ++(Vector3D value) => throw new NotSupportedException();

        public static Vector3D MultiplicativeIdentity => throw new NotSupportedException();

        static Vector3D IMultiplyOperators<Vector3D, Vector3D, Vector3D>.operator *(Vector3D left, Vector3D right) => throw new NotSupportedException();

        static bool INumberBase<Vector3D>.TryConvertFromChecked<TOther>(TOther value, out Vector3D result)
        {
            if (value is Vector3D other)
            {
                result = other;
                return true;
            }
            result = default;
            return false;
        }


        static bool INumberBase<Vector3D>.TryConvertFromSaturating<TOther>(TOther value, out Vector3D result)
        {
            if (value is Vector3D other)
            {
                result = other;
                return true;
            }
            result = default;
            return false;
        }
        static bool INumberBase<Vector3D>.TryConvertFromTruncating<TOther>(TOther value, out Vector3D result)
        {
            if (value is Vector3D other)
            {
                result = other;
                return true;
            }
            result = Vector3D.Zero;
            return false;

        }
        static bool INumberBase<Vector3D>.TryConvertToChecked<TOther>(Vector3D value, out TOther result)
        {
            if (typeof(TOther) == typeof(Vector3D))
            {
                result = (TOther)(object)value;
                return true;
            }
            result = TOther.Zero;
            return false;
        }
        static bool INumberBase<Vector3D>.TryConvertToSaturating<TOther>(Vector3D value, out TOther result)
        {
            if (typeof(TOther) == typeof(Vector3D))
            {
                result = (TOther)(object)value;
                return true;
            }
            result = TOther.Zero;
            return false;

        }
        static bool INumberBase<Vector3D>.TryConvertToTruncating<TOther>(Vector3D value, out TOther result)
        {
            if (typeof(TOther) == typeof(Vector3D))
            {
                result = (TOther)(object)value;
                return true;
            }
            result = TOther.Zero;
            return false;

        }

    }
    public readonly struct ComplexVector3D : IEquatable<ComplexVector3D>
    {
        public static readonly ComplexVector3D Zero = new ComplexVector3D(0, 0, 0);

        public Complex X { get; }
        public Complex Y { get; }
        public Complex Z { get; }

        public Vector3D Real => new Vector3D(X.Real, Y.Real, Z.Real);
        public Vector3D Imaginary => new Vector3D(X.Imaginary, Y.Imaginary, Z.Imaginary);


        public ComplexVector3D(Complex x, Complex y, Complex z)
            => (X, Y, Z) = (x, y, z);

        public static implicit operator ComplexVector3D(in Vector3D v) => new ComplexVector3D(v.X, v.Y, v.Z);

        public override string ToString() => $"({X}, {Y}, {Z})";

        public override bool Equals(object? obj) => obj is ComplexVector3D v && Equals(v);

        public override int GetHashCode() => HashCode.Combine(X, Y, Z);

        public bool Equals(ComplexVector3D a) => a.X == X && a.Y == Y && a.Z == Z;

        #region Static operators

        public static ComplexVector3D operator -(in ComplexVector3D a) => new ComplexVector3D(-a.X, -a.Y, -a.Z);

        public static ComplexVector3D operator +(in ComplexVector3D a) => a;

        public static ComplexVector3D operator *(Complex a, in ComplexVector3D b) => new ComplexVector3D(a * b.X, a * b.Y, a * b.Z);

        public static ComplexVector3D operator *(in ComplexVector3D b, Complex a) => new ComplexVector3D(a * b.X, a * b.Y, a * b.Z);

        public static ComplexVector3D operator /(in ComplexVector3D a, Complex v) => new ComplexVector3D(a.X / v, a.Y / v, a.Z / v);

        public static ComplexVector3D operator +(in ComplexVector3D a, in ComplexVector3D b) => new ComplexVector3D(a.X + b.X, a.Y + b.Y, a.Z + b.Z);

        public static ComplexVector3D operator -(in ComplexVector3D a, in ComplexVector3D b) => new ComplexVector3D(a.X - b.X, a.Y - b.Y, a.Z - b.Z);

        public static Complex operator *(in ComplexVector3D a, in ComplexVector3D b)
            => a.X * Complex.Conjugate(b.X) + a.Y * Complex.Conjugate(b.Y) + a.Z * Complex.Conjugate(b.Z);

        public static ComplexVector3D Cross(Vector3D v1, Vector3D v2) =>
            new ComplexVector3D(v1.Y * v2.Z - v2.Y * v1.Z, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);

        public static ComplexVector3D Cross(ComplexVector3D v1, Vector3D v2) =>
            new ComplexVector3D(v1.Y * v2.Z - v2.Y * v1.Z, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);

        public static ComplexVector3D Cross(Vector3D v1, ComplexVector3D v2) =>
            new ComplexVector3D(v1.Y * v2.Z - v2.Y * v1.Z, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);

        public static ComplexVector3D Cross(ComplexVector3D v1, ComplexVector3D v2) =>
            new ComplexVector3D(v1.Y * v2.Z - v2.Y * v1.Z, v1.Z * v2.X - v1.X * v2.Z, v1.X * v2.Y - v1.Y * v2.X);


        #endregion
    }
    public static class VectorExtensions
    {
        public static Vector3D Sum(this IEnumerable<Vector3D> vectors)
        {
            try
            {
                return vectors.Aggregate(Vector3D.Sum);
            }
            catch (InvalidOperationException)
            {
                return default;
            }
        }
        public static Vector3D WeightedSum(this IEnumerable<Vector3D> vectors, IEnumerable<double> weights)
        {
            return vectors.Zip(weights, (a, b) => a * b).Aggregate(Vector3D.Sum);
        }

        public static Vector3D WeightedSum(this Vector3D[] vectors, double[] weights)
        {
            double x = 0, y = 0, z = 0;
            for (int i = 0; i < vectors.Length; i++)
            {
                x += vectors[i].X * weights[i];
                y += vectors[i].Y * weights[i];
                z += vectors[i].Z * weights[i];
            }
            return new Vector3D(x, y, z);
        }
        public static Vector3D WeightedSum(this ReadOnlySpan<Vector3D> vectors, ReadOnlySpan<double> weights)
        {
            double x = 0, y = 0, z = 0;
            for (int i = 0; i < vectors.Length; i++)
            {
                x += vectors[i].X * weights[i];
                y += vectors[i].Y * weights[i];
                z += vectors[i].Z * weights[i];
            }
            return new Vector3D(x, y, z);
        }

        public static Vector3D CenterMass(this IEnumerable<Vector3D> vectors)
        {
            return vectors.Aggregate(Vector3D.Sum) / vectors.Count();
        }

        public static Vector2D Center(this IEnumerable<Vector2D> vectors)
        {
            return (vectors.Aggregate(Vector2D.Sum) / vectors.Count());
        }

        public static Vector3D CenterBox(this IEnumerable<Vector3D> vectors)
        {
            var min = vectors.Aggregate(vectors.First(), Vector3D.Min);
            var max = vectors.Aggregate(vectors.First(), Vector3D.Max);
            return (min + max) / 2;
        }

        public static Vector3D Center(this IEnumerable<Vector3D> vectors)
        {
            return vectors.CenterMass();
        }

        public static double Radius(this IEnumerable<Vector3D> vectors)
        {
            var c = vectors.Center();
            return vectors.Max(v => (v - c).Norm);
        }

        public static Func<Vector3D, int> SplitBox(this IEnumerable<Vector3D> vectors)
        {
            var min = vectors.Aggregate(Vector3D.Min);
            var max = vectors.Aggregate(Vector3D.Max);
            var center = (min + max) / 2;
            double diameter = Vector3D.Axes.Max(a => (max - min) * a);
            var v = Vector3D.Axes.Where(a => (max - min) * a > diameter / 2).ToArray();
            return p => v.Select((x, i) => (p - center) * x >= 0 ? 1 << i : 0).Sum();
        }

        public static Func<Vector3D, int> SplitSphere(this IEnumerable<Vector3D> vectors)
        {
            var min = vectors.Aggregate(vectors.First(), Vector3D.Min);
            var max = vectors.Aggregate(vectors.First(), Vector3D.Max);
            var center = (min + max) / 2;

            var e = new Vector3D[3];
            e[0] = vectors.Aggregate(vectors.First() - center,
                                     (acc, curr) => (curr - center).Norm > acc.Norm ? curr - center : acc);

            if (e[0].Norm > 1e-10)
            {
                Vector3D pr1(Vector3D x) => x - x * e[0] / (e[0] * e[0]) * e[0];
                e[1] = vectors.Aggregate(pr1(vectors.First() - center),
                                           (acc, curr) => pr1(curr - center).Norm > acc.Norm ? pr1(curr - center) : acc);
                if (e[1].Norm > 1e-10)
                {
                    Vector3D pr2(Vector3D x) => x - x * e[0] / (e[0] * e[0]) * e[0] - x * e[1] / (e[1] * e[1]) * e[1];
                    e[2] = vectors.Aggregate(pr2(vectors.First() - center),
                                               (acc, curr) => pr2(curr - center).Norm > acc.Norm ? pr2(curr - center) : acc);
                }
            }

            var v = e.Where(a => 2 * a.Norm > e[0].Norm).ToArray();
            return p => v.Select((x, i) => (p - center) * x >= 0 ? 1 << i : 0).Sum();
        }

        public static double[,] JacobiEigenValues(double[,] a)
        {
            int n = a.GetLength(0);
            var u = new double[3, 3];
            for (int i = 0; i < n; i++)
                u[i, i] = 1;

            do
            {
                int maxi = 1, maxj = 0;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < i; j++)
                        if (Math.Abs(a[i, j]) > Math.Abs(a[maxi, maxj]))
                        {
                            maxi = i;
                            maxj = j;
                        }
                if (Math.Abs(a[maxi, maxj]) < 1e-3) break;

                double phi = 0.5 * Math.Atan2(2 * a[maxi, maxj], a[maxi, maxi] - a[maxj, maxj]);
                double cos = Math.Cos(phi);
                double sin = Math.Sin(phi);

                for (int i = 0; i < n; i++)
                {
                    double x = a[i, maxi];
                    double y = a[i, maxj];
                    a[i, maxi] = cos * x - sin * y;
                    a[i, maxj] = sin * x + cos * y;
                }
                for (int i = 0; i < n; i++)
                {
                    double x = a[maxi, i];
                    double y = a[maxj, i];
                    a[maxi, i] = cos * x - sin * y;
                    a[maxj, i] = sin * x + cos * y;
                }

                for (int i = 0; i < n; i++)
                {
                    double x = u[maxi, i];
                    double y = u[maxj, i];
                    u[maxi, i] = cos * x - sin * y;
                    u[maxj, i] = sin * x + cos * y;
                }
            } while (true);

            return u;
        }

        public static Func<Vector3D, int> Principial(this IEnumerable<Vector3D> vectors)
        {
            var mean = vectors.Center();
            int n = vectors.Count();
            var C = new double[3, 3];
            double a = 1.0 / (n - 1);
            foreach (var v in vectors.Select(v => (v - mean)))
            {
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        C[i, j] += a * (v * Vector3D.Axes[i]) * (v * Vector3D.Axes[j]);
            }
            double trace = 0;
            for (int i = 0; i < 3; i++) trace += C[i, i];
            var u = JacobiEigenValues(C);

            var eigenvalues = Enumerable.Range(0, 3).Select(
                i => new { lambda = C[i, i], x = new Vector3D(u[i, 0], u[i, 1], u[i, 2]) }).
                OrderByDescending(e => e.lambda).TakeWhile(next => next.lambda > trace / 3).Select(r => r.x).ToArray();

            return p => eigenvalues.Select((x, i) => (p - mean) * x >= 0 ? 1 << i : 0).Sum();


            //Random r = new Random(12345);
            //Vector3D x = new Vector3D(r.NextDouble(), r.NextDouble(), r.NextDouble()).Normalize();
            //Vector3D dx;
            //do
            //{
            //    dx.X = C[0, 0] * x.X + C[0, 1] * x.Y + C[0, 2] * x.Z;
            //    dx.Y = C[1, 0] * x.X + C[1, 1] * x.Y + C[1, 2] * x.Z;
            //    dx.Z = C[2, 0] * x.X + C[2, 1] * x.Y + C[2, 2] * x.Z;
            //    dx /= dx.Norm;
            //    dx -= x;
            //    x += dx;
            //} while (dx.Norm > 1e-2);
            //return x;
        }
    }
}
