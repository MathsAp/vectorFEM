using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/* Локальные матрицы одномерных элементов без домножения на шаг и коэффициент */

namespace Core;

public static class LinearLocalMatrices
{
    public static readonly double[,] G = { {  1,  -1 },

                                           { -1,   1 } };

    public static readonly double[,] M = { { 1 / 3d,  1 / 6d },

                                           { 1 / 6d,  1 / 3d } };
}

public static class QuadraticLocalMatrices
{
    public static readonly double[,] G = { {  7 / 3d,  -8 / 3d,   1 / 3d },

                                           { -8 / 3d,  16 / 3d,  -8 / 3d },

                                           {  1 / 3d,  -8 / 3d,   7 / 3d } };


    public static readonly double[,] M = { {  2 / 15d,  1 / 15d,  -1 / 30d },

                                           {  1 / 15d,  8 / 15d,   1 / 15d },

                                           { -1 / 30d,  1 / 15d,   2 / 15d } };

}

public static class CubicLocalMatrices
{
    public static readonly double[,] G = { {   37 / 10d,  -189 / 40d,    27 / 20d,   -13 / 40d },

                                           { -189 / 40d,     54 / 5d,  -297 / 40d,    27 / 20d },

                                           {   27 / 20d,  -297 / 40d,     54 / 5d,  -189 / 40d },

                                           {  -13 / 40d,    27 / 20d,  -189 / 40d,    37 / 10d } };


    public static readonly double[,] M = { {   8 / 105d,    33 / 560d,   -3 / 140d,  19 / 1680d },

                                           {  33 / 560d,     27 / 70d,  -27 / 560d,   -3 / 140d },

                                           {  -3 / 140d,   -27 / 560d,    27 / 70d,   33 / 560d },

                                           { 19 / 1680d,    -3 / 140d,   33 / 560d,    8 / 105d } };
}

public static class RectangleRZHierarhicalQuadraticLocalMatrices
{
    public static readonly double[,] G = { { 1d, -1d, 0d },

                                           { -1d, 1d, 0d },

                                           { 0d, 0d, 1d / 3d } };

    public static readonly double[,] M = { { 1d / 3d, 1d / 6d, 1d / 12d },

                                           { 1d / 6d, 1d / 3d, 1d / 12d },

                                           { 1d / 12d, 1d / 12d, 1d / 30d } };

    public static readonly double[,] Grz = { { 1 / 2d, -1 / 2d, 1 / 6d },

                                             { -1 / 2d, 1 / 2d, -1 / 6d },

                                             { 1 / 6d, -1 / 6d, 1 / 6d } };

    public static readonly double[,] Mrz = { { 1 / 12d, 1 / 12d, 1 / 30d },

                                             { 1 / 12d, 1 / 4d, 1 / 20d },

                                             { 1 / 30d, 1 / 20d, 1 / 60d } };
}

public static class RectangleHierarhicalQuadraticLocalMatrices
{
    public static readonly double[,] G = { { 1d, -1d, 0d },

                                           { -1d, 1d, 0d },
        
                                           { 0d, 0d, 1d / 3d } };

    public static readonly double[,] M = { { 1d / 3d, 1d / 6d, 1d / 12d },

                                           { 1d / 6d, 1d / 3d, 1d / 12d },

                                           { 1d / 12d, 1d / 12d, 1d / 30d } };
}
