using System;
using System.Numerics;

namespace Kursovaya 
{
    internal class Program
    {
        class SparseMatrix
        {
            public int N;
            public double[] di;
            public int[] ig;
            public int[] jg;
            public double[] ggl;
            public double[] ggu;
            // Конструктор матрицы
            public SparseMatrix(double[] di, int[] ig, int[] jg, double[] ggl, double[] ggu)
            {
                this.N = di.Length;
                this.di = di;
                this.ig = ig;
                this.jg = jg;
                this.ggl = ggl;
                this.ggu = ggu;
            }
            // Умножение матрицы на вектор x
            public double[] VectorMultiply(double[] x)
            {
                double[] y = new double[N];

                for (int i = 0; i < N; i++)
                {
                    y[i] = x[i] * di[i];

                    for (int j = ig[i]; j < ig[i + 1]; j++)
                    {
                        y[i] += ggl[j - 1] * x[jg[j - 1] - 1];
                        y[jg[j - 1] - 1] += ggu[j - 1] * x[i];
                    }
                }

                return y;
            }
        }
        static void Main(string[] args)
        {
            double[] di = new double[] { 1, 5, 8, 12, 15, 18 };
            int[] ig = new int[] { 1, 1, 2, 2, 4, 5, 7 };
            int[] jg = new int[] { 1, 1, 2, 4, 2, 3 };
            double[] ggl = new double[] { 4, 10, 11, 14, 16, 17 };
            double[] ggu = new double[] { 2, 3, 6, 13, 7, 9 };

            SparseMatrix test_matrix = new SparseMatrix(di, ig, jg, ggl, ggu);

            double[] x = new double[] { 1, 2, 3, 4, 5, 6 };

            double[] y = test_matrix.VectorMultiply(x);

            for (int i = 0; i < di.Length; i++)
            {
                Console.WriteLine(y[i]);
            }
        }
    }
}