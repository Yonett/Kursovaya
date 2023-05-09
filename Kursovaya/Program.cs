using System;
using System.Numerics;

namespace Kursovaya 
{
    internal class Program
    {
        // Скалярное произведение векторов (x, y)
        static double ScalarMultiply(double[] x, double[] y)
        {
            double result = 0;
            for (int i = 0; i < x.Length; i++)
            {
                result += x[i] * y[i];
            }
            return result;
        }
        // Сумма векторов (a + b)
        static double[] VectorSum(double[] a, double[] b)
        {
            double[] result = new double[a.Length];
            for (int i = 0; i <= a.Length; i++)
            {
                result[i] = a[i] + b[i];
            }
            return result;
        }
        // Разность векторов (a - b)
        static double[] VectorDiff(double[] a, double[] b)
        {
            double[] result = new double[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                result[i] = a[i] - b[i];
            }
            return result;
        }
        // Разреженая матрица в строчно-столбцовом формате
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
        // Класс-решатель СЛАУ методом ЛОС
        class SLAESolver
        {
            public double[] LOS(SparseMatrix A, double[] b)
            {
                int maxIter = 1000;
                double[] x = new double[A.N];
                for (int i = 0; i < A.N; i++)
                {
                    x[i] = 1;
                }

                double[] r, z, p, temp;

                double alpha, beta, normSqr, nev = 0;
                double eps = 1e-7;

                r = VectorDiff(b, A.VectorMultiply(x));
                z = r;
                p = A.VectorMultiply(r);

                normSqr = ScalarMultiply(r, r);

                for (int i = 0; i < maxIter && normSqr > eps && normSqr != nev; i++)
                {
                    nev = normSqr;

                    alpha =
                        ScalarMultiply(p, r)
                    / //--------------------
                        ScalarMultiply(p, p);

                    for (int j = 0; j < A.N; j ++)
                    {
                        x[j] += alpha * z[j];
                        r[j] -= alpha * p[j];
                    }

                    temp = A.VectorMultiply(r);

                    beta =
                        -1 * ScalarMultiply(p, temp)
                    / //-------------------------
                          ScalarMultiply(p, p);

                    normSqr -= alpha * alpha * ScalarMultiply(p, p);

                    if (i % 2 == 0)
                    {
                        r = VectorDiff(b, A.VectorMultiply(x));
                        z = r;
                        p = A.VectorMultiply(r);
                    }
                    else
                    {
                        for (int j = 0; j < A.N; j++)
                        {
                            z[j] = r[j] + beta * z[j];
                            p[j] = temp[j] + beta * p[j];
                        }
                    }

                    Console.WriteLine("iter: {0}\t nev: {1}", i, Math.Sqrt(normSqr));
                }

                return x;
            }
        }

        static void Main(string[] args)
        {
            //double[] di = new double[] { 1, 5, 8, 12, 15, 18 };
            //int[] ig = new int[] { 1, 1, 2, 2, 4, 5, 7 };
            //int[] jg = new int[] { 1, 1, 2, 4, 2, 3 };
            //double[] ggl = new double[] { 4, 10, 11, 14, 16, 17 };
            //double[] ggu = new double[] { 2, 3, 6, 13, 7, 9 };

            double[] di = new double[] { 1, 5, 8 };
            int[] ig = new int[] { 1, 1, 2, 3, 3 };
            int[] jg = new int[] { 1, 2 };
            double[] ggl = new double[] { 4, 1 };
            double[] ggu = new double[] { 2, 1 };

            SparseMatrix test_matrix = new(di, ig, jg, ggl, ggu);

            //double[] x = new double[] { 1, 1, 1, 1, 1, 1 };
            double[] x = new double[] { 1, 2, 3 };

            double[] y = test_matrix.VectorMultiply(x);

            SLAESolver solver = new();
            double[] res = solver.LOS(test_matrix, y);

            for (int i = 0; i < res.Length; i++)
            {
                Console.WriteLine("{0} : {1}", res[i], y[i]);
            }
        }
    }
}