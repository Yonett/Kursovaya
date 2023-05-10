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
            // Локально-оптимальная схема
            public double[] LOS(SparseMatrix A, double[] b)
            {
                int maxIter = 100000;
                double[] x = new double[A.N];

                for (int i = 0; i < A.N; i++)
                {
                    x[i] = 1;
                }

                double[] r, z, p, temp;

                double alpha, beta, pp, nev = 0;
                double eps = 1e-16;

                r = new double[A.N];
                z = new double[A.N];

                temp = A.VectorMultiply(x);

                for (int i = 0; i < A.N; i++)
                {
                    r[i] = b[i] - temp[i];
                    z[i] = r[i];
                }
                p = A.VectorMultiply(r);

                nev = ScalarMultiply(r, r);

                for (int i = 0; i < maxIter && Math.Abs(nev) > eps; i++)
                {
                    pp = ScalarMultiply(p, p);

                    alpha = ScalarMultiply(p, r) / pp;

                    for (int j = 0; j < A.N; j++)
                    {
                        x[j] += alpha * z[j];
                        r[j] -= alpha * p[j];
                    }

                    temp = A.VectorMultiply(r);

                    beta = -ScalarMultiply(p, temp) / pp;

                    for (int j = 0; j < A.N; j++)
                    {
                        z[j] = r[j] + beta * z[j];
                        p[j] = temp[j] + beta * p[j];
                    }

                    nev = ScalarMultiply(r, r);
                }

                return x;
            }
            
            // Локально-оптимальная схема с LU-факторизацией
            //public double[] LOS_LU(SparseMatrix A, double[] b)
            //{
                
            //}
        }


        static void Main(string[] args)
        {
            double[] di = new double[] { 1, 5, 8, 12, 15, 18 };
            int[] ig = new int[] { 1, 1, 2, 2, 4, 5, 7 };
            int[] jg = new int[] { 1, 1, 2, 4, 2, 3 };
            double[] ggl = new double[] { 4, 10, 11, 14, 16, 17 };
            double[] ggu = new double[] { 2, 3, 6, 13, 7, 9 };

            //double[] di = new double[] { 1, 5, 8 };
            //int[] ig = new int[] { 1, 1, 2, 3 };
            //int[] jg = new int[] { 1, 2 };
            //double[] ggl = new double[] { 4, 1 };
            //double[] ggu = new double[] { 2, 1 };

            SparseMatrix test_matrix = new(di, ig, jg, ggl, ggu);

            double[] x = new double[] { 1, 2, 3, 4, 5, 6 };
            //double[] x = new double[] { 1, 2, 3 };

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