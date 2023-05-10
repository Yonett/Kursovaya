using System;
using System.Numerics;

namespace Kursovaya 
{
    internal class Program
    {
        // Класс содержащий данные о задаче
        public class Data
        {
            public int nodes;    // Кол-во узлов
            public int cells;    // Кол-во элементов
            public int maxIter;  // Максимальное кол-во итерации для решателя СЛАУ
            public double eps;   // Точность решения

            public int[] ig;     // Массив ig разреженой матрицы (кол-во элементов в строке-столбце)
            public int[] jg;     // Массив jg разреженой матрицы (номера столбцов-строк элементов матрицы)

            public double[] di;  // Массив di разреженой матрицы (диагональ)
            public double[] ggl; // Массив ggl разреженой матрицы (нижний треугольник)
            public double[] ggu; // Массив ggu разреженой матрицы (верхний треугольник)

            public double[] d;   // Массив d LU-разложения матрицы (диагональ)
            public double[] l;   // Массив l LU-разложения матрицы (нижний треугольник)
            public double[] u;   // Массив u LU-разложения матрицы (верхний треугольник)

            public double[] r;   // Массив r используемый в ЛОС
            public double[] z;   // Массив z используемый в ЛОС
            public double[] p;   // Массив p используемый в ЛОС

            public double[] b;   // Массив-вектор правой части
            public double[] x;   // Массив-вектор решения

            public double[] temp1, temp2;   // Вспомогательные массивы

            // Конструктор класса данных
            public Data(int nodes, int cells, int maxIter, double eps)
            {
                this.nodes = nodes;
                this.cells = cells;
                this.maxIter = maxIter;
                this.eps = eps;

                ig = new int[nodes + 1];

                di = new double[nodes];
                d = new double[nodes];

                b = new double[nodes];
                x = new double[nodes];

                temp1 = new double[nodes];
                temp2 = new double[nodes];

                r = new double[nodes];
                z = new double[nodes];
                p = new double[nodes];
            }
        }

        // Класс-решатель СЛАУ
        class SLAESolver
        {
            // LU факторизация
            static void CalcLU(Data data)
            {
                data.l = new double[data.ggl.Length];
                data.u = new double[data.ggu.Length];

                for (int i = 0; i < data.di.Length; i++)
                    data.d[i] = data.di[i];
                for (int i = 0; i < data.ggl.Length; i++)
                    data.l[i] = data.ggl[i];
                for (int i = 0; i < data.ggu.Length; i++)
                    data.u[i] = data.ggu[i];

                double sumU, sumL, sumD;
                int N = data.nodes;

                for (int i = 0; i < N; i++)
                {
                    sumD = 0;

                    int begI = data.ig[i];
                    int endI = data.ig[i + 1];
                    for (int igi = begI; igi < endI; igi++)
                    {
                        sumU = 0;
                        sumL = 0;

                        int Jindex = data.jg[igi];

                        for (int igj = begI; igj < igi; igj++)
                        {
                            int begJ = data.ig[Jindex];
                            int endJ = data.ig[Jindex + 1];

                            for (int jgi = begJ; jgi < endJ; jgi++)
                            {
                                if (data.jg[igj] == data.jg[jgi])
                                {
                                    sumL += data.l[igj] * data.u[jgi];
                                    sumU += data.l[jgi] * data.u[igj];
                                }
                            }
                        }
                        data.l[igi] -= sumL;
                        data.u[igi] -= sumU;
                        data.u[igi] /= data.d[Jindex];
                        sumD += data.l[igi] * data.u[igi];
                    }

                    data.d[i] -= sumD;
                }
            }

            // Прямой ход Ly = F
            static void CalcDir(Data data, double[] y, double[] F)
            {
                double sum, buf;
                int N = data.nodes;

                for (int i = 0; i < N; i++)
                {
                    y[i] = F[i];
                }

                for (int i = 0; i < N; i++)
                {
                    sum = 0;

                    int begI = data.ig[i];
                    int endI = data.ig[i + 1];

                    for (int igi = begI; igi < endI; igi++)
                    {
                        sum += y[data.jg[igi]] * data.l[igi];
                    }

                    buf = y[i] - sum;
                    y[i] = buf / data.d[i];
                }
            }

            // Обратный ход Ux = y
            static void CalcRev(Data data, double[] X, double[] y)
            {
                int N = data.nodes;

                for (int i = 0; i < N; i++)
                {
                    X[i] = y[i];
                }

                for (int i = N - 1; i >= 0; i--)
                {
                    int begI = data.ig[i];
                    int endI = data.ig[i + 1];

                    for (int igi = begI; igi < endI; igi++)
                    {
                        X[data.jg[igi]] -= X[i] * data.u[igi];
                    }
                }
            }

            // Локально-оптимальная схема с LU-факторизацией
            public void LOS_LU(Data data, double[] b)
            {
                double alpha, beta, norm, temp_nev = 0;

                int N = data.nodes;

                for (int i = 0; i < N; i++)
                {
                    data.x[i] = 1;
                }

                CalcLU(data);
                // A * x0
                data.temp1 = VectorMultiply(data, data.x);

                // f - A * x0
                for (int i = 0; i < N; i++)
                {
                    data.temp1[i] = data.b[i] - data.temp1[i];
                }

                // L * r0 = f - A * x0
                CalcDir(data, data.r, data.temp1);

                // U * z0 = r0
                CalcRev(data, data.z, data.r);

                // A * z0
                data.temp1 = VectorMultiply(data, data.z);

                // L * p0 = A * z0
                CalcDir(data, data.p, data.temp1);

                norm = ScalarMultiply(data.x, data.x);

                int k;

                for (k = 0; k < data.maxIter && Math.Abs(norm) > data.eps && norm != temp_nev; k++)
                {
                    // если невязка не изменилась, то выходим из итерационного процесса
                    temp_nev = norm;

                    alpha = ScalarMultiply(data.p, data.x) / ScalarMultiply(data.p, data.p);

                    for (int i = 0; i < N; i++)
                    {
                        data.x[i] += alpha * data.z[i];
                        data.r[i] -= alpha * data.p[i];
                    }

                    // U * temp = x
                    CalcRev(data, data.temp1, data.x);

                    // A * U-1 * x = temp0
                    data.temp2 = VectorMultiply(data, data.temp1);

                    // L * temp = A * U-1 * x 
                    CalcDir(data, data.temp1, data.temp2);

                    double Scal = ScalarMultiply(data.p, data.p);
                    if (Math.Abs(Scal) < 1e-60) break;

                    beta = -1 * ScalarMultiply(data.p, data.temp1) / Scal;

                    // U * temp0 = x
                    CalcRev(data, data.temp2, data.x);

                    norm = norm - alpha * alpha * ScalarMultiply(data.p, data.p);

                    for (int i = 0; i < N; i++)
                    {
                        data.z[i] = data.temp2[i] + beta * data.z[i];
                        data.p[i] = data.temp1[i] + beta * data.p[i];
                    }
                }
            }

            // Локально-оптимальная схема
            public void LOS(Data data, double[] b)
            {
                int N = data.nodes;

                for (int i = 0; i < N; i++)
                {
                    data.x[i] = 1; // Начальное приближение
                }

                double alpha, beta, nev = 0;

                data.temp1 = VectorMultiply(data, data.x);

                for (int i = 0; i < data.nodes; i++)
                {
                    data.r[i] = b[i] - data.temp1[i];
                    data.z[i] = data.r[i];
                }
                data.p = VectorMultiply(data, data.r);

                nev = ScalarMultiply(data.r, data.r);

                for (int i = 0; i < data.maxIter && Math.Abs(nev) > data.eps; i++)
                {

                    alpha = ScalarMultiply(data.p, data.r)
                        / //------------------------------
                            ScalarMultiply(data.p, data.p);

                    for (int j = 0; j < data.nodes; j++)
                    {
                        data.x[j] += alpha * data.z[j];
                        data.r[j] -= alpha * data.p[j];
                    }

                    data.temp1 = VectorMultiply(data, data.r);

                    beta = (-1) * ScalarMultiply(data.p, data.temp1)
                       / //----------------------------------------
                           ScalarMultiply(data.p, data.p);

                    for (int j = 0; j < data.nodes; j++)
                    {
                        data.z[j] = data.r[j] + beta * data.z[j];
                        data.p[j] = data.temp1[j] + beta * data.p[j];
                    }

                    nev = ScalarMultiply(data.r, data.r);
                }
            }
        }

        // Скалярное произведение векторов (x, y)
        public static double ScalarMultiply(double[] x, double[] y)
        {
            double result = 0;
            for (int i = 0; i < x.Length; i++)
            {
                result += x[i] * y[i];
            }
            return result;
        }
        
        // Умножение разреженой матрицы на вектор
        public static double[] VectorMultiply(Data data, double[] x)
        {
            double[] y = new double[data.nodes];

            for (int i = 0; i < data.nodes; i++)
            {
                y[i] = x[i] * data.di[i];

                for (int j = data.ig[i]; j < data.ig[i + 1]; j++)
                {
                    y[i] += data.ggl[j] * x[data.jg[j]];
                    y[data.jg[j]] += data.ggu[j] * x[i];
                }
            }

            return y;
        }
        
        static void Main(string[] args)
        {
            //int N = 6;
            //Data data = new(N, 0, 1000, 1e-16);

            //data.ig = new int[] {0, 0, 1, 1, 3, 4, 6 };
            //data.jg = new int[] {0, 0, 1, 3, 1, 2 };

            //data.di = new double[] { 1, 5, 8, 12, 15, 18 };
            //data.ggl = new double[] { 4, 10, 11, 14, 16, 17 };
            //data.ggu = new double[] { 2, 3, 6, 13, 7, 9 };

            //double[] x = new double[] { 1, 1, 1, 1, 1, 1 };

            int N = 3;
            Data data = new(N, 0, 1000, 1e-16)
            {
                ig = new int[] { 0, 0, 1, 2 },
                jg = new int[] { 0, 1 },

                di = new double[] { 1, 5, 8 },
                ggl = new double[] { 4, 1 },
                ggu = new double[] { 2, 1 }
            };

            double[] x = new double[] { 1, 2, 3 };

            double[] res = VectorMultiply(data, x);

            for (int i = 0; i < N; i++)
            {
                Console.WriteLine(res[i]);
            }

            SLAESolver solver = new();

            solver.LOS(data, res);

            for (int i = 0; i < N; i++)
            {
                Console.WriteLine(data.x[i]);
            }

            solver.LOS_LU(data, res);

            for (int i = 0; i < N; i++)
            {
                Console.WriteLine(data.x[i]);
            }
        }
    }
}