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


        // Класс-решатель СЛАУ методом ЛОС
        class SLAESolver
        {
            // Локально-оптимальная схема с LU-факторизацией
            public void LOS_LU(Data data, double[] b)
            {

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

        
        static void Main(string[] args)
        {
            //int N = 6;
            //Data data = new(N, 0, 1000, 1e-7);

            //data.ig = new int[] {0, 0, 1, 1, 3, 4, 6 };
            //data.jg = new int[] {0, 0, 1, 3, 1, 2 };

            //data.di = new double[] { 1, 5, 8, 12, 15, 18 };
            //data.ggl = new double[] { 4, 10, 11, 14, 16, 17 };
            //data.ggu = new double[] { 2, 3, 6, 13, 7, 9 };

            //double[] x = new double[] { 1, 1, 1, 1, 1, 1 };

            int N = 3;
            Data data = new(N, 0, 1000, 1e-16);

            data.ig = new int[] { 0, 0, 1, 2 };
            data.jg = new int[] { 0, 1 };

            data.di = new double[] { 1, 5, 8 };
            data.ggl = new double[] { 4, 1};
            data.ggu = new double[] { 2, 1 };

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
        }
    }
}