using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static Kursovaya.Program;

namespace Kursovaya
{
    // Класс-решатель СЛАУ
    internal class SLAESolver
    {
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

        // Вычисление невязки
        public double CalcDiscrepancy(Data data)
        {
            double sum1 = 0, sum2 = 0;

            data.temp1 = VectorMultiply(data, data.x);

            for (int i = 0; i < data.nodes; i++)
            {
                sum1 += (data.b[i] - data.temp1[i]) * (data.b[i] - data.temp1[i]);
                sum2 += data.b[i] * data.b[i];
            }

            return Math.Sqrt(sum1 / sum2);
        }

        // Разложение LUsq
        public static void LU_sq(Data data)
        {
            for (int i = 0; i < data.l.Length; i++)
            {
                data.l[i] = data.ggl[i];
                data.u[i] = data.ggu[i];
            }

            for (int i = 0; i < data.d.Length; i++)
                data.d[i] = data.di[i];

            for (int i = 0; i < data.d.Length; i++)
            {
                double sumd = 0;
                int i0 = data.ig[i];
                int i1 = data.ig[i + 1];
                for (int k = i0; k < i1; k++)
                {
                    int j = data.jg[k];
                    double sl = 0, su = 0;
                    int j0 = data.ig[j];
                    int j1 = data.ig[j + 1];
                    int ki = i0;
                    int kj = j0;
                    for (; ki < k && kj < j1;)
                    {
                        int jl = data.jg[ki];
                        int ju = data.jg[kj];
                        if (jl == ju)
                        {
                            sl += data.u[kj] * data.l[ki];
                            su += data.l[kj] * data.u[ki];
                            ki++; kj++;
                        }
                        else if (jl < ju) ki++;
                        else kj++;
                    }
                    data.u[k] = (data.u[k] - su) / data.d[j];
                    data.l[k] = (data.l[k] - sl) / data.d[j];
                    sumd += data.u[k] * data.l[k];
                }
                data.d[i] = Math.Sqrt(Math.Abs(data.d[i] - sumd));
            }
        }

        // Прямой ход (стр. 875)
        public static void Straight(Data data, double[] a, double[] c)
        {
            for (int i = 0; i < a.Length; i++)
            {
                double sum = 0;
                int i0 = data.ig[i];
                int i1 = data.ig[i + 1];
                for (int k = i0; k < i1; k++)
                {
                    int j = data.jg[k];
                    sum += a[j] * data.l[k];
                }
                a[i] = (c[i] - sum) / data.d[i];
            }
        }

        // Обратный ход (стр. 876)
        public static void Reverse(Data data, double[] a, double[] c)
        {
            int n = a.Length;
            for (int i = 0; i < n; i++)
                a[i] = c[i];
            for (int i = n - 1; i >= 0; i--)
            {
                int i0 = data.ig[i];
                int i1 = data.ig[i + 1];
                a[i] /= data.d[i];
                for (int k = i1 - 1; k >= i0; k--)
                {
                    int j = data.jg[k];
                    a[j] -= a[i] * data.u[k];
                }
            }
        }

        // Локально-оптимальная схема с неполной факторизацией
        public void LOS_LUsq(Data data)
        {
            int n = data.di.Length;
            double scalar1 = 0;
            double scalar2 = 0;

            int iters = 0;

            LU_sq(data);

            data.temp1 = VectorMultiply(data, data.x);
            for (int i = 0; i < n; i++)
            {
                data.temp2[i] = data.b[i] - data.temp1[i];
            }

            Straight(data, data.r, data.temp2);

            Reverse(data, data.z, data.r);

            data.temp1 = VectorMultiply(data, data.z);

            Straight(data, data.p, data.temp1);

            double nev = ScalarMultiply(data.r, data.r);
            for (int k = 0; k < data.maxIter && nev > data.eps; k++)
            {
                iters++;
                scalar1 = ScalarMultiply(data.p, data.r);
                scalar2 = ScalarMultiply(data.p, data.p);
                double alpha = scalar1 / scalar2;
                for (int i = 0; i < n; i++)
                {
                    data.x[i] += alpha * data.z[i];
                    data.r[i] -= alpha * data.p[i];
                }

                Reverse(data, data.temp1, data.r);

                data.temp2 = VectorMultiply(data, data.temp1);

                Straight(data,data.temp1, data.temp2);

                scalar1 = ScalarMultiply(data.p, data.temp1);

                double beta = -scalar1 / scalar2;

                Reverse(data, data.temp2, data.r);

                for (int i = 0; i < n; i++)
                {
                    data.z[i] = data.temp2[i] + beta * data.z[i];
                    data.p[i] = data.temp1[i] + beta * data.p[i];
                }
                nev = ScalarMultiply(data.r, data.r);
            }
        }

        // Локально-оптимальная схема
        public void LOS(Data data)
        {
            int N = data.nodes;

            for (int i = 0; i < N; i++)
            {
                data.x[i] = 0; // Начальное приближение
            }

            double alpha, beta, nev;

            data.temp1 = VectorMultiply(data, data.x);

            for (int i = 0; i < data.nodes; i++)
            {
                data.r[i] = data.b[i] - data.temp1[i];
                data.z[i] = data.r[i];
            }
            data.p = VectorMultiply(data, data.r);

            nev = ScalarMultiply(data.r, data.r);

            for (int i = 0; i < data.maxIter && Math.Abs(nev) > data.eps; i++)
            {

                alpha = ScalarMultiply(data.p, data.r)
                    / //-------------------------------
                        ScalarMultiply(data.p, data.p);

                for (int j = 0; j < data.nodes; j++)
                {
                    data.x[j] += alpha * data.z[j];
                    data.r[j] -= alpha * data.p[j];
                }

                data.temp1 = VectorMultiply(data, data.r);

                beta = (-1) * ScalarMultiply(data.p, data.temp1)
                   / //-----------------------------------------
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
}
