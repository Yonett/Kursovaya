using System;
using System.Net;
using System.Numerics;
using System.Text.Json.Serialization;

namespace Kursovaya
{
    internal class Program
    {
        // Функция f правой части уравнения
        public static double Target(double x, double y, double t, int area)
        {
            double result = 0;

            switch (area)
            {
                case 1:
                    result = 2;
                    break;
                case 2:
                    result = 2 * x;
                    break;
                case 3:
                    result = 2 * x * y;
                    break;
                case 4:
                    result = 2 * x * t + 2 * x;
                    break;
                case 5:
                    result = 2 * x * x * t + 2 * x * x - 1;
                    break;
                case 6:
                    result = 2 * x * x * y * y * t + 2 * x * x * y * y - x * x - y * y;
                    break;
                default:
                    break;
            }

            return result;
        }

        // Параметр лямбда
        public static double Lambda(int area)
        {
            double result = 0;

            switch (area)
            {
                case 1:
                    result = 1;
                    break;
                case 2:
                    result = 1;
                    break;
                case 3:
                    result = 1;
                    break;
                case 4:
                    result = 1;
                    break;
                case 5:
                    result = 1;
                    break;
                case 6:
                    result = 1;
                    break;
                default:
                    break;
            }

            return result;
        }

        // Параметр гамма
        public static double Gamma(int area)
        {
            double result = 0;

            switch (area)
            {
                case 1:
                    result = 3;
                    break;
                case 2:
                    result = 2;
                    break;
                case 3:
                    result = 2;
                    break;
                case 4:
                    result = 2;
                    break;
                default:
                    break;
            }

            return result;
        }

        // Параметр хи
        public static double Hi(int area)
        {
            double result = 0;

            switch (area)
            {
                case 1:
                    result = 2;
                    break;
                case 2:
                    result = 2;
                    break;
                case 3:
                    result = 2;
                    break;
                case 4:
                    result = 2;
                    break;
                case 5:
                    result = 2;
                    break;
                case 6:
                    result = 2;
                    break;
                default:
                    break;
            }

            return result;
        }

        // Параметр сигма
        public static double Sigma(int area)
        {
            double result = 0;

            switch (area)
            {
                case 1:
                    result = 2;
                    break;
                case 2:
                    result = 2;
                    break;
                case 3:
                    result = 2;
                    break;
                case 4:
                    result = 2;
                    break;
                case 5:
                    result = 2;
                    break;
                case 6:
                    result = 2;
                    break;
                default:
                    break;
            }

            return result;
        }

        // Функция u истинная
        public static double Actual(double x, double y, double t, int type)
        {
            double result = 0;

            switch(type)
            {
                case 1:
                    result = t;
                    break;
                case 2:
                    result = x * t;
                    break;
                case 3:
                    result = x * y * t;
                    break;
                case 4:
                    result = x * t * t;
                    break;
                case 5:
                    result = x * x * t * t;
                    break;
                case 6:
                    result = x * x * y * y * t * t;
                    break;
                default:
                    break;
            }

            return result;
        }

        // Рассчёт detD для конечного элемента (5.74)
        public static void CalcDetD(Cell cell, List<Node> nodes)
        {
            cell.detD += (nodes[cell.v[1]].x - nodes[cell.v[0]].x) * (nodes[cell.v[2]].y - nodes[cell.v[0]].y);
            cell.detD -= (nodes[cell.v[2]].x - nodes[cell.v[0]].x) * (nodes[cell.v[1]].y - nodes[cell.v[0]].y);
        }

        // Рассчёт коэффициентов alpha для конечного элемента (5.69, 5.75)
        public static void CalcAlphas(Cell cell, List<Node> nodes)
        {
            if (cell.detD != 0)
            {
                cell.alpha[0, 0] = (nodes[cell.v[1]].x * nodes[cell.v[2]].y - nodes[cell.v[2]].x * nodes[cell.v[1]].y) / cell.detD;
                cell.alpha[0, 1] = (nodes[cell.v[1]].y - nodes[cell.v[2]].y) / cell.detD;
                cell.alpha[0, 2] = (nodes[cell.v[2]].x - nodes[cell.v[1]].x) / cell.detD;

                cell.alpha[1, 0] = (nodes[cell.v[2]].x * nodes[cell.v[0]].y - nodes[cell.v[0]].x * nodes[cell.v[2]].y) / cell.detD;
                cell.alpha[1, 1] = (nodes[cell.v[2]].y - nodes[cell.v[0]].y) / cell.detD;
                cell.alpha[1, 2] = (nodes[cell.v[0]].x - nodes[cell.v[2]].x) / cell.detD;

                cell.alpha[2, 0] = (nodes[cell.v[0]].x * nodes[cell.v[1]].y - nodes[cell.v[1]].x * nodes[cell.v[0]].y) / cell.detD;
                cell.alpha[2, 1] = (nodes[cell.v[0]].y - nodes[cell.v[1]].y) / cell.detD;
                cell.alpha[2, 2] = (nodes[cell.v[1]].x - nodes[cell.v[0]].x) / cell.detD;
            }
        }

        // Рассчёт локальной матрицы жёсткости
        public static void CalcG(Cell cell, List<Node> nodes, double[,] G)
        {
            double lambda = Lambda(cell.area);

            // Рассчёт коэффициентов alpha для конечного элемента (5.69, 5.75)
            cell.alpha[0, 0] = (nodes[cell.v[1]].x * nodes[cell.v[2]].y - nodes[cell.v[2]].x * nodes[cell.v[1]].y) / cell.detD;
            cell.alpha[0, 1] = (nodes[cell.v[1]].y - nodes[cell.v[2]].y) / cell.detD;
            cell.alpha[0, 2] = (nodes[cell.v[2]].x - nodes[cell.v[1]].x) / cell.detD;

            cell.alpha[1, 0] = (nodes[cell.v[2]].x * nodes[cell.v[0]].y - nodes[cell.v[0]].x * nodes[cell.v[2]].y) / cell.detD;
            cell.alpha[1, 1] = (nodes[cell.v[2]].y - nodes[cell.v[0]].y) / cell.detD;
            cell.alpha[1, 2] = (nodes[cell.v[0]].x - nodes[cell.v[2]].x) / cell.detD;

            cell.alpha[2, 0] = (nodes[cell.v[0]].x * nodes[cell.v[1]].y - nodes[cell.v[1]].x * nodes[cell.v[0]].y) / cell.detD;
            cell.alpha[2, 1] = (nodes[cell.v[0]].y - nodes[cell.v[1]].y) / cell.detD;
            cell.alpha[2, 2] = (nodes[cell.v[1]].x - nodes[cell.v[0]].x) / cell.detD;

            double detD = Math.Abs(cell.detD);

            // Рассчёт значений компонент матрицы жесткости (5.80)
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    G[i, j] = lambda * detD * (cell.alpha[i, 1] * cell.alpha[j, 1] + cell.alpha[i, 2] * cell.alpha[j, 2]) / 2;
                }
            }
        }

        // Рассчёт локальной матрицы массы
        public static void CalcM(Cell cell, double[,] M)
        {
            double gamma = Gamma(cell.area);

            double detD = Math.Abs(cell.detD);

            double[,] values = new double[3, 3]
            {
                { 2, 1, 1 },
                { 1, 2, 1 },
                { 1, 1, 2 }
            };

            // Рассчёт значений компонент матрицы масс (5.81)
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    M[i, j] = gamma * detD * values[i, j] / 24;
                }
            }
        }

        public static void CalcMHi(Cell cell, double[,] M)
        {
            double hi = Hi(cell.area);

            double detD = Math.Abs(cell.detD);

            double[,] values = new double[3, 3]
            {
                { 2, 1, 1 },
                { 1, 2, 1 },
                { 1, 1, 2 }
            };

            // Рассчёт значений компонент матрицы масс (5.81)
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    M[i, j] = hi * detD * values[i, j] / 24;
                }
            }
        }

        public static void CalcMSigma(Cell cell, double[,] M)
        {
            double sigma = Sigma(cell.area);

            double detD = Math.Abs(cell.detD);

            double[,] values = new double[3, 3]
            {
                { 2, 1, 1 },
                { 1, 2, 1 },
                { 1, 1, 2 }
            };

            // Рассчёт значений компонент матрицы масс (5.81)
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    M[i, j] = sigma * detD * values[i, j] / 24;
                }
            }
        }

        // Рассчёт локального вектора правой части
        public static void CalcB(Cell cell, List<Node> nodes, double[] b, double t)
        {
            double detD = Math.Abs(cell.detD);
            double[] f = new double[3];
            for (int i = 0; i < 3; i++)
                f[i] = Target(nodes[cell.v[i]].x, nodes[cell.v[i]].y, t, cell.area);

            b[0] = detD * (f[0] / 12 + f[1] / 24 + f[2] / 24);
            b[1] = detD * (f[0] / 24 + f[1] / 12 + f[2] / 24);
            b[2] = detD * (f[0] / 24 + f[1] / 24 + f[2] / 12);
        }

        // Генерация портрета матрицы
        public static void GenerateSparseGlobal(Data data)
        {
            int size = -1;

            for (int i = 0; i < data.nodes; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    if (i == j)
                        data.di[i] = data.global[i, j];
                    else
                    {
                        if (data.global[i, j] != 0)
                        {
                            size++;
                            data.ggl[size] = data.global[i, j];
                            data.ggu[size] = data.global[i, j];
                            data.ig[i + 1] = size + 1;
                            data.jg[size] = j;
                        }
                        else
                        {
                            data.ig[i + 1] = size + 1;
                        }
                    }
                }
            }
        }

        // Учёт первых краевых условий
        public static void Consider1(List<Node> nodes, Data data, double t)
        {
            for (int i = 0; i < data.nodes; i++)
            {
                // Если для узла задано первое краевое условие
                if (nodes[i].condition1 > 0)
                {
                    for (int k = 0; k < data.nodes; k++)
                        if (k == i)
                            data.global[k, k] = 1;
                        else
                            data.global[i, k] = 0;

                    // Ставим на i-ом (глобальном) элементе диагонали единицу
                    data.di[i] = 1;

                    // Обнуляем внедиагональные элементы i-ой строки в ggl
                    for (int j = data.ig[i]; j < data.ig[i + 1]; j++)
                        data.ggl[j] = 0;

                    // Обнуляем внедиагональные элементы i-ой строки в ggu
                    for (int j = 0; j < data.ig[data.nodes]; j++)
                        if (data.jg[j] == i)
                            data.ggu[j] = 0;

                    // В правой части замещаем значение на значение функции первого краевого условия
                    data.b[i] = Actual(nodes[i].x, nodes[i].y, t, nodes[i].condition1);
                }
            }
        }

        public static double[] MatrixVector(double[,] matrix, double[] vector)
        {
            int n = vector.Length;
            double[] result = new double[n];

            for (int i = 0; i < n; i++)
                for (int k = 0; k < n; k++)
                    result[i] += matrix[i, k] * vector[k];

            return result;
        }

        static void Main(string[] args)
        {
            int test = 1;

            List<Node> nodes = new(); // Узлы сетки
            List<Cell> cells = new(); // Конечные элементы
            List<double[]> solutions = new(); // Решения по времени
            List<double> times = new(); // Временные точки

            double deltaT, deltaT1, deltaT0;
            double actual_value, t, error, d;

            string path = @"C:\Users\User\Documents\kurs_test";

            string[] files = { "nodes.txt", "cells.txt", "condition1.txt", "t.txt" };

            string? s;
            StreamReader reader;
            StreamWriter writer;

            // Заполнение списка узлов расчётной обалсти
            using (reader = new(Path.Combine(path, files[0])))
            {
                while ((s = reader.ReadLine()) != null)
                {
                    string[] coords = s.Split('\t');

                    nodes.Add(new Node(
                        double.Parse(coords[0]), // x
                        double.Parse(coords[1])  // y
                        ));
                }
                reader.Close();
            }

            // Заполнение списка конечных элементов
            using (reader = new(Path.Combine(path, files[1])))
            {
                while ((s = reader.ReadLine()) != null)
                {
                    string[] values = s.Split('\t');

                    cells.Add(new Cell(
                        int.Parse(values[0]) - 1, // 1 вершина
                        int.Parse(values[1]) - 1, // 2 вершина
                        int.Parse(values[2]) - 1, // 3 вершина
                        /*int.Parse(values[3])*/ test      // Область
                        ));
                }
                reader.Close();
            }

            // Считывание данных о первых краевых условиях
            using (reader = new(Path.Combine(path, files[2])))
            {
                int i = 0;
                while ((s = reader.ReadLine()) != null)
                {
                    string[] values = s.Split('\t');
                    i = int.Parse(values[0]) - 1;
                    //nodes[i].condition1 = int.Parse(values[1]);
                    nodes[i].condition1 = test;
                }
                reader.Close();
            }

            // Считывание данных о временных промежутках
            using (reader = new(Path.Combine(path, files[3])))
            {
                while ((s = reader.ReadLine()) != null)
                    times.Add(double.Parse(s));
                reader.Close();
            }

            // Инициализация класса данных
            Data data = new(nodes.Count, cells.Count, 10000, 1e-30);

            double[] temp1 = new double[nodes.Count];
            double[] temp2 = new double[nodes.Count];
            double[] temp3 = new double[nodes.Count];
            double[] temp4 = new double[nodes.Count];

            for (int i = 0; i < 2; i++)
            {
                solutions.Add(new double[nodes.Count]);
                for (int j = 0; j < data.nodes; j++)
                {
                    solutions[i][j] = Actual(nodes[j].x, nodes[j].y, times[i], test);
                }
            }

            double[,] G = new double[3, 3]; // Локальная масса жёсткости
            double[,] MHi = new double[3, 3]; // Локальная масса масс хи
            double[,] MSigma = new double[3, 3]; // Локальная масса масс сигма
            double[] b = new double[3];     // Локальный вектор правой части

            foreach (Cell cell in cells)
            {
                CalcDetD(cell, nodes);

                CalcG(cell, nodes, G);
                CalcMHi(cell, MHi);
                CalcMSigma(cell, MSigma);

                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        data.G[cell.v[i], cell.v[j]] += G[i, j];
                        data.MHi[cell.v[i], cell.v[j]] += MHi[i, j];
                        data.MSigma[cell.v[i], cell.v[j]] += MSigma[i, j];
                    }
                }
            }

            for (int k = 2; k < times.Count; k++)
            {
                error = 0;
                d = 0;
                solutions.Add(new double[nodes.Count]);

                deltaT = times[k] - times[k - 2];
                deltaT1 = times[k - 1] - times[k - 2];
                deltaT0 = times[k] - times[k - 1];

                temp1 = MatrixVector(data.MHi, solutions[k - 2]);
                temp2 = MatrixVector(data.MHi, solutions[k - 1]);
                temp3 = MatrixVector(data.MSigma, solutions[k - 2]);
                temp4 = MatrixVector(data.MSigma, solutions[k - 1]);

                // Сборка глобальной матрицы
                foreach (Cell cell in cells)
                {
                    CalcB(cell, nodes, b, times[k]);

                    for (int i = 0; i < 3; i++)
                    {
                        data.b[cell.v[i]] += b[i];
                    }
                }

                for (int i = 0; i < data.nodes; i++)
                {
                    data.b[i] -= temp1[i] * 2 / (deltaT1 * deltaT);
                    data.b[i] += temp2[i] * 2 / (deltaT1 * deltaT0);

                    data.b[i] -= temp3[i] * deltaT0 / (deltaT1 * deltaT);
                    data.b[i] += temp4[i] * deltaT / (deltaT1 * deltaT0);

                    for (int j = 0; j < data.nodes; j++)
                        data.global[i, j] +=
                            data.MHi[i, j] * 2 / (deltaT * deltaT0) +
                            data.MSigma[i, j] * (deltaT + deltaT0) / (deltaT * deltaT0) -
                            data.G[i, j];
                }

                GenerateSparseGlobal(data);
                Consider1(nodes, data, times[k]);

                SLAESolver solver = new();

                solver.LOS_LUsq(data);

                Console.WriteLine("{0}", times[k]);
                Console.WriteLine();
                for (int i = 0; i < data.nodes; i++)
                {
                    Console.WriteLine(data.x[i]);
                    solutions[k][i] = data.x[i];
                }
                Console.WriteLine();

                using (writer = new(Path.Combine(path, test.ToString(), times[k].ToString() + ".csv")))
                {
                    for (int i = 0; i < data.x.Length; i++)
                    {
                        actual_value = Actual(nodes[i].x, nodes[i].y, times[k], cells[0].area);

                        t = Math.Abs(actual_value - data.x[i]);

                        error += t * t;
                        d += actual_value * actual_value;

                        writer.WriteLine("{0:E}; {1:E}; {2:E}; {3:E}; {4:E}",
                            nodes[i].x, nodes[i].y, actual_value, data.x[i], t);

                    }
                    writer.WriteLine("Error = {0:E}", Math.Sqrt(error / d));
                }

                for (int i = 0; i < data.nodes; i++)
                {
                    data.b[i] = 0;

                    for (int j = 0; j < data.nodes; j++)
                        data.global[i, j] = 0;
                }
            }
        }
    }
}