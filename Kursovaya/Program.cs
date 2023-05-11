using System;
using System.Numerics;

namespace Kursovaya 
{
    internal class Program
    {        
        // Класс содержащий данные о СЛАУ задачи
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

        static void Main(string[] args)
        {
            {
                List<Node> nodes = new();
                List<Cell> cells = new();

                string path = @"C:\Users\User\Documents\kurs_test";

                string[] files = { "nodes.txt", "cells.txt", "" };

                string? s;
                StreamReader reader;

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
                        string[] vertexes = s.Split('\t');

                        cells.Add(new Cell(
                            int.Parse(vertexes[0]) - 1, // 1 вершина
                            int.Parse(vertexes[1]) - 1, // 2 вершина
                            int.Parse(vertexes[2]) - 1  // 3 вершина
                            ));
                    }
                    reader.Close();
                }
            }


            //int N = 6;
            //Data data = new(N, 0, 1000, 1e-16);

            //data.ig = new int[] {0, 0, 1, 1, 3, 4, 6 };
            //data.jg = new int[] {0, 0, 1, 3, 1, 2 };

            //data.di = new double[] { 1, 5, 8, 12, 15, 18 };
            //data.ggl = new double[] { 4, 10, 11, 14, 16, 17 };
            //data.ggu = new double[] { 2, 3, 6, 13, 7, 9 };

            //double[] x = new double[] { 1, 1, 1, 1, 1, 1 };

            //int N = 3;
            //Data data = new(N, 0, 1000, 1e-16)
            //{
            //    ig = new int[] { 0, 0, 1, 2 },
            //    jg = new int[] { 0, 1 },

            //    di = new double[] { 1, 5, 8 },
            //    ggl = new double[] { 4, 1 },
            //    ggu = new double[] { 2, 1 }
            //};

            int N = 3;
            Data data = new(N, 0, 10000, 1e-16)
            {
                ig = new int[] { 0, 0, 1, 2 },
                jg = new int[] { 0, 0 },

                di = new double[] { 4, 3, 17 },
                ggl = new double[] { 8, 16 },
                ggu = new double[] { 1, 2 }
            };

            double[] x = new double[] { 1, 1, 1 };

            double[] res = VectorMultiply(data, x);

            for (int i = 0; i < N; i++)
            {
                Console.WriteLine(res[i]);
            }

            SLAESolver solver = new();

            //solver.LOS(data, res);

            //for (int i = 0; i < N; i++)
            //{
            //    Console.WriteLine(data.x[i]);
            //}

            solver.LOS_LU(data, res);

            for (int i = 0; i < N; i++)
            {
                Console.WriteLine(data.x[i]);
            }
        }
    }
}