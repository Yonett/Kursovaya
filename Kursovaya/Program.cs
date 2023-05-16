using System;
using System.Net;
using System.Numerics;
using System.Text.Json.Serialization;

namespace Kursovaya
{
    internal class Program
    {
        // Функция f правой части уравнения
        public static double Target(double x, double y, int area)
        {
            double result = 0;

            switch (area)
            {
                case 1:
                    result = 3;
                    break;
                case 2:
                    result = 2 * x;
                    break;
                case 3:
                    result = 2 * x * y;
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
                    result = 2;
                    break;
                case 2:
                    result = 1;
                    break;
                case 3:
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
                    result = 3;
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
                    result = 3;
                    break;
                default:
                    break;
            }

            return result;
        }

        // Функция u истинная
        public static double Actual(double x, double y, int type)
        {
            double result = 0;

            switch(type)
            {
                case 1:
                    result = 1;
                    break;
                case 2:
                    result = x;
                    break;
                case 3:
                    result = x * y;
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

            //for (int i = 0; i < 3; i++)
            //{
            //    double sum = 0;
            //    for (int j = 0; j < 3; j++)
            //    {
            //        sum += G[i, j];
            //    }
            //    Console.WriteLine("[{0}]\tSum in matrix G = {1}", i + 1, sum);
            //}
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

        // Рассчёт локального вектора правой части
        public static void CalcB(Cell cell, List<Node> nodes, double[] b)
        {
            double koef = Math.Abs(cell.detD) / 6;
            double[] f = new double[3];
            for (int i = 0; i < 3; i++)
            {
                f[i] = Target(nodes[cell.v[i]].x, nodes[cell.v[i]].y, cell.area);
                b[i] = koef * f[i];
            }
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
        public static void Consider1(List<Node> nodes, Data data)
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
                    data.b[i] = Actual(nodes[i].x, nodes[i].y, nodes[i].condition1);
                }
            }
        }

        static void Main(string[] args)
        {
            int test = 3;
            
            List<Node> nodes = new();
            List<Cell> cells = new();

            string path = @"C:\Users\User\Documents\kurs_test";

            string[] files = { "nodes.txt", "cells.txt", "condition1.txt" };

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

            // Инициализация класса данных
            Data data = new(nodes.Count, cells.Count, 10000, 1e-30);

            double[] u0 = new double[nodes.Count];
            double[] u1 = new double[nodes.Count];

            for (int i = 0; i < nodes.Count; i++)
            {
                u0[i] = 1;
                u1[i] = 2;
            }

            double[,] G = new double[3, 3]; // Локальная масса жёсткости
            double[,] M = new double[3, 3]; // Локальная масса масс
            double[] b = new double[3];     // Локальный вектор правой части

            // Сборка глобальной матрицы
            foreach (Cell cell in cells)
            {
                CalcDetD(cell, nodes);

                CalcG(cell, nodes, G);
                CalcM(cell, M);

                CalcB(cell, nodes, b);

                for (int i = 0; i < 3; i++)
                {
                    data.b[cell.v[i]] += b[i];
                    for (int j = 0; j < 3; j++)
                        data.global[cell.v[i], cell.v[j]] += G[i, j] + M[i, j];
                }
            }

            for (int i = 0; i < data.nodes; i++)
            {
                for (int j = 0; j < data.nodes - 1; j++)
                    Console.Write("{0} ", data.global[i, j]);
                Console.WriteLine(data.global[i, data.nodes - 1]);
            }

            GenerateSparseGlobal(data);
            Consider1(nodes, data);

            for (int i = 0; i < data.nodes; i++)
            {
                for (int j = 0; j < data.nodes - 1; j++)
                    Console.Write("{0} ", data.global[i, j]);
                Console.WriteLine(data.global[i, data.nodes - 1]);
            }

            for (int i = 0; i < data.nodes; i++)
            {
                Console.WriteLine(data.b[i]);
            }

            SLAESolver solver = new();

            //solver.LOS(data);

            //for (int i = 0; i < data.nodes; i++)
            //    Console.WriteLine(data.x[i]);

            solver.LOS_LUsq(data);

            for (int i = 0; i < data.nodes; i++)
            {
                Console.WriteLine(data.x[i]);
            }
        }
    }
}