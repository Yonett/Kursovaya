﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kursovaya
{
    // Класс содержащий данные о СЛАУ задачи
    internal class Data
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

        public double[,] G; // Глобальная матрица жесткости
        public double[,] MHi; // Глобальная матрица массы компоненты хи
        public double[,] MSigma; // Глобальная матрица массы компоненты сигма
        public double[,] global; // Глобальная матрица A

        public double[] temp1, temp2;   // Вспомогательные массивы

        // Конструктор класса данных
        public Data(int nodes, int cells, int maxIter, double eps)
        {
            this.nodes = nodes;
            this.cells = cells;
            this.maxIter = maxIter;
            this.eps = eps;

            int arrSize = (nodes * (nodes - 1)) / 2;

            ig = new int[nodes + 1];
            jg = new int[arrSize];

            di = new double[nodes];
            ggl = new double[arrSize];
            ggu = new double[arrSize];

            d = new double[nodes];
            l = new double[arrSize];
            u = new double[arrSize];

            b = new double[nodes];
            x = new double[nodes];

            G = new double[nodes, nodes]; 
            MHi = new double[nodes, nodes];
            MSigma = new double[nodes, nodes];
            global = new double[nodes, nodes];

            temp1 = new double[nodes];
            temp2 = new double[nodes];

            r = new double[nodes];
            z = new double[nodes];
            p = new double[nodes];
        }
    }
}
