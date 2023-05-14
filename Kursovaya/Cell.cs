using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kursovaya
{
    // Класс конечного элемента расчётной области
    internal class Cell
    {
        public int[] v = new int[3]; // Вершины элемента
        public double[,] alpha = new double[3, 3]; // Коэффициенты alpha
        public double detD = 0; // Определитель матрицы D на данном элементе
        public int area = 0;

        // Конструктор класса
        public Cell(int v0, int v1, int v2, int area)
        {
            v[0] = v0;
            v[1] = v1;
            v[2] = v2;
            this.area = area;
        }
    }
}
