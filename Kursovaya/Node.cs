using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Kursovaya
{
    // Класс узла расчётной области
    internal class Node
    {
        public double x;
        public double y;
        // Конструктор класса
        public Node(double x, double y)
        {
            this.x = x;
            this.y = y;
        }
    }
}
