#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>


struct Gauss
{
public:
   Gauss()
   {
      this->size_ = 0;
      this->G_=std::vector <std::vector<double>>();
      this->b_ = std::vector<double>();
      this->x_ = std::vector<double>();
   }
   Gauss(size_t size,
      std::vector <std::vector<double>> G,
      std::vector<double> b,
      std::vector<double> x
   )
   {
      size_ = size;
      G_ = G;
      b_ = b;
      x_ = x;
   }

   int size() const
   {
      return size_;
   }

   void load_data()
   {
      std::ifstream input;
      input.open("input.txt");
      size_t n = 0;
      input >> n;
      size_ = n;

      G_ = std::vector <std::vector<double>>();
      G_.resize(n);
      for (size_t i = 0; i < size(); i++)
      {
         G_[i].resize(size_);
      }
      b_.reserve(size_);
      x_.reserve(size_);

      for (size_t i = 0; i < size(); ++i)
      {
         for (size_t j = 0; j < size(); ++j)
         {
            double value = 0.0;
            input >> value;
            G_[i][j]= value;
         }
      }

      for (size_t i = 0; i < size(); ++i)
      {
         int value = 0;
         input >> value;
         b_.push_back(value);
      }

      for (size_t i = 0; i < size(); ++i)
      {
         double value = 0.0;
         input >> value;
         x_.push_back(value);
      }
      input.close();
   }

   void transform(int k)
   {
      if (k > 0)
      {
         G_[0][0] = G_[0][0] - (1 - pow(10, -k));
         b_[0] = b_[0] - (1 - pow(10, -k));
      }
   }

   void gauss()
   {
      double zero = 1e-15;
      for (int k = 0; k < size() - 1; k++)
      {
         //находим максимум
         int max = abs(G_[k][k]);
         int index = k;
         for (int j = k + 1; j < size(); j++)
         {
            if (abs(G_[j][k]) > max)
            {
               max = abs(G_[j][k]);
               index = j;
            }
         }
         if (max > zero)
         {
            //меняем столбцы
            if (index != k)
            {
               double temp;
               for (int m = k; m < size(); m++)
               {
                  temp = G_[k][m];
                  G_[k][m] = G_[index][m];
                  G_[index][m] = temp;
               }
               temp = b_[k];
               b_[k] = b_[index];
               b_[index] = temp;
            }
         }
         //вычитаем строки
         for (int j = k + 1; j < size(); j++)
         {
            double coefficient = G_[j][k] / G_[k][k];
            G_[j][k] = 0;
            b_[j] -= coefficient * b_[k];
            for (int i = k + 1; i < size(); i++)
               G_[j][i] -= coefficient * G_[k][i];
         }
      }
      //решаем системы верхнетреугольного вида
      for (int k = size() - 1; k >= 0; k--)
      {
         double sum = 0;
         for (int j = k + 1; j < size(); j++)
            sum += G_[k][j] * b_[j];
         b_[k] = (b_[k] - sum) / G_[k][k];
      }
   };

      void print()
      {
         std::ofstream output;
         output.open("out.txt", std::ios_base::app);
         output << "\nxk:";
         for (auto x : b_)
         {
            output << "\n" << x << " ";
         }
         output << "\n";

         output << "\nx*-xk:";
         for (int i = 0; i < size(); i++)
         {
            output << "\n" << x_[i] - b_[i] << " ";
         }
         output << "\n";
         output.close();
      };

private:
   size_t size_;
   std::vector <std::vector<double>> G_;
   std::vector<double> b_;
   std::vector<double> x_;
};


int main()
{
   Gauss Mx;

   Mx.load_data();
   Mx.gauss();
   Mx.print();

   return 0;
}
