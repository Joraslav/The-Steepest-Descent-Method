#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;
using type = long double;

template <class T>
vector<type> operator+(vector<T> const& l, vector<T> const& r)
{
  if (l.size() != r.size()){throw std::invalid_argument("wrong sizes");}  //Проверка размерности

  auto Rez{l};
  auto const nRows{Rez.size()};

  for (auto i{0u}; i < nRows; ++i)
  {
    Rez[i] += r[i];
  }
  
  return Rez;
}

template <class T>
vector<type> operator-(vector<T> const& l, vector<T> const& r)
{
  if (l.size() != r.size()){throw std::invalid_argument("wrong sizes");}  //Проверка размерности

  auto Rez{l};
  auto const nRows{Rez.size()};

  for (auto i{0u}; i < nRows; ++i)
  {
    Rez[i] -= r[i];
  }
  
  return Rez;
}

template<class T>
type operator*(vector<T> const& l, vector<T> const& r)
{
  if (l.size() != r.size()){throw std::invalid_argument("wrong sizes");}  //Проверка размерности

  type Rez = 0.0;
  auto const nRows{l.size()};

  for (auto i{0u}; i < nRows; ++i)
  {
    Rez += l[i] * r[i];
  }
  
  return Rez;
}

template<class T>
ostream& operator<<(ostream& out, vector<T> const& v)
{
  for (auto const &i : v)
  {
    out << i << '\t' << '\t';
  }
  return out;
}

type Norm(vector<type> const& v)
{
  auto Rez = 0.0;
  auto ScalarProd = v*v;

  Rez = sqrt(ScalarProd);

  return Rez;
}

vector<type> Grad_f(vector<type> const x)
{
    vector<type> Rez;

    Rez.push_back(4.*pow(x[0]-2,3) + 2.*(x[0]-2.*x[1]));
    Rez.push_back(-4.*(x[0]-2.*x[1]));

    return Rez;
}

type f(vector<type> const x)
{
    type Rez;
    Rez = pow(x[0]-2.,4) + pow(x[0]-2.*x[1],2);
    return Rez;
}

vector<type> X_Next(vector<type> v, type p)
{
    vector<type> Rez;
    Rez.push_back(v[0] - p*Grad_f(v)[0]);
    Rez.push_back(v[1] - p*Grad_f(v)[1]);
    return Rez;
}

type FindLyam(vector<type> x_k,type const e)
{
    type a=0, b=1, c;
    while (abs(a-b) > e)
    {
        c = (a+b)/2;
        vector<type> x_ki_c = X_Next(x_k,c);
        // vector<type> x_ki_b = X_Next(x_k,b);
        if (f(x_k) * f(x_ki_c) <= 0)
        {
            a = c;
        }
        else
        {
            b = c;
        }
    }
    return c;
}

vector<type> Spusk(vector<type> x_k, type const eps)
{
    vector<type> x_next;
    type lyam = FindLyam(x_k,eps);
    x_next = X_Next(x_k,lyam);
    unsigned short int iter = 0;
    while (Norm(x_next-x_k) > eps)
    {
        if (iter%100 == 0)
        {
            cout << "Iteration\t" << iter << endl;
        }
        x_k.clear();
        x_k = x_next;
        x_next.clear();
        type lyam_next = FindLyam(x_k,eps);
        x_next = X_Next(x_k,lyam_next);
        iter++;
    }
    return x_next;
}

int main()
{
    vector<type> x_f{2,1};
    vector<type> x_0{2.2,0.5};
    vector<type> grad_test = Grad_f(x_0);
    cout << "Grad in x_0\t" << grad_test << endl;

    type cen = 0.5, st = 0;
    type f0 = f(x_f);
    cout << "Solve\t" << f0 << endl;

    type lyam = FindLyam(x_0,0.01);
    cout << "Test Luambda\t" << lyam << endl;

    vector<type> ans = Spusk(x_0,0.0001);
    cout << "x = " << ans << endl;

    return 0;
}
