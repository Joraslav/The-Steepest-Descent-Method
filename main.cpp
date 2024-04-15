#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>

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
vector<T> operator*(type const& c, vector<T> const& v)
{
  vector<type> Rez{v};
  auto const nRows{Rez.size()};
  for (auto i{0u}; i < nRows; ++i)
  {
    Rez[i] = c*Rez[i];
  }
  return Rez;
}

template<class T>
ostream& operator<<(ostream& out, vector<T> const& v)
{
  for (auto const &i : v)
  {
    out << i << '\t';
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

type f(vector<type> const& x)
{
  return pow(x[0]-2,4)+pow(x[0]-2*x[1],2);
}

vector<type> Grad(vector<type> const& x)
{
  vector<type> Rez;
  Rez.push_back(4.*pow(x[0]-2,3)+2.*(x[0]-2.*x[1]));
  Rez.push_back(-4.*(x[0]-2.*x[1]));
  return Rez;
}

vector<type> X_Next(type const& lyam, vector<type> const& x, function<vector<type>(vector<type> const&)> func)
{
  vector<type> Rez;
  for (auto i{0u}; i < x.size(); ++i)
  {
    Rez.push_back(x[i]-lyam*func(x)[i]);
  }
  return Rez;
}

type df_lyam(vector<type> const& x, type const& lyam)
{
  vector<type> x_next = x - lyam*Grad(x);
  return 4*pow(x_next[0]-2,3)*Grad(x)[0]
        +2*(x_next[0]-2*x_next[1])*(Grad(x)[1]-2*Grad(x)[1]);
}

type FindLyam(vector<type> const& x, type const& eps)
{
  type a=0, b=1, c;
  while (abs(a-b)/2 > eps)
  {
    c = (a+b)/2;
    if (df_lyam(x,a) * df_lyam(x,c) <= 0)     {b = c;}
    else      {a = c;}
  }
  return c;  
}

vector<type> X_Solve(vector<type> &x, type const& eps)
{
  type lyam = FindLyam(x,eps);
  vector<type> x_next = x - lyam*Grad(x);
  unsigned short int iter = 1;
  while (Norm(Grad(x)) > eps)
  {
    cout << "Iteration\t" << iter << endl;
    x.clear();
    x = x_next;
    x_next.clear();
    lyam = FindLyam(x,eps);
    x_next = x - lyam*Grad(x);
    iter++;
  }
  return x_next;
}

int main()
{
  vector<type> x_0{5,5};
  type eps = 0.0001;

  vector<type> x = X_Solve(x_0,eps);
  cout << "Ans is\t" << x << endl;

  return 0;
}
