#ifndef VECT2D_H
#define VECT2D_H

#if  defined(_MSC_VER)
	#include <type_traits>
#elif defined(__APPLE__)
	#include <type_traits>
#else
	#include <bits/cpp_type_traits.h>
#endif
//#include <g++-3/std/cpp_type_traits.h>
#include <sstream>
//#include "defs.h"

#define SQR(x) ((x)*(x))
using namespace std;
//class ostream;

template <typename T>
class Vect2d
{
   public:
      Vect2d (T xc=0, T yc=0) : x (xc), y (yc) { } 
      template<typename U> Vect2d(const Vect2d<U>&);

      T xcomp () const;
      T ycomp () const;
      void setX(T x);
      void setY(T y);
      void set(T X, T Y);

      template<typename U> T distsq(const Vect2d<U>&);

      template<typename U> Vect2d<T>& operator=(const Vect2d<U>&);
      template<typename U> Vect2d<T>& operator+=(const Vect2d<U>&);
      template<typename U> Vect2d<T>& operator-=(const Vect2d<U>&);
      template<typename U> Vect2d<T>& operator^=(const Vect2d<U>&);
      template<typename U> Vect2d<T>& operator*=(const U&);
      template<typename U> Vect2d<T>& operator/=(const U&);

   private:
      T x,y;
};

template <class T> inline T
Vect2d<T>::xcomp() const
{
   return x;
}

template <class T> inline T
Vect2d<T>::ycomp() const
{
   return y;
}


template <class T> 
template<typename U>
   inline 
Vect2d<T>::Vect2d(const Vect2d<U>& __v)
   : x(__v.x()), y(__v.y()) { }


   template<typename T>
   template<typename U>
   Vect2d<T>&
Vect2d<T>::operator*=(const U& __t)
{
   x *= __t;
   y *= __t;
   return *this;
}

template<typename T>
template<typename U>
   Vect2d<T>&
Vect2d<T>::operator/=(const U& __t)
{
   x /= __t;
   y /= __t;
   return *this;
}

template<typename T>
template<typename U>
   Vect2d<T>&
Vect2d<T>::operator=(const Vect2d<U>& __r)
{
   x = __r.xcomp();
   y = __r.ycomp();
   return *this;
}

template<typename T>
template<typename U>
   Vect2d<T>&
Vect2d<T>::operator+=(const Vect2d<U>& __r)
{
   x += __r.xcomp();
   y += __r.ycomp();
   return *this;
}


template<typename T>
template<typename U>

T Vect2d<T>::distsq(const Vect2d<U>& __r)
{
   return (SQR(x-__r.xcomp()) + SQR(y-__r.ycomp()));
}


template<typename T>
template<typename U>
   Vect2d<T>&
Vect2d<T>::operator-=(const Vect2d<U>& __r)
{
   x -= __r.xcomp();
   y -= __r.ycomp();
   return *this;
}

   template<typename T>
inline void Vect2d<T>::setX(T X)
{
   x = X;
}

   template<typename T>
inline void Vect2d<T>::setY(T Y)
{
   y = Y;
}

   template<typename T>
inline void Vect2d<T>::set(T X, T Y)
{
   x = X;
   y = Y;
}

/*
template<typename T>
template<typename U>
   Vect2d<T>&
Vect2d<T>::operator^=(const Vect2d<U>& __r)
{
   x = y*__r.zcomp() - z* __r.ycomp();
   y = z*__r.xcomp() - x* __r.zcomp();
   z =  x*__r.ycomp() - y* __r.xcomp();
   return *this;
}
*/

// Operators:
template<typename T>
   inline Vect2d<T>
operator+(const Vect2d<T>& __x, const Vect2d<T>& __y)
{ return Vect2d<T> (__x) += __y; }

template<typename T>
   inline Vect2d<T>
operator-(const Vect2d<T>& __x, const Vect2d<T>& __y)
{ return Vect2d<T> (__x) -= __y; }

template<typename T>
   inline Vect2d<T>
operator*(const Vect2d<T>& __x, const T& __y)
{ return Vect2d<T> (__x) *= __y; }

template<typename T>
   inline Vect2d<T>
operator*(const T& __x, const Vect2d<T>& __y)
{ return Vect2d<T> (__y) *= __x; }

/*
template<typename T>
   inline Vect2d<T>
operator^(const Vect2d<T>& u, const Vect2d<T>& v)
{ 
   return Vect2d<T>
      (u.ycomp()*v.zcomp()- u.zcomp()*v.ycomp(), 
       u.zcomp()*v.xcomp()- u.xcomp()*v.zcomp(), 
       u.xcomp()*v.ycomp()- u.ycomp()*v.xcomp()); 
}
*/
template<typename T>
   inline Vect2d<T>
operator/(const Vect2d<T>& __x, const T& __y)
{ return Vect2d<T> (__x) /= __y; }

template<typename T>
   inline bool
operator==(const Vect2d<T>& u, const Vect2d<T>& v)
{ return u.xcomp() == v.xcomp() && u.ycomp() == v.ycomp();}

#if 1
template<typename _Tp, typename _CharT, class _Traits>
   basic_ostream<_CharT, _Traits>&
operator<<(basic_ostream<_CharT, _Traits>& __os, const Vect2d<_Tp>& __v)
{
   basic_ostringstream<_CharT, _Traits> __s;
   __s.flags(__os.flags());
   __s.imbue(__os.getloc());
   __s.precision(__os.precision());
   __s << '[' << __v.xcomp() << "," << __v.ycomp() << ']';
   return __os << __s.str();
}
#else
template <class FLOAT> ostream&
operator << (ostream& os, const Vect2d<FLOAT>& v)
{
     return os << '(' << v.xcomp()  << ',' << v.ycomp() <<  ')';
}
#endif


template<typename T>
   inline T
operator*(const Vect2d<T>& u, const Vect2d<T>& v)
{return (u.xcomp()*v.xcomp() +u.ycomp()*v.ycomp());}


typedef Vect2d<float> vect_f;
typedef Vect2d<double> vect_d;
typedef Vect2d<int> vect_i;

static Vect2d<double> x_hat = Vect2d<double>(1.0,0.0);
static Vect2d<double> y_hat = Vect2d<double>(0.0,1.0);
static Vect2d<float> i_hat = Vect2d<float>(1.0,0.0);
static Vect2d<float> j_hat = Vect2d<float>(0.0,1.0);


#endif//VECT2D_H
