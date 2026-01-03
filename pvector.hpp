#ifndef _PVECTOR_
#define _PVECTOR_
// EXERCISE create an efficient vector class in C++ with methods and overloaded operators which you can find in the TODO LIST.
// TODO LIST:
// * = TODO, X = done, E = exercise
// [X] standard constructor/destructor
// [X] constructor (with overloading->initializer_list)
// [X] show()
// [X] overloading of = assignment 
// [X] overloading of + operator for addition of vectors
// [X] overloading of - operator for substraction of vectors
// [X] sum() method to sum two vectors
// [X] get( ) (get i-th element)
// [X] set( , ) (set i-th element)
// [X] overloading of += and -= operators with vectors
// [X] overloading of * operator for scalar product (vector times vector)
// [X] norm() norm of vector, implemented by using overloaded operator "*"
// [X] overloading of * operator for "vector times scalar" product
// [X] overloading of / operator for "vector divided by scalar" product
// [X] overloading of *= and /= (with scalars)
// [X] overloading of operator == (check whether two vectors are equal)
// [X] overloading of ^ operator for cross product 
// [X] overloading of (...) operator () to use vector as C arrays, e.g. v(0)=1
// [X] "scalar times vector" through friend function
// [X] rint(), i.e. if V is a vector rint(V) is the vector with components rounded to nearest integer
// [X] mulcw, multiplication element wise: C=mulcw(A,B), i.e. C(i) = A(i)*B(i)
// [X] divcw, division element wise: C=divcw(A,B), i.e. C(i) = A(i)*B(i)
// [X] overloading of << for output (friend function)
//
// Discuss random number generator class with Mersenne-Twister
// [X] random(L), random vector inside a box of length L
// [X] random_orient() random unit on unit sphere (marsaglia) 

#include<iostream> // input/output
#include<cmath>
#include<string>
#include "./randnumgen.hpp"
// STRATEGIES TO MAKE CLASS GENERIC
//#define NT 3
// constexpr tells compiler that NT can be calculated at compile time
//constexpr int NT=3;

//typedef float ntype; // as in C
//using ntype=double; // use "using" to define an alias for the keyword "double"

template <typename ntype, int NT>
class pvector
{
  ntype v[NT]; // private member
public:
  // since assign operator returns pvector objects by reference:
  // Note this method changes data of object hence it is *not* declared "const"
  pvector& assign( const pvector& P1)
    {
      for (int i=0; i < NT; i++)
        v[i] = P1.v[i];
      return (*this);
    }

  pvector() // void constructor
    {
      // initialize the vector
      for (int i=0; i < NT; i++)
        v[i]=0;
      //std::cout << "void constructor\n";
    }

  ~pvector() // destructor
    {
      // nothing to do
      //std::cout << "void destructor\n";
    }

  // curly brace initialization, i.e.
  // pvector V={1.0,2.0,3.0};
  pvector(std::initializer_list<ntype> L)  // overloading constructor
    {
      int c=0;
      for (ntype el: L) // loop over all element in the list L
                        // hence "el" will each of list elements  
        {
          if (c < NT)
            {
              v[c] = el;
            }
          c++;
        }
      for ( ;c < NT; c++) // gli elementi sono in numero di NT allora inizializzo 0
        {
          v[c]=0.0;
        }
      //std::cout << "curly brace constructor\n";
    }

  // Note this method does not change data of object hence it is declared "const"
  void show(std::string s="") const // argomento di default 
    {
      std::cout << s << "("; 
      for (int i=0; i < NT; i++)
        {
          std::cout << v[i];
          if (i < NT-1)
            std::cout << ",";
        }
      std::cout << ")\n";
    }

  pvector& operator=(const pvector& V2)
    {
      for (int i=0; i < NT; i++)
        v[i] = V2.v[i];
      return (*this);

    }

  // method to sum two vectors
  pvector sum(const pvector &V2) const
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        {
          VT.v[i] = v[i] + V2.v[i];
        }
      return VT;
    }

  // addition operator overloading
  // Note this method does not change data of object hence it is declared "const"
  pvector operator+(const pvector& V2) const
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        {
          VT.v[i] = v[i] + V2.v[i];
        }
      return VT;
    } 

  // subtraction operator overloading
  // Note this method does not change data of object hence it is declared "const"
  pvector operator-(const pvector& V2) const
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        {
          VT.v[i] = v[i] - V2.v[i];
          // equivalently:
          // vs.v[i] = (*this).v[i] - v2.v[i];
        }
      return VT;
    } 

  // get value of a vector element
  // Note this method does not change data of object hence it is declared "const"
  ntype get(int i) const
    {
      return v[i];
    }

  // set value of vector element
  ntype set(int i, ntype val) 
    {
      return v[i]=val;
    }

  // add vector and assign
  pvector& operator+=(const pvector& V2)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] += V2.v[i];
        }
      return (*this);
    }

  // subtract vector and assign
  pvector& operator-=(const pvector& V2)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] -= V2.v[i];
        }
      return (*this);
    }

  ntype operator*(const pvector& V2) const
    {
      ntype sp=0;
      for (int i=0; i < NT; i++)
        sp += v[i]*V2.v[i];
      return sp;
    }

  // norm of vector
  ntype norm(void) const
    {
      return sqrt((*this)*(*this));
    }

  // vector times scalar 
  // A*s <=> A.operator*(s)
  pvector operator*(ntype s) const
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        // note that inside the class is possible to use private member of objects
        // belonging to the class (such as member v of VT in this case)
        VT.v[i] = v[i]*s;
      return VT;
    }

  // vector divided by scalar 
  pvector operator/(ntype s) const
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        // note that inside the class is possible to use private member of objects
        // belonging to the class (such as vt in this case)
        VT.v[i] = v[i]/s;
      return VT;
    }

  // multiply by scalar and assign result to vector
  pvector& operator*=(ntype s)
    {
      for (int i=0; i < NT; i++)
        v[i] *= s;
      return (*this);
    }

  // divide by scalar and assign result to vector
  pvector& operator /=(ntype s)
    {
      for (int i=0; i < NT; i++)
        v[i] /= s;
      return (*this);
    }

  // compare two vectors
  bool operator==(const pvector& V2) const
    {
      for (int i=0; i < NT; i++)
        {
          if (v[i] != V2.v[i])
            return 0;
        }
      return 1;
    }

  // only for 3D vectors 
  pvector operator^(const pvector& V2) const
    {
      if (NT==3)
        {
          pvector VT;
          VT.v[0] = v[1]*V2.v[2]-v[2]*V2.v[1];
          VT.v[1] = v[2]*V2.v[0]-v[0]*V2.v[2];
          VT.v[2] = v[0]*V2.v[1]-v[1]*V2.v[0];
          return VT;
        }
      else
        {
          std::cout << "Cross product not defined\n";
          exit(1);
        }
    }
  // >>>>>>>>> 28/10/2025 <<<<<<<
  // This is needed to access elements of const vector objects
  // this is like having following arguments: operator()(const pvector *this, int i) 
  // where (implicit, i.e. (*this) ) first argument is the "calling" object
#if 1 
  // V(i) = 1.0;
  // cout << V(0);
  // V         (   i          )
  // operando      operando
  ntype operator()(int idx) const
    {
      return v[idx]; 
    }
#endif
  // this is like having following arguments: operator()(pvector *this, int i) 
  // where first (implicit, i.e. (*this) ) argument is the "calling" object
  // pvector A;
  // double val;
  // std::cout << A(1);
  // A(1)=val;
  ntype& operator()(int idx)
    {
      return v[idx]; 
    }

  // scalar time vector 
  // (*this) is not implicitly passed to operator* in this case, 
  // being a friend operator of pvector class
  // pvector A;
  // double s;
  // s*A
  friend pvector operator*(ntype s, const pvector& V2)
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        {
          VT.v[i] = s*V2.v[i];
        }
      return VT;
    }

  // rint() method
  // pvector<double,3> A, B={1.1,2.2,3.3};
  // A = rint(B);
  friend pvector rint(const pvector& V2)
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        {
          VT.v[i] = rint(V2.v[i]);
        }
      return VT;
    }

  // mulcw: component wise multiplication, i.e. C=mulcw(A,B), i.e. C(i) = A(i)*B(i)
  friend pvector mulcw(const pvector& V1, const pvector& V2)
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        {
          VT(i) = V1(i)*V2(i);
        }
      return VT;
    }

  // divcw: component wise division, i.e. C=divcw(A,B), i.e. C(i) = A(i)/B(i)
  friend pvector divcw(const pvector& V1, const pvector& V2)
    {
      pvector VT;
      for (int i=0; i < NT; i++)
        {
          VT(i) = V1(i)/V2(i);
        }
      return VT;
    }

  // pvector V;
  // (std::cout  <<  V ) << "\n";
  //     os      ,     vec
  friend std::ostream& operator<<(std::ostream& os, const pvector& vec)
    {
      os << "(";
      for (int i=0; i < NT; i++)
        {
          os << vec.v[i];
          if (i < NT-1)
           os << ","; 
        }
      os << ")";
      return os;
    }
  
  // ============== RANDOM VECTORS
  // generate a random vector inside a box and assign it to (*this) object
  void random(const ntype& L)
    {
      for (int i=0; i < NT; i++)
        {
          v[i] = (rng.ranf()-0.5)*L; // assign a random value in [-L/2,L/2]
        }
    }

  // generate a random vector on unit sphere and assign it to (*this) object
  void random_orient(void)
    {
      ntype rS, S, V1, V2;
      if (NT==3)
        {
          do
            {
              V1 = 2.0*rng.ranf()-1.0;
              V2 = 2.0*rng.ranf()-1.0;
              S = V1*V1+V2*V2;
            }
          while (S >= 1.0);
          rS = sqrt(1.0-S);
#if 0
          v[0] = 2.0*rS*V1;
          v[1] = 2.0*rS*V2;
          v[2] = 1.0-2.0*S;
#else
          (*this) = {2.0*rS*V1, 2.0*rS*V2, 1.0-2.0*S};
#endif
        }
      else
        {
          std::cout << "[random_orient] Only 3D vectors are supported\n";
        }
    }
};

#endif
