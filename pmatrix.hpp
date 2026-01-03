// X = done
// * = to discuss
// E = exercise for students
// N = don't do
#ifndef _PMATRIX_
#define _PMATRIX_ 
// TODO implement following methods in this lab #5 (04/11)
// pmatrix<double,3> M;
// M = {1.,2.,3.,4.,5.,6.,7.,8.,9.};
// 
// [X] constructor with initializer_list
// [X] overloading of () operator (m(i,j))
// [X] overloading of << operator to print a matrix
// [X] matrix addition and subtraction
// 
// =================================================
//
// during next lab #6 on 11/11
// [X] += and -= operators
// [X] matrix times scalar (also *=)
// [X] matrix divided by scalar (also /=) 
// [X] scalar times matrix <--- with friend function
// [X] matrix times vector 
// [X] vector times matrix <--- with friend function
// [X] matrix times matrix (also *=)
// [X] transpose
// [X] operator == to pmatrix class
// 
#include "./pvector.hpp"
template <typename ntype, int NT>
class pmatrixq
{
  ntype m[NT][NT];
public:
  ~pmatrixq()
    {

    }
  pmatrixq()
    {
      for (int i=0; i < NT; i++)
        for (int j=0; j < NT; j++)
          m[i][j]=0;
    }

  // to initialize a matrix as follows: pmatrix<double,2>={1,2,3,4};
  pmatrixq(std::initializer_list<ntype> list)
    {
      int cc=0;
      for (auto el: list)
        {
          m[cc/NT][cc%NT] = el;
          cc++;  
        }
    } 

  // to initialize a matrix as follows: pmatrix<double,2>={{1,2},{3,4}};
  pmatrixq(std::initializer_list<std::initializer_list<ntype>> list)
    {
      int i=0, j=0;
      for (auto r: list)
        {
          j=0;
          for (auto c: r)
            {
              if (i < NT && j < NT)
                m[i][j] = c;
              j++;
            }
          i++;
        }
    } 

  // ============================================================
  // overloading of () operator to get/set (i,j)-th element
  // m(i,j)=2.0;
  // std::cout << m(i, j);
  ntype& operator()(int i, int j)
    {
      return m[i][j];
     }
  
  ntype operator()(int i, int j) const
    {
      return m[i][j];
    }

  // output (std:cout << M) << "\n";
  friend std::ostream& operator<<(std::ostream& os, const pmatrixq& M1)
    {
      int i, j;
      os << "{";
      for (i=0; i < NT; i++)
        {
          os << "{";
          for (j=0; j < NT; j++)
            {
              os << M1(i,j);
              if (j < NT-1)
                os << ",";
            }
          os << "}";
          if (i < NT-1)
            os << ",";
        }
      os << "}";
      return os;
    }

  // addition 
  pmatrixq operator+(const pmatrixq &M2) const
    {
      // MT = M1 + M2
      pmatrixq MT;
      for (auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              MT(i,j) = (*this)(i,j) + M2(i,j);
            }
        }
      return MT;
    }
  // subtraction
  pmatrixq operator-(const pmatrixq &M2) const
    {
      // MT = M1 - M2
      pmatrixq MT;
      for (auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              MT(i,j) = (*this)(i,j) - M2(i,j);
            }
        }
      return MT;
    }

  pmatrixq& operator+=(const pmatrixq &M2)
    {
      return (*this = *this + M2);
     }
  
  pmatrixq<ntype,NT>& operator-=(const pmatrixq &M2)
    {
      return (*this = *this - M2);
    }

  // matrix /= scalar
  pmatrixq& operator/=(const ntype& s) 
    {
      for(auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              m[i][j] /= s;
            }
        }
      return (*this);
    } 
  
  // matrix *= scalar
  pmatrixq& operator*=(const ntype& s) 
    {
      for(auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              m[i][j] *= s;
            }
        }
      return (*this);
    } 

  // matrix times scalar
  pmatrixq operator*(const ntype& s) const
    {
      // write your code here
      pmatrixq MT;
      for(auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              MT(i,j) = m[i][j]*s;
            }
        }
      return MT;
    }   
 
  // matrix divided by scalar
  pmatrixq operator/(const ntype& s) 
    {
      pmatrixq MT;
      for(auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              MT(i,j) = m[i][j]/s;
            }
        }
      return MT;
    }
 
  // scalar times matrix 
  friend pmatrixq operator*(const ntype& s, const pmatrixq &M2)
    {
      pmatrixq MT;
      for(auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              MT(i,j) = s*M2(i,j);
            }
        }
      return MT;
    }   
 
  // matrix times vector
  pvector<ntype,NT> operator*(const pvector<ntype,NT> &V2) const 
    {
      pvector<ntype,NT> VT;
      for (int i=0; i < NT; i++)
        {
          VT(i) = 0.0;
           for (int j=0; j < NT; j++)
            VT(i) += m[i][j]*V2(j); 
         }
      return VT;
    }
 
  // vector times matrix
  friend pvector<ntype,NT> operator*(const pvector<ntype,NT>& V1, const pmatrixq& M2)
    {
      pvector<ntype,NT> VT;
      for (int j=0; j < NT; j++)
        {
          VT(j) = 0.0;
          for (int i=0; i < NT; i++)
            VT(j) += V1(i)*M2(i,j); 
        }
      return VT;
    }
 
  // matrix times matrix
  // m1 è *this
  // m1_{il}*m2_{lj} 
  pmatrixq operator*(const pmatrixq &M2) const
    {
      pmatrixq MT;
      for (int i=0; i < NT; i++)
        for (int j=0; j < NT; j++)
          { 
            MT(i,j) = 0.0;
            for (int k=0; k < NT; k++)
              MT(i,j)+=m[i][k]*M2(k,j);
          }
      return MT;
    }
  
  pmatrixq transpose(void)
    {
      pmatrixq MT;
      for (int i=0; i < NT; i++)
        for (int j=0; j < NT; j++)
          MT(i,j) = m[j][i];
      return MT;
    }

  bool operator==(const pmatrixq& M2) const
    {
      for (auto i=0; i < NT; i++)
        {
          for (auto j=0; j < NT; j++)
            {
              if (M2(i,j)!=m[i][j])
                {
                  return 0;
                }
            }
        }
      return 1;
    }

  // =================================== up here on 11/11/25
};
#endif
