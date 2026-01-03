#include <iostream>
#include <cmath> // C math function pow, log
#include <string>
class customtype
{
//public:
  int a; // private members
  void prfunc(int i)
    { 
      a = a+i;
    }
protected:
  int b;
public: 
  int c; // public member

  void set_a(int i)
    {
      prfunc(i);
      a = i;
    }
  void set_b(int i)
    {
      prfunc(i);
      b = i;
    }
   void byref(int& a)
     {
       a=2;
     } 
   int get_a(void)
    {
      return a;
    }
   int get_b(void)
    {
      return b;
    }
  
  customtype()//constructor
    {
      a = 1;
      b = 2;
      c = 3;
      std::cout << "Constructor\n";
    }
  ~customtype()//destructor
    {
      std::cout << "Destructor\n";
    }
};

class derivedtype: public customtype
{
//
public:
 void show(void)
   {
     std::cout << "b=" << b << "\n";
   }
 void show(std::string txt, std::string txt2="!") const
   {
     //(*this).b = 3; 
     std::cout << txt << (*this).b << txt2 << "\n";
   }

};

