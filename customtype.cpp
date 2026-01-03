#include "./customtype.hpp"
int main(void)
{
  int A=1;
  customtype T; 
  derivedtype TD;
  T.set_a(A);
  T.byref(A);
  std::cout << "A=" << A << "\n";
  TD.show();
  TD.show("Il valore di b è ", "!!!");
}

