#ifndef TRABAJADORES_H
#define TRABAJADORES_H

#include <vector>

class trabajadores{
public:

  /* Defino qué tipo de persona es
   * 0:sus, 1:susT,
   * 2:exp, 3:expT, 4:expA
   * 5:pre, 6:preT, 7:preA, 8:preM
   * 9:lev, 10:levT, 11:levA, 12:levM
   * 13:infA,
   * 14:recI, 15:recT, 16:recA
  */
  unsigned int kind;

  //Defino la variable que me contabiliza el tiempo que dura testeado
  double time;

  //Defino la variable que me dice el tiempo máximo que dura testeado antes de aislarse
  double tmax;

  //Defino la variable que me dice el tiempo que permaence en el estado actual
  double tstate;

  //Guardo a las personas que infecto si soy infeccioso
  std::vector<int> my_inf;

  //Guardo la identidad de la persona que me infectó
  int DF;

  //Función de iniciación
  void init();

  //Función de campio de estado
  void change(unsigned int now, unsigned int past);
};

#endif
