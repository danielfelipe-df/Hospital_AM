#ifndef TRABAJADORES_H
#define TRABAJADORES_H

#include <vector>

class trabajadores{
public:

  /* Defino qué tipo de persona es
   * 0:sus, 1:susT,
   * 2:exp, 3:expT, 4:expA
   * 5:pre, 6:preT, 7:preTA
   * 8:lev, 9:levT, 10:levTA,
   * 11:infA,
   * 12:recI, 13:recT, 14:recA
  */
  unsigned int kind;

  //Defino la variable que me contabiliza el tiempo que dura testeado
  double time;

  //Defino la variable que me dice el tiempo máximo que dura testeado antes de aislarse
  double tmax;

  //Defino la variable que me dice el tiempo que permaence en el estado actual
  double tstate;

  //Defino la variable que me dice el tiempo que el leve iota lleva desde que se testeó
  double tlev;

  //Defino la variable que me dice cuánto tiempo tiene que esperar el leve iota para saber el resultado del test
  double tlevmax;

  //Defino el tiempo que se demoran en entregarme la prueba
  double myTt;

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
