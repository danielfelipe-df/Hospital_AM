/**
 * @file workers.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing Workers class
 */


#ifndef WORKERS_H
#define WORKERS_H

#include <vector>

class Workers{
public:

  /* Defino qué tipo de persona es
   * 0:sus, 1:susT, 2:susA
   * 3:exp, 4:expT, 5:expA
   * 6:pre, 7:preT, 8:preTA
   * 9:lev, 10:levT, 11:levTA,
   * 12:infA,
   * 13:recI, 14:recT, 15:recA
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

  //Defino la variable que me dice si está en rastreo
  bool btrace;

  //Defino la variable que me dice cuánto tiempo lleva rastreado
  double ttrace;

  //Guardo a las personas que infecto si soy infeccioso
  std::vector<int> my_inf;

  //Guardo la identidad de la persona que me infectó
  int DF;

  //Función de iniciación
  void init();

  //Función de campio de estado
  void change(unsigned int now, unsigned int past);
};

#endif /* WORKERS_H */
