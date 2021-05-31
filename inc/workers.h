/**
 * @file workers.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing Workers class and Stages enum.
 */


#ifndef WORKERS_H
#define WORKERS_H


// Enumeration of the stages in the model
enum Stages
  {
   SUS = 0, // Susceptible
   SUST = 1, // Susceptible Tested
   SUSA = 2, // Susceptible Isolated
   EXP = 3, // Exposed
   EXPT = 4, // Exposed Tested
   EXPA = 5, // Exposed Isolated
   PRE = 6, // Pre-symptomatic
   PRET = 7, // Pre-symptomatic Tested
   PREA = 8, // Pre-symptomatic Isolated
   MSYM = 9, // Mild-symptomatic
   MSYMT = 10, // Mild-symptomatic Tested
   MSYMA = 11, // Mild-symptomatic Isolated
   SSYMA = 12, // Severe-symptomatic Isolated
   RECI = 13, // Recovered No-detected
   RECT = 14, // Recovered Tested
   RECA = 15 // Recovered Detected
};


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
  Stages kind;

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
  void change(Stages now, Stages past);
};

#endif /* WORKERS_H */
