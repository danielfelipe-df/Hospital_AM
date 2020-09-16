#ifndef TRABAJADORES_H
#define TRABAJADORES_H

class trabajadores{
public:

  /* Defino qué tipo de persona es
   * 0:sus, 1:exp, 2:pre, 3:preT, 4:preTA,
   * 5:lev, 6:levT, 7:levTA, 8:levA,
   * 9:infA, 10:rec;
  */
  bool kind[11];

  //Defino la variable que me contabiliza el tiempo en el estado
  double time;

  //Defino la variable que me dice el tiempo máximo que pasa en el estado
  double tmax;

  //Función de iniciación
  void init();

  //Función de campio de estado
  void change(int now, int past);
};

#endif
