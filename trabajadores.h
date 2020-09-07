#ifndef TRABAJADORES_H
#define TRABAJADORES_H

class trabajadores{
public:

  /* Defino qué tipo de persona es
   * 0:sus, 1:exp, 2:pre, 3:preTA, 4:lev
   * 5:levTA, 6:levA, 7:infA, 8:rec;
  */
  bool kind[9];

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
