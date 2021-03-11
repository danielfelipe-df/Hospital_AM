#include <iostream>
#include <stdlib.h>
#include <trabajadores.h>

void trabajadores::init(){
  kind = 0;
  time = 0;
  tmax = 0;
  tstate = 0;
  tlev = 0;
  tlevmax = 0;
  myTt = 0;
  btrace = false;
  ttrace = 0;
  my_inf.clear();
  DF = -2;
}


void trabajadores::change(unsigned int now, unsigned int past){
  if(kind == past){kind = now;}
  else{
    std::cerr << "Error: Tipo de trabajador no coincide.\nSe debio ingresar " << kind << " pero se ingresa " << past << std::endl;
    std::exit(EXIT_FAILURE);}
}
