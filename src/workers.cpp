#include <iostream>
#include <stdlib.h>
#include <workers.h>

void Workers::init(){
  kind = SUS;
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


void Workers::change(Stages now, Stages past){
  if(this->kind == past){this->kind = now;}
  else{
    std::cerr << "Error: Tipo de trabajador no coincide.\nSe debio ingresar " << this->kind << " pero se ingresa " << past << std::endl;
    std::exit(EXIT_FAILURE);}
}
