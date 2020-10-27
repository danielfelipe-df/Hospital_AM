#include "trabajadores.h"


void trabajadores::init(){
  kind[0] = true;
  for(unsigned int i=1; i<11; i++){kind[i] = false;}
  time = 0;
  tmax = 0;
  tstate = 0;
  my_inf.clear();
  DF = -2;
}


void trabajadores::change(int now, int past){
  kind[past] = false;  kind[now] = true;
}
