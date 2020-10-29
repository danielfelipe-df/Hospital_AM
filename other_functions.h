#ifndef OTHER_FUNCTIONS_H
#define OTHER_FUNCTIONS_H


#include "bases.h"
#include "trabajadores.h"


/* Con esta función actualizo los tiempos de estado de ese grupo */
void update_times(grupo &G, trabajadores *family, double time);


/* Con esta función actualizo los tiempos de los testeados masivamente */
void update_massive(grupo &G, trabajadores *family, double time);


#endif
