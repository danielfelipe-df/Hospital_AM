#ifndef OTHER_FUNCTIONS_H
#define OTHER_FUNCTIONS_H

#include <string>
#include "bases.h"
#include "trabajadores.h"


/* Con esta función actualizo los tiempos de estado de ese grupo */
void update_times(grupo &G, trabajadores *family, double time);


/* Con esta función actualizo los tiempos de los testeados masivamente */
void update_massive(grupo &G, trabajadores *family, double time);


/* Con esta función actualizo los tiempos de estado de todos los grupos */
void update_times_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time);


/* Con esta función actualizo los tiempos de todos los testeados masivamente */
void update_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time);


/* Con esta función imprimo cada uno de los vectores */
void print_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, double time, std::string name);


/* Con esta función imprimo la cantidad que hayan por cajita */
void print_types(std::vector<grupo> &Val, std::vector<grupo> &Vba, double time, std::string name);


/* Con esta función imprimo todos los infectados (total) divididos por sus tipo (alto y bajo) y el total (alto + bajo)*/
void print_inf(std::vector<grupo> &Val, std::vector<grupo> &Vba, double time, std::string name);

#endif
