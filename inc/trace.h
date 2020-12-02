#ifndef TRACE_H
#define TRACE_H

#include <Random64.h>
#include <bases.h>
#include <trabajadores.h>


/* Esta es la función para hacer rastreo */
void main_trace(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran);


/* Esta es una función auxiliar para el 'main_trace' */
void aux_main(int num, std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, Crandom &ran, double cons1, double cons2, bool type, unsigned int index, bool pre);


/* Esta es la función para hacer la reaccion de rastreo */
int reaction_trace(std::vector<grupo> &V, trabajadores *family, Crandom &ran, int index);


/* Esta es una función auxiliar para hallar a la persona en 'reaction_trace' */
void aux_trace(grupo &G, grupo &T, trabajadores *family, int typeout, int typein, int index, bool normal);


/* Eliminar repetidos en el vector */
void eliminar_repetidos(std::vector<int> &y);

#endif
