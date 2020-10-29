#ifndef TRACE_H
#define TRACE_H

#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"


/* Esta es la función para hacer rastreo */
void main_trace(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran);


/* Esta es una función auxiliar para el 'main_trace' */
void aux_main(int num, grupo &G, grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &RIa, grupo &RIb, grupo &RTa, grupo &RTb, trabajadores *altos, trabajadores *bajos, Crandom &ran, double cons1, double cons2, bool type);


/* Esta es la función para hacer la reaccion de rastreo */
int reaction_trace(int index, grupo &S, grupo &ST, grupo &E, grupo &ET, grupo &P, grupo &PT, grupo &L, grupo &LT, grupo &RI, grupo &RT, trabajadores *family);


/* Esta es una función auxiliar para hallar a la persona en 'reaction_trace' */
int aux_trace(grupo &G, grupo &T, trabajadores *family, int typeout, int typein, int index);


#endif
