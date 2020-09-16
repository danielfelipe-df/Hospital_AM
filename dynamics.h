#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"


/* Esta función me dice cuánto tiempo se demora en hacerse una reacción
 * y qué reacción es.
 */
std::vector<double> contagio(double Sa, double Sb, double Ea, double Eb, double Pa, double Pb, double PTa, double PTb, double PTAa, double PTAb, double La, double Lb, double LTa, double LTb, double LTAa, double LTAb, double LAa, double LAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran, double t);


/* Esta función es la función plantilla para hacer las reacciones */
void mother_reaction(grupo &Out, grupo &In, Crandom &ran, trabajadores *family, int typeout, int typein);


/* Esta función es la que genera la reacción en el testeo masivo */
void massive_reaction(grupo &S, grupo &E, grupo &P, grupo &PT, grupo &L, grupo &LT, grupo &R, Crandom &ran, trabajadores *family);


/* Esta función es la que genera la reacción en el testeo continuo */
void continue_reaction(grupo &L, grupo &LT, trabajadores *family, Crandom &ran);


/* Esta función implementa el método de bisección para hallar el tiempo */
double biseccion(double A, double prom, double sigma, double t, double B, double ranr);


/* Es la función integrada de la prevalencia externa */
double function(double A, double prom, double sigma, double t, double tau);


/* Esta función me hace la reacción de contagio para un susceptible alto */
void reaction0(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de contagio para un susceptible bajo */
void reaction1(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un expuesto alto */
void reaction2(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un expuesto bajo */
void reaction3(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve alto */
void reaction4(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve bajo */
void reaction5(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado alto */
void reaction6(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado bajo */
void reaction7(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a recuperado alto*/
void reaction8(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de tránsito de un pre-sintomático a recuperado bajo*/
void reaction9(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de leve a recuperado alto */
void reaction10(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de leve a recuperado bajo */
void reaction11(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de infeccioso grave a recuperado alto */
void reaction12(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


/* Esta función me hace la reacción de infeccioso grave a recuperado bajo */
void reaction13(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);


#endif
