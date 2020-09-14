#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <vector>
#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"


/* Esta función me dice cuánto tiempo se demora en hacerse una reacción
 * y qué reacción es.
 */
std::vector<double> contagio(double Sa, double Sb, double Ea, double Eb, double Pa, double Pb, double PTa, double PTb, double La, double Lb, double LTa, double LTb, double Ia, double Ib, double ITa, double ITb, double Na, double Nb, double prev, Crandom &ran, double t);


/* Esta función es la función plantilla para hacer las reacciones */
void mother_reaction(grupo &Out, grupo &In, Crandom &ran, trabajadores *family, int typeout, int typein);


/* Esta función es la que genera la reacción en el testeo masivo */
void massive_reaction(grupo &S, grupo &E, grupo &P, grupo &PTA, grupo &L, grupo &LTA, grupo &R, Crandom &ran, trabajadores *family);


/* Esta función implementa el método de bisección para hallar el tiempo */
double biseccion(double A, double prom, double sigma, double t, double B, double ranr);


/* Es la función integrada de la prevalencia externa */
double function(double A, double prom, double sigma, double t, double tau);


/* Esta función me hace la reacción de contagio para un susceptible alto */
inline void reaction0(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Sa, Ea, ran, altos, 0, 1);
}


/* Esta función me hace la reacción de contagio para un susceptible bajo */
inline void reaction1(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Sb, Eb, ran, bajos, 0, 1);
}


/* Esta función me hace la reacción de tránsito de un expuesto alto */
inline void reaction2(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Ea, Pa, ran, altos, 1, 2);
}


/* Esta función me hace la reacción de tránsito de un expuesto bajo */
inline void reaction3(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Eb, Pb, ran, bajos, 1, 2);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve alto */
inline void reaction4(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pa, La, ran, altos, 2, 4);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a leve bajo */
inline void reaction5(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pb, Lb, ran, bajos, 2, 4);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado alto */
inline void reaction6(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pa, IAa, ran, altos, 2, 7);
}


/* Esta función me hace la reacción de tránsito de un pre-sintomático a infeccioso grave aislado bajo */
inline void reaction7(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pb, IAb, ran, bajos, 2, 7);
}


/* Esta función me hace la reacción de tránsito de un leve a aislado leve alto */
inline void reaction8(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(La, LAa, ran, altos, 4, 6);
}


/* Esta función me hace la reacción de tránsito de un leve a aislado leve bajo */
inline void reaction9(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Lb, LAb, ran, bajos, 4, 6);
}


/* Esta función me hace la reacción de pre-sintomático a recuperado alto */
inline void reaction10(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pa, Ra, ran, altos, 2, 8);
}


/* Esta función me hace la reacción de pre-sintomático a recuperado bajo */
inline void reaction11(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pb, Rb, ran, bajos, 2, 8);
}


/* Esta función me hace la reacción de leve a recuperado alto */
inline void reaction12(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(La, Ra, ran, altos, 4, 8);
}


/* Esta función me hace la reacción de leve a recuperado bajo */
inline void reaction13(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Lb, Rb, ran, bajos, 4, 8);
}


/* Esta función me hace la reacción de pre-sintomático testeado y aislado a recuperado alto */
inline void reaction14(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(PTAa, Ra, ran, altos, 3, 8);
}


/* Esta función me hace la reacción de pre-sintomático testeado y aislado a recuperado bajo */
inline void reaction15(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(PTAb, Rb, ran, bajos, 3, 8);
}


/* Esta función me hace la reacción de leve aislado a recuperado alto */
inline void reaction16(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(LAa, Ra, ran, altos, 6, 8);
}


/* Esta función me hace la reacción de leve aislado a recuperado bajo */
inline void reaction17(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(LAb, Rb, ran, bajos, 6, 8);
}


/* Esta función me hace la reacción de leve aislado y testeado a recuperado alto */
inline void reaction18(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(LTAa, Ra, ran, altos, 5, 8);
}


/* Esta función me hace la reacción de leve aislado y testeado a recuperado bajo */
inline void reaction19(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(LTAb, Rb, ran, bajos, 5, 8);
}


/* Esta función me hace la reacción de tránsito de infeccioso grave a recuperado alto */
inline void reaction20(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(IAa, Ra, ran, altos, 7, 8);
}


/* Esta función me hace la reacción de tránsito de infeccioso grave a recuperado bajo*/
inline void reaction21(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(IAb, Rb, ran, bajos, 7, 8);
}



#endif
