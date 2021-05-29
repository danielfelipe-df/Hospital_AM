/**
 * @file trace.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing the tracing functions and reactions
 */


#ifndef TRACE_H
#define TRACE_H

#include <Random64.h>
#include <bases.h>
#include <workers.h>


/* Esta es la función para hacer la reacción del rastreo */
void trace_reaction(grupo &Out, grupo &In, int index, Workers *family, int typeout, int typein, double time, bool is_traced);


/* Esta es la función para hacer rastreo */
void main_trace(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos, double time, Crandom &ran);


/* Esta es una función auxiliar para el 'main_trace' */
void aux_main(int num, std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos, Crandom &ran, double cons1, double cons2, bool type, unsigned int index, bool pre);


/* Esta es la función para hacer la reaccion de rastreo */
int reaction_trace(std::vector<grupo> &V, Workers *family, Crandom &ran, int index);


/* Esta es una función auxiliar para hallar a la persona en 'reaction_trace' */
void aux_trace(grupo &G, grupo &T, Workers *family, int typeout, int typein, int index, bool normal);


/* Eliminar repetidos en el vector */
void eliminar_repetidos(std::vector<int> &y);


/* Esta función mira cuáles rastreados ya cumplieron el tiempo, y los devuelve */
void trace_massive(grupo &R, grupo &G, Workers *family, double time, int typeout, int typein);


/* Esta función cambia a todos los rastreados que ya cumplieron el tiempo */
void trace_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos, double time);


#endif /* TRACE_H */
