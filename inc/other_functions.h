/**
 * @file other_functions.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing the functions to update times and print prevalences
 */


#ifndef OTHER_FUNCTIONS_H
#define OTHER_FUNCTIONS_H

#include <map>
#include <bases.h>
#include <workers.h>


/* Con esta función actualizo los tiempos de estado de ese grupo */
void update_times(group &G, Workers *family, double time);


/* Con esta función actualizo los tiempos de los testeados masivamente */
void update_massive(group &G, Workers *family, double time);


/* Con esta función actualizo los tiempos de estado de todos los grupos */
void update_times_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, double time);


/* Con esta función actualizo los tiempos de todos los testeados masivamente */
void update_massive_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, double time);


/* Con esta función imprimo cada uno de los vectores */
void print_all(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, double time, std::string name);


/* Con esta función imprimo la cantidad que hayan por cajita */
void print_types(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, double time, std::string name);


/* Con esta función imprimo la cantidad de infecciosos por tipo (alto y bajo) y el total */
void print_inf(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, double time, std::string name);


/* Con esta función imprimo los vínculos de la red de contagio */
void print_net(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Workers *altos, Workers *bajos, std::string name);


#endif /* OTHER_FUNCTIONS_H */
