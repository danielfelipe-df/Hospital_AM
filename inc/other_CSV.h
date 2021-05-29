/**
 * @file other_CSV.h
 * @author Daniel Felipe
 * @date 2020
 * @brief Header containing the functions to read files
 */


#ifndef OTHER_CSV_H
#define OTHER_CSV_H

#include <vector>
#include <string>
#include <algorithm>
#include <functional>

//Lee el archivo CSV y los retorna en una matriz de std::wstring
//Además, retorna todos los caracteres en mayúscula.
void csv_ws(std::vector<std::vector<std::wstring> > &w, std::string archivo, wint_t wc = L'\t');

//Lee el archivo CSV y los retorna en una matriz de long long
void csv_ll(std::vector<std::vector<long long> > &w, std::string archivo, wint_t wc = L'\t');

//Lee el archivo CSV y los retorna en una matriz de double
void csv_d(std::vector<std::vector<double> > &w, std::string archivo, wint_t wc = L'\t');


//Lee el archivo CSV y los retorna en una matriz de int
void csv_i(std::vector<std::vector<int> > &w, std::string archivo, wint_t wc = L'\t');


//Lee el archivo CSV y los retorna en una matriz de float
void csv_f(std::vector<std::vector<float> > &w, std::string archivo, wint_t wc = L'\t');


//Esta función quita los espacios de más en el std::wstring al inicio y al final
void elimina_espacios(std::wstring & ws);


//Lee el archivo CSV y los retorna en una matriz de std::wstring. Pero utlizando la función elimina espacios.
//Además, retorna todos los caracteres en mayúscula.
void csv_ws_without_spaces(std::vector<std::vector<std::wstring> > &w, std::string archivo, wint_t wc = L'\t');

//Esta función me da los archivos que hay en la carpeta name
void read_directory(const std::string& name, std::vector<std::string> &v);


//Esta función me elimina los elementos repetidos en el vector
template <typename T>
void eliminar_repetidos(std::vector<T> &y)
{
  auto end = y.end();
  for(auto it = y.begin(); it != end; it++){
    end = std::remove(it + 1, end, *it);
  }
  y.erase(end, y.end());
}


//Esta función ordena los elementos del vector de menor a mayor
template <typename T>
void menor_a_mayor(std::vector<T> &y)
{
  std::sort(y.begin(), y.end());
}


//Esta función ordena los elementos del vector de mayor a menor
template <typename T>
void mayor_a_menor(std::vector<T> &y)
{
  std::sort(y.begin(), y.end(), std::greater<T>());
}


#endif /* OTHER_CSV_H */
