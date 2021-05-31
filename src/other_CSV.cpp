#include <fstream>
#include <sys/types.h>
#include <dirent.h>
#include <locale.h>
#include <other_CSV.h>

void csv_ws(std::vector<std::vector<std::wstring> > &w, std::string archivo, wint_t wc){
  char *locale = setlocale(LC_ALL, "");
  wint_t c;

  FILE * pFile2 = fopen(archivo.c_str(), "r");
  std::wstring ws = L"";

  int k=0;
  w.push_back( std::vector<std::wstring>() );
  while ((c = fgetwc(pFile2)) != WEOF){
    if(c == L'\n'){
      w[k].push_back(ws);
      w.push_back( std::vector<std::wstring>() );
      k++;
      ws.clear();
    }
    else if(c == wc){
      w[k].push_back(ws);
      ws.clear();
    }
    else{
      ws.push_back(towupper(c));
    }
  }
  fclose(pFile2);
  w.pop_back();
  *(locale) += 'a';
}


void csv_ll(std::vector<std::vector<long long> > &w, std::string archivo, wint_t wc){
  char *locale = setlocale(LC_ALL, "");
  wint_t c;

  FILE * pFile2 = fopen(archivo.c_str(), "r");
  std::wstring ws = L"";

  int k=0;
  w.push_back( std::vector<long long>() );
  while ((c = fgetwc(pFile2)) != WEOF){
    if(c == L'\n'){
      if(ws.length() != 0){w[k].push_back(std::stoll(ws));}
      w.push_back( std::vector<long long>() );
      k++;
      ws.clear();
    }
    else if(c == wc){
      if(ws.length() != 0){w[k].push_back(std::stoll(ws));}
      ws.clear();
    }
    else{
      ws.push_back(c);
    }
  }
  fclose(pFile2);
  w.pop_back();
  *(locale) += 'a';
}

void csv_d(std::vector<std::vector<double> > &w, std::string archivo, wint_t wc){
  char *locale = setlocale(LC_ALL, "");
  wint_t c;

  FILE * pFile2 = fopen(archivo.c_str(), "r");
  std::wstring ws = L"";

  int k=0;
  w.push_back( std::vector<double>() );
  while ((c = fgetwc(pFile2)) != WEOF){
    if(c == L'\n'){
      if(ws.length() != 0){w[k].push_back(std::stod(ws));}
      w.push_back( std::vector<double>() );
      k++;
      ws.clear();
    }
    else if(c == wc){
      if(ws.length() != 0){w[k].push_back(std::stod(ws));}
      ws.clear();
    }
    else{
      ws.push_back(c);
    }
  }
  fclose(pFile2);
  w.pop_back();
  *(locale) += 'a';
}


void csv_i(std::vector<std::vector<int> > &w, std::string archivo, wint_t wc){
  char *locale = setlocale(LC_ALL, "");
  wint_t c;

  FILE * pFile2 = fopen(archivo.c_str(), "r");
  std::wstring ws = L"";

  int k=0;
  w.push_back( std::vector<int>() );
  while ((c = fgetwc(pFile2)) != WEOF){
    if(c == L'\n'){
      if(ws.length() != 0){w[k].push_back(std::stoi(ws));}
      w.push_back( std::vector<int>() );
      k++;
      ws.clear();
    }
    else if(c == wc){
      if(ws.length() != 0){w[k].push_back(std::stoi(ws));}
      ws.clear();
    }
    else{
      ws.push_back(c);
    }
  }
  fclose(pFile2);
  w.pop_back();
  *(locale) += 'a';
}


void csv_f(std::vector<std::vector<float> > &w, std::string archivo, wint_t wc){
  char *locale = setlocale(LC_ALL, "");
  wint_t c;

  FILE * pFile2 = fopen(archivo.c_str(), "r");
  std::wstring ws = L"";

  int k=0;
  w.push_back( std::vector<float>() );
  while ((c = fgetwc(pFile2)) != WEOF){
    if(c == L'\n'){
      if(ws.length() != 0){w[k].push_back(std::stof(ws));}
      w.push_back( std::vector<float>() );
      k++;
      ws.clear();
    }
    else if(c == wc){
      if(ws.length() != 0){w[k].push_back(std::stof(ws));}
      ws.clear();
    }
    else{
      ws.push_back(c);
    }
  }
  fclose(pFile2);
  w.pop_back();
  *(locale) += 'a';
}


void elimina_espacios(std::wstring & ws){
  while(iswspace(ws.at(0))){ws.erase(ws.begin() + 0);}
  while(iswspace(ws.at(ws.length() - 1))){ws.erase(ws.end() - 1);}
}


void csv_ws_without_spaces(std::vector<std::vector<std::wstring> > &w, std::string archivo, wint_t wc){
  char *locale = setlocale(LC_ALL, "");
  wint_t c;

  FILE * pFile2 = fopen(archivo.c_str(), "r");
  std::wstring ws = L"";

  int k=0;
  w.push_back( std::vector<std::wstring>() );
  while ((c = fgetwc(pFile2)) != WEOF){
    if(c == L'\n'){
      if(ws.length() != 0){elimina_espacios(ws);}
      w[k].push_back(ws);
      w.push_back( std::vector<std::wstring>() );
      k++;
      ws.clear();
    }
    else if(c == wc){
      if(ws.length() != 0){elimina_espacios(ws);}
      w[k].push_back(ws);
      ws.clear();
    }
    else{
      ws.push_back(towupper(c));
    }
  }
  fclose(pFile2);
  w.pop_back();
  *(locale) += 'a';
}


void read_directory(const std::string &name, std::vector<std::string> &v){
  DIR* dirp = opendir(name.c_str());
  struct dirent * dp;
  while((dp = readdir(dirp)) != NULL){
    v.push_back(dp->d_name);
  }
  closedir(dirp);
}
