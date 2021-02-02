#include <iostream>
#include <algorithm>
#include <other_CSV.h>

int main(void)
{
  std::vector<std::string> names;
  read_directory("Data/", names);
  std::sort(names.begin(), names.end());
  names.erase(names.begin() + 0);
  names.erase(names.begin() + 0);
  names.erase(names.begin() + 0);
  
  std::vector<std::vector<int> > data;
  std::vector<int> my_data, mycopy;
  int value, index;
  std::vector<int>::iterator it;

  for(size_t j = 0; j<names.size(); j++){
    csv_i(data, "Data/" + names[j]);

    my_data.resize(data.size()), mycopy.resize(data.size());
    for(size_t i=0; i<data.size(); i++){my_data[i] = data[i][0];    mycopy[i] = data[i][0];}
    
    std::sort(mycopy.begin(), mycopy.end());

    value = mycopy[(int)(mycopy.size()*0.5)];
    it = std::find(my_data.begin(), my_data.end(), value);
    index = std::distance(my_data.begin(), it);

    std::cout << names[j] << '\t' << value << '\t' << index << std::endl;

    data.clear();
    my_data.clear();
    mycopy.clear();
  }
    
  return 0;
}
