#include <iostream>
#include <fstream>
using namespace std;

int main() {
  ifstream fichero;
  std::ifstream file("/export/home/shared/Projects/Microbiome/Input_Data/micro_dataset.txt", ios::in);
  std::string linea;            //variable para almacenar linea
  std::string contenido_archivo;//variable para almacenar contenido del fichero
  
  while (std::getline(file, linea))
    {
      contenido_archivo += linea;
      contenido_archivo.push_back('\n');
      cout << linea;
      cout << '\n';
     }
  fichero.close();
  return 0;
    }
