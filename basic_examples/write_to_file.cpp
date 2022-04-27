#include <iostream>
#include <fstream>
using namespace std;

int main() {
  ofstream myfile;
  myfile.open("/export/home/shared/Projects/Microbiome/Output_Data/test.txt");
  myfile << "Prueba para escribir a archivo.\n";
    myfile.close();
  return 0;
}
