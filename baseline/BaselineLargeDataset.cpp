//23/07/2018. Este codigo es un modelo Baseline que hace una prediccion de 0 para todas las parejas
//paciente-microbio.

//0. LIBRARIES
//+++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++
#include <iostream>
#include <fstream>
#include "armadillo"

using namespace std;
using namespace arma;
//+++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++

//1. MAIN
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char *argv[])
{
  //1.1. Read TestFile
  //=========================================================================

  //1.1.1. Choose Dataset
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  string FoldNumber="0";                       //Choose input fold
  string Dataset="5Fold_LargeData/Quartiles/"; //Choose input Dataset "5Fold_LargeData/Categories/"
  string InputPath="/export/home/shared/Projects/Microbiome/Input_Data/" + Dataset;
  string InputFile=InputPath + FoldNumber + "Test.txt";
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.1.2. Read
  //:::::::::::::::::::::::
  mat Test;
  Test.load(InputFile);
  int Testsize=Test.n_rows;
  //:::::::::::::::::::::::
  
  //=========================================================================
 
  
  //1.2. Write scores to file
  //=========================================================================
  
  //1.2.1. Choose Dataset
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  static std::vector<std::string> Outputs;
  Outputs={"Categories/", "Quartiles/"};
  string Output=Outputs[1];
  string OutputPath="/export/home/shared/Projects/Microbiome/Output_Data/LargeDataset/" + Output;
  string OutputFile;
 
  OutputFile= OutputPath +"Baseline_"+ FoldNumber + "_scores.txt";
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.2.2. Write to file
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ofstream scores;
  scores.open(OutputFile);

  rowvec score={1,0,0,0,0};
  int Links=5;
  mat Scores(Testsize,Links,fill::zeros);
  for (int index=0; index<Testsize;index++){
    scores<< Test(index,0) <<" "<<Test(index,1)<<" "<<score<<endl;
 }
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  //=========================================================================

return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
