//28/04/2022. This code finds latent groups of hosts and microbes

//K: Groups of hosts. L: Groups of microbes


#include <iostream>
#include <fstream>
#include "armadillo"

using namespace std;
using namespace arma;

int People;
int MBiome;

//0. VARIABLES                                                      
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct MatricesMemberships{
  mat theta;    //Person Memberships 
  mat eta;      //Microbe Memberships
  cube P;       //P Matrices
};


float nan_saver=1e-10;//This number prevents singularities in numerical calculations
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//1. FUNCTIONS
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//1.1. Get Initial \eta and \theta and ps    
//=========================================================================================
MatricesMemberships GetMatrices\
(int People,int Groups_People, int Microbes, int Groups_Microbes, int Number_of_Links) {
  MatricesMemberships Init_Outputs;
  //Generate Thetas matrix
  //:::::::::::::::::::::::::::::::::::::::::::::                           
  Init_Outputs.theta=randu(People, Groups_People);    	  //Generate at random 
  Init_Outputs.theta=normalise(Init_Outputs.theta, 1, 1 );//Normalise
  //Generate Etas matrix
  //:::::::::::::::::::::::::::::::::::::::::::::
  Init_Outputs.eta  =randu(Microbes, Groups_Microbes);   //Generate at random
  Init_Outputs.eta  =normalise( Init_Outputs.eta, 1, 1 );//Normalise
  //Generate p matrices
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  Init_Outputs.P=zeros(Groups_People,Groups_Microbes,Number_of_Links);
   for (int i=0; i<Groups_People; i++){
   //Create a matrix of dimensions 'number of links x number of group of microbes
    mat M=randu(Number_of_Links,Groups_Microbes);    
    M=normalise(M,1,0);//Normalise columns
    
    for (int k=0; k<Number_of_Links; k++){
      //For each tensor slice, fill each of the k rows with all p(k,:) probabilities.
      //Iterate over rows
      Init_Outputs.P.slice(k).row(i)=M.row(k);	
      }}    
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  return Init_Outputs;
}
//=========================================================================================


//1.2. Find Omega matrix
//=========================================================================================
void GetOmega\
(int People,int Groups_People, int Microbes, int Groups_Microbes, int Number_of_Links,\
 MatricesMemberships InitMatrices, mat Data, cube ListAbund, vec Tope, field<mat>* Omega) {

  //1.2.0. Compute omega numerators
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for (int k=0; k<Groups_People; k++){
    for (int l=0; l<Groups_Microbes; l++){
      for (int abund=0; abund<Number_of_Links; abund++){	  
  	for (int index=0; index<Tope(abund);index++){
  	  (*Omega)(k,l)( ListAbund(index,1,abund) , ListAbund(index,0,abund) ) = \
  	    InitMatrices.eta( ListAbund(index,0,abund) ,l)*InitMatrices.theta( ListAbund(index,1,abund) ,k) \
  	    *InitMatrices.P(k,l,abund);
	}}}}
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.2.1. Compute normalizing constants
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  mat Normal(People, Microbes,fill::zeros);
  for (int k=0; k<Groups_People; k++){
    for (int l=0; l<Groups_Microbes; l++){
      Normal+=(*Omega)(k,l);
    }}
  Normal+=nan_saver; //Avoid singularities
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.2.2. Compute final \omegas
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for (int k=0; k<Groups_People; k++)
    for (int l=0; l<Groups_Microbes; l++)
      for (int person=0; person<People; person++)
  	for (int microbe=0; microbe<Microbes; microbe++)
  	  (*Omega)(k,l)(person, microbe) /= Normal(person, microbe);
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
  return;
 }
//=========================================================================================

//1.3. Process Etas,Thetas and Ps
//================================================================================================
MatricesMemberships UpdateParameters\
(field<mat> Omega, int People, int Groups_People, int Microbes, int Groups_Microbes,\
 int Number_of_Links, mat Data, cube ListAbund, vec Tope,vec dm, vec dp ){
  MatricesMemberships Outputs;
  
  //1.3.1. Define Variables
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //theta
  mat theta_raw(People,Groups_People,fill::zeros);
  Outputs.theta=zeros(People,Groups_People);
  //eta
  mat eta_raw(Groups_Microbes,Microbes,fill::zeros);
  Outputs.eta=zeros(Microbes,Groups_Microbes);
  //pkl
  Outputs.P=zeros(Groups_People,Groups_Microbes,Number_of_Links);
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.3.2. Compute theta and eta
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //1.3.2.1. Compute non-normalized version
  //-----------------------------------------
  for (int k=0;k<Groups_People; k++){
    for (int l=0;l<Groups_Microbes; l++){      
	theta_raw.col(k)+=sum(Omega(k,l),1);
	eta_raw.row(l)+=sum(Omega(k,l),0);
      }}
  //-----------------------------------------
  
  //1.3.2.2. Transpose eta matrix
  //-----------------------------------------
  eta_raw=eta_raw.t();
  //-----------------------------------------
  
  //1.3.2.3. Normalize (divide by degree)
  //-----------------------------------------
  //Normalize theta
  for (int person=0;person<People;person++){
    Outputs.theta.row(person)=theta_raw.row(person)/dp(person);
  }
  
  //Normalize eta
  for (int microbe=0;microbe<Microbes;microbe++){
    Outputs.eta.row(microbe)=eta_raw.row(microbe)/dm(microbe);
  }
  //-----------------------------------------
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.3.3. Compute pkl matrices
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for (int abund=0; abund<Number_of_Links; abund++){
    for (int k=0; k<Groups_People; k++){
      for (int l=0; l<Groups_Microbes; l++){
  	//.....................................................................
	float Numerator=0;
  	for (int index=0; index<Tope(abund);index++){
  	  Numerator+= Omega(k,l)( ListAbund(index,1,abund) , ListAbund(index,0,abund) );
	}
  	//.....................................................................
  	Outputs.P(k,l,abund)=Numerator/accu(Omega(k,l));
        }}}
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  return Outputs;
    }
//================================================================================================

//1.4. Compute LogLikelihood
//================================================================================================
double ComputeLogLikelihood\
(int K,int L,int Links,vec number_of_links,cube ListAbund,MatricesMemberships Update_Matrices){

  double LogLikelihoodFun=0;
  mat Likelihood(MBiome,People,fill::zeros);
  
  for (int abund=0; abund<Links; abund++){	  
    for (int index=0; index<number_of_links(abund);index++){
      int person_L =ListAbund(index,1,abund);
      int microbe_L=ListAbund(index,0,abund);
      for (int k=0;k<K;k++){
	for (int l=0;l<L;l++){
	  Likelihood(microbe_L,person_L)+=\
  	    Update_Matrices.theta(person_L,k)*Update_Matrices.eta(microbe_L,l)*Update_Matrices.P(k,l,abund);
  	      }}
      LogLikelihoodFun+=log(Likelihood(microbe_L,person_L));

    }}

  return LogLikelihoodFun;
}
//================================================================================================

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


//MAIN
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char *argv[])
{
  wall_clock timer;

  //1.0. INPUT VARIABLES AND INPUT/OUTPUT FILES
  //=============================================================

  //0.1. Input Variables
  //:::::::::::::::::::::::::::::::::::::::::::::::::
  
  //0.1.1. Command line arguments
  //--------------------------------------------
  //Seed
  std::stringstream InputSeed(argv[1]);int Seed;
  if (!(InputSeed >> Seed)) // do the conversion
    cout<<"No Seed introduced"<<endl;
  arma_rng::set_seed(Seed);

  //K
  std::stringstream Input_K(argv[2]);int K;
  if (!(Input_K >> K)) // do the conversion
      cout<<"No K introduced"<<endl;

  //L
  std::stringstream Input_L(argv[3]);int L;
  if (!(Input_L >> L)) // do the conversion
      cout<<"No L introduced"<<endl;

  //Fold Number
  std::stringstream Input_FoldNumber(argv[4]);string FoldNumber;
  if (!(Input_FoldNumber >> FoldNumber)) // do the conversion
      cout<<"No Fold Number introduced"<<endl;
  //--------------------------------------------

  //0.1.2. Static Arguments
  //---------------------------------------
  int Limit=7000;    //Iterations                         
  float NanRate=0.2;//Number of nans
  //---------------------------------------

  //:::::::::::::::::::::::::::::::::::::::::::::::::

  //0.2. Input/Output Files
  //:::::::::::::::::::::::::::::::::::::::::::::::::

  //Choose Dataset
  //--------------------------------------------------------
  static std::vector<std::string> Datasets;
  Datasets={"Leave_One_Out_S-8_Tot/", "Leave_One_Out_V-10_Tot/", "Leave_One_Out_V-22_Tot/", "Leave_One_Out_V-23_24_Tot/", "Leave_One_Out_V-25_Tot"};

  static std::vector<std::string> Outputs;
  Outputs={"Leave_One_Out_S-8_Tot/", "Leave_One_Out_V-10_Tot/", "Leave_One_Out_V-22_Tot/", "Leave_One_Out_V-23_24_Tot/", "Leave_One_Out_V-25_Tot"};

  int Index=0; 	 //0: S-8, 1: V-10, 2: V-22, 3: V-23_24, 4: V-25
  string Dataset=Datasets[Index];
  string Output=Outputs[Index];
  //--------------------------------------------------------

  //Choose whether you are executing from cluster or locally
  //--------------------------------------------------------

  string Exec="local";
  string InputPath;string OutputPath;string InputPathTest;

  if (Exec=="local"){
    
    InputPath="/home/sergio/work/Tarragona/Microbiome/Input_Data/" + Dataset;
    OutputPath="/home/sergio/work/Tarragona/Microbiome/Output_Data/" + Output;
  }
  else {cout<<"you have to choose an execution mode"<<endl;

    }
    //--------------------------------------------------------
  
  //=============================================================
  

  //1.1. Read Data            
  //===========================================
  string InputFile;
  InputFile= InputPath + "P_51_F_" +  FoldNumber+ "_Train.txt";
  
  mat Data;
  Data.load(InputFile);
  People=Data.n_cols;MBiome=Data.n_rows; //People= NACEs. MBiome=LURs
  //===========================================

  
  //1.2. Make lists of abundance-person-microbe and lists of Nans
  //=============================================================================================
  
  //1.2.1.Compute number of slices of cube (hypermatrix). This is the number of types of links
  //------------------------------------------------------------------------
  int Links= Data.max() - Data.min() +1;
  //------------------------------------------------------------------------
  
  //1.2.2. Compute most frequent link and how many elements does it contain
  //------------------------------------------------------
  vec number_of_links=zeros<vec>(Links);
  
  for (int abundance=0;abundance<Links;abundance++){
    uvec Census=find(Data==abundance);      //Vector of indices of instances of an abundance
    int  Population=Census.n_rows;          //Vector dimension (number of links)
    number_of_links(abundance)=Population;  //Fill vector of abundances
  }
 
  int Top_Abundance=number_of_links.max(); //Maximum of vector
  //------------------------------------------------------


  //1.2.3. Save list of links in a cube. Save list of nans as a matrix
  //-----------------------------------------------------------------------------------

  //1.2.3.1. Define cube and counters
  //................................................................
  cube ListAbund(Top_Abundance,2,Links,fill::zeros); //Define cube
  vec contador=zeros<vec>( Links );		     //Create vector of counters
  //................................................................
  
  //1.2.3.2. Define nan list and counter
  //..........................................
  uvec RecordOfNans=find_nonfinite(Data); //Vector with indices of nans in the train matrix
  int NanRows=RecordOfNans.n_rows;        //Number of nans
    
  mat ListNans(NanRows,2,fill::zeros);//NANS LIST 
  int count_nans=0;
  //..........................................


  //1.2.3.3. Fill Lists
  //............................................................
  for(int row=0;row<MBiome;row++){
    for(int col=0;col<People;col++){
      
      //Nans: Skipt them and safe list
      //...............................
      if ( isnan( Data(row,col)) ){
	ListNans(count_nans,0)=row; //Row where nan lives
	ListNans(count_nans,1)=col; //Column where nan lives
	count_nans+=1;
	  continue;
	}
      //...............................

      ListAbund( contador( Data(row,col) ) ,0, Data(row,col) )=row; //First columna (Microbe)
      ListAbund( contador( Data(row,col) ) ,1, Data(row,col) )=col; //Second column (Host)
      contador( Data(row,col) )+=1; //Activate corresponding counter
    }}
  //............................................................


  //------------------------------------------------------
  //=============================================================================================

  
  //1.3. Make list of neighbors of each node
  //=============================================================================================

  //1.3.1. Create list as a vector and fill it automatically
  //----------------------------------------------------
  vec Nodes_Microbe(MBiome);Nodes_Microbe.fill(People);
  vec Nodes_People(People) ;Nodes_People.fill(MBiome);
  //----------------------------------------------------

  //1.3.2. Complete list. We remove nans from each corresponding element
  //----------------------------------------------------
  for (int i=0; i<NanRows;i+=1)
    { Nodes_Microbe( ListNans(i,0) )-=1;
      Nodes_People(  ListNans(i,1) )-=1;
    }
  //----------------------------------------------------
  //=============================================================================================



  //1.4. Declare matrices
  //===========================================

  //1.4.1. Matrices
  //:::::::::::::::::::::::::::::::::::::::::::::::::::
  MatricesMemberships Update_Matrices;
  Update_Matrices=GetMatrices(People,K,MBiome,L,Links);
  //:::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.4.2. Omega
  //:::::::::::::::::::::::::::::::::::::::::::::::::::
  field<mat> Omega(K,L);
  for (int i=0; i<K; i++)
    for (int j=0; j<L; j++)
      Omega(i,j) = zeros(People, MBiome);
  //:::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.4.3. Memory for computing Likelihood
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  double LMemory=ComputeLogLikelihood(K,L,Links, number_of_links, ListAbund, Update_Matrices);
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //===========================================

  //1.5. Expectation-Maximization Loop
  //================================================================================================

  //1.5.1. Declare variables
  //::::::::::::::::::::::::::::::::::::::::::::::::::::
  double LogLikelihood=0.0;
  
  int LikelihoodCounter=0;int BreakCounter=0;

  double DeltaLikelihood=0.0;
  //::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.5.2. Variables to export likelihood to file
  //::::::::::::::::::::::::::::::::::::::::::::::::::::
  std::string strSeed = std::to_string(Seed);
  std::string strK = std::to_string(K);
  std::string strL = std::to_string(L);

  ofstream OutputLogLikelihood;
  OutputLogLikelihood.open(OutputPath + FoldNumber + "_" +\
  "Seed_"+ strSeed + "_" + "K_"+strK+ "_" + "L_"+ strL +"_LogLikelihood_Test.txt");

  OutputLogLikelihood<< LMemory <<endl;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.5.3. Loop
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  timer.tic();
  for (int Iterations=0;Iterations<Limit;Iterations++){


    GetOmega(People,K,MBiome,L,Links,Update_Matrices,Data,ListAbund,number_of_links,&Omega);
    
    Update_Matrices = UpdateParameters				\
      (Omega,People,K,MBiome,L,Links,Data,ListAbund,number_of_links,Nodes_Microbe,Nodes_People);


    //1.5.2.1. Compute and save Likelihood + control convergence
    //---------------------------------------------------------------------------------------
    if (LikelihoodCounter==10){
      
      cout<<Iterations<<endl;

      //Compute Likelihood and broadcast to file
      //...............................................................................
      LogLikelihood=\
      	  ComputeLogLikelihood(K,L,Links, number_of_links, ListAbund, Update_Matrices);

      OutputLogLikelihood<< LogLikelihood <<endl;

      DeltaLikelihood = LogLikelihood - LMemory;
      //...............................................................................

      
      //Convergence Control
      //...............................
      if (abs(DeltaLikelihood)<0.5){
      	BreakCounter+=1;
      }
      else {
      	BreakCounter=0;
      }

      if (BreakCounter>5) {
	cout<<"Rompo la ejecucion"<<endl;
      	break;
      }
      //...............................
      
      LMemory=LogLikelihood;
      LikelihoodCounter=0;

    }
    //---------------------------------------------------------------------------------------

    LikelihoodCounter+=1;
    
  }
  //::::::::::::::::::::::::::::::::::::::::::::::::::::
  double n = timer.toc();
  cout << "number of seconds: " << n << endl;


  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.6. Decide whether to prompt Outputs during execution
  //==========================================================
  int outputs=0;
  if (outputs==1){
  	cout<<"\033[1;31m List Abundances: \033[0m \n"<<endl;
  	cout<<ListAbund<<"\n"<<endl;

	cout<<"\033[1;31m List of nans: \033[0m \n"<<endl;
	cout<<ListNans<<endl;
	
	cout<<"\033[1;31m Theta Matrix: \033[0m \n"<<endl;
	cout<<Update_Matrices.theta<<endl;

	cout<<"\033[1;31m Eta Matrix: \033[0m \n"<<endl;
	cout<<Update_Matrices.eta<<endl;

	cout<<"\033[1;31m P Matrices: \033[0m \n"<<endl;
	cout<<Update_Matrices.P<<endl;
  }
  //==========================================================

  //1.7. Write model to file
  //=========================================================================
  string ModelParameters;
  ModelParameters=OutputPath + FoldNumber + "_" +\
  "Seed_" + strSeed + "_K_" + strK + "_L_" + strL + "_Parameters_Test.txt";

  ofstream Parameters;
  Parameters.open(ModelParameters);

  Parameters << Update_Matrices.theta << "\n" << Update_Matrices.eta << "\n" << Update_Matrices.P << endl;
  //=========================================================================
  
  //1.8. Write scores to file
  //=========================================================================
  string OutputFile;
  OutputFile= OutputPath + FoldNumber + "_" + \
  "Seed_" + strSeed + "_" + "K_" + strK + "L_" + strL + "_scores_Test.txt";
    
  ofstream scores;
  scores.open(OutputFile);

  mat Scores(NanRows,Links,fill::zeros);
  for (int index=0; index<NanRows;index++){
    for (int type_link=0; type_link<Links; type_link++){
      for (int k=0;k<K;k++){
	for (int l=0;l<L;l++){
	  Scores(index,type_link)+=\
	    Update_Matrices.theta( ListNans(index,1), k)*Update_Matrices.eta( ListNans(index,0), l)*\
	    Update_Matrices.P(k,l,type_link);
	}}}
    scores<< ListNans(index,0) <<" "<<ListNans(index,1)<<" "<<Scores.row(index)<<endl;
  }

  //=========================================================================
  
  return 0;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
