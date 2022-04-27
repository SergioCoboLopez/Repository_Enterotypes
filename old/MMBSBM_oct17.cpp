//19/03/2018. MixMembership Stochastic Block Model
//Incluimos una estructura para samplear la likelihood; de cada
//x iteraciones, calculamos la likelihood solo en una. 
//En la matriz \omega: filas=personas, columnas=microbios
//En el dataset:       filas=microbios, columnas=personas

//Toy Dataset: 1220 Microbes, 370 people


#include <iostream>
#include <fstream>
#include "armadillo"

using namespace std;
using namespace arma;

int People;
int MBiome;

//0. VARIABLES                                                      
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct MatricesMemberships{
  mat theta;    //Person Memberships 
  mat eta;      //Microbe Memberships
  cube P;       //P Matrices
};


float nan_saver=1e-10;//--->numero para prevenir singularidades en los calculos
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//1. FUNCTIONS
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//1.1. Get Initial \eta and \theta and ps                                                               
//================================================================================================
MatricesMemberships GetMatrices				\
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
   //Crear una matriz de dimensiones (Numero de links,Numero de grupos de microbios) 
    mat M=randu(Number_of_Links,Groups_Microbes);    
    M=normalise(M,1,0);//Normalizar sus columnas
    
    for (int k=0; k<Number_of_Links; k++){
      //Para cada slice del tensor, rellenar cada una de las k filas con las probabilidades p(k,:)
      //Iterar sobre filas
      Init_Outputs.P.slice(k).row(i)=M.row(k);	
      }}    
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  return Init_Outputs;
}
//================================================================================================

//1.2. Get \Omegas Marzo 2018
//================================================================================================
void GetOmegaMarch\
(int People,int Groups_People, int Microbes, int Groups_Microbes, int Number_of_Links,\
 MatricesMemberships InitMatrices, mat Data, cube ListAbund, vec Tope, field<mat>& Omega) {

  //1.2.0. Declare tensors and subtensors (matrices)
  //:::::::::::::::::::::::::::::::::
  //mat KL = zeros(People,Microbes);
  //field<mat> Omega_raw(Groups_People,Groups_Microbes);
  //:::::::::::::::::::::::::::::::::

  //1.2.1. Fill tensor with zero matrices
  //:::::::::::::::::::::::::::::::::
  //Omega_raw.fill(zeros(People,Microbes));
  // for (int i=0; i<Groups_People; i++){
  //   for (int j=0; j<Groups_Microbes; j++){
  // 	Omega_raw(i,j) = zeros(People,Microbes);
  //     }}
  //:::::::::::::::::::::::::::::::::

  //1.2.2. Compute \omega numerators
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for (int k=0; k<Groups_People; k++){
    for (int l=0; l<Groups_Microbes; l++){
      for (int abund=0; abund<Number_of_Links; abund++){	  
  	for (int index=0; index<Tope(abund);index++){
  	  Omega(k,l)( ListAbund(index,1,abund) , ListAbund(index,0,abund) ) = \
  	    InitMatrices.eta( ListAbund(index,0,abund) ,l)*InitMatrices.theta( ListAbund(index,1,abund) ,k) \
  	    *InitMatrices.P(k,l,abund);
	}}}}

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.2.3. Compute normalizing constants
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  mat Normal(People, Microbes,fill::zeros);
  for (int k=0; k<Groups_People; k++){
    for (int l=0; l<Groups_Microbes; l++){
      Normal+=Omega(k,l);
    }}
  Normal+=nan_saver; //Evitamos singularidades
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.2.4. Compute final \omegas
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // field<mat> OmegaTest(Groups_People,Groups_Microbes);
  for (int k=0; k<Groups_People; k++)
    for (int l=0; l<Groups_Microbes; l++)
      for (int person=0; person<People; person++)
  	for (int microbe=0; microbe<Microbes; microbe++)
  	  Omega(k,l)(person, microbe) /= Normal(person, microbe);
  // for (int k=0; k<Groups_People; k++)
  //   for (int l=0; l<Groups_Microbes; l++)
  // 	  Omega(k,l) = Omega_raw(k,l) /Normal;

  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  
  return;
 }
//================================================================================================

//1.2. Get \Omegas
//================================================================================================
field<mat> GetOmega\
(int People,int Groups_People, int Microbes, int Groups_Microbes, int Number_of_Links,\
 MatricesMemberships InitMatrices, mat Data, cube ListAbund, vec Tope) {
  
  //1.2.0. Declare tensors and subtensors (matrices)
  //:::::::::::::::::::::::::::::::::
  mat KL = zeros(People,Microbes);
  field<mat> Omega_raw(Groups_People,Groups_Microbes);
  //:::::::::::::::::::::::::::::::::

  //1.2.1. Fill tensor with zero matrices
  //:::::::::::::::::::::::::::::::::
  for (int i=0; i<Groups_People; i++){
    for (int j=0; j<Groups_Microbes; j++){
  	Omega_raw(i,j) = KL;
      }}
  //:::::::::::::::::::::::::::::::::

  //1.2.2. Compute \omega numerators
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  for (int k=0; k<Groups_People; k++){
    for (int l=0; l<Groups_Microbes; l++){
      for (int abund=0; abund<Number_of_Links; abund++){	  
  	for (int index=0; index<Tope(abund);index++){
  	  Omega_raw(k,l)( ListAbund(index,1,abund) , ListAbund(index,0,abund) ) = \
  	    InitMatrices.eta( ListAbund(index,0,abund) ,l)*InitMatrices.theta( ListAbund(index,1,abund) ,k) \
  	    *InitMatrices.P(k,l,abund);
  	    }}}}
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.2.3. Compute normalizing constants
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  mat Normal(People, Microbes,fill::zeros);
    for (int k=0; k<Groups_People; k++){
      for (int l=0; l<Groups_Microbes; l++){
  	Normal+=Omega_raw(k,l);
      }}
    Normal+=nan_saver; //Evitamos singularidades
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.2.4. Compute final \omegas
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  field<mat> OmegaTest(Groups_People,Groups_Microbes);
    for (int k=0; k<Groups_People; k++){
      for (int l=0; l<Groups_Microbes; l++){
  	OmegaTest(k,l)=Omega_raw(k,l)/Normal;
      }}
  //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  return OmegaTest;
}
//================================================================================================


//1.3. Process Etas,Thetas and Ps (Marzo 2018)
//================================================================================================
MatricesMemberships UpdateParametersMarch\
(field<mat> Omega, int People, int Groups_People, int Microbes, int Groups_Microbes,\
 int Number_of_Links, mat Data, cube ListAbund, vec Tope,vec dm, vec dp, MatricesMemberships Outputs){
  
  //1.3.1. Define Variables
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //theta
  //----------------------------------------------
  mat theta_raw(People,Groups_People,fill::zeros);
  //----------------------------------------------
  //eta
  //----------------------------------------------
  mat eta_raw(Groups_Microbes,Microbes,fill::zeros);
  //----------------------------------------------
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

//1.3. Process Etas,Thetas and Ps
//================================================================================================
MatricesMemberships UpdateParameters\
(field<mat> Omega, int People, int Groups_People, int Microbes, int Groups_Microbes,\
 int Number_of_Links, mat Data, cube ListAbund, vec Tope,vec dm, vec dp ){
  MatricesMemberships Outputs;
  
  //1.3.1. Define Variables
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //theta
  //----------------------------------------------
  mat theta_raw(People,Groups_People,fill::zeros);
  Outputs.theta=zeros(People,Groups_People);
  //----------------------------------------------
  //eta
  //----------------------------------------------
  mat eta_raw(Groups_Microbes,Microbes,fill::zeros);
  Outputs.eta=zeros(Microbes,Groups_Microbes);
  //----------------------------------------------
  //pkl
  //----------------------------------------------
  Outputs.P=zeros(Groups_People,Groups_Microbes,Number_of_Links);
  //----------------------------------------------
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
float ComputeLogLikelihood\
(int K,int L,int Links,vec number_of_links,cube ListAbund,MatricesMemberships Update_Matrices){

  float LogLikelihood=0;
  mat Likelihood(MBiome,People,fill::zeros);
  
  for (int abund=0; abund<Links; abund++){	  
    for (int index=0; index<number_of_links(abund);index++){
      int person_L =ListAbund(index,1,abund);
      int microbe_L=ListAbund(index,0,abund);
      for (int k=0;k<K;k++){
	for (int l=0;l<L;l++){
	  Likelihood(microbe_L,person_L)+=					\
  	    Update_Matrices.theta(person_L,k)*Update_Matrices.eta(microbe_L,l)*Update_Matrices.P(k,l,abund);
  	      }}
      LogLikelihood+=log(Likelihood(microbe_L,person_L));

    }}

  return LogLikelihood;
}
//================================================================================================

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//MAIN
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int main(int argc, char *argv[])
{
  wall_clock timer;

  //0. INPUT VARIABLES AND INPUT/OUTPUT FILES
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

  string  FoldNumber="0";
  //--------------------------------------------

  //0.1.2. Static Arguments
  //---------------------------------------
  int Limit=400;    //Iterations                         
  float NanRate=0.2;//Number of nans
  //---------------------------------------

  //:::::::::::::::::::::::::::::::::::::::::::::::::

  //0.2. Input/Output Files
  //:::::::::::::::::::::::::::::::::::::::::::::::::

  //Choose whether you are executing from cluster or locally
  //--------------------------------------------------------
  string Exec="local";
  string InputPath;string OutputPath;string InputPathTest;
  
  if (Exec=="cluster"){
    InputPath="/home/sees/cobo/5Fold_ToyData/";
    OutputPath="/home/sees/cobo/Output_Data/";}
  else if( Exec=="local" ) {
    InputPath="/export/home/shared/Projects/Microbiome/Input_Data/5Fold_ToyData/";
    InputPathTest="/export/home/shared/Projects/Microbiome/Input_Data/Micro5Fold/";
    OutputPath="/export/home/shared/Projects/Microbiome/Output_Data/";}
  else {cout<<"you have to choose an execution mode"<<endl;
    }
    //--------------------------------------------------------
  
  //=============================================================
  
  
  //1. Read Data            
  //===========================================
  string InputFile;
  InputFile= InputPath + FoldNumber+ "Train.txt";
  
  mat Data;
  Data.load(InputFile);
  People=Data.n_cols;MBiome=Data.n_rows;
  //===========================================

  
  //1.2. Make lists of abundance-person-microbe and lists of Nans
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  //1.2.1.Calculamos el numero de slices que tiene el cubo (numero de links)
  //------------------------------------------------------------------------
  int Links= Data.max() - Data.min() +1;
  //------------------------------------------------------------------------
  
  //1.2.2.Calculamos que tipo de link tiene el mayor numero de elementos y cuantos elementos tiene
  //------------------------------------------------------
  vec number_of_links=zeros<vec>(Links);
  
  for (int abundance=0;abundance<Links;abundance++){
    uvec Census=find(Data==abundance);      //Vector de indices de las ocurrencias de una abundancia
    int  Population=Census.n_rows;          //Dimension del vector (numero de links)
    number_of_links(abundance)=Population;  //Rellenar el vector de abundancias
  }
 
  int Top_Abundance=number_of_links.max(); //Coger el maximo del vector
  //------------------------------------------------------


  //1.2.3. Guardamos la lista de links en forma de cubo y la de nans en forma de matriz 
  //-----------------------------------------------------------------------------------

  //1.2.3.1. Definir cubo y contadores
  //................................................................
  cube ListAbund(Top_Abundance,2,Links,fill::zeros);//Crear el cubo
  vec contador=zeros<vec>( Links );		     //Crear vector de contadores
  //................................................................
  
  //1.2.3.2. Definir lista de nans y contador
  //..........................................
  int NanRows=round(NanRate*People*MBiome);
  mat ListNans(NanRows,2,fill::zeros);//LISTA DE NANS. 
  int count_nans=0;
  //..........................................


  //1.2.3.3. Rellenar listas
  //............................................................
  for(int row=0;row<MBiome;row++){
    for(int col=0;col<People;col++){
      
      //Nans: Saltarlos y guardar lista
      //...............................
      if ( isnan( Data(row,col)) ){
	ListNans(count_nans,0)=row; //Fila donde vive el nan
	ListNans(count_nans,1)=col; //Columna donde vive el nan
	count_nans+=1;
	  continue;
	}
      //...............................
      
      ListAbund( contador( Data(row,col) ) ,0, Data(row,col) )=row; //Primera columna (Microbio)
      ListAbund( contador( Data(row,col) ) ,1, Data(row,col) )=col; //Segunda columna (Persona)

      contador( Data(row,col) )+=1; //Activar contador que toque
    }}
  //............................................................
  
  //------------------------------------------------------

  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  //1.3. Make list of neighbors of each node
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  //1.3.1. Crear lista como vector y rellenarlo por defecto
  //----------------------------------------------------
  vec Nodes_Microbe(MBiome);Nodes_Microbe.fill(People);
  vec Nodes_People(People) ;Nodes_People.fill(MBiome);
  //----------------------------------------------------

  //1.3.2. Completar lista. A cada elemento le quitamos la cantidad de nans detectados
  //----------------------------------------------------
  for (int i=0; i<NanRows;i+=1)
    { Nodes_Microbe( ListNans(i,0) )-=1;
      Nodes_People(  ListNans(i,1) )-=1;
    }
  //----------------------------------------------------
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



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
  float LMemory=ComputeLogLikelihood(K,L,Links, number_of_links, ListAbund, Update_Matrices);
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  //===========================================


  //1.5. Expectation-Maximization Loop
  //=============================================================================================================

  //1.5.1. Declare variables
  //::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  float LSampling=0.1;                  //Porcentaje de veces que calculas la likelihood
  int LSize= round(LSampling*Limit);    //Dimension del vector que almacena la likelihood
  int LStep= round(Limit/LSize);        //Iteraciones en las que hay que samplear
  vec LogLikelihoodVector(LSize,fill::zeros); //Vector de Likelihood
  
  int LikelihoodCounter=0;int BreakCounter=0;
  int indexL=0;

  float DeltaLikelihood=0;
  //::::::::::::::::::::::::::::::::::::::::::::::::::::
  

  //1.5.2. Loop
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  timer.tic();
  for (int Iterations=0;Iterations<Limit;Iterations++){

    //Funciones Viejas
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Omega=GetOmega(People,K,MBiome,L,Links,Update_Matrices,Data,ListAbund,number_of_links);

    // Update_Matrices=							\
    //   UpdateParameters(Omega,People,K,MBiome,L,Links,Data,ListAbund,number_of_links,Nodes_Microbe,Nodes_People);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    //Funciones Nuevas
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GetOmegaMarch(People,K,MBiome,L,Links,Update_Matrices,Data,ListAbund,number_of_links,Omega);
    //  cout << Omega(0, 0)(0, 0) << endl;
    
    Update_Matrices = UpdateParametersMarch				\
      (Omega,People,K,MBiome,L,Links,Data,ListAbund,number_of_links,Nodes_Microbe,Nodes_People,Update_Matrices);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    //1.5.2.1. Compute and save Likelihood + control convergence
    //---------------------------------------------------------------------------------------
    if (LikelihoodCounter==10){
      // indexL=Iterations/LStep;

      cout<<Iterations<<endl;

      
      LogLikelihoodVector(indexL)=\
      	  ComputeLogLikelihood(K,L,Links, number_of_links, ListAbund, Update_Matrices);

      DeltaLikelihood = LogLikelihoodVector(indexL) - LMemory;
      
      //Convergence Control
      //...............................
      if (abs(DeltaLikelihood)<0.5){
      	BreakCounter+=1;
      }
      else {
      	BreakCounter=0;
      }

      if (BreakCounter>5) {
      	break;
      }
      //...............................
      
      LMemory=LogLikelihoodVector(indexL);
      LikelihoodCounter=0;
      indexL+=1;
    }
    //---------------------------------------------------------------------------------------

    LikelihoodCounter+=1;
    
  }
  //::::::::::::::::::::::::::::::::::::::::::::::::::::
  double n = timer.toc();
  cout << "number of seconds: " << n << endl;


  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  

  
  //MODO WHILE
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  // int WhileCounter=0;
  // timer.tic();
  // while (BreakCounter<=3){
    
  //   Omega=GetOmega(People,K,MBiome,L,Links,Update_Matrices,Data,ListAbund,number_of_links);

  //   // Update_Matrices=							\
  //   //   UpdateParameters(Omega,People,K,MBiome,L,Links,Data,ListAbund,number_of_links,Nodes_Microbe,Nodes_People);


  //   Update_Matrices=UpdateParametersMarch\
  //     (Omega,People,K,MBiome,L,Links,Data,ListAbund,number_of_links,Nodes_Microbe,Nodes_People,Update_Matrices);

  //   //Bucle para calcular y comparar Likelihood
  //   //-----------------------------------------------------------------------------------
  //   if (WhileCounter==10){

  //     //Calcular Likelihood
  //     //...............................................................................
  //     LogLikelihoodVector(indexL)=\
  //     	  ComputeLogLikelihood(K,L,Links, number_of_links, ListAbund, Update_Matrices);
  //     //...............................................................................

  //     //Calcular diferencia de Likelihoods
  //     //......................................................
  //     DeltaLikelihood = LogLikelihoodVector(indexL) - LMemory;
  //     //......................................................

  //     //...............................................
  //     if (abs(DeltaLikelihood)<0.5){
  // 	BreakCounter+=1;
  // 	cout<<"Contador "<<BreakCounter<<endl;}
  //     else {
  // 	BreakCounter=0;
  //     }
  //     //...............................................
      
  //     //Actualizar Memoria
  //     //..................................
  //     LMemory=LogLikelihoodVector(indexL);
  //     //..................................
      
  //     WhileCounter=0;
  //     indexL+=1;
  //   }

  //   WhileCounter+=1;  
  // }
  // double n = timer.toc();
  // cout << "number of seconds: " << n << endl;
  //-----------------------------------------------------------------------------------

  
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  

  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  //=============================================================================================================


  //1.6. Decide whether to print Outputs?
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

  
  //1.7. Write scores to file
  //=========================================================================
  string OutputFile;
  std::string strSeed = std::to_string(Seed);
  std::string strK = std::to_string(K);
  std::string strL = std::to_string(L);
  OutputFile= OutputPath + FoldNumber + "_" + "Seed_" + strSeed + "_" + "K_" + strK + "L_" + strL + "_scores.txt";
    
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



  //1.8. Export LogLikelihood
  //=========================================================================
  ofstream OutputLogLikelihood;
  OutputLogLikelihood.open(OutputPath + FoldNumber + "_" +"Seed_"+ strSeed + "_" + "K_"+strK+ "_" + "L_"+ strL \
			   +"_LogLikelihood.txt");
  
  for (const auto& j : LogLikelihoodVector){
    OutputLogLikelihood<< j <<endl;
  }

  OutputLogLikelihood<<Update_Matrices.theta<<'\n'<<Update_Matrices.eta<<'\n'<<Update_Matrices.P<<'\n'<<Omega<<endl;
  //=========================================================================
  
  return 0;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
