#include "argUtils.h"
#include <string.h>
#include <sys/stat.h>
//#include <stdio.h>
#include <iostream>

using namespace std;

//// helper functions

/* string nervousSystemName = "NervousSystem";
string nervousSystemNameForSim = "nmlNervousSystem";
string nervousSystemNameForEvol = "NervousSystem";
string output_dir_name = "";
bool randomInit = 0;
int pop_size = 96;
bool simRandomInit = 0;
bool do_evol = 1;
bool do_nml = 0;
int traceDuration = 24; */


//SuppliedArgs supArgs1;  

/* string rename_file(const string & file_name){
  if (supArgs1.output_dir_name != "") return supArgs1.output_dir_name + "/" + file_name;
  return file_name;
} */

/* bool checkNervousSystemForJson(){
return (strcmp(supArgs1.nervousSystemName.c_str(), "NervousSystem") == 0);
} */

SuppliedArgs::SuppliedArgs()
{
  nervousSystemNameForSim = "NervousSystem";
  //nervousSystemNameForEvol = "NervousSystem";
  output_dir_name = "";
  randomInit = 0;
  pop_size = 96;
  simRandomInit = 0;
  do_evol = 1;
  do_nml = 0;
  traceDuration = 24;
  doOrigNS = 1;
  max_gens = 10;
}

/* void SuppliedArgs::setDefaultArgs()
{
  nervousSystemNameForSim = "NervousSystem";
  //nervousSystemNameForEvol = "NervousSystem";
  output_dir_name = "";
  randomInit = 0;
  pop_size = 96;
  simRandomInit = 0;
  do_evol = 1;
  do_nml = 0;
  traceDuration = 24;
  doOrigNS = 1;

} */

SuppliedArgs2018::SuppliedArgs2018():SuppliedArgs()
{
  do_evol = 1;
  pop_size = 96;
  max_gens = 1000;
  printToFile = 1;
  evo_seed = 0;
  output = 0;
  speedoutput = 0;
}


void SuppliedArgs::setSimRandomInit()
{
  randomInit = simRandomInit;
}

string SuppliedArgs::rename_file(const string & file_name)
{
  if (output_dir_name != "") return output_dir_name + "/" + file_name;
  return file_name;
}

void SuppliedArgs::writeMessage()
{
  cout << "Run evaluation with seed: " << randomseed << ", pop size: " << pop_size << endl;
}

bool SuppliedArgs::setArgs(int argc, const char* argv[], const long & randomseed1)
{

if (((argc-1) % 2) != 0)
     {cout << "The arguments are not configured correctly." << endl;return 0;}
    
    bool seed_flag = 1;
    
    for (int arg = 1; arg<argc; arg+=2)
    { 
    if (strcmp(argv[arg],"--maxgens")==0) max_gens = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"--doevol")==0) do_evol = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"--dorandinit")==0) simRandomInit = atoi(argv[arg+1]);
    //if (strcmp(argv[arg],"--skipOrigSim")==0) skipOrigSim = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"--donml")==0) do_nml = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"--folder")==0) {
      output_dir_name = argv[arg+1];
      struct stat sb;
      if (stat(output_dir_name.c_str(), &sb) != 0) 
      {cout << "Directory doesn't exist." << endl;return 0;}
    }

    if (seed_flag){ 
    if (strcmp(argv[arg],"-R")==0) randomseed = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"-r")==0) randomseed = randomseed1 + atoi(argv[arg+1]);
    seed_flag = 0;
    }

    if (strcmp(argv[arg],"-p")==0) pop_size = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"-d")==0) traceDuration = atoi(argv[arg+1]);
    if (strcmp(argv[arg],"--nervous")==0) 
    {
      nervousSystemNameForSim = argv[arg+1];
    }
    }
  return 1;

  }

