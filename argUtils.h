#pragma once
#include <string>
//#include <vector>
//#include "VectorMatrix.h"


using std::string;
//using std::vector;


/* extern string nervousSystemName; 
extern string nervousSystemNameForSim; 
extern string nervousSystemNameForEvol; 
extern string output_dir_name; 
extern bool randomInit; // if true random nervous initial state

extern int pop_size;
extern bool simRandomInit;
extern bool do_evol;
extern bool do_nml;
extern int traceDuration; */

class SuppliedArgs
{

public:
SuppliedArgs();

bool setArgs(int argc, const char* argv[], const long & randomseed1);
void writeMessage();
string rename_file(const string & file_name);
void setSimRandomInit();
//void setDefaultArgs();


//string nervousSystemName;
string nervousSystemNameForSim;
//string nervousSystemNameForEvol;
string output_dir_name;
bool randomInit;
int pop_size;
bool simRandomInit;
bool do_evol;
bool do_nml;
int traceDuration;
bool doOrigNS;
long randomseed;
int max_gens;


};
class SuppliedArgs2018 : public SuppliedArgs
{

public:
SuppliedArgs2018();

bool printToFile;
bool evo_seed;
bool output;
bool speedoutput;

};

class SuppliedArgs2021 : public SuppliedArgs
{

public:
SuppliedArgs2021();
bool printToFile;
bool output;

};








