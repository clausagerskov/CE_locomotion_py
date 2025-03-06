#pragma once
#include <string>
#include <vector>
#include "VectorMatrix.h"


using std::string;
using std::vector;


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

string nervousSystemName;
string nervousSystemNameForSim;
string nervousSystemNameForEvol;
string output_dir_name;
bool randomInit;
int pop_size;
bool simRandomInit;
bool do_evol;
bool do_nml;
int traceDuration;
bool doOrigNS;
long randomseed;

};


extern SuppliedArgs supArgs1;  

//bool checkNervousSystemForJson();
//bool setArgs(int argc, const char* argv[], long & randomseed);

//string rename_file(const string & file_name);


template<class T>
vector<T> & append(vector<T> & v1, const vector<T> & v2)
{
v1.insert(v1.end(), v2.begin(), v2.end());
return v1;
}    

template<class T> 
vector<T> getVector(TVector<T> & vec, int size)
{ 
vector<T> retvec;    
for (int i = 1; i <= size; i++)
        retvec.push_back(vec[i]);   
return retvec;    
}

template<class T> 
TVector<T> getTVector(vector<T> & vec)
{ 
TVector<T> retvec;
retvec.SetBounds(1,vec.size());    
for (int i = 0; i < vec.size(); i++) retvec[i+1]=vec[i];
return retvec;    
}




