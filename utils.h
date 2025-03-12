#pragma once
#include <string>
#include <vector>


using std::string;
using std::vector;

template <class T>
struct Params {
Params(){}    
vector<string> names;
vector<T> vals;
vector<int> messages_inds;
vector<string> messages;
};

template <class T>
struct ParamsHead : Params<T> {
ParamsHead(string head_val, Params<T> par_val):Params<T>(par_val){head=head_val;}
ParamsHead():Params<T>(){head = "NULL";}
string head;
};

struct doubIntParamsHead
{
ParamsHead<double> parDoub;
ParamsHead<long> parInt;
};





