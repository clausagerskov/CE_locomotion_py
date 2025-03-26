//
//  Worm.hpp
//  one
//
//  Created by Eduardo Izquierdo on 9/25/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//

//#include "VectorMatrix.h"
//#include "random.h"
//#include "WormBody21.h"
//#include "NervousSystem.h"
//#include "Muscles.h"
#include "Worm2D21.h"

#include <cmath>

#define PI 3.14159265

using namespace std;

// Parameters
//const int N_muscles = 24;           // Number of muscles alongside the body
//const int N_units = 7;              // Number of neural units
//const int N_neuronsperunit = 7;     // Number of neurons in a neural unit

//const double T_muscle = 0.1;        // Muscle time constant



// Body segment name conventions
//const int Head = 1;
//const int Tail = N_segments;

class Worm21 : public Worm2D21 {
public:
    
    Worm21(TVector<double> &v);
    
    void InitializeState(RandomState &rs);
    void DumpParams(ofstream &ofs);
    
    
    NervousSystem & n;

    protected:
  
    
    void addParsToJson(json & j){
        string nsHead = "Nervous system";
        appendAllNSJson(j[nsHead], n);
        Worm2D21::addParsToJson(j);}

    
};
