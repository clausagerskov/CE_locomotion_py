//
//  Worm.hpp
//  one
//
//  Created by Eduardo Izquierdo on 9/25/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//
#pragma once
//#include "VectorMatrix.h"
//#include "random.h"
//#include "WormBody.h"
//#include "NervousSystem.h"
//#include "Muscles.h"
#include "StretchReceptor.h"
#include "Worm2D.h"

#include <cmath>

#define PI 3.14159265

using namespace std;

struct RS18Macros{
bool headsr;
bool vncsr;
};

class Worm18 : public WormIzq {
public:
    
    Worm18(TVector<double> &v, double output);
   
    void InitializeState(RandomState &rs);
    void HeadStep(double StepSize, double output);
    void Step(double StepSize, double output); 
    
    void DumpActState(ofstream &ofs, int skips);
    void DumpVoltage(ofstream &ofs, int skips);
    void DumpParams(ofstream &ofs);
    
    
    private:
    
    void addExtraParsToJson(json & j);
    const vector<string> & getCellNames() {return {"DB", "DD", "VBA", "VDA", "VBP", "VDP"};}
    Params<double> getWormParams();

    //WormBody b;
    //Muscles m;
    //NervousSystem n;
    StretchReceptor sr;
    NervousSystem h;
    
    
    
    // Neuromuscular junctions
    double NMJ_DB, NMJ_VBa, NMJ_VBp, NMJ_DD, NMJ_VDa, NMJ_VDp;
    double NMJ_SMDD, NMJ_RMDD, NMJ_SMDV, NMJ_RMDV;
    double NMJ_Gain_Map;
    
    TVector<double> NMJ_Gain;
    
    // Head oscillator
    //double dorsalinput1, ventralinput1, dorsalinput2, ventralinput2;
    //double headFreq, headDelay, headGain, headBias;

// Parameters
     //int N_muscles = 24;           // Number of muscles alongside the body
     //int N_units = 6;              // Number of neural units
     //int N_neuronsperunit = 6;     // Number of neurons in a neural unit

    const int N_stretchrec = 6;         // N_units + 1 // Number of stretch receptors
    //
     //double T_muscle = 0.1;        // Muscle time ant

    const int HeadMotorNeuronMuscles = 6;  // Head motorneurons innervate first 8 muscles (temporarily first 6)
    const int VNCMuscleStart = 7;           // VNC motorneurons innervate starting from 7th muscle
    const int NmusclePerNU = 3;             // All the way down to 24, in groups of 3 per unit

    // Neuron name conventions
    const int DB = 1;
    const int DD = 2;
    const int VBA = 3;
    const int VDA = 4;
    const int VBP = 5;
    const int VDP = 6;

    // Neuron name conventions
    const int SMDD = 1;
    const int RMDD = 2;
    const int SMDV = 3;
    const int RMDV = 4;

    //const int Head = 1;
    //const int Tail = N_segments;

    const RS18Macros rS18Macros;
    RS18Macros setMacros();

};
