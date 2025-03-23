#pragma once
//#include "random.h"

class NSForW2D {
    public:
    
virtual double NeuronOutput(int i) = 0;
virtual void SetNeuronExternalInput(int i, double value) = 0;
virtual void EulerStep(double stepsize) = 0;
virtual ~NSForW2D(){};

};

