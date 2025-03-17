//
//  Worm21.cpp
//  one
//
//  Created by Eduardo Izquierdo on 9/25/15.
//  Copyright © 2015 Eduardo Izquierdo. All rights reserved.
//

#include "Worm21.h"
#include "../argUtils.h"


extern SuppliedArgs2021 supArgs1;

// The constructor
Worm21::Worm21(TVector<double> &v)
:WormIzq({7,24,0.1,7}, new NervousSystem),n(static_cast<NervousSystem&>(*n_ptr))
{
    // Muscles
    //m.SetMuscleParams(N_muscles, T_muscle);
    
    // Nervous system // Ventral cord
    n_ptr->SetCircuitSize(par1.N_units*par1.N_neuronsperunit, 9, 6);
    
    int as, da, db, dd, vd, vb, va;
    int asNext, dbNext, ddNext, vdNext, vbNext, vaNext ;
    
    for (int u = 1; u <= par1.N_units; u++){
        as = nn(AS, u);
        da = nn(DA, u);
        db = nn(DB, u);
        dd = nn(DD, u);
        vd = nn(VD, u);
        vb = nn(VB, u);
        va = nn(VA, u);

        asNext = nn(AS, u+1);
        dbNext = nn(DB, u+1);
        ddNext = nn(DD, u+1);
        vdNext = nn(VD, u+1);
        vbNext = nn(VB, u+1);
        vaNext = nn(VA, u+1);
        
        // Bias, Time Constant and Self Connections
        n_ptr->SetNeuronBias(as, v(1));
        n_ptr->SetNeuronBias(da, v(2));
        n_ptr->SetNeuronBias(db, v(3));
        n_ptr->SetNeuronBias(dd, v(4));
        n_ptr->SetNeuronBias(vd, v(5));
        n_ptr->SetNeuronBias(vb, v(6));
        n_ptr->SetNeuronBias(va, v(7));

        n_ptr->SetNeuronTimeConstant(as, v(8));
        n_ptr->SetNeuronTimeConstant(da, v(9));
        n_ptr->SetNeuronTimeConstant(db, v(10));
        n_ptr->SetNeuronTimeConstant(dd, v(11));
        n_ptr->SetNeuronTimeConstant(vd, v(12));
        n_ptr->SetNeuronTimeConstant(vb, v(13));
        n_ptr->SetNeuronTimeConstant(va, v(14));
        
        n_ptr->SetChemicalSynapseWeight(as, as, v(15));
        n_ptr->SetChemicalSynapseWeight(da, da, v(16));
        n_ptr->SetChemicalSynapseWeight(db, db, v(17));
        n_ptr->SetChemicalSynapseWeight(dd, dd, v(18));
        n_ptr->SetChemicalSynapseWeight(vd, vd, v(19));
        n_ptr->SetChemicalSynapseWeight(vb, vb, v(20));
        n_ptr->SetChemicalSynapseWeight(va, va, v(21));
        
        // --------
        // Chemical Synapses minimal network
        n_ptr->SetChemicalSynapseWeight(as, da, v(22));
        n_ptr->SetChemicalSynapseWeight(as, vd, v(23));
        n_ptr->SetChemicalSynapseWeight(da, db, v(24));
        n_ptr->SetChemicalSynapseWeight(db, as, v(25));
        n_ptr->SetChemicalSynapseWeight(vd, va, v(26));
        n_ptr->SetChemicalSynapseWeight(vd, vb, v(27));

        n_ptr->SetChemicalSynapseWeight(da, dd, v(28));
        n_ptr->SetChemicalSynapseWeight(vb, dd, v(29));
        n_ptr->SetChemicalSynapseWeight(va, dd, v(30));

        // Electrical Synapse minimal network
        n_ptr->SetElectricalSynapseWeight(vd, dd, v(31));

//        // Intersegment connections
//        // Chemicals
        if (u < par1.N_units){
            n_ptr->SetChemicalSynapseWeight(db, ddNext, v(40));
            n_ptr->SetChemicalSynapseWeight(vaNext, dd, v(41));
        }
//        // Electricals
        if (u < par1.N_units){
//        // Interclasses
            n_ptr->SetElectricalSynapseWeight(as, vaNext, v(42));
            n_ptr->SetElectricalSynapseWeight(da, asNext, v(43));
            n_ptr->SetElectricalSynapseWeight(vb, dbNext, v(44));
        // Intraclasses
//            n_ptr->SetElectricalSynapseWeight(db, dbNext, v(32));
//            n_ptr->SetElectricalSynapseWeight(vb, vbNext, v(32));
//            n_ptr->SetElectricalSynapseWeight(vd, vdNext, v(32));
//            n_ptr->SetElectricalSynapseWeight(dd, ddNext, v(32));
        }
    }
    
    // Interneuron inputs (AVB)
    wAVB_DB = 1;
    wAVB_VB = 1;
    // Interneuron inputs (AVB)
    wAVA_DA = 1;
    wAVA_VA = 1;
    
    
    // NMJ Weight
    NMJ_AS = v(32);
    NMJ_DA = v(33);
    NMJ_DB = v(34);
    NMJ_DD = v(35);
    NMJ_VD = v(36);
    NMJ_VB = v(37);
    NMJ_VA = v(38);
    
    // NMJ Gain XXX
    NMJ_Gain_Map = v(39);
    NMJ_Gain.SetBounds(1, par1.N_muscles);
    for (int i=1; i<=par1.N_muscles; i++)
    {
        NMJ_Gain(i) = 0.7*(1.0 - (((i-1)*NMJ_Gain_Map)/par1.N_muscles));
    }
}

void Worm21::InitializeState(RandomState &rs)
{    
    WormIzq::InitializeState(rs);
    n_ptr->RandomizeCircuitOutput(0.5, 0.5, rs); //fix this error!! adam (should be -0.5?)
    
}

void Worm21::Step(double StepSize)
{
    int mi;
    double dorsalInput, ventralInput;
    
    // Update Body
    b.StepBody(StepSize);
    
    
    // Update Nervous System
    n_ptr->EulerStep(StepSize);
    
    // Interneuron input  //////////////////////
    for (int i = 1; i <= par1.N_units; i++){
        n_ptr->SetNeuronExternalInput(nn(DB, i), wAVB_DB * AVB);
        n_ptr->SetNeuronExternalInput(nn(VB, i), wAVB_VB * AVB);
        n_ptr->SetNeuronExternalInput(nn(DA, i), wAVA_DA * AVA);
        n_ptr->SetNeuronExternalInput(nn(VA, i), wAVA_VA * AVA);
    }
    
    // Set input to Muscles
    // Head: 4 muscles one neural unit  //////////////////////
    mi = 1;
    dorsalInput  = NMJ_AS*n_ptr->NeuronOutput(nn(AS,mi)) + NMJ_DA*n_ptr->NeuronOutput(nn(DA,mi)) + NMJ_DB*n_ptr->NeuronOutput(nn(DB,mi)) + NMJ_DD*n_ptr->NeuronOutput(nn(DD,mi));
    ventralInput = NMJ_VD*n_ptr->NeuronOutput(nn(VD,mi)) + NMJ_VA*n_ptr->NeuronOutput(nn(VA,mi)) + NMJ_VB*n_ptr->NeuronOutput(nn(VB,mi));
    for (int i = 1; i < 5; i++){
                m.SetVentralMuscleInput(i, NMJ_Gain(i)*ventralInput);
                m.SetDorsalMuscleInput(i, NMJ_Gain(i)*dorsalInput);
    }
    
    //  Body Anterior: 4 segments, 3 muscles each  //////////////////////
    for (int mi = 2; mi <= 5; mi++){
        dorsalInput  = NMJ_AS*n_ptr->NeuronOutput(nn(AS,mi)) + NMJ_DA*n_ptr->NeuronOutput(nn(DA,mi)) + NMJ_DB*n_ptr->NeuronOutput(nn(DB,mi)) + NMJ_DD*n_ptr->NeuronOutput(nn(DD,mi));
        ventralInput = NMJ_VD*n_ptr->NeuronOutput(nn(VD,mi)) + NMJ_VA*n_ptr->NeuronOutput(nn(VA,mi)) + NMJ_VB*n_ptr->NeuronOutput(nn(VB,mi));
        for (int i = 5 + 3*(mi-2); i < 5 + 3*(mi-1); i++){
            m.SetVentralMuscleInput(i, NMJ_Gain(i)*ventralInput);
            m.SetDorsalMuscleInput(i, NMJ_Gain(i)*dorsalInput);
        }
    }
    
    //  Posterior Body: 2 segments, 4 muscles each  //////////////////////
    for (int mi = 6; mi <= 7; mi++){
        dorsalInput  = NMJ_AS*n_ptr->NeuronOutput(nn(AS,mi)) + NMJ_DA*n_ptr->NeuronOutput(nn(DA,mi)) + NMJ_DB*n_ptr->NeuronOutput(nn(DB,mi)) + NMJ_DD*n_ptr->NeuronOutput(nn(DD,mi));
        ventralInput = NMJ_VD*n_ptr->NeuronOutput(nn(VD,mi)) + NMJ_VA*n_ptr->NeuronOutput(nn(VA,mi)) + NMJ_VB*n_ptr->NeuronOutput(nn(VB,mi));
        for (int i = 17 + 4*(mi-6); i < 17 + 4*(mi-5); i++){
            m.SetVentralMuscleInput(i, NMJ_Gain(i)*ventralInput);
            m.SetDorsalMuscleInput(i, NMJ_Gain(i)*dorsalInput);
        }
    }
    
    // Update Muscle activation
    m.EulerStep(StepSize);
    
    // Set input to Body
    //  First two segments receive special treatment because they are only affected by a single muscle
    b.SetDorsalSegmentActivation(1, m.DorsalMuscleOutput(1)/2);
    b.SetVentralSegmentActivation(1, m.VentralMuscleOutput(1)/2);
    b.SetDorsalSegmentActivation(2, m.DorsalMuscleOutput(1)/2);
    b.SetVentralSegmentActivation(2, m.VentralMuscleOutput(1)/2);
    
    //  All other segments receive force from two muscles
    for (int i = 3; i <= N_segments-2; i++)
    {
        mi = (int) ((i-1)/2);
        b.SetDorsalSegmentActivation(i, (m.DorsalMuscleOutput(mi) + m.DorsalMuscleOutput(mi+1))/2);
        b.SetVentralSegmentActivation(i, (m.VentralMuscleOutput(mi) + m.VentralMuscleOutput(mi+1))/2);
    }
    
    //  Last two segments receive special treatment because they are only affected by a single muscle
    b.SetDorsalSegmentActivation(N_segments-1, m.DorsalMuscleOutput(par1.N_muscles)/2);
    b.SetVentralSegmentActivation(N_segments-1, m.VentralMuscleOutput(par1.N_muscles)/2);
    b.SetDorsalSegmentActivation(N_segments, m.DorsalMuscleOutput(par1.N_muscles)/2);
    b.SetVentralSegmentActivation(N_segments, m.VentralMuscleOutput(par1.N_muscles)/2);
    
    // Time
    t += StepSize;
}

vector<doubIntParamsHead> Worm21::getWormParams(){

    vector<doubIntParamsHead> parvec;
    doubIntParamsHead var1;
    
    var1.parDoub.head = "Worm";
    var1.parDoub.names = {"NMJ_AS", "NMJ_DA", "NMJ_DB", "NMJ_VD", "NMJ_VB", "NMJ_VA", "NMJ_DD"};
    var1.parDoub.vals = {NMJ_AS, NMJ_DA, NMJ_DB, NMJ_VD, NMJ_VB, NMJ_VA, NMJ_DD};

    append<string>(var1.parDoub.names,{"wAVA_DA", "wAVA_VA", "wAVB_DB", "wAVB_VB", "AVA, AVB"});
    append<double>(var1.parDoub.vals, {wAVA_DA, wAVA_VA, wAVB_DB, wAVB_VB, AVA, AVB});
   
  
    var1.parInt.head = "Worm";
    var1.parInt.vals = {startingMuscleA,NmusclePerNUA, startingMuscleB,NmusclePerNUB};
    var1.parInt.names = {"startingMuscleA","NmusclePerNUA", "startingMuscleB","NmusclePerNUB"};
    
  
    parvec.push_back(var1);
    return parvec;
  
}




void Worm21::DumpActState(ofstream &ofs, int skips)
{
    static int tt = skips;
    
    if (++tt >= skips) {
        tt = 0;
        //time
        ofs << t;

        // Ventral Cord Motor Neurons
        //ofs << "\nV: ";
        for (int i = 1; i <= par1.N_units; i++) {
            for (int j = 1; j <= par1.N_neuronsperunit; j++) {
                ofs <<  " " << n_ptr->NeuronOutput(nn(j,i));
            }
        }
        // Muscles
        //ofs << "\nM: ";
        for (int i = 1; i <= par1.N_muscles; i++) {
            ofs <<  " " << m.DorsalMuscleOutput(i) << " " << m.VentralMuscleOutput(i);
        }
        ofs << "\n";
    }
}

void Worm21::DumpCurvature(ofstream &ofs, int skips)
{
    
    double dx1,dy1,dx2,dy2,a,a1,a2,seg;
    static int tt = skips;
    
    if (++tt >= skips) {
        tt = 0;
        //time
        ofs << t;
        
        for (int i = 3; i < N_segments-1; i+=2)
        {
            dx1 = b.X(i) - b.X(i-2);
            dy1 = b.Y(i) - b.Y(i-2);
            dx2 = b.X(i+2) - b.X(i);
            dy2 = b.Y(i+2) - b.Y(i);
            
            a1 = atan2(dy1,dx1);
            a2 = atan2(dy2,dx2);
            
            if (a1 > PI/2 and a2 < -PI/2)
            a = (a1 - 2*PI) - a2;
            else
            if (a1 < -PI/2 and a2 > PI/2)
            a = a1 - (a2 - 2*PI);
            else
            a = a1-a2;
            
            seg = sqrt(pow(b.X(i-2)-b.X(i+2),2) + pow(b.Y(i-2)-b.Y(i+2),2));
            ofs <<  " " << (2*sin(a)/seg)/1000;
        }
        ofs << "\n";
    }
}


void Worm21::DumpParams(ofstream &ofs)
{ofs << "Biases: \n DB: " << n_ptr->NeuronBias(DB) << "\n VB/P: " << n_ptr->NeuronBias(VB) << " / " << n_ptr->NeuronBias(VB)  << "\n VDA/P: " << n_ptr->NeuronBias(VD) <<  " / " << n_ptr->NeuronBias(VD) << endl;
}
