#include "Worm.h"
#include <nlohmann/json.hpp>
#include <vector>
#include "jsonUtils.h"
#include <iomanip>
#include "../argUtils.h"

extern SuppliedArgs2018 supArgs1;

using std::vector;
using json = nlohmann::json;

struct toFromWeight{
    
    toFromWeight(weightentry w_val, int to_val){w=w_val;to=to_val;}
    toFromWeight(){}
    weightentry w;
    int to;
};


void to_json(json & j, const weightentry & w)
{
  j = json{{"from", w.from}, {"weight", w.weight}};
}

void to_json(json & j, const toFromWeight & w)
{
  j = json{{"to", w.to}, {"from", w.w.from}, {"weight", w.w.weight}};
}

void from_json(const json& j, toFromWeight & w) 
{
        j.at("to").get_to(w.to);
        j.at("from").get_to(w.w.from);
        j.at("weight").get_to(w.w.weight);
}

template<class T>
void appendToJson(json & j, const Params<T> & par)
{
    size_t mess_ind = 0;
    for (size_t i=0;i<par.names.size(); i++) {
        if (par.messages_inds.size()>mess_ind && par.messages_inds[mess_ind]==static_cast<int>(i)) 
        {j[par.names[i]]["message"] = par.messages[i];mess_ind++;}
        j[par.names[i]]["value"] = par.vals[i];
        }
               
}

void appendMatrixToJson(json & j, TMatrix<weightentry> & vec, TVector<int> & sizes, int tot_size)
{    
    vector<toFromWeight> newvec;
    for (int i=1; i<=tot_size; i++){    
        for (int j=1; j<=sizes[i]; j++) { 
            toFromWeight tv(vec[i][j], i);
            newvec.push_back(tv);}        
    }
    j["value"] = newvec;

}
void appendNSToJson(json & j, NervousSystem& c)
{
    j["Chemical weights"]["message"] = "chemical weights in sparse format";
    appendMatrixToJson(j["Chemical weights"], c.chemicalweights, c.NumChemicalConns, c.size);
    appendMatrixToJson(j["Electrical weights"], c.electricalweights, c.NumElectricalConns, c.size);
    j["Electrical weights"]["message"] = "electrical weights in sparse format";
}

Params< vector<double> > getNervousSysParamsDoubleNH(NervousSystem& c)
{
Params< vector<double> > par;
par.names = {"taus", "biases", "gains", "outputs", "states", "paststates", "Rtaus", "externalinputs"};
par.vals = {
getVector<double>(c.taus, c.size), 
getVector<double>(c.biases, c.size), 
getVector<double>(c.gains, c.size),
getVector<double>(c.outputs, c.size),
getVector<double>(c.states, c.size),
getVector<double>(c.paststates, c.size),
getVector<double>(c.Rtaus, c.size),
getVector<double>(c.externalinputs, c.size),
};
return par;
}




// Neuron name conventions



Params< vector<string> > getNervousSysCellNames(vector<string> & cell_names, int neurons_per_unit)
{
Params< vector<string> > par;
par.names = {"Cell name"};
vector<string> cell_names_all;
for (int i=0;i<neurons_per_unit;i++) cell_names_all.insert(cell_names_all.end(),cell_names.begin(),cell_names.end());
par.vals = {cell_names_all};
return par;
}




Params<int> getNervousSysParamsIntNH(NervousSystem& c)
{
Params<int> par;    
par.names = {"size", "maxchemcons", "maxelecconns"};
par.vals = {c.size, c.maxchemconns, c.maxelecconns};
return par;
}

Params< vector<int> > getNervousSysVecInt(NervousSystem& c)
{
Params< vector<int> > par;
par.names = {"NumChemicalConns", "NumElectricalConns"};
par.vals = {
getVector<int>(c.NumChemicalConns, c.size), 
getVector<int>(c.NumElectricalConns, c.size), 
};
return par;
}


void getNSJson(NervousSystem & n, json & j)
{

{Params<vector<double> > parvec = getNervousSysParamsDoubleNH(n);
appendToJson<vector<double> >(j,parvec);}
        
{Params<int> parvec = getNervousSysParamsIntNH(n);
appendToJson<int>(j,parvec);}

{Params< vector<int> > parvec = getNervousSysVecInt(n);
appendToJson<vector<int> >(j,parvec);}

appendNSToJson(j, n);

}

void addNSJson(Worm & w, json & j)
{

{string nsHead = "Nervous system";
getNSJson(w.n, j[nsHead]);
vector<string> cell_names = {"DB", "DD", "VBA", "VDA", "VBP", "VDP"};
{Params< vector<string> > parvec = getNervousSysCellNames(cell_names, 6);
appendToJson<vector<string> >(j[nsHead],parvec);}}

{string nsHead = "Head Nervous system";
getNSJson(w.h, j[nsHead]);
vector<string> cell_names = {"SMDD", "RMDD", "SMDV", "RMDV"};
{Params< vector<string> > parvec = getNervousSysCellNames(cell_names, 1);
appendToJson<vector<string> >(j[nsHead],parvec);}}   

}

Params<double> getStretchReceptorParams(StretchReceptor& s)
{
Params<double> par;
par.names = {"NSR", "NSEGS", "NSEGSSR", "VNC_gain", "Head_gain"};
par.vals = {s.NSR, s.NSEGS, s.NSEGSSR, s.SRvncgain, s.SRheadgain};
par.messages = {"Number of stretch receptor",
                                "Number of segments in the body",
                                "Number of segments sensed by each stretch receptor"};
par.messages_inds = {0,1,2}; //must be ordered
return par;
}

void addSRJson(Worm & w, json & j)
{
    Params<double> par = getStretchReceptorParams(w.sr);
    appendToJson<double>(j["Stretch receptor"],par);
}


Params<double> getWormParams(Worm & w)
{
Params<double> par;
par.vals = {w.NMJ_DB, w.NMJ_VBa, w.NMJ_VBp, w.NMJ_DD, w.NMJ_VDa, w.NMJ_VDp, 
    w.NMJ_SMDD, w.NMJ_SMDV, w.NMJ_RMDD, w.NMJ_RMDV, w.NMJ_Gain_Map};
par.names = {"NMJ_DB", "NMJ_VBa", "NMJ_VBp", "NMJ_DD", "NMJ_VDa", "NMJ_VDp",
    "NMJ_SMDD", "NMJ_SMDV", "NMJ_RMDD", "NMJ_RMDV", "NMJ Gain"};
return par;
}

void addWormJson(Worm & w, json & j)
{
    Params<double> par = getWormParams(w);
    appendToJson<double>(j["Worm"],par);
   
}

void writeParsToJson(Worm & w)
{
    if (supArgs1.doOrigNS){

        json j;

        cout << "making json" << endl;
        vector<doubIntParamsHead> parvec;
        //doubIntParamsHead var1;
        doubIntParamsHead var1 = supArgs1.getParams();
        /* var1.parInt.head = "Evolutionary Optimization Parameters";
        var1.parInt.names = {"pop_size", "randomseed"};
        var1.parInt.vals = {supArgs1.pop_size, supArgs1.randomseed};
        var1.parInt.messages ={"population size", "seed"};
        var1.parInt.messages_inds = {0,1}; */
        parvec.push_back(var1);

        for (size_t i=0;i<parvec.size(); i++) {
            if (strcmp(parvec[i].parDoub.head.c_str(),"NULL")!=0)
            appendToJson<double>(j[parvec[i].parDoub.head],parvec[i].parDoub);
            if (strcmp(parvec[i].parInt.head.c_str(),"NULL")!=0)
            appendToJson<long>(j[parvec[i].parInt.head],parvec[i].parInt);
            }

       // writeParsToJson(w, "worm_data.json", evolutionParams);
        addNSJson(w,j);
        addWormJson(w,j);
        addSRJson(w,j);

        ofstream json_out(supArgs1.rename_file("worm_data.json"));
        json_out << std::setw(4) << j << std::endl;
        json_out.close();


    }
}


/* 

void DumpParamsOrig(ofstream &ofs)
{
    ofs << "Time-constants: \n DB: " << n.NeuronTimeConstant(DB) << "\n VBA/P: " << n.NeuronTimeConstant(VBA) << " / " << n.NeuronTimeConstant(VBP) << "\n DD: " << n.NeuronTimeConstant(DD) << "\n VDA/P: " << n.NeuronTimeConstant(VDA) << " / " << n.NeuronTimeConstant(VDP) << endl;
    ofs << "Biases: \n DB: " << n.NeuronBias(DB) << "\n VBA/P: " << n.NeuronBias(VBA) << " / " << n.NeuronBias(VBP)  <<  "\n DD: " << n.NeuronBias(DD) << "\n VDA/P: " << n.NeuronBias(VDA) <<  " / " << n.NeuronBias(VDP) << endl;
    ofs << "Self conns: \n DB: " << n.ChemicalSynapseWeight(DB, DB) << "\n VBA/P: " << n.ChemicalSynapseWeight(VBA, VBA) << " / " << n.ChemicalSynapseWeight(VBP, VBP) << "\n DD: " << n.ChemicalSynapseWeight(DD, DD) <<  "\n VDA/P: " << n.ChemicalSynapseWeight(VDA, VDA) <<  " / " << n.ChemicalSynapseWeight(VDP, VDP) << endl;
    ofs << "Chem Conns: \n DB->DD: " << n.ChemicalSynapseWeight(DB, DD) <<  "\n DB->VDA/VDP: " << n.ChemicalSynapseWeight(DB, VDA) << " / " << n.ChemicalSynapseWeight(DB, VDP) << "\n VBA/P->DD: " << n.ChemicalSynapseWeight(VBA, DD) << " / " << n.ChemicalSynapseWeight(VBP, DD) << "\n VBA/P->VDA/P: " << n.ChemicalSynapseWeight(VBA, VDA) << " / " << n.ChemicalSynapseWeight(VBP, VDP) << "\n VDA/P->VBA/P: " << n.ChemicalSynapseWeight(VDA, VBA) << " / " << n.ChemicalSynapseWeight(VDP, VBP) << "\n DD->VDA: " << n.ChemicalSynapseWeight(DD, VDA) <<endl;
    ofs << "Gap Juncs: \n DB-DB+1: " << n.ElectricalSynapseWeight(DB, DB+N_neuronsperunit) << "\n VBA-VBP / VBP-VBP+1: " << n.ElectricalSynapseWeight(VBA, VBP) << " / " << n.ElectricalSynapseWeight(VBP, VBA+N_neuronsperunit) << "\n VBP-DB+1: " << n.ElectricalSynapseWeight(VBP, DB+N_neuronsperunit) << "\n DD-VDA/P: " << n.ElectricalSynapseWeight(DD, VDA) << " / " << n.ElectricalSynapseWeight(DD, VDP) << "\n DD-DD+1: " << n.ElectricalSynapseWeight(DD, DD+N_neuronsperunit) << "\n VDA-VDP / VDP-VDP+1: " << n.ElectricalSynapseWeight(VDA, VDP) << " / " << n.ElectricalSynapseWeight(VDP, VDA+N_neuronsperunit) <<  endl;
    ofs << "SR Gain (VNC and Head): " << sr.SRvncgain << " " << sr.SRheadgain << endl;
    ofs << "NMJ weights: \n B: " << NMJ_DB << " " << NMJ_VBa << " " << NMJ_VBp << "\n D: " <<  NMJ_DD << " " << NMJ_VDa << " " << NMJ_VDp << endl;
    ofs << "Head: \nBiases: \n SMD(D/V): " << h.NeuronBias(SMDD) << " " << h.NeuronBias(SMDV) << "\n RMD(D/V): "<< h.NeuronBias(RMDD) << " "<< h.NeuronBias(RMDV) << endl;
    ofs << "Time-constants: \n SMD(D/V): " << h.NeuronTimeConstant(SMDD) << " " << h.NeuronTimeConstant(SMDV) << "\n RMD(D/V): " << h.NeuronTimeConstant(RMDD) << " " << h.NeuronTimeConstant(RMDV) << endl;
    ofs << "Self conns: \n SMD(D/V): " <<h.ChemicalSynapseWeight(SMDD,SMDD) << " " << h.ChemicalSynapseWeight(SMDV,SMDV) << "\n RMD(D/V): " << h.ChemicalSynapseWeight(RMDD,RMDD) << " "<< h.ChemicalSynapseWeight(RMDV,RMDV) << endl;
    ofs << "Chem conns: " << "\n SMDD->RMDD: " << h.ChemicalSynapseWeight(SMDD, RMDD) << "\n SMDD->SMDV: " << h.ChemicalSynapseWeight(SMDD, SMDV) << "\n SMDD->RMDV: " << h.ChemicalSynapseWeight(SMDD, RMDV) << "\n RMDD->SMDD: " << h.ChemicalSynapseWeight(RMDD, SMDD) << "\n RMDD->SMDV: " << h.ChemicalSynapseWeight(RMDD, SMDV) << "\n RMDD->RMDV: " << h.ChemicalSynapseWeight(RMDD, RMDV) << "\n SMDV->SMDD: " << h.ChemicalSynapseWeight(SMDV, SMDD) << "\n SMDV->RMDD: " << h.ChemicalSynapseWeight(SMDV, RMDD) << "\n SMDV->RMDV: " << h.ChemicalSynapseWeight(SMDV, RMDV) << "\n RMDV->SMDD: " << h.ChemicalSynapseWeight(RMDV, SMDD) << "\n RMDV->RMDD: " << h.ChemicalSynapseWeight(RMDV, RMDD) << "\n RMDV->SMDV: " << h.ChemicalSynapseWeight(RMDV, SMDV) << endl;
    ofs << "Gap Juncs: " << "\n SMD-RMD: " << h.ElectricalSynapseWeight(SMDD, RMDD)
                         << "\n RMD-RMD: " << h.ElectricalSynapseWeight(RMDD, RMDV) << endl;
    ofs << "NMJ weights: \n SMD(D/V): " << NMJ_SMDD << " " << NMJ_SMDV << "\n RMD(D/V): " <<  NMJ_RMDD << " " << NMJ_RMDV << endl;
    ofs << "NMJ Gain: " << NMJ_Gain_Map << endl;
}
 */