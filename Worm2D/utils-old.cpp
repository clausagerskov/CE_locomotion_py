#include "Worm.h"
#include <nlohmann/json.hpp>
#include <vector>

using std::vector;
using json = nlohmann::json;

void dumpParamsAsJson(Worm & w, ofstream &ofs)
{
json j;
json jns = j["Nervous System"];
vector<string> names = {"DB", "VBA", "VBP", "DD", "VDA", "VDP"};
{json jtcs = jns["Time constants"];
    vector<double> vals = {w.n.NeuronTimeConstant(DB), w.n.NeuronTimeConstant(VBA),
    w.n.NeuronTimeConstant(VBP), w.n.NeuronTimeConstant(DD),
    w.n.NeuronTimeConstant(VDA), w.n.NeuronTimeConstant(VDP)
};
for (int i=0;i<names.size();i++) jtcs[names[i]]["value"] = vals[i];
}
{json jtcs = jns["Bias"];
    vector<double> vals = {w.n.NeuronBias(DB), w.n.NeuronBias(VBA),
    w.n.NeuronBias(VBP), w.n.NeuronBias(DD),
    w.n.NeuronBias(VDA), w.n.NeuronBias(VDP)
};
for (int i=0;i<names.size();i++) jtcs[names[i]]["value"] = vals[i];
}
{json jtcs = jns["Chemical connections"];
    vector<double> vals = {w.n.NeuronBias(DB), w.n.NeuronBias(VBA),
    w.n.NeuronBias(VBP), w.n.NeuronBias(DD),
    w.n.NeuronBias(VDA), w.n.NeuronBias(VDP)
};
for (int i=0;i<names.size();i++) jtcs[names[i]]["value"] = vals[i];
}


}


/* jtcs['DB']['value'] = w.n.NeuronTimeConstant(DB);
jtcs['VBA']['value'] = w.n.NeuronTimeConstant(VBA);
jtcs['VBP']['value'] = w.n.NeuronTimeConstant(VBP);
jtcs['DD']['value'] = w.n.NeuronTimeConstant(DD);
jtcs['VDA']['value'] = w.n.NeuronTimeConstant(VDA);
jtcs['VDP']['value'] = w.n.NeuronTimeConstant(VDP); */




void dumpLine(ofstream &ofs, const string & str1, const double & val)
{
ofs << str1.c_str() << " " << val << endl;
}


void DumpParams(Worm & w, ofstream &ofs)
{
    dumpLine(ofs, "DB", w.n.NeuronTimeConstant(DB));
    dumpLine(ofs, "VBA", w.n.NeuronTimeConstant(VBA));
    dumpLine(ofs, "VBP", w.n.NeuronTimeConstant(VBP));
    dumpLine(ofs, "DD", w.n.NeuronTimeConstant(DD));
    dumpLine(ofs, "VDA", w.n.NeuronTimeConstant(VDA));
    dumpLine(ofs, "VDP", w.n.NeuronTimeConstant(VDP));


    ofs << w.n.NeuronTimeConstant(DB) << endl; 
    ofs << n.NeuronTimeConstant(VBA) << endl;
    << n.NeuronTimeConstant(VBP) << "\n DD: " 
    << n.NeuronTimeConstant(DD) << "\n VDA/P: " 
    << n.NeuronTimeConstant(VDA) << " / " 
    << n.NeuronTimeConstant(VDP) << endl;
    
    ofs << "Biases: \n DB: " << 
    n.NeuronBias(DB) << "\n VBA/P: " 
    << n.NeuronBias(VBA) << " / " 
    << n.NeuronBias(VBP)  <<  "\n DD: " 
    << n.NeuronBias(DD) << "\n VDA/P: " 
    << n.NeuronBias(VDA) <<  " / " 
    << n.NeuronBias(VDP) << endl;


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