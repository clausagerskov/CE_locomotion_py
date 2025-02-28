import sys
import os
import json
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns


sys.path.append("..")

from run_main import make_directory


def getCellNames(network_json_data):
    return network_json_data["Nervous system"]["Cell name"]["value"]


def getNervousSystemVal(network_json_data, val):
    return network_json_data["Nervous system"][val]["value"]


def getPopNamesUnsorted(network_json_data, pop_num=6):
    cell_names = getCellNames(network_json_data)
    return cell_names[:pop_num]


def check_equal(list):
    return all(i == list[0] for i in list)


def getWeight(weights, pre_cell_ind, post_cell_ind):
    return [
        weight["weight"]
        for weight in weights
        if ((weight["from"] == pre_cell_ind) & (weight["to"] == post_cell_ind))
    ][0]


def getWeightsDict(weights, cell_names):
    new_weights = {}
    for weight in weights:
        kk = cell_names[weight["from"] - 1] + cell_names[weight["to"] - 1]
        if kk in new_weights:
            if new_weights[kk] != weight["weight"]:
                print("Weights not equal!")
                sys.exit(1)
        else:
            new_weights[kk] = weight["weight"]
    return new_weights

def getUniques(list1):
    result1 = defaultdict(list)
    for d in list1:
        for key, value in d.items():
            result1[key].append(value)
    return result1

path_list = []
outFolderBases = ["varyEvolSeeds", "varyEvolSeeds1"]
current = os.path.dirname(os.path.realpath(__file__))  # location of this file!
for outFolderBase in outFolderBases:
    path = current + "/" + outFolderBase
    dir_list = sorted(os.listdir(path))
    path_list += [path + "/" + dir for dir in dir_list]

print(len(path_list))
# sys.exit(1)

""" con_list = [
DA->VD: 15.0622
 DB->VD: 6.38036
 VD->VA: -4.39552
 VD->VB: -10.1964
 VA->DD: 15.0622
 VA->VD: 0
 VB->DD: 6.38036
 VB->VD: 0
] """

fitness_keys = ['Best', 'Average', 'Variance']
weights_list = []
fitness_list = []
for dir in path_list:
    json_file = dir + "/worm_data.json"
    with open(json_file, "r") as file:
        worm_data = json.load(file)
    cell_names = getCellNames(worm_data)
    chem_weights = getNervousSystemVal(worm_data, "Chemical weights")
    weights_list.append(getWeightsDict(chem_weights, cell_names))
    ol1 = None
    with open(dir + '/fitness.dat', "r") as file:
        while True:
            l1=file.readline()
            if l1:
                ol1 = l1
            else:  
                break
            #print(l1)
        file.close()
    if ol1 is not None:
        d1 = {}
        for key, val in zip(fitness_keys, ol1.split()):
            d1[key]=float(val)
        fitness_list.append(d1)
    else:
        print("No fitness")
        sys.exit(1)


""" print(fitness_list)
sys.exit(1)
print(weights_list)
sys.exit(1)
 """
all_weights = getUniques(weights_list)
all_fitnesses = getUniques(fitness_list)

print(all_weights)
print(all_fitnesses)


#sys.exit(1)
make_directory("results")
fit_level = 0.95

for key, val in all_weights.items():
    sns.displot(val, bins=10, kde=True)
    plt.savefig("results/" + key + ".png")
    plt.close()
    best_fit = [val1 for val1, fitness in zip(val, all_fitnesses['Best']) if fitness > fit_level]
    sns.displot(best_fit, bins=10, kde=True)
    plt.savefig("results/best_fit_" + key + ".png")
    plt.close()
