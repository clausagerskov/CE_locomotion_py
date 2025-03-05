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


def getParsDict(vals):
    new_vals = {}
    for key, val in vals.items():
        new_vals[key] = val["value"]
    return new_vals


def getValsDict(vals, cell_names):
    new_vals = {}
    for ind, val in enumerate(vals):
        kk = cell_names[ind]
        if kk in new_vals:
            if new_vals[kk] != val:
                print("Vals not equal!")
                sys.exit(1)
        else:
            new_vals[kk] = val
    return new_vals


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
outFolderBases = ["varyEvolSeeds", "varyEvolSeeds1", "varyEvolSeeds2", "varyEvolSeeds3"]
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

fitness_keys = ["Best", "Average", "Variance"]
weights_list = []
elec_weights_list = []
fitness_list = []
biases_list = []
worm_vals_list = []
for dir in path_list:
    json_file = dir + "/worm_data.json"
    if not os.path.isfile(json_file):
        break
    if not os.path.isfile(dir + "/fitness.dat"):
        break
    with open(json_file, "r") as file:
        worm_data = json.load(file)
    cell_names = getCellNames(worm_data)
    chem_weights = getNervousSystemVal(worm_data, "Chemical weights")
    weights_list.append(getWeightsDict(chem_weights, cell_names))
    elec_weights = getNervousSystemVal(worm_data, "Electrical weights")
    elec_weights_list.append(getWeightsDict(elec_weights, cell_names))
    biases = getNervousSystemVal(worm_data, "biases")
    biases_list.append(getValsDict(biases, cell_names))
    worm_vals = getParsDict(worm_data["Worm"])
    worm_vals["SR_A_gain"] = worm_data["Stretch receptor"]["SR_A_gain"]["value"]
    worm_vals["SR_B_gain"] = worm_data["Stretch receptor"]["SR_B_gain"]["value"]
    worm_vals_list.append(worm_vals)
    ol1 = None
    with open(dir + "/fitness.dat", "r") as file:
        while True:
            l1 = file.readline()
            if l1:
                ol1 = l1
            else:
                break
            # print(l1)
        file.close()
    if ol1 is not None:
        d1 = {}
        for key, val in zip(fitness_keys, ol1.split()):
            d1[key] = float(val)
        fitness_list.append(d1)
    else:
        print("No fitness")
        sys.exit(1)


""" print(fitness_list)
sys.exit(1)
print(weights_list)
sys.exit(1)
 """
# print(worm_vals_list)
# sys.exit(1)

# all_weights = getUniques(weights_list)
all_fitnesses = getUniques(fitness_list)

# print(all_weights)
print(all_fitnesses)


# sys.exit(1)

if not make_directory("results", overwrite=True):
    sys.exit(1)
fit_level = 0.95

results_titles = ["ChemWei", "ElectWei", "Bias", "Worm"]
data_results_list = [weights_list, elec_weights_list, biases_list, worm_vals_list]

for title, data_result in zip(results_titles, data_results_list):
    res1 = getUniques(data_result)
    for key, val in res1.items():
        sns.displot(val, bins=10, kde=True)
        title_str = title + "_" + key + ".png"
        plt.savefig("results/" + title_str)
        plt.close()
        best_fit = [
            val1
            for val1, fitness in zip(val, all_fitnesses["Best"])
            if fitness > fit_level
        ]
        sns.displot(best_fit, bins=10, kde=True)
        plt.savefig("results/best_fit_" + title_str)
        plt.close()
