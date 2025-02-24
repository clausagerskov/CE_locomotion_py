import subprocess
import argparse
import os
import sys
import helper_funcs as hf
# from pyneuroml.utils.cli import build_namespace
import random
from datetime import datetime
import json

DEFAULTS = {
    "popSize": None, #96,
    "duration": None, #24,
    #"simsep": False,
    "RandSeed": None,
    "folderName": None,
    "doEvol": False,
    "overwrite": False,
    "doNML": False,
    #"nervousSystemName" : 'nmlNervousSystem',
    #"nmlOutputFolderName": None,
}


def process_args():
    """Parse command-line arguments.

    :returns: None
    """
    parser = argparse.ArgumentParser(
        description=("A script for supplying arguments to execute Worm2D")
    )

    parser.add_argument(
        "-f",
        "--folderName",
        type=str,
        metavar="<folder name>",
        default=DEFAULTS["folderName"],
        help=(
            "Name of directory for output. This must be supplied."
            "If the directory does not exist it will be created and"
            "the evolution and simulation will be performed"
            "with supplied parameters, or defaults and a random seed if not supplied."
            "If the directory exists and overwrite is set true its contents will"
            "be modified. If doEvol is True the evolution and simulation will be"
            "executed with supplied parameters, or the values  existing in the folder if not supplied."
            "Simulation results are placed in the simulation subfolder."
            "If doEvol is False, only the simulation will be performed using the"
            "supplied parameters, or the values existing in the simulation subfolder if not supplied."
        ),
    )

    """     parser.add_argument(
        "-i",
        "--nmlOutputFolderName",
        type=str,
        metavar="<nml outpul folder name>",
        default=DEFAULTS["nmlOutputFolderName"],
        help=(
            "Name of directory for output from neuroml simulation.\n"
            "If not supplied neuroml simulation will not be peformed."
        ),
    )

    parser.add_argument(
        "-n",
        "--nervousSystemName",
        type=str,
        metavar="<nervous system name>",
        default=DEFAULTS["nervousSystemName"],
        help=(
            "Name of nervous system for neuroml simulation" 
        ),
    ) """


    parser.add_argument(
        "-N",
        "--doNML",
        action="store_true",
        default=DEFAULTS["doNML"],
        help=("Run the equivalent neuroML simulation too if true. Results will be"
              "placed in the nml_simulation subfolder. The simulation will use" 
              "the parameters supplied or if not those in the folder if it exists, or"
              "defaults otherwise."
              ),
    )

    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        default=DEFAULTS["overwrite"],
        help=("Overwrite the contents of the folder. If doEvol is set True" 
              "all contents will be overwritten. If doEvol is False"
              "only the simulation subfolder results will be overwritten."
              ),
    )

    parser.add_argument(
        "-E",
        "--doEvol",
        action="store_true",
        default=DEFAULTS["doEvol"],
        help=(
            "If used the evolutionary algorithm is executed, the best worm simulation performed,"
            "and the phenotype results are deposited in the folder. If parameters are not supplied"
            "parameters already existing in the folder will be used, or defaults if not."
            "Simulation results are placed in a subfolder."
            "If not used the only the simulation"
            "in the subfolder is executed and results deposited in it, again using parameters supplied,"
            "existing in the subfolder, or default if not."
        ),
    )

    """     parser.add_argument(
        "-S",
        "--simsep",
        action="store_true",
        default=DEFAULTS["simsep"],
        help=("If used, user input of the directory name is interactively requested."),
    ) """

    parser.add_argument(
        "-d",
        "--duration",
        type=float,
        metavar="<duration>",
        default=DEFAULTS["duration"],
        help="Duration of simulation for evolution and best worm in ms, default: %sms"
        % DEFAULTS["duration"],
    )

    parser.add_argument(
        "-p",
        "--popSize",
        type=int,
        metavar="<pop size>",
        default=DEFAULTS["popSize"],
        help="Population size for evolutionary algorithm, default: %s"
        % DEFAULTS["popSize"],
    )

    parser.add_argument(
        "-R",
        "--RandSeed",
        type=int,
        metavar="<rand seed>",
        default=DEFAULTS["RandSeed"],
        help="Seed value for evolution and simulation, or just simulation if doEvol is False." 
             "If not set the seed in the appropriate evolution or simulation directory will be used." 
             "If there is no such seed, a random seed will be generated."
        % DEFAULTS["RandSeed"],
    )

    return parser.parse_args()


def make_directory(directory_name, overwrite, str1 = 'the contents'):
    try:
        os.mkdir(directory_name)
        print(f"Directory '{directory_name}' created successfully.")
        return True
    except FileExistsError:
        if overwrite:
            print(
                f"Directory '{directory_name}' already exists and " + str1 + " will be overwritten."
            )
            return True
        else:
            print(f"Directory '{directory_name}' already exists and overwrite is false.")
            return False
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory_name}'.")
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)


def run_main(args=None):
    if args is None:
        args = process_args()
    run(a=args)


def build_namespace(DEFAULTS={}, a=None, **kwargs):
    if a is None:
        a = argparse.Namespace()

    # Add arguments passed in by keyword.
    for key, value in kwargs.items():
        setattr(a, key, value)

    # Add defaults for arguments not provided.
    for key, value in DEFAULTS.items():
        if not hasattr(a, key):
            setattr(a, key, value)

    return a

def setDict(dictval, keyval, parval, default_val):
        if parval is not None:
            dictval[keyval] = parval
            return
        if keyval not in dictval:
            dictval[keyval] = default_val
        return

        

def run(a=None, **kwargs):
    a = build_namespace(DEFAULTS, a, **kwargs)
    
    #folder_name = ""
    #do_evol = 1
    #nml_folder_name = None
    #do_nml = True
    if a.doEvol:
        do_evol = 1
    else:
        do_evol = 0
    if not os.path.isdir(a.folderName):
        do_evol = 1

    if do_evol:
           str1 = 'all the contents'
    else:
           str1 = 'the simulation subfolder(s) contents'    
    nml_sim_folder_name = None
    if a.folderName: 
        if not make_directory(a.folderName, a.overwrite, str1):
            print(f"Please change phenotype directory name or set overwrite to True\n"
                  "and doEvol to True to overwrite the evolution results.\n"
                  "or set overwrite to True and doEvol to False (the default) if you want\n"
                   "to just overwrite the simulation."
                  )
            sys.exit(1)
        sim_folder_name = a.folderName + '/sim_results'
        if not make_directory(sim_folder_name, a.overwrite):
                sys.exit(1)
        if a.doNML:
            nml_sim_folder_name = a.folderName + '/nml_sim_results'
            if not make_directory(nml_sim_folder_name, a.overwrite):
                sys.exit(1)
    else:
        print("You need to supply a folder name.")
        sys.exit(1)

    sim_par_file = sim_folder_name + '/pars.json'
    if os.path.isfile(sim_par_file):
        with open(sim_par_file) as f:
            sim_data = json.load(f)
    else:
        sim_data = {}

    evol_par_file = a.folderName + '/pars.json'
    if os.path.isfile(evol_par_file):
        with open(evol_par_file) as f:
            evol_data = json.load(f)
    else:
        evol_data = {}    
   
    random.seed(datetime.now().timestamp())
    random_seed = random.randint(1, 1000000)
    
    
    setDict(sim_data, 'seed', a.RandSeed, random_seed)
    setDict(sim_data, 'duration', a.duration, 24)
    if do_evol:
        setDict(evol_data, 'seed', a.RandSeed, random_seed)
        setDict(evol_data, 'popSize', a.popSize, 96)
        setDict(evol_data, 'duration', a.duration, 24)

    with open(evol_par_file, 'w', encoding='utf-8') as f:
        json.dump(evol_data, f, ensure_ascii=False, indent=4)
    with open(sim_par_file, 'w', encoding='utf-8') as f:
        json.dump(sim_data, f, ensure_ascii=False, indent=4)
    if a.doNML:
        with open(nml_sim_folder_name + '/pars.json', 'w', encoding='utf-8') as f:
            json.dump(sim_data, f, ensure_ascii=False, indent=4)  

    
    """     if a.doNML:
        do_nml = 1
    else:
        do_nml = 0 """


    cmd = ["./main",]
    cmd += ["-R", str(evol_data['seed'])]
    cmd += ["-sr", str(sim_data['seed'])]
    cmd += ["-p", str(evol_data['popSize'])]
    cmd += ["-d", str(evol_data['duration'])]
    cmd += ["-sd", str(sim_data['duration'])]
    cmd += ["--doevol", str(do_evol)]
    if a.doNML:
        cmd += ["--nmlfolder", str(nml_sim_folder_name)]
    cmd += ["--simfolder", str(sim_folder_name)]
    #cmd += ["--donml", str(do_nml)]
    cmd += ["--folder", str(a.folderName)]

    """ if a.RandSeed is not None:
        cmd = ["./main", "-R", str(a.RandSeed)]
    else:
        cmd = ["./main", "-r", str(a.randSeed)]

    cmd += [
        "-p",
        str(a.popSize),
        "-d",
        str(a.duration),
        "--doevol",
        str(do_evol),
        "--folder",
        folder_name,
        "--nervous",
        a.nervousSystemName,
    ]

    if nml_folder_name is not None:
        cmd += [
        "--nmlfolder",
        nml_folder_name
        ] """


    # Run the C++
    result = subprocess.run(cmd, capture_output=True, text=True)
    # result = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # output, errors = result.communicate()

    # p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # print(output)

    if result.stdout:
        print(result.stdout)

    if result.stderr:
        print("Error:")
        print(result.stderr)

    
    if nml_sim_folder_name is not None:
        hf.dir_name = nml_sim_folder_name
        from load_data import reload_single_run
        reload_single_run(show_plot=False)    
    if sim_folder_name is not None:
        hf.dir_name = sim_folder_name
        from load_data import reload_single_run
        reload_single_run(show_plot=False)


if __name__ == "__main__":
    run_main()
