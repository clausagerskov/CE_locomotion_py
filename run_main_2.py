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
    "randSeed": None,
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
            "Name of directory for output.\n"
            "If not supplied, both evolutionary algorithm and simulation of best worm are performed,\n"
            "and results placed in current directory."
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
        help=("Run the equivalent neuroML simulation too if true."),
    )

    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        default=DEFAULTS["overwrite"],
        help=("Overwrite the contents of the specified simulation output directory."),
    )

    parser.add_argument(
        "-E",
        "--doEvol",
        action="store_true",
        default=DEFAULTS["doEvol"],
        help=(
            "If used the evolutionary algorithm is executed, the best worm simulation performed,"
            "and results are deposited in the directory."
            "If not used the simulation"
            "in the directory is executed and results deposited in it."
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
        "-r",
        "--randSeed",
        type=int,
        metavar="<rand seed>",
        default=DEFAULTS["randSeed"],
        help="Seed value for evolution and simulation, or just simulation if doEvol is False." 
             "If not set the seed in the directory will be used. If there is no such seed"
             "a random seed will be generated."
        % DEFAULTS["randSeed"],
    )

    return parser.parse_args()


def make_directory(directory_name, overwrite):
    try:
        os.mkdir(directory_name)
        print(f"Directory '{directory_name}' created successfully.")
        return True
    except FileExistsError:
        if overwrite:
            print(
                f"Directory '{directory_name}' already exists and contents will be overwritten."
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

    nml_sim_folder_name = None
    if a.folderName:
        if not make_directory(a.folderName, a.overwrite):
            print(f"Please change phenotype directory name or set overwrite to True.")
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
    
    
    setDict(sim_data, 'seed', a.randSeed, random_seed)
    setDict(sim_data, 'duration', a.duration, 24)
    if a.doEvol:
        setDict(evol_data, 'seed', a.randSeed, random_seed)
        setDict(evol_data, 'popSize', a.popSize, 96)
        setDict(evol_data, 'duration', a.duration, 24)

    with open(evol_par_file, 'w', encoding='utf-8') as f:
        json.dump(evol_data, f, ensure_ascii=False, indent=4)
    with open(sim_par_file, 'w', encoding='utf-8') as f:
        json.dump(sim_data, f, ensure_ascii=False, indent=4)
    if a.doNML:
        with open(nml_sim_folder_name + '/pars.json', 'w', encoding='utf-8') as f:
            json.dump(sim_data, f, ensure_ascii=False, indent=4)  

    if a.doEvol:
        do_evol = 1
    else:
        do_evol = 0
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
