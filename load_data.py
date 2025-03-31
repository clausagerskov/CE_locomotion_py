#####
# Script to load generated data on worm motion/cell activity & generate graphical output
#####


import numpy as np
from matplotlib import pyplot as plt
import sys
import argparse

sys.path.append("..")

# import random
import helper_funcs as hf

plot_formats = {}
plot_formats["RS18"] = {}
plot_formats["RS18"]["fig_titles"] = [
    "Stretch receptors",
    "Head Neurons",
    "Body Neurons",
    "Muscles",
]
plot_formats["RS18"]["data_sizes"] = [20, 4, 36, 48]
plot_formats["RS18"]["fig_labels"] = ["SR", "Neu", "Neu", "Mu"]


plot_formats["Net21"] = {}
plot_formats["Net21"]["fig_titles"] = ["Neurons", "Muscles"]
plot_formats["Net21"]["data_sizes"] = [49, 48]
plot_formats["Net21"]["fig_labels"] = ["Neu", "Mu"]
plot_formats["Net21"]["cell_names"] = ["AS", "DA", "DB", "DD", "VD", "VB", "VA"]

plot_formats["CE"] = {}
plot_formats["CE"]["fig_titles"] = ["Stretch receptors", "Neurons", "Muscles"]
plot_formats["CE"]["data_sizes"] = [40, 60, 48]
plot_formats["CE"]["fig_labels"] = ["SR", "Neu", "Mu"]
plot_formats["CE"]["cell_names"] = ["DA", "DB", "DD", "VD", "VA", "VB"]

DEFAULTS = {"modelName": None, "showPlot": True, "folderName": None, "verbose": False}


def process_args():
    """Parse command-line arguments.

    :returns: None
    """
    parser = argparse.ArgumentParser(
        description=("A script for supplying arguments to execute Worm2D")
    )

    parser.add_argument(
        "-m",
        "--modelName",
        type=str,
        metavar="<model name>",
        default=DEFAULTS["modelName"],
        help=(
            "Name of model is required.\nOptions include: RS18, CE, Net21."
            # "Default is: %s" % DEFAULTS["modelName"]
        ),
    )

    parser.add_argument(
        "-s",
        "--showPlot",
        action="store_true",
        # metavar="<run NML>",
        default=DEFAULTS["showPlot"],
        help=("Show plot."),
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        # metavar="<run NML>",
        default=DEFAULTS["verbose"],
        help=("Verbose."),
    )

    parser.add_argument(
        "-f",
        "--folderName",
        type=str,
        metavar="<folder name>",
        default=DEFAULTS["folderName"],
        help=("Required name of data folder."),
    )


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


def run_main(args=None):
    if args is None:
        args = process_args()
    reload_single_run(a=args)


# def reload_single_run(show_plot=True, verbose=False, plot_format_name=None):
def reload_single_run(a=None, **kwargs):
    a = build_namespace(DEFAULTS, a, **kwargs)

    if a.modelName is None:
        print("plot_format is required to make figure.")
        return

    plot_format = plot_formats[a.modelName]

    if a.folderName is None:
        print("Folder name is required for data.")
        return

    hf.dir_name = a.folderName

    # N_muscles_perside = 24  # Number of muscles alongside the body
    # N_muscles = N_muscles_perside * 2
    # N_units = 10  # Number of neural units in VNC
    # N_neuronsperunit = 6  # Number of neurons in a VNC neural unit (6 neurons)
    # N_stretchrec_units = 10  # Number of stretch receptors
    # N_stretchrec = N_stretchrec_units * 4  # Number of stretch receptors

    # N_stretchrec = 2 + 6 * 3  # number of streatch receptors
    # N_hneurons = 4
    # N_vneurons = 36
    # N_muscles = 24 * 2

    # N_neurons = N_neuronsperunit * N_units

    act_data = np.loadtxt(hf.rename_file("act.dat")).T

    def makeFigure(data_offset, data_size, title, label, plot_num):
        axs[plot_num, 0].set_title(title, fontsize=title_font_size)
        axs[plot_num, 1].set_title(title, fontsize=title_font_size)

        for i in range(data_offset, data_size + data_offset):
            axs[plot_num, 0].plot(
                act_data[0][data_seg],
                act_data[i][data_seg],
                label=label + " %i" % (i - offset),
                linewidth=0.5,
            )
            # axs[plot_num, 0].xaxis.set_ticklabels([])
        # plt.legend()

        data_list = act_data[data_offset : data_size + data_offset, data_seg]
        axs[plot_num, 1].imshow(data_list, aspect="auto", interpolation="nearest")
        # axs[plot_num, 1].xaxis.set_ticklabels([])

    t_data = act_data[0]

    fig, axs = plt.subplots(4, 1, figsize=(16, 8))
    t_start = 60
    t_end = 70
    data_seg = (t_data >= t_start) & (t_data < t_end)

    fig, axs = plt.subplots(len(plot_format["fig_titles"]) + 1, 2, figsize=(16, 8))

    title_font_size = 10

    ###  Worm neuron/muscle activation

    t_start = 0
    t_end = 10000
    data_seg = (t_data >= t_start) & (t_data < t_end)

    offset = 1
    count_num = 0
    for val in zip(
        plot_format["data_sizes"], plot_format["fig_titles"], plot_format["fig_labels"]
    ):
        makeFigure(offset, *val, count_num)
        offset += val[0]
        count_num += 1

    ###  Worm body curvature

    curv_data = np.loadtxt(hf.rename_file("curv.dat"))
    curv_data_less_time = curv_data.T[1:, :]

    axs[count_num, 1].set_title("Body curvature", fontsize=title_font_size)
    axs[count_num, 1].imshow(curv_data_less_time, aspect="auto")

    ###  Body position

    body_data = np.loadtxt(hf.rename_file("body.dat")).T

    tmax = 1520
    if tmax >= body_data.shape[1]:
        tmax = body_data.shape[1]
    num = 60.0

    axs[count_num, 0].set_title("2D worm motion", fontsize=title_font_size)

    for t in range(1, tmax, int(tmax / num)):
        f = float(t) / tmax

        color = "#%02x%02x00" % (int(0xFF * (f)), int(0xFF * (1 - f) * 0.8))
        # color2 = "#%06x" % random.randint(0, 0xFFFFFF)

        point_start = 1
        for i in range(point_start, 50):
            x = body_data[i * 3 + 1][t]
            y = body_data[i * 3 + 2][t]
            # y1 = body_data[i * 3 + 2][t]
            if i == 1 and a.verbose:
                print(
                    "%s + Plotting %i at t=%s (%s,%s), %s"
                    % ("\n" if i == point_start else "", i, t, x, y, color)
                )

            axs[count_num, 0].plot(
                [x], [y], ".", color=color, markersize=3 if t == 1 else 0.4
            )

            # print("%s - Plotting %i at t=%s (%s,%s), %s"%('\n' if i==point_start else '', i, t,x,y1, color))
            # plt.plot([x],[y1],'.',color=color)

    axs[count_num, 0].set_aspect("equal")

    filename = hf.rename_file("ExampleActivity.png")
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    print("Saved plot image to: %s" % filename)

    if a.showPlot:
        plt.show()
    plt.close()

    from F2_fig_behavior import make_fig

    make_fig(plot_format=plot_format)


if __name__ == "__main__":
    import sys

    reload_single_run(showPlot=False, modelName="Net21", folderName=sys.argv[1])
