#####
# Script to load generated data on worm motion/cell activity & generate graphical output
#####


import numpy as np
from matplotlib import pyplot as plt
import sys

sys.path.append("..")

# import random
import helper_funcs as hf


def reload_single_run(show_plot=True, verbose=False, plot_format=None):
    if plot_format is None:
        print("plot_format is required to make figure.")
        return

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

    fig, axs = plt.subplots(len(plot_format["fig_titles"]) + 1, 2, figsize=(16, 8))

    title_font_size = 10

    ###  Worm neuron/muscle activation

    act_data = np.loadtxt(hf.rename_file("act.dat")).T

    def makeFigure(data_offset, data_size, title, label, plot_num):
        data_list = act_data[data_offset : data_size + data_offset]
        axs[plot_num, 0].set_title(title, fontsize=title_font_size)
        axs[plot_num, 1].set_title(title, fontsize=title_font_size)

        for i in range(data_offset, data_size + data_offset):
            axs[plot_num, 0].plot(
                act_data[0],
                act_data[i],
                label=label + " %i" % (i - offset),
                linewidth=0.5,
            )
            axs[plot_num, 0].xaxis.set_ticklabels([])
        # plt.legend()

        axs[plot_num, 1].imshow(data_list, aspect="auto", interpolation="nearest")
        axs[plot_num, 1].xaxis.set_ticklabels([])

    # fig_titles = ["Stretch receptors", "Head Neurons", "Body Neurons", "Muscles"]
    # data_sizes = [N_stretchrec, N_hneurons, N_vneurons, N_muscles]
    # fig_labels = ["SR", "Neu", "Neu", "Mu"]

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
            if i == 1 and verbose:
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

    if show_plot:
        plt.show()
    plt.close()


if __name__ == "__main__":
    import sys

    reload_single_run(show_plot="-nogui" not in sys.argv)
