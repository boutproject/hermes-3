
# Data consists of
#    (path, delta*, setting changes)
#
# These are arranged in approximate order of importance, so
# the list can be truncated
dataset = [('delta-1', 1.0, [('nout', 15)])
           ,('delta-0.5', 0.5, [('nout', 10),
                                ('mesh:Lrad', 0.025),
                                ('mesh:Lpol', 0.025)])
           ,('delta-0.25', 0.25, [('nout', 10),
                                  ('mesh:Lrad', 0.0125),
                                  ('mesh:Lpol', 0.0125)])
           ,('delta-8', 8.0, [('nout', 20),
                              ('mesh:Lrad', 0.2),
                              ('mesh:Lpol', 0.2),
                              ('Ne:width', 0.1)])
           ,('delta-16', 16.0, [('nout', 20),
                                ('mesh:Lrad', 0.4),
                                ('mesh:Lpol', 0.4),
                                ('Ne:width', 0.1)])
           ,('delta-0.125', 0.125, [('nout', 10),
                                    ('timestep', 25),
                                    ('mesh:Lrad', 0.00625),
                                    ('mesh:Lpol', 0.00625)])
           ,('delta-0.35', 0.35, [('nout', 10),
                                  ('mesh:Lrad', 0.0175),
                                  ('mesh:Lpol', 0.0175)])
           ,('delta-2', 2.0, [('nout', 20),
                              ('mesh:Lrad', 0.1),
                              ('mesh:Lpol', 0.1)])
           ,('delta-2.8', 2.8, [('nout', 20),
                                ('mesh:Lrad', 0.14),
                                ('mesh:Lpol', 0.14)])
           ,('delta-4', 4.0, [('nout', 15),
                              ('mesh:Lrad', 0.1),
                              ('mesh:Lpol', 0.1),
                              ('Ne:width', 0.1)])
           ,('delta-11.2', 11.2, [('nout', 20),
                                  ('mesh:Lrad', 0.2),
                                  ('mesh:Lpol', 0.2),
                                  ('Ne:width', 0.14)])
           ,('delta-5.6', 5.6, [('nout', 15),
                                ('mesh:Lrad', 0.1),
                                ('mesh:Lpol', 0.1),
                                ('Ne:width', 0.14)])
           ]

from boutdata.data import BoutOptionsFile
import os


def create_inputs(origin, destination, datasets):
    """
    Parameters
    ----------
    origin : str
      The directory with the starting BOUT.inp file
    destination : str
      Directory where new inputs will be created
    datasets : list
      Structured data specifying the directories and settings
      [('path', [('setting', value),...]), ...]
    """

    print("\nGenerating inputs")

    for path, _, changes in datasets:
        newdir = os.path.join(destination, path)
        print("  Input directory {}".format(newdir))

        # Read the original options
        options = BoutOptionsFile(os.path.join(origin, "BOUT.inp"))

        # Make the modifications
        for setting, value in changes:
            options[setting] = value

        # Save to new location
        try:
            os.mkdir(newdir)
            print("   -> Created new directory")
            options.write(os.path.join(newdir, "BOUT.inp"))
        except FileExistsError:
            print("   -> Already exists. Skipping")


def run_cases(executable, destination, dataset):
    """
    Parameters
    ----------

    executable : str
      The path to the executable (including name of executable)
    destination : str
      Directory containing the run inputs as subdirectories
    dataset : list
    """

    from boututils.run_wrapper import launch

    print("\nRunning simulations")

    for path, _, _ in dataset:
        rundir = os.path.join(destination, path)

        print("  Input directory {}".format(rundir))

        if os.path.exists(os.path.join(rundir, "BOUT.dmp.0.nc")):
            print("   -> Already exists. Skipping")
            continue

        print("   -> Running executable {}".format(executable))
        launch(
            executable + " -d " + rundir,
            nproc=1,
            output=os.path.join(rundir, "run.log"),
        )


def analyse_cases(destination, dataset):
    """
    Analyse the output from previous runs

    Parameters
    ----------

    destination : str
      Directory containing the runs as subdirectories
    dataset : list
      List of tuples with first element as the subdirectory
    """

    import numpy as np
    from boutdata.collect import collect
    import blob_velocity

    deltas = []
    maxvelocity = []

    print("\nAnalysing data")

    for path, delta, _ in dataset:
        datadir = os.path.join(destination, path)

        print("  Input directory {}".format(datadir))

        vx, vz, _, _ = blob_velocity.calculate_SI(datadir)

        deltas.append(delta)
        maxvelocity.append(np.amax(vx))

    return np.array(deltas), np.array(maxvelocity)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--exec", help="The Hermes-3 executable to run", default="../../hermes-3"
    )
    parser.add_argument(
        "--input",
        help="The directory containing the BOUT.inp file to start from",
        default=".",
    )
    parser.add_argument(
        "--workdir",
        help="The directory to use. Subdirectories will be generated for each run.",
        default=".",
    )
    parser.add_argument(
        "--nruns",
        help="The number of runs to perform",
        type=int,
        default=len(dataset),
        choices=range(2, len(dataset) + 1),
    )
    args = parser.parse_args()

    # Select the number of runs needed. Later runs tend to fill in gaps between earlier runs
    mydataset = dataset[: args.nruns]

    create_inputs(args.input, args.workdir, mydataset)

    run_cases(args.exec, args.workdir, mydataset)

    deltas, maxv = analyse_cases(args.workdir, mydataset)

    # Estimate exponents from gradients
    # v = v0 * d^a
    # log(v) = log(v0) + a * log(d)
    # a = log(v2 / v1) / log(d2 / d1)

    import numpy as np

    sortind = np.argsort(deltas)
    deltas = deltas[sortind]
    maxv = maxv[sortind]

    alow = np.log(maxv[1] / maxv[0]) / np.log(deltas[1] / deltas[0])
    ahigh = np.log(maxv[-1] / maxv[-2]) / np.log(deltas[-1] / deltas[-2])

    print("Exponents: {}, {}".format(alow, ahigh))

    try:
        import matplotlib.pyplot as plt

        plt.plot(deltas, maxv, "ko")
        plt.yscale("log")
        plt.xscale("log")

        plt.plot(
            [deltas[0], deltas[-1]],
            [maxv[0], maxv[0] * (deltas[-1] / deltas[0]) ** alow],
            "--k",
        )
        plt.text(deltas[1], maxv[1] * 1.4, "{:.2f}".format(alow))

        plt.plot(
            [deltas[0], deltas[-1]],
            [maxv[-1] * (deltas[0] / deltas[-1]) ** ahigh, maxv[-1]],
            "--k",
        )
        plt.text(deltas[-2], maxv[-2] * 1.2, "{:.2f}".format(ahigh))

        plt.xlabel(r"Blob size $\delta^*$")
        plt.ylabel(r"Radial velocity [m/s]")

        plt.ylim(top=1000)

        plt.savefig("blob_velocity.png")
        plt.savefig("blob_velocity.pdf")

        plt.show()
    except ImportError:
        pass
