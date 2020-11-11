#!usr/local/bin/python3

"""
Program to test and visualize different phasing depths in primers for NGS on Illumina platforms as well as the
impact of different amounts of PhiX library in the sequencing run.

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2020
"""

import os
import random
import dill
import pysam
import argparse
import numpy
import collections
import pandas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from distutils.util import strtobool
from matplotlib.lines import Line2D
from Valkyries import Tool_Box, Version_Dependencies as VersionDependencies, Sequence_Magic

__author__ = 'Dennis A. Simpson'
__version__ = '0.5.0'
__package__ = 'PrimerPhasing'


def error_checking(args):
    """
    Check parameter file for errors.
    :param args:
    """

    if not os.path.exists(args.WorkingFolder):
        print("\033[1;31mERROR:\n\tWorking Folder Path: {} Not Found.  Check Options File."
              .format(args.WorkingFolder))
        raise SystemExit(1)

    if not getattr(args, "RefSeq", False):
        if not os.path.isfile(args.RefSeq):
            print("\033[1;31mERROR:\n\tRefSeq File: {} Not Found.  Check Options File."
                  .format(args.RefSeq))
            raise SystemExit(1)


def string_to_boolean(parser):
    """
    Converts strings to boolean.  Done to keep the eval() function out of the code.
    :param parser:
    :return:
    """

    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()
    options_parser.set_defaults(PairedEnd=bool(strtobool(args.PairedEnd)))
    options_parser.set_defaults(Build_PhiX_DataFrame=bool(strtobool(args.Build_PhiX_DataFrame)))

    args = options_parser.parse_args()

    return args, options_parser


def phix_processing(args):
    """
    Create and pickle dataframe of PhiX reads.
    :param args:
    """
    r1_list = []
    r2_list = []

    # Read in PhiX sequence.  Place each read in a list.
    with open("{}PhiX_R1.txt".format(args.WorkingFolder), 'r') as f:
        for line in f:
            r1_list.append(line.strip("\n"))
    with open("{}PhiX_R2.txt".format(args.WorkingFolder), 'r') as f:
        for line in f:
            r2_list.append(line.strip("\n"))

    # create Pandas row and column indices
    idx = []
    column_label = []

    for i in range(1, 51):
        idx.append(i)
    for i in range(150):
        column_label.append(i)

    # Group PhiX nucleotides by cycle number for each percentage spike
    max_phix_fraction = 50
    phix_fraction = 1
    read_count = len(r1_list)
    read1_phix_dict = collections.defaultdict(list)
    read2_phix_dict = collections.defaultdict(list)

    while phix_fraction <= max_phix_fraction:
        list_count = 0
        r1_row_data = []
        r2_row_data = []
        for i in range(150):
            r1_row_data.append([])
            r2_row_data.append([])

        while list_count < phix_fraction*1000:
            x = random.randint(0, read_count - 1)
            read1 = r1_list[x]
            read2 = r2_list[x]

            for i in range(150):
                r1_row_data[i].append(read1[i])
                r2_row_data[i].append(read2[i])

            list_count += 1

        read1_phix_dict[phix_fraction] = r1_row_data
        read2_phix_dict[phix_fraction] = r2_row_data
        if phix_fraction % 10 == 0:
            print("PhiX Fraction Processed: {}".format(phix_fraction))

        phix_fraction += 1

    r1_dataframe = pandas.DataFrame.from_dict(read1_phix_dict, orient='index')
    r2_dataframe = pandas.DataFrame.from_dict(read2_phix_dict, orient='index')

    print("Begin pickling read 1")
    with open('PhiX_150_R1.df', 'wb') as file:
        dill.dump(r1_dataframe, file, protocol=-1)
    file.close()
    print("Read 1 pickled, begin pickling read 1")
    with open('PhiX_150_R2.df', 'wb') as file:
        dill.dump(r2_dataframe, file, protocol=-1)
    file.close()


def main():
    VersionDependencies.python_check()

    # run_start = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    parser = \
        argparse.ArgumentParser(description="A package to visualize primer phasing.\n {0} v{1}"
                                .format(__package__, __version__),
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    args, options_parser = string_to_boolean(parser)

    # Check options file for errors.
    error_checking(args)

    if args.Build_PhiX_DataFrame:
        phix_processing(args)

    if ":" in args.Amplicon:
        refseq = pysam.FastaFile(args.RefSeq)
        chrm = args.Amplicon.split(":")[0].split("chr")[1]
        start = int(args.Amplicon.split(":")[1].split("-")[0])
        stop = int(args.Amplicon.split(":")[1].split("-")[1])

        target_region = refseq.fetch(chrm, start, stop)
        # Tool_Box.debug_messenger(target_region)

    else:
        target_region = args.Amplicon

    read1_seq = target_region[:int(args.SeqLength)]
    read2_seq = Sequence_Magic.rcomp(target_region[-int(args.SeqLength):])
    read1_phasing_dict = collections.defaultdict(list)
    read2_phasing_dict = collections.defaultdict(list)
    read1_phase_results_dict = collections.defaultdict(list)
    read2_phase_results_dict = collections.defaultdict(list)
    phase_dict = collections.defaultdict(lambda: collections.defaultdict(str))

    # Calculate the total number of reads required
    if int(args.PhiX_Fraction) > 0:
        # Load the pickled dataframes containing the PhiX data
        try:
            with open('{}PhiX_150_R1.df'.format(args.PickleFileFolder), 'rb') as file:
                r1_dataframe = dill.load(file)
            with open('{}PhiX_150_R2.df'.format(args.PickleFileFolder), 'rb') as file:
                r2_dataframe = dill.load(file)

        except FileNotFoundError:
            raise SystemExit("PhiX data not found with --PhiX_Fraction > 0.  Do pickled PhiX dataframes exist?")

        r1_phix = r1_dataframe.loc[int(args.PhiX_Fraction), :].tolist()
        r2_phix = r2_dataframe.loc[int(args.PhiX_Fraction), :].tolist()
        read_number = (100-int(args.PhiX_Fraction))*1000
    else:
        read_number = 1

    # Keep the original copy of the reads.
    r1_seq = read1_seq
    r2_seq = read2_seq

    # Build phasing dictionary.
    for ix in range(len(args.PhasingSeq)):
        read1_phasing_dict.clear()
        read2_phasing_dict.clear()

        # Determine the number of reads for each phase.
        if read_number > 1:
            reads = int(read_number/(ix+2))
        else:
            reads = 1

        # All phasing plots include Phase 0.  This is added here.
        for i in range(150):
            read1_phasing_dict[i].extend([r1_seq[i]] * reads)
            read2_phasing_dict[i].extend([r2_seq[i]] * reads)
            if int(args.PhiX_Fraction) > 0:
                read1_phasing_dict[i].extend(r1_phix[i])
                read2_phasing_dict[i].extend(r2_phix[i])

        # Add the phased data to the dictionaries
        for i in range(int(args.SeqLength)-1):
            if i > 0:
                read1_phasing_dict[i+1].extend([read1_seq[i]]*reads)
                read2_phasing_dict[i+1].extend([read2_seq[i]]*reads)
            else:
                read1_phasing_dict[i].extend([args.PhasingSeq[ix]]*reads)
                read1_phasing_dict[i+1].extend([r1_seq[i]]*reads)

                read2_phasing_dict[i].extend([args.PhasingSeq[ix]]*reads)
                read2_phasing_dict[i+1].extend([r2_seq[i]]*reads)

        # Add the previous phased data to the dictionaries
        if ix > 0:
            for phase in phase_dict:
                phased_read1 = phase_dict[phase]['read1']
                phased_read2 = phase_dict[phase]['read2']
                for cycle in range(int(args.SeqLength)-1):
                    read1_phasing_dict[cycle].extend([phased_read1[cycle]]*reads)
                    read2_phasing_dict[cycle].extend([phased_read2[cycle]]*reads)

        read1_seq = "{}{}".format(args.PhasingSeq[ix], read1_seq)
        read2_seq = "{}{}".format(args.PhasingSeq[ix], read2_seq)

        # Need to keep the previous phased reads.
        phase_dict[ix]['read1'] = read1_seq
        phase_dict[ix]['read2'] = read2_seq

        # Each phase is processed into the final data dictionaries.
        read1_phase_results_dict, read2_phase_results_dict = \
            process_data(read1_phasing_dict, read2_phasing_dict, ix+1, read1_phase_results_dict,
                         read2_phase_results_dict)

    plot_data(read1_phase_results_dict, read2_phase_results_dict, args)


def plot_data(read1_phase_results_dict, read2_phase_results_dict, args):
    """
    Plot the data and save to a pdf file.
    :param read1_phase_results_dict:
    :param read2_phase_results_dict:
    :param args:
    :return:
    """
    def build_axes(phase_results_dict, read1):
        print("Plotting Data")
        fig, axs = plt.subplots(nrows=4, ncols=2, sharex='col')
        for phase in phase_results_dict:
            g = []
            a = []
            t = []
            c = []

            for data in phase_results_dict[phase]:
                g.append(data[0][1] * 100)
                a.append(data[1][1] * 100)
                t.append(data[2][1] * 100)
                c.append(data[3][1] * 100)

            numpy.array(g)
            numpy.array(a)
            numpy.array(t)
            numpy.array(c)

            # Position the plots on the page
            if phase <= 2:
                ax = axs[0, phase - 1]
            elif phase <= 4:
                ax = axs[1, phase - 3]
            elif phase <= 6:
                ax = axs[2, phase - 5]
            else:
                ax = axs[3, phase - 7]

            # Set the scale and font size of each axis.
            ax.set_ylim(0, 100)
            ax.set_xlim(0, 150)
            ax.tick_params(axis='both', which='major', labelsize=10)

            # Label and draw each plot
            ax.set_title("Phase {}".format(phase), fontsize=10)
            # ax.plot(g, lw=0.5, color='blue', label="G")
            ax.plot(g, lw=0.5, color='blue')
            # ax.plot(a, lw=0.5, color='red', label="A")
            ax.plot(a, lw=0.5, color='red')
            # ax.plot(t, lw=0.5, color='black', label="T")
            ax.plot(t, lw=0.5, color='black')
            # ax.plot(t, lw=0.5, color='green', label="C")
            ax.plot(c, lw=0.5, color='green')

        # Define and draw a common figure legend centered at top of page.
        custom = [Line2D([0], [0], color='blue', lw=1),
                  Line2D([0], [0], color='red', lw=1),
                  Line2D([0], [0], color='black', lw=1),
                  Line2D([0], [0], color='green', lw=1)]

        fig.legend(custom, ["G", "A", "T", "C"],
                   fontsize=10,
                   bbox_to_anchor=(0.1, 0.84, 0.80, 0.1),
                   ncol=4,
                   loc='upper center',
                   shadow=True)

        # Write common X-axes label
        fig.text(0.5, 0.07, 'Cycle Number', ha='center', fontsize=12)

        # Write common Y-axes label
        fig.text(0.06, 0.5, 'Read Fraction',
                 horizontalalignment='right',
                 verticalalignment='center',
                 rotation='vertical',
                 fontsize=12)

        if read1:
            fig.suptitle("{} Read 1".format(args.ChartTitle))
        else:
            fig.suptitle("{} Read 2".format(args.ChartTitle))

        fig.set_size_inches(8.5, 11)
        return fig

    outfile = "{}{}.pdf".format(args.WorkingFolder, args.OutFileName)
    # output a multipage PDF.
    with PdfPages(outfile) as pdf:
        pdf.savefig(build_axes(read1_phase_results_dict, read1=True))
        # plt.show()
        plt.close()

        if args.PairedEnd:
            pdf.savefig(build_axes(read2_phase_results_dict, read1=False))
            # plt.show()
            plt.close()


def process_data(read1_phasing_dict, read2_phasing_dict, phase, read1_phase_results_dict, read2_phase_results_dict):
    """
    Get the data into a format that can be easily plotted.
    :param read1_phasing_dict:
    :param read2_phasing_dict:
    :param phase:
    :param read1_phase_results_dict:
    :param read2_phase_results_dict:
    :return:
    """
    for (key1, value1), (key2, value2) in zip(read1_phasing_dict.items(), read2_phasing_dict.items()):
        depth1 = len(value1)
        depth2 = len(value2)

        read1_phase_results_dict[phase].append([("G", value1.count('G')/depth1), ("A", value1.count('A')/depth1),
                                                ("T", value1.count('T')/depth1), ("C", value1.count('C')/depth1)])

        read2_phase_results_dict[phase].append([("G", value2.count('G')/depth2), ("A", value2.count('A')/depth2),
                                                ("T", value2.count('T')/depth2), ("C", value2.count('C')/depth2)])
    return read1_phase_results_dict, read2_phase_results_dict


if __name__ == '__main__':
    main()
