

import configparser
import argparse


def cnf_file_read(cnf_file):
    # Reads the configuration file
    config = configparser.ConfigParser()
    config.read(cnf_file)  # "./data_specs.cnf"

    # Determines the list of available datasets

    num_of_datasets = int(config.get('Available_datasets', 'number'))

    lst_datasets = []
    for n_set in range(0, num_of_datasets):
        set_n = "set_{0:d}".format(n_set)
        lst_datasets.append(config.get('Available_datasets', set_n))

    # Determines the list of available scenarios
    num_of_sc = int(config.get('Scenarios', 'number'))

    lst_scenarios = []
    for n_sc in range(0, num_of_sc):
        sc_n = "scenario_{0:d}".format(n_sc)
        lst_scenarios.append(config.get('Scenarios', sc_n))


    conf_data = {"num_of_datasets": num_of_datasets,
                 "list_of_datasets": lst_datasets,
                 "number_of_scenarios": num_of_sc,
                 "list_of_scenarios": lst_scenarios}
    return (conf_data)


def cnf_file_scenario_select(cnf_file):
    config = configparser.ConfigParser()
    config.read(cnf_file)  # "./analysis.cnf"



    # Read a path to a home folder
    path_home_folder = config.get('Paths', 'home_dir')
    # Read a path to a data folder
    path_data_folder = path_home_folder + config.get('Paths', 'data_dir')

    path_to_data = config.get('Available_datasets', 'set_1') + '/'

    filename = config.get('Scenario_0', 'file_00')

    file_path = path_data_folder + path_to_data + filename

    return (file_path)


def parse_CMDLine(cnf_file):

    conf_data = cnf_file_read(cnf_file)
    # Parses a set of input arguments coming from a command line
    parser = argparse.ArgumentParser(
        description='''
                            Script downloads the data
                            prepared in a dedicated folder according to a
                            pre-defined scenario. Parameters are specified
                            in a configuration file. Scenario has to be
                            selected by an argument.''')
    #      Read command line arguments to get a scenario
    parser.add_argument("-s", "--scenario", help='''Sets an analysis to a given
                                                  scenario. The scenario has to
                                                  be one from an existing ones.''')
    #      Select the radar to process
    parser.add_argument("-l", "--list", action="store_true",
                        help="Prints a list of available scenarios")
    #      Output folder
    parser.add_argument("-o", "--output",
                        help="Sets path to the folder where output files will be stored.")
    #      Select a scenario
    argv = parser.parse_args()


    return (conf_data_out)
