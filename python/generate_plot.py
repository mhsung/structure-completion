#!/usr/bin/python

##########
#
#  Generate HTML table
#
##########

import lib_result_report as librr
import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import os.path
import shutil
import sys
from collections import namedtuple
from PIL import Image
from enum import Enum


max_x_value = 0.05
input_path_postfix = '/output/'
output_path_root = '/home/mhsung/app/cuboid-prediction/report/'
output_dir_prefix = ''


def load_instances(input_filepath, output_filepath):
    dirnames = glob.glob(input_filepath + '/*')

    all_accu_values = [[], [], [], []]
    all_comp_values = [[], [], [], []]
    x_values = []
    count_instances = 0

    for dirname in dirnames:
        if not os.path.isdir(dirname):
            continue

        prefix = os.path.basename(dirname)
        #print prefix

        is_loaded = True

        candidate_index = librr.find_best_candidate(dirname, prefix)
        #print('Candidate index: ' + str(candidate_index))

        relative_image_filepath = []
        image_filenames = []
        image_filenames.append(prefix + '_view.png')
        image_filenames.append(prefix + '_input.png')
        image_filenames.append(prefix + '_' + str(candidate_index) + '.png')
        image_filenames.append(prefix + '_' + str(candidate_index) + '_symmetry_accuracy.png')
        image_filenames.append(prefix + '_' + str(candidate_index) + '_symmetry_completeness.png')
        image_filenames.append(prefix + '_' + str(candidate_index) + '_database_accuracy.png')
        image_filenames.append(prefix + '_' + str(candidate_index) + '_database_completeness.png')
        image_filenames.append(prefix + '_' + str(candidate_index) + '_fusion_accuracy.png')
        image_filenames.append(prefix + '_' + str(candidate_index) + '_fusion_completeness.png')

        for image_filename in image_filenames:
            if not os.path.exists(dirname + '/' + image_filename):
                print 'Warning: File does not exist: "' + (dirname + '/' + image_filename) + '"'
                is_loaded = False
                break

        if not is_loaded:
            continue

        csv_filename_postfixes = []
        csv_filename_postfixes.append(prefix + '_' + str(candidate_index) + '_symmetry')
        csv_filename_postfixes.append(prefix + '_' + str(candidate_index) + '_database')
        csv_filename_postfixes.append(prefix + '_' + str(candidate_index) + '_fusion')
        csv_filename_postfixes.append('/../../part_assembly/' + prefix + '/' + prefix +'_assembly')
        

        for i in range(len(csv_filename_postfixes)):
            csv_filename_postfix = csv_filename_postfixes[i]
            csv_filename = dirname + '/' + csv_filename_postfix + '.csv'
            (accu_values, comp_values, x_values) = librr.get_csv_all_value(csv_filename)

            if not all_accu_values[i]:
                all_accu_values[i] = accu_values
            if not all_comp_values[i]:
                all_comp_values[i] = comp_values
                
            all_accu_values[i] = [x+y for x,y in zip(all_accu_values[i], accu_values)]
            all_comp_values[i] = [x+y for x,y in zip(all_comp_values[i], comp_values)]

        count_instances += 1

    if count_instances > 0:
        for i in range(len(all_accu_values)):
            all_accu_values[i] = [x / count_instances for x in all_accu_values[i]]
        for i in range(len(all_comp_values)):
            all_comp_values[i] = [x / count_instances for x in all_comp_values[i]]

    return (all_accu_values, all_comp_values, x_values)


def main():
    input_path, output_path, dataset_name, symmetry_part_names = librr.parse_arguments(
            input_path_postfix, output_path_root, output_dir_prefix)
    
    data_dirname = os.path.basename(os.path.normpath(sys.argv[1]))
    print(" >> " + data_dirname)

    (all_accu_values, all_comp_values, x_values) = load_instances(input_path, output_path)

    title = 'accu'
    plt.plot(x_values, all_accu_values[0], label='Symmetry')
    plt.plot(x_values, all_accu_values[1], label='Database')
    plt.plot(x_values, all_accu_values[2], label='Fusion')
    plt.plot(x_values, all_accu_values[3], label='Part Assembly')
    plt.legend(loc=4)
    plt.xlim(0, max_x_value)
    plt.ylim(0, 1)
    plt.xlabel('Neighbor Distance')
    plt.savefig(output_path_root + '/' + data_dirname + "_" + title)
    plt.clf()
    
    title = 'comp'
    plt.plot(x_values, all_comp_values[0], label='Symmetry')
    plt.plot(x_values, all_comp_values[1], label='Database')
    plt.plot(x_values, all_comp_values[2], label='Fusion')
    plt.plot(x_values, all_comp_values[3], label='Part Assembly')
    plt.legend(loc=4)
    plt.xlim(0, max_x_value)
    plt.ylim(0, 1)
    plt.xlabel('Neighbor Distance')
    plt.savefig(output_path_root + '/' + data_dirname + "_" + title)

main()
