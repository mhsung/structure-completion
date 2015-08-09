#!/usr/bin/python

##########
#
#  Generate HTML table
#
##########

import lib_result_report as librr
import csv
import glob
import math
import numpy as np
import os.path
import shutil
import sys
from collections import namedtuple
from PIL import Image
from enum import Enum


input_path_postfix = '/output/'
output_path_root = '/home/mhsung/app/cuboid-prediction/output/'
output_dir_prefix = ''

gt_n = [1, 0, 0]
gt_t = 0


def load_instances(input_filepath, output_filepath, symemtry_part_index):

    output_file_prefix = os.path.basename(os.path.normpath(input_filepath + '/../'))

    dirnames = glob.glob(input_filepath + '/*')

    ours_file = open(output_file_prefix + '_ours.csv', 'w')
    podolak_file = open(output_file_prefix + '_podolak.csv', 'w')
    
    for dirname in dirnames:
        if not os.path.isdir(dirname):
            continue

        prefix = os.path.basename(dirname)
        print prefix

        is_loaded = True

        candidate_index = librr.find_best_candidate(dirname, prefix)
        print('Candidate index: ' + str(candidate_index))

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

            # Get relative file path.
            relative_image_filepath.append('./' + image_filename)

        if not is_loaded:
            continue

        ours_symm_filename_postfix = '_' + str(candidate_index) + '_symmetry_info.txt'
        ours_symm_filename = dirname + '/' + prefix + ours_symm_filename_postfix

        if os.path.exists(ours_symm_filename):
            with open(ours_symm_filename, 'r') as csv_file:
                data = csv.reader(csv_file, delimiter=',')
                x_data = data.next()
                n = np.array([float(x_data[1]), float(x_data[2]), float(x_data[3])])
                t = float(x_data[4])
                n = n / np.linalg.norm(n)
                ours_file.write(str(n[0]) + ',' + str(n[1]) + ',' + str(n[2]) + ',' + str(t) + '\n')

        podolak_symm_filename_postfix = '_symmetry_info.txt'
        podolak_symm_filename = dirname + '/../../symmetry_detection/' + prefix + '/' + prefix + podolak_symm_filename_postfix
        print(podolak_symm_filename)

        if os.path.exists(podolak_symm_filename):
            with open(podolak_symm_filename, 'r') as csv_file:
                data = csv.reader(csv_file, delimiter=',')
                x_data = data.next()
                n = np.array([float(x_data[1]), float(x_data[2]), float(x_data[3])])
                t = float(x_data[4])
                n = n / np.linalg.norm(n)
                podolak_file.write(str(n[0]) + ',' + str(n[1]) + ',' + str(n[2]) + ',' + str(t) + '\n')

    ours_file.close()
    podolak_file.close()


def main():
    input_path, output_path, dataset_name, symmetry_part_names = librr.parse_arguments(
            input_path_postfix, output_path_root, output_dir_prefix)
    
    load_instances(input_path, output_path, -1)

main()
