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
import os.path
import shutil
import sys
from collections import namedtuple
from PIL import Image
from enum import Enum


input_path_postfix = '/output/'
output_path_root = '/home/mhsung/app/cuboid-prediction/report/'
output_dir_prefix = ''

attr_names = ['Name', 'View_Image',
              'Input_Image', 'Structure_Reconstruction_Image',
              'Symmetry_Accuracy_Image', 'Symmetry_Completeness_Image',
              'Database_Accuracy_Image', 'Database_Completeness_Image',
              'Fusion_Accuracy_Image', 'Fusion_Completeness_Image',
              'Symmetry_Accuracy', 'Symmetry_Completeness',
              'Database_Accuracy', 'Database_Completeness',
              'Fusion_Accuracy', 'Fusion_Completeness']

attr_types = [librr.AttrType.text, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.number, librr.AttrType.number,
              librr.AttrType.number, librr.AttrType.number,
              librr.AttrType.number, librr.AttrType.number]

OutputInstance = namedtuple('OutputInstance', attr_names)


def load_instances(input_filepath, output_filepath, symemtry_part_index):
    dirnames = glob.glob(input_filepath + '/*')

    instances = []

    for dirname in dirnames:
        if not os.path.isdir(dirname):
            continue

        prefix = os.path.basename(dirname)
        print prefix

        abs_dirname = os.path.abspath(dirname)

        is_loaded = True

        candidate_index = librr.find_best_candidate(dirname, prefix)
        print('Candidate index: ' + str(candidate_index))

        absolute_image_filepath = []
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

            if not os.path.exists(output_filepath + '/' + image_filename):
                '''
                # Copy the image.
                shutil.copy(dirname + '/' + image_filename, output_filepath)

                # Create a thumbnail.
                librr.create_thumbnails(output_filepath + '/' + image_filename, librr.thumbname_width)
                '''

            # Get relative file path.
            absolute_image_filepath.append(abs_dirname + '/' + image_filename)

        if not is_loaded:
            continue

        accuracy_values = []
        completeness_values = []
        csv_filename_postfixes = []
        csv_filename_postfixes.append('_' + str(candidate_index) + '_symmetry')
        csv_filename_postfixes.append('_' + str(candidate_index) + '_database')
        csv_filename_postfixes.append('_' + str(candidate_index) + '_fusion')

        for csv_filename_postfix in csv_filename_postfixes:
            csv_filename = dirname + '/' + prefix + csv_filename_postfix + '.csv'
            if symemtry_part_index >= 0:
                # Read per-part files.
                csv_filename = dirname + '/' + prefix + csv_filename_postfix\
                               + '_' + str(symemtry_part_index) + '.csv'

            all_values = librr.get_csv_value(csv_filename, librr.threshold)

            if not all_values:
                accuracy_values.append(float("NaN"))
                completeness_values.append(float("NaN"))
            else:
                accuracy_values.append(all_values[0])
                completeness_values.append(all_values[1])

        instance = OutputInstance(prefix, absolute_image_filepath[0],
                                  absolute_image_filepath[1], absolute_image_filepath[2],
                                  absolute_image_filepath[3], absolute_image_filepath[4],
                                  absolute_image_filepath[5], absolute_image_filepath[6],
                                  absolute_image_filepath[7], absolute_image_filepath[8],
                                  accuracy_values[0], completeness_values[0],
                                  accuracy_values[1], completeness_values[1],
                                  accuracy_values[2], completeness_values[2])

        instances.append(instance)

    return instances


def main():
    input_path, output_path, dataset_name, symmetry_part_names = librr.parse_arguments(
            input_path_postfix, output_path_root, output_dir_prefix)

    instances = load_instances(input_path, output_path, -1)
    html_filename = output_path + '/output.tex'
    librr.write_latex_table(instances, attr_names, attr_types, dataset_name + ' (All)', html_filename)

    # For each part
    '''
    for i in range(len(symmetry_part_names)):
        instances = load_instances(input_path, output_path, i)
        html_filename = output_path + '/' + symmetry_part_names[i] + '.html'
        librr.write_html_table(instances, attr_names, attr_types,
                dataset_name + ' (' + symmetry_part_names[i].title() + ')', html_filename)
    '''

main()
