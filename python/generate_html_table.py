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

attr_names = ['Name', 'View_Image',
              'Input_Image', 'Structure_1Ours2_Image',
              'Symmetry3only_1Ours2_Accuracy_Image', 'Symmetry3only_1Ours2_Completeness_Image',
              'Database3only_1Ours2_Accuracy_Image', 'Database3only_1Ours2_Completeness_Image',
              'Complete_1Ours2_Accuracy_Image', 'Complete_1Ours2_Completeness_Image',
              'Symmetry3only_1Ours2_Accuracy', 'Symmetry3only_1Ours2_Completeness',
              'Databass3only_1Ours2_Accuracy', 'Database3only_1Ours2_Completeness',
              'Complete_1Ours2_Accuracy', 'Complete_1Ours2_Completeness',
              'Per3point_Labeling_Accuracy', 'Cuboid_Distance_to_Ground_Truth']

attr_types = [librr.AttrType.text, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.number, librr.AttrType.number,
              librr.AttrType.number, librr.AttrType.number,
              librr.AttrType.number, librr.AttrType.number,
              librr.AttrType.number, librr.AttrType.number]

OutputInstance = namedtuple('OutputInstance', attr_names)


def read_value_from_csv_file(csv_filename, symemtry_part_index, value_index):
    ret = float("NaN")
    with open(csv_filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if symemtry_part_index >= 0:
                if row[0] == str(symemtry_part_index):
                    if row[value_index] == 'none':
                        ret = float("NaN")
                    else:
                        ret = float(row[value_index])
            else:
                if row[0] == 'all':
                    if row[value_index] == 'none':
                        ret = float("NaN")
                    else:
                        ret = float(row[value_index])
    return ret

def load_instances(input_filepath, output_filepath, symemtry_part_index):
    dirnames = glob.glob(input_filepath + '/*')

    instances = []

    for dirname in dirnames:
        if not os.path.isdir(dirname):
            continue

        prefix = os.path.basename(dirname)
        print prefix

        is_loaded = True

        candidate_index = librr.find_best_candidate(dirname, prefix)
        print('Candidate index: ' + str(candidate_index))

        # Read images.
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

            if not os.path.exists(output_filepath + '/' + image_filename):
                # Copy the image.
                shutil.copy(dirname + '/' + image_filename, output_filepath)

                # Create a thumbnail.
                librr.create_thumbnails(output_filepath + '/' + image_filename, librr.thumbname_width)

            # Get relative file path.
            relative_image_filepath.append('./' + image_filename)

        if not is_loaded:
            continue

        # Read stats.
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

            if not os.path.exists(csv_filename):
                print 'Warning: File does not exist: "' + csv_filename + '"'
                is_loaded = False
                break
            else:
                all_values = librr.get_csv_value(csv_filename, librr.threshold)
                if not all_values:
                    accuracy_values.append(float("NaN"))
                    completeness_values.append(float("NaN"))
                else:
                    accuracy_values.append(all_values[0])
                    completeness_values.append(all_values[1])

        if not is_loaded:
            continue

        # Read per-point labeling accuracy.
        csv_filename = dirname + '/' + prefix + '_' + str(candidate_index) + '_labeling_stats.csv'
        if not os.path.exists(csv_filename):
            print 'Warning: File does not exist: "' + csv_filename + '"'
            continue
        else:
            per_point_labeling_accuracy = read_value_from_csv_file(csv_filename, symemtry_part_index, 3)

        # Read cuboid distance to ground truth.
        csv_filename = dirname + '/' + prefix + '_' + str(candidate_index) + '_cuboid_distance.csv'
        if not os.path.exists(csv_filename):
            print 'Warning: File does not exist: "' + csv_filename + '"'
            continue
        else:
            cuboid_distance_to_ground_truth = read_value_from_csv_file(csv_filename, symemtry_part_index, 1)

        instance = OutputInstance(prefix, relative_image_filepath[0],
                                  relative_image_filepath[1], relative_image_filepath[2],
                                  relative_image_filepath[3], relative_image_filepath[4],
                                  relative_image_filepath[5], relative_image_filepath[6],
                                  relative_image_filepath[7], relative_image_filepath[8],
                                  accuracy_values[0], completeness_values[0],
                                  accuracy_values[1], completeness_values[1],
                                  accuracy_values[2], completeness_values[2],
                                  per_point_labeling_accuracy, cuboid_distance_to_ground_truth)

        instances.append(instance)

    return instances


def main():
    input_path, output_path, dataset_name, symmetry_part_names = librr.parse_arguments(
            input_path_postfix, output_path_root, output_dir_prefix)

    instances = load_instances(input_path, output_path, -1)
    html_filename = output_path + '/index.html'
    librr.write_html_table(instances, symmetry_part_names, attr_names, attr_types,
            dataset_name + ' (All)', html_filename)

    # For each part
    for i in range(len(symmetry_part_names)):
        instances = load_instances(input_path, output_path, i)
        html_filename = output_path + '/' + symmetry_part_names[i] + '.html'
        librr.write_html_table(instances, symmetry_part_names, attr_names, attr_types,
                dataset_name + ' (' + symmetry_part_names[i].title() + ')', html_filename)

main()
