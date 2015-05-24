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


input_path_postfix = '/part_assembly/'
output_path_root = '/home/mhsung/app/cuboid-prediction/output/'
output_dir_prefix = 'assembly_'

attr_names = ['Name', 'View_Image',  'Input_Image',
              'Assembly_Accuracy_Image', 'Assembly_Completeness_Image',
              'Assembly_Accuracy', 'Assembly_Completeness']

attr_types = [librr.AttrType.text, librr.AttrType.image, librr.AttrType.image,
              librr.AttrType.image, librr.AttrType.image,
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

        is_loaded = True

        relative_image_filepath = []
        image_filenames = []
        image_filenames.append(prefix + '_view.png')
        image_filenames.append(prefix + '_input.png')
        image_filenames.append(prefix + '_assembly_accuracy.png')
        image_filenames.append(prefix + '_assembly_completeness.png')

        # FIXME
        if not os.path.exists(dirname + '/../../output/' + prefix + '/' + prefix + '_input.png'):
            continue

        shutil.copy(dirname + '/../../output/' + prefix + '/' + prefix + '_input.png',
                dirname + '/' + prefix + '_input.png')


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

        accuracy_values = []
        completeness_values = []
        csv_filename_postfixes = []
        csv_filename_postfixes.append('_assembly')

        for csv_filename_postfix in csv_filename_postfixes:
            csv_filename = dirname + '/' + prefix + csv_filename_postfix + '.csv'
            if symemtry_part_index >= 0:
            # Real per-part files.
                csv_filename = dirname + '/' + prefix + csv_filename_postfix\
                               + '_' + str(symemtry_part_index) + '.csv'

            all_values = librr.get_csv_value(csv_filename, librr.threshold)

            if not all_values:
                accuracy_values.append(float("NaN"))
                completeness_values.append(float("NaN"))
            else:
                accuracy_values.append(all_values[0])
                completeness_values.append(all_values[1])

        instance = OutputInstance(prefix, relative_image_filepath[0], relative_image_filepath[1],
                                  relative_image_filepath[2], relative_image_filepath[3],
                                  accuracy_values[0], completeness_values[0])

        instances.append(instance)

    return instances


def main():
    input_path, output_path, dataset_name, symmetry_part_names = librr.parse_arguments(
            input_path_postfix, output_path_root, output_dir_prefix)

    instances = load_instances(input_path, output_path, -1)
    html_filename = output_path + '/index.html'
    librr.write_html_table(instances, attr_names, attr_types, dataset_name + ' (All)', html_filename)

    # For each part
    for i in range(len(symmetry_part_names)):
        instances = load_instances(input_path, output_path, i)
        html_filename = output_path + '/' + symmetry_part_names[i] + '.html'
        librr.write_html_table(instances, attr_names, attr_types,
                dataset_name + ' (' + symmetry_part_names[i].title() + ')', html_filename)

main()
