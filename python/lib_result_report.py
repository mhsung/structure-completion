#!/usr/bin/python

##########
#
#  Generate HTML table
#
##########

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


# Parameters
threshold = 0.02
thumbname_width = 150


##
dataset_name = ''
symmetry_part_names = ''

dataset_name_set = [
        'assembly_chairs',
        'assembly_airplanes',
        'assembly_bicycles',
        'coseg_chairs']
symmetry_part_names_set = [
        ['seat', 'back', 'legs', 'wheels', 'leg_column', 'armrests'],
        ['body', 'wings', 'tail_wings', 'fuselages', 'body'],
        ['wheel', 'handle', 'paddle', 'front_column', 'rear_body', 'front_body', 'seat', 'chain'],
        ['seat', 'back', 'legs', 'wheels', 'leg_column', 'armrests']]
##

class AttrType(Enum):
    text = 1
    image = 2
    number = 3


# Functions
def get_csv_value(input_filepath, threshold):
    values = []

    if os.path.exists(input_filepath):
        with open(input_filepath, 'r') as csv_file:
            data = csv.reader(csv_file, delimiter=',')

            x_data = data.next()

            threshold_index = 0
            while threshold_index < len(x_data) and float(x_data[threshold_index]) < threshold:
                threshold_index += 1

            for y_data in data:
                values.append(float(y_data[threshold_index]))

    return values


def get_csv_all_value(input_filepath):
    values = []

    if os.path.exists(input_filepath):
        with open(input_filepath, 'r') as csv_file:
            data = csv.reader(csv_file, delimiter=',')

            x_values = map(float, data.next())
            accu_values = map(float, data.next())
            comp_values = map(float, data.next())

    return (accu_values, comp_values, x_values)


def create_thumbnails(image_filepath, basewidth):
    img = Image.open(image_filepath)
    wpercent = (basewidth/float(img.size[0]))
    hsize = int((float(img.size[1])*float(wpercent)))
    img = img.resize((basewidth, hsize), Image.ANTIALIAS)

    filename, fileext = os.path.splitext(image_filepath)
    thumbname_filename = filename + '_thumb' + fileext
    img.save(thumbname_filename)


def find_best_candidate(dirname, prefix):
    candidate_index = 0

    max_accuracy_value = 0
    best_candidate_index = 0

    while True:
        csv_filename_postfix = '_' + str(candidate_index) + '_symmetry'
        csv_filename = dirname + '/' + prefix + csv_filename_postfix + '.csv'
        if not os.path.exists(csv_filename):
            break

        all_values = get_csv_value(csv_filename, threshold)
        accuracy_value = all_values[0]

        if accuracy_value > max_accuracy_value:
            max_accuracy_value = accuracy_value
            best_candidate_index = candidate_index

        candidate_index += 1

    return best_candidate_index


def write_html_header(file):
    file.write('<head>\n')
    file.write('<script src="//code.jquery.com/jquery-1.11.2.min.js"></script>\n')
    file.write('<script src="//code.jquery.com/jquery-migrate-1.2.1.min.js"></script>\n')
    file.write('<script src="//cdn.datatables.net/1.10.6/js/jquery.dataTables.min.js"></script>\n')
    file.write('<script src="//datatables.net/release-datatables/extensions/ColVis/js/dataTables.colVis.js"></script>\n')
    file.write('<link rel="stylesheet" type="text/css" href="//cdn.datatables.net/1.10.6/css/jquery.dataTables.css">\n')
    file.write('<link rel="stylesheet" type="text/css" href="//datatables.net/release-datatables/extensions/ColVis/css/dataTables.colVis.css">\n')
    file.write('\n')
    file.write('<script type="text/javascript">\n')
    file.write('$(document).ready(function() {\n')

    file.write('\tvar table = $(\'\#example\').DataTable( {\n')
    file.write('\t\tdom: \'C<"clear">lfrtip\'\n')
    file.write('\t} );\n')
    file.write('\n')
    file.write('\t$(\'a.toggle-vis\').on( \'click\', function (e) {\n')
    file.write('\t\te.preventDefault();\n')
    file.write('\n')
    file.write('\t// Get the column API object\n')
    file.write('\tvar column = table.column( $(this).attr(\'data-column\') );\n')
    file.write('\n')
    file.write('\t// Toggle the visibility\n')
    file.write('\tcolumn.visible( ! column.visible() );\n')
    file.write('\t} );\n')
    file.write('} );\n')
    file.write('</script>\n')
    file.write('</head>\n')
    file.write('\n')


def write_html_image(file, image_filepath):
    filename, fileext = os.path.splitext(image_filepath)
    thumbname_filename = filename + '_thumb' + fileext

    file.write('<a href="' + image_filepath + '" target="_blank">')
    file.write('<img src="' + thumbname_filename + '" width="' + str(thumbname_width) + '">')
    file.write('</a>')


def write_html_part_links(file):
    file.write('Per-part Statistics: \n')
    file.write('<a href="./index.html"> All </a>\n')
    for symmetry_part_name in symmetry_part_names:
        file.write(' - <a href="./' + symmetry_part_name + '.html"> ' + symmetry_part_name.title() + ' </a>\n')
    file.write('<p></p>\n')


def write_html_overall_stats(file, instances, attr_names, attr_types):
    assert len(attr_names) == len(attr_types)
    if len(instances) == 0:
        return

    file.write('<b>Distance threshold</b>: ' + str(threshold) + '<p></p>\n')

    num_attrs = len(attr_names)
    num_instances = len(instances)
    sum_attr_values = [0] * num_attrs

    file.write('<table border="1" cellpadding="4" cellspacing="0">\n')

    file.write('<thead>\n')
    file.write('<tr>\n')
    file.write('<th></th> ')
    file.write('<th>Mean</th> ')
    file.write('<th>Median</th> ')
    file.write('<th>Stdev</th>')
    file.write('</tr>\n')
    file.write('</thead>\n')

    file.write('<tbody>\n')
    for i in range(num_attrs):
        if attr_types[i] == AttrType.number:
            all_attr_values = []
            for instance in instances:
                if not math.isnan(instance[i]):
                    all_attr_values.append(float(instance[i]))

            mean_value = np.average(np.array(all_attr_values))
            median_value = np.median(np.array(all_attr_values))
            stdev_value = np.std(np.array(all_attr_values))

            attr_name = attr_names[i].replace('_', ' ')
            file.write('<td>' + attr_name + '</td>\n')
            file.write('<td>' + "{0:0.3f}".format(mean_value) + '</td>\n')
            file.write('<td>' + "{0:0.3f}".format(median_value) + '</td>\n')
            file.write('<td>' + "{0:0.3f}".format(stdev_value) + '</td>\n')
            file.write('</tr>\n')
    file.write('</tbody>\n')

    file.write('</table>\n')
    file.write('<p></p>\n')


def write_html_table(instances, attr_names, attr_types, title, filename):
    assert len(attr_names) == len(attr_types)
    num_attrs = len(attr_names)

    print "Saving the file..."
    file = open(filename, 'w')

    file.write('<html>\n')
    write_html_header(file)

    file.write('<body>\n')
    file.write('<h2>' + title + '</h2>\n')
    write_html_part_links(file)
    write_html_overall_stats(file, instances, attr_names, attr_types)


    file.write('<table id="example" class="cell-border" cellspacing="0" width="100%">\n')

    file.write('<thead>\n')
    file.write('<tr>\n')
    for attr_name in attr_names:
        attr_name = attr_name.replace('_', '<br>')
        file.write('<th>'); file.write(attr_name); file.write('</th>\n')
    file.write('</tr>\n')
    file.write('</thead>\n')

    file.write('<tfoot>\n')
    file.write('<tr>\n')
    for attr_name in attr_names:
        attr_name = attr_name.replace('_', '<br>')
        file.write('<th>'); file.write(attr_name); file.write('</th>\n')
    file.write('</tr>\n')
    file.write('</tfoot>\n')

    file.write('<tbody>\n')
    for instance in instances:
        file.write('<tr>\n')
        for i in range(num_attrs):
            file.write('<td align="center">')
            if attr_types[i] == AttrType.text or attr_types[i] == AttrType.number:
                file.write(str(instance[i]))
            elif attr_types[i] == AttrType.image:
                write_html_image(file, instance[i])
            file.write('</td>\n')
        file.write('</tr>\n')
    file.write('</tbody>\n')

    file.write('</table>\n')
    file.write('</body>\n')


def write_latex_header(file):
    file.write('\\documentclass[]{article}\n')
    file.write('\\usepackage{graphicx}\n')
    file.write('\\usepackage[margin=0.10in]{geometry}\n')
    file.write('\n')


def write_latex_image(file, image_size, image_filepath):
    filename, fileext = os.path.splitext(image_filepath)
    #thumbname_filename = filename + '_thumb' + fileext
    file.write('\\begin{minipage}{' + str(image_size) + '\\textwidth} \\includegraphics[width=\linewidth, clip, trim={50 -10 50 -10}]{' + image_filepath + '} \\end{minipage} \n')


def open_latex_table(filename, attr_names, attr_types):
    assert len(attr_names) == len(attr_types)
    num_attrs = len(attr_names)

    num_image_attrs = 0
    for i in range(num_attrs):
        if attr_types[i] == AttrType.image:
            num_image_attrs += 1

    print "Saving the file..."
    file = open(filename, 'w')

    write_latex_header(file)
    file.write('\\begin{document}\n')
    file.write('\\begin{table}[h]\n')
    file.write('\\begin{tabular}{')
    for count_image_attrs in range(num_image_attrs):
        file.write('c')
        if (count_image_attrs + 1) < num_image_attrs:
            file.write('|')
    file.write('}\n')

    '''
    count_image_attrs = 0
    for i in range(num_attrs):
        if attr_types[i] == AttrType.image:
            attr_name = attr_names[i].replace('_', ' ')
            file.write(attr_name)
            if (count_image_attrs + 1) < num_image_attrs:
                file.write(' &\n')
                count_image_attrs += 1
    file.write('\\\\ \\hline\n')
    '''

    return file;


def write_latex_table(file, instances, attr_names, attr_types):
    assert len(attr_names) == len(attr_types)
    num_attrs = len(attr_names)

    num_image_attrs = 0
    for i in range(num_attrs):
        if attr_types[i] == AttrType.image:
            num_image_attrs += 1

    image_size = 1 / (1.12 * num_image_attrs)

    num_instances = len(instances)
    for instance_id in range(num_instances):
        instance = instances[instance_id]

        count_image_attrs = 0
        for i in range(num_attrs):
            #if attr_types[i] == AttrType.text or attr_types[i] == AttrType.number:
            #    file.write(str(instance[i]))
            if attr_types[i] == AttrType.image:
                write_latex_image(file, image_size, instance[i])
                if (count_image_attrs + 1) < num_image_attrs:
                    file.write(' &\n')
                count_image_attrs += 1

        if (instance_id + 1) < num_instances:
            file.write('\\\\ \\hline\n')


def close_latex_table(file):
    file.write('\\end{tabular}\n')
    file.write('\\end{table}\n')
    file.write('\\end{document}\n')


def parse_arguments(input_path_postfix, output_path_root, output_dir_prefix):
    ## Parse arguments
    if len(sys.argv) < 2:
        print("Usage: " + sys.argv[0] + " experiment_path")
        exit()
    else:
        data_path = sys.argv[1]
        data_dirname = os.path.basename(os.path.normpath(data_path))
        input_path = data_path + input_path_postfix
        output_path = output_path_root + output_dir_prefix + os.path.basename(data_dirname)
        print(input_path)
        print(output_path)
        
        for i in range(len(dataset_name_set)):
            if dataset_name_set[i] in data_dirname:
                dataset_name = dataset_name_set[i]
                symmetry_part_names = symmetry_part_names_set[i]
                break

        print(dataset_name)
        print(symmetry_part_names)
        if dataset_name == '':
            print('Error: No dataset name exist.')
            exit()

    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    return [input_path, output_path, dataset_name, symmetry_part_names]

