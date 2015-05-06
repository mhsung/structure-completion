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
from collections import namedtuple
from PIL import Image
from enum import Enum


# Parameters
threshold = 0.001
thumbname_width = 600


#
dataset_name = 'Assembly Chairs'
symmetry_part_names = ['seat', 'back', 'legs', 'wheels', 'leg_column', 'armrests']
input_path = '/home/mhsung/app/cuboid-prediction/experiments/exp8_assembly_chairs/output'
output_path = '/home/mhsung/Dropbox/Public/web/exp8_assembly_chairs'
name_filters = ['4049']

#dataset_name = 'Assembly Airplanes'
#symmetry_part_names = ['body', 'wings', 'tail_wings', 'fuselages', 'body']
#input_path = '/home/mhsung/app/cuboid-prediction/experiments/exp1_assembly_airplanes/output'
#output_path = '/home/mhsung/Dropbox/Public/web/exp1_assembly_airplanes'

#
class AttrType(Enum):
    text = 1
    image = 2
    number = 3

attr_names = ['Name', 'Input_Image', 'Structure_Reconstruction_Image',
              'Symmetry_Accuracy_Image', 'Symmetry_Completeness_Image',
              'Database_Accuracy_Image', 'Database_Completeness_Image',
              'Fusion_Accuracy_Image', 'Fusion_Completeness_Image',
              'Symmetry_Accuracy', 'Symmetry_Completeness',
              'Database_Accuracy', 'Database_Completeness',
              'Fusion_Accuracy', 'Fusion_Completeness']

attr_types = [AttrType.text, AttrType.image, AttrType.image,
              AttrType.image, AttrType.image,
              AttrType.image, AttrType.image,
              AttrType.image, AttrType.image,
              AttrType.number, AttrType.number,
              AttrType.number, AttrType.number,
              AttrType.number, AttrType.number]

OutputInstance = namedtuple('OutputInstance', attr_names)


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


def create_thumbnails(image_filepath, basewidth):
    img = Image.open(image_filepath)
    wpercent = (basewidth/float(img.size[0]))
    hsize = int((float(img.size[1])*float(wpercent)))
    img = img.resize((basewidth, hsize), Image.ANTIALIAS)

    filename, fileext = os.path.splitext(image_filepath)
    thumbname_filename = filename + '_thumb' + fileext
    img.save(thumbname_filename)


def load_instances(input_filepath, output_filepath, symemtry_part_index):
    dirnames = glob.glob(input_filepath + '/*')

    instances = []

    for dirname in dirnames:
        if not os.path.isdir(dirname):
            continue

        prefix = os.path.basename(dirname)
        print prefix

        if not prefix in name_filters:
            continue

        is_loaded = True

        relative_image_filepath = []
        image_filenames = []
        image_filenames.append(prefix + '_view.png')
        image_filenames.append(prefix + '_0.png')
        image_filenames.append(prefix + '_0_symmetry_accuracy.png')
        image_filenames.append(prefix + '_0_symmetry_completeness.png')
        image_filenames.append(prefix + '_0_database_accuracy.png')
        image_filenames.append(prefix + '_0_database_completeness.png')
        image_filenames.append(prefix + '_0_fusion_accuracy.png')
        image_filenames.append(prefix + '_0_fusion_completeness.png')

        for image_filename in image_filenames:
            if not os.path.exists(dirname + '/' + image_filename):
                print 'Warning: File does not exist: "' + (dirname + '/' + image_filename) + '"'
                is_loaded = False
                break

            if not os.path.exists(output_filepath + '/' + image_filename):
                # Copy the image.
                shutil.copy(dirname + '/' + image_filename, output_filepath)

                # Create a thumbnail.
                create_thumbnails(output_filepath + '/' + image_filename, thumbname_width)

            # Get relative file path.
            relative_image_filepath.append('./' + image_filename)

        if not is_loaded:
            continue

        accuracy_values = []
        completeness_values = []
        csv_filename_postfixes = []
        csv_filename_postfixes.append('_0_symmetry')
        csv_filename_postfixes.append('_0_database')
        csv_filename_postfixes.append('_0_fusion')

        for csv_filename_postfix in csv_filename_postfixes:
            csv_filename = dirname + '/' + prefix + csv_filename_postfix + '.csv'
            if symemtry_part_index >= 0:
            # Real per-part files.
                csv_filename = dirname + '/' + prefix + csv_filename_postfix\
                               + '_' + str(symemtry_part_index) + '.csv'

            all_values = get_csv_value(csv_filename, threshold)

            if not all_values:
                accuracy_values.append(float("NaN"))
                completeness_values.append(float("NaN"))
            else:
                accuracy_values.append(all_values[0])
                completeness_values.append(all_values[1])

        instance = OutputInstance(prefix, relative_image_filepath[0], relative_image_filepath[1],
                                  relative_image_filepath[2], relative_image_filepath[3],
                                  relative_image_filepath[4], relative_image_filepath[5],
                                  relative_image_filepath[6], relative_image_filepath[7],
                                  accuracy_values[0], completeness_values[0],
                                  accuracy_values[1], completeness_values[1],
                                  accuracy_values[2], completeness_values[2])

        instances.append(instance)

    return instances


def write_header(file):
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

def write_header_figure(file):
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
    file.write('\t$(\'#example\').dataTable( {')
    file.write('\t\t"paging":   false,')
    file.write('\t\t"ordering": false,')
    file.write('\t\t"info":     false,')
    file.write('\t\t"columnDefs": [')
    file.write('\t\t{')
    file.write('\t\t\t"targets": [9, 10, 11, 12, 13, 14],')
    file.write('\t\t\t"visible": false')
    file.write('\t\t}')
    file.write('\t\t]')
    file.write('\t} );')
    file.write('} );')
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


def write_html_overall_stats(file, instances):
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


def write_html_table(instances, title, filename):
    assert len(attr_names) == len(attr_types)
    num_attrs = len(attr_names)

    print "Saving the file..."
    file = open(filename, 'w')

    file.write('<html>\n')
    write_header(file)

    file.write('<body>\n')
    #file.write('<h2>' + title + '</h2>\n')
    #write_html_part_links(file)
    #write_html_overall_stats(file, instances)


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


def main():
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    instances = load_instances(input_path, output_path, -1)
    html_filename = output_path + '/index.html'
    write_html_table(instances, dataset_name + ' (All)', html_filename)

    # For each part
    for i in range(len(symmetry_part_names)):
        instances = load_instances(input_path, output_path, i)
        html_filename = output_path + '/' + symmetry_part_names[i] + '.html'
        write_html_table(instances, dataset_name + ' (' + symmetry_part_names[i].title() + ')', html_filename)


main()