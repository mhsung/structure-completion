#!/usr/bin/python

##########
#
#  Get Statistics
#
##########

import os, string, glob, sys, csv
import numpy as np
import matplotlib.pyplot as plt


##########
def getAttribute(fileList, attributeName):
    data = []
    count = 0

    for fileName in fileList:
        ++count
        with open(fileName, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            for row in csvreader:
                #print(row[0])
                if row[0] == attributeName:
                    data = np.append(data, float(row[1]))

    return data


##########
def plotAccuracy(data, x_range_1, x_range_2, title):
    x_axis = np.linspace(x_range_1, x_range_2, num=30)
    y_axis = []

    for x in np.nditer(x_axis):
        if x_range_1 < x_range_2:
            cumulative_count = np.count_nonzero(data < x)
        else:
            cumulative_count = np.count_nonzero(data > x)
        y_axis = np.append(y_axis, (float)(cumulative_count) / len(data))

    plt.plot(x_axis.tolist(), y_axis.tolist())
    plt.xlim(x_range_1, x_range_2)
    plt.xlabel(title)
    plt.ylabel("Percentage of Objects")

    plt.axvline(data.mean(), color='r')
    plt.annotate("Average = " + str(data.mean()), xy=(data.mean(), 0.5),
                 xytext=(data.mean() + (x_range_2 - x_range_1)/30, 0.2))
    #plt.show()
    plt.savefig("_" + title)


##########
dataDir = ""

## Parse arguments
if len(sys.argv) < 1:
    print("Usage: " + sys.argv[0] + " dataset_path")
    exit()
else:
    dataDir = sys.argv[1]

fileList = [os.path.basename(x) for x in glob.glob(dataDir + "/*.csv")];

## Batch jobs
os.chdir(dataDir)
count = 0

data = getAttribute(fileList, "correct_sample_point_label_ratio")
plotAccuracy(data, 1.0, 0.7, "Point Segmentation Accuracy")

data = getAttribute(fileList, "max_cuboid_distance")
plotAccuracy(data, 0.0, 0.4, "Max Cuboid Distance")

