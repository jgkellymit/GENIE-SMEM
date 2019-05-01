import csv
import numpy as np
import matplotlib.pyplot as plt

models = [
    '[10, 10]',
    '[10, 100]',
    '[10, 1000]',
    '[100, 100]',
    '[100, 1000]',
    '[100]',
    '[1000, 1000]',
    '[1000, 10000]',
    '[1000]',
]

for model in models:
    plt.figure()

    x = []
    skip = True
    with open('results_' + model + '.csv') as csvfile:
        csvreader = csv.reader(csvfile)
        for line in csvreader:
            if skip:
                skip = False
                continue
            x.append(int(line[0]))

    x_99_percentile = np.percentile(x, 99)
    x_filtered = [x_p for x_p in x if x_p <= x_99_percentile]

    num_bins = 100
    n, bins, patches = plt.hist(x_filtered, num_bins, facecolor='blue', alpha=0.5)

    plt.title('Distribution of Error Distances for ' + model)
    plt.xlabel('Error Distance')
    plt.ylabel('Count')

    plt.savefig('plot_' + model + '.png')