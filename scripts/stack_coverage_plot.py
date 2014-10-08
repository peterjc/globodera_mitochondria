#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg') # Want to use without X windows
from matplotlib import pyplot as plt
from optparse import OptionParser

VERSION = "0.0.2"

USAGE = """usage: %prog -o OUTPUT.png INPUT.cov [INPUT2.cov ...]

Default output is one PNG per input coverage file, appending .png
to the input filename.

Coverage files are expected in a FASTA QUAL like format, with
a header line ">identifier" followed by multiple (long) lines
of space separated coverage numbers for each base of the
specificed sequence. Each line represents a different coverage
category.
"""

parser = OptionParser(USAGE,
                      version="%prog "+VERSION)

parser.add_option("-o", "--output", dest="output",
                  type="string", metavar="FILE",
                  help="PNG file to create.")

(options, args) = parser.parse_args()

if not args:
    sys.stderr.write("At least one coverage file required.")
    sys.exit(1)

def load(filename):
    h = open(filename)
    line = h.readline()
    assert line.startswith(">")
    while line and line[0] == ">":
        name = line[1:].split(None,1)[0]
        values = []
        while line:
            line = h.readline()
            if not line or line[0] == ">":
                break
            values.append([float(v) for v in line.rstrip("\n").split("\t")])
        yield name, np.array(values, np.float)
    h.close()

def make_colors(start, end, steps):
    delta = (end - start) / float(steps-1)
    return ["#%02x%02x%02x" % tuple(start + i*delta) for i in range(steps)]

def stack(data, filename, colors=None):
    total = len(data)
    max_value = 0
    for names, values in data:
        max_value = max(max_value, values.sum(axis=0).max())
    plt.ylim([0, max_value])

    fig = plt.figure(figsize=(12,2*total))
    if not colors:
        #Assumes all the examples have same number of colors:
        if data[0][1].shape[0] == 3:
            colors = ["#CC6666", "#1DACD6", "#6E5160"]
        elif data[0][1].shape[0] == 5:
            colors = ["#CDCDC1", "#8B8B83", "#FF6A6A", "#F0E68C", "#CDC673"]
        else:
            colors = make_colors(np.array([0xCC, 0x66, 0x66]),
                                 #np.array([0x6E, 0x51, 0x60]),
                                 np.array([0x90, 0x41, 0x50]),
                                 #np.array([0x20, 0xF0, 0x60]),
                                 data[0][1].shape[0])
        print colors
    for i, (name, values) in enumerate(data):
        x = range(values.shape[1])
        print i, name, values.shape, "coverage:"
        print "\t".join("%0.1f" % v for v in values.sum(axis=1))
        y_stack = np.cumsum(values, axis=0)
        ax1 = fig.add_subplot(total, 1, i+1)
        ax1.set_autoscaley_on(False)
        ax1.set_ylim([0, max_value])
        ax1.set_title(name.split(None,1)[0], fontsize="xx-small")
        ax1.fill_between(x, 0, y_stack[0,:], facecolor=colors[0], alpha=.7)
        for i in range(0, values.shape[0]-1):
            ax1.fill_between(x, y_stack[i,:], y_stack[i+1,:], facecolor=colors[i+1], alpha=.7)
    #fig.tight_layout()
    plt.show()
    plt.savefig(filename)


if options.output:
    print("Drawing %i coverage files to %s" % (len(args), options.output))
    #One combined PNG
    data = []
    for filename in args:
        data.extend(load(filename))
    stack(data, options.output)
else:
    #One PNG per coverage file
    for filename in args:
        if not filename.endswith(".cov"):
            continue
        print "-"*60
        print filename
        print "-"*60
        data = list(load(filename))
        stack(data, filename+".png")
print "Done"
