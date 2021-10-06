#!/usr/bin/env python3
"""
Author: Joris Louwen
Script to plot convergence of LDA model.
"""

import re
import matplotlib.pyplot as plt
from sys import argv

if __name__ == '__main__':
    logfile= argv[1]
    outfile = logfile.split('.txt')[0]+'_convergence_likelihood.pdf'
    p = re.compile("(-*\d+\.\d+) per-word .* (\d+\.\d+) perplexity")
    matches = [p.findall(l) for l in open(logfile)]
    matches = [m for m in matches if len(m) > 0]
    tuples = [t[0] for t in matches]
    perplexity = [float(t[1]) for t in tuples]
    liklihood = [float(t[0]) for t in tuples]
    cores=10 #eval_every * cores used
    iters = list(range(0,len(tuples)*10,cores))
    plt.plot(iters,liklihood,c="black")
    plt.ylabel("log likelihood")
    plt.xlabel("iteration")
    plt.title("Topic Model Convergence")
    plt.grid()
    plt.savefig(outfile)
    # plt.show
    plt.close()
