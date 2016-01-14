#! /usr/bin/env python3
# FitMix.py
# Gao Wang (c) 2015
from OmicsBMA.mix_opt import mixIP
import numpy as np

if __name__ == '__main__':
    res = mixIP(np.loadtxt("data.txt"), np.loadtxt("prior.txt"))
    print('\t'.join(list(map(str, res[0]))))
    print(res[1])
