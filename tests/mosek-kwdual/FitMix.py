#! /usr/bin/env python3
# FitMix.py
# Gao Wang (c) 2015
from OmicsBMA.mix_opt import mixIP
import numpy as np

if __name__ == '__main__':
    data = np.loadtxt("data.txt")
    res = mixIP(data, np.loadtxt("prior.txt"))
    print('\t'.join(list(map(str, res[0]))))
    print(res[1])
    res = mixIP(data / np.max(data), np.loadtxt("prior.txt"))
    print('\t'.join(list(map(str, res[0]))))
    print(res[1])
