from OmicsBMA.utils import openFile, get_value_type
from pandas import DataFrame
import tables as tb
import warnings
warnings.simplefilter(action = "ignore", category = tb.NaturalNameWarning)
from OmicsBMA.deepdishio import load
import os

COLNAMES = ['gen.' + str(i + 1) for i in range(25)] + \
  ['gen-fix.' + str(i + 1) for i in range(25)] + \
  ['gen-maxh.' + str(i + 1) for i in range(25)]
for j in range(44):
    COLNAMES.extend(['cfg.{}.{}'.format(j + 1, i + 1) for i in range(10)])

def load_eqtlbma_bf(filename):
    with openFile(filename) as f:
        # skip header
        f.readline()
        # load all data
        data = [[float(x) if get_value_type(x) in ['float', 'int'] else x.decode("utf-8") for x in line.strip().split()] for line in f.readlines()]
    res = {}
    for line in data:
        if not line[0] in res:
            res[line[0]] = {}
        if not line[1] in res[line[0]]:
            res[line[0]][line[1]] = []
        res[line[0]][line[1]].extend([x for x in line[3:] if x == x])
    for key in list(res.keys()):
        res[key] = DataFrame.from_dict(res[key], orient = 'index')
        res[key].columns = COLNAMES
    # print(res)
    return res

def load_omicsbma_bf(filename):
    res = load(filename, "/")
    # print(res)
    return res
