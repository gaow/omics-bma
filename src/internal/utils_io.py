#! /usr/bin/env python3
# utils_io.py
# Gao Wang (c) 2015
import numpy as np
import deepdish as dd
import pandas as pd
import re
from .utils import env, rename_stamp

def str2list(value):
    return [x.strip() for x in re.split(" |\+|,", value) if x.strip()]

class InputChecker:
    def __init__(self, procedure):
        if procedure == 'test_association':
            self.apply = self.check_association_input
        elif procedure == 'fit_hm':
            self.apply = self.check_fit_hm_input
        elif procedure == 'posterior_inference':
            self.apply = self.check_posterior_inference_input
        else:
            env.error("Undefined procedure for input check %s" % procedure, exit = True)

    def check_association_input(self, params):
        #
        if "bfs" in params:
            params["bfs"] = str2list(params["bfs"])
        else:
            params["bfs"] = ["gen"]
        #
        params["sbgrp"] = str2list(params["sbgrp"])
        #
        if params["sumstats"] is not None and params["sumstats"] != 'None':
            params["sumstats"] = str2list(params["sumstats"])
            if len(params["sumstats"]) != 2:
                env.error("Parameter 'sumstats' should have 2 elements: filename, data_path", exit = True)
            if not params["sumstats"][1].startswith("/"):
                params["sumstats"][1] = "/" + params["sumstats"][1]
        #
        if "wrtsize" not in params:
            params["wrtsize"] = 10
        try:
            params["output"] = rename_stamp(re.match(r'stamp\((.*)\)', params["output"]).group(1))
        except AttributeError:
            pass
        return dict2map(params)

    def check_fit_hm_input(self, params):
        return params

    def check_posterior_inference_input(self, params):
        return params

def map2pandas(data, to_obj, rownames = None, colnames = None):
    return Map2DataFrame(to_obj).convert(data, rownames, colnames)

def load_ddm(filename, data_path):
    res = dd.io.load(filename, data_path)
    for k1, v1 in list(res.items()):
        for k2, v2 in list(v1.items()):
            res[k1][k2] = v2.transpose().to_dict()
    return res

class Map2DataFrame(object):
    '''
    Convert data from SWIG map objects to Python Pandas objects
    '''
    def __init__(self, data_type):
        if data_type == "ddm":
            self.convert = self.get_dict_x2_obj
        elif data_type == "dm":
            self.convert = self.get_dict_obj
        elif data_type == "m":
            self.convert = self.get_matrix_obj
        else:
            env.error("Unknown data type {}".format(data_type), exit = True)

    def get_dict_x2_obj(self, value, rownames, colnames):
        if colnames is None and "colnames" in rownames:
            colnames = rownames["colnames"]
        res = {}
        for k1, val1 in list(dict(value).items()):
            res[k1] = {}
            for k2, val2 in list(dict(val1).items()):
                res[k1][k2] = pd.DataFrame(np.matrix(val2), index = rownames[k1][k2],
                                           columns = colnames)
        return res

    def get_dict_obj(self, value, rownames, colnames):
        if colnames is None and "colnames" in rownames:
            colnames = rownames["colnames"]
        res = {}
        for k1, val1 in list(dict(value).items()):
            res[k1] =  pd.DataFrame(np.matrix(val1), index = rownames[k1],
                                    columns = colnames)
        return res

    def get_matrix_obj(self, value, rownames, colnames):
        return pd.DataFrame(np.matrix(value), index = rownames, columns = colnames)

def dict2map(value):
    '''
    Convert data from Python heterogeneous dict to SWIG homogeneous map objects
    '''
    params = {}
    for item in ["string", "float", "int", "vectors", "vectorf", "vectori", "dict"]:
        params[item] = {}
    params["None"] = []
    for k, val in list(value.items()):
        if type(val) == str:
            if val != "None":
                params["string"][k] = val
            else:
                params["None"].append(k)
        elif type(val) == float:
            params["float"][k] = val
        elif type(val) == int:
            params["int"][k] = val
        elif type(val) == list:
            if type(val[0]) == int:
                params['vectori'][k] = val
            if type(val[0]) == float:
                params['vectorf'][k] = val
            if type(val[0]) == str:
                params['vectors'][k] = val
        elif type(val) == dict:
            params["dict"][k] = self.get_dict(val)
        elif val is None:
            params["None"].append(k)
        else:
            pass
    return params
