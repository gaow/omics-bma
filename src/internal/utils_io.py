#! /usr/bin/env python3
# utils_io.py
# Gao Wang (c) 2015
import numpy as np
import deepdish as dd
import pandas as pd
import re
from .utils import env, rename_stamp
import tables as tb
from .pyeqtlbma import get_eqtlbma_configurations, dict_x2_vectors

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
        try:
            params["output"] = rename_stamp(re.match(r'stamp\((.*)\)', params["output"]).group(1))
        except AttributeError:
            pass
        if not 'optimizer_control' in params:
            params['optimizer_control'] = {}
        if "thread" in params:
            params['optimizer_control']['iparam.num_threads'] = params['thread']
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
            self.convert = self.get_dict_x2_matrix_obj
        elif data_type == "dm":
            self.convert = self.get_dict_matrix_obj
        elif data_type == "m":
            self.convert = self.get_matrix_obj
        elif data_type == "dv":
            self.convert = self.get_dict_vector_obj
        else:
            env.error("Unknown data type {}".format(data_type), exit = True)

    def get_dict_x2_matrix_obj(self, value, rownames, colnames):
        if colnames is None and "colnames" in rownames:
            colnames = rownames["colnames"]
        res = {}
        for k1, val1 in list(dict(value).items()):
            res[k1] = {}
            for k2, val2 in list(dict(val1).items()):
                res[k1][k2] = pd.DataFrame(np.matrix(val2),
                                           index = rownames[k1][k2] if type(rownames) is dict_x2_vectors else rownames,
                                           columns = colnames)
        return res

    def get_dict_matrix_obj(self, value, rownames, colnames):
        if colnames is None and "colnames" in rownames:
            colnames = rownames["colnames"]
        res = {}
        for k1, val1 in list(dict(value).items()):
            res[k1] =  pd.DataFrame(np.matrix(val1), index = rownames[k1],
                                    columns = colnames)
        return res

    def get_dict_vector_obj(self, value, rownames = None, colnames = None):
        return dict(value)

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
            # FIXME: should write some recursive function so that the yaml can have multiple layers
            params["dict"][k] = val
        elif val is None:
            params["None"].append(k)
        else:
            pass
    return params


def get_tb_groups(filenames, group_name = None, verbose = True):
    if verbose:
        env.log('Collecting group names from input files ...')
    names = set()
    for filename in filenames:
        if verbose:
            env.log(filename, flush = True)
        with tb.open_file(filename) as f:
            names.update([node._v_name for node in \
                          (f.root if group_name is None else getattr(f.root, '{}'.format(group_name)))])
        if verbose:
            env.log('%s unique groups identified from %s files\n' \
                    % (len(names), len(filenames)), flush = True)
    return sorted(names)

def load_priors(prior_path, dim, config):
    '''
    Input
    -----
    prior_path: (file, table)
    dim: integer of number of groups
    config: str
      choose from 'gen', 'sin', 'all'

    Output
    ------
    A dictionary of prior matrices with keys matching ID's for Bayes Factors and weights
    '''
    assert config in ['gen', 'sin', 'all']
    data = dd.io.load(prior_path[0], prior_path[1])
    res = {}
    # FIXME: check prior input data make sure they are good
    # Add null model
    res['null'] = np.zeros((dim, dim))
    # expand customized matrices
    for k, value in data.items():
        if k in ['gridL', 'gridS', 'gridM']:
            continue
        for idx, item in enumerate(data['gridM'][0]):
            res["{}.{}".format(k, idx + 1)] = value * item
    # expand "consistent configuration"
    for idx, item in enumerate(data['gridL']):
        # "gen"
        res['gen.{}'.format(idx + 1)] = np.full((dim, dim), item[1])
        np.fill_diagonal(res['gen.{}'.format(idx + 1)], np.sum(item))
        # "gen-fix"
        res['gen-fix.{}'.format(idx + 1)] = np.full((dim, dim), np.sum(item))
        # "gen-maxh"
        res['gen-maxh.{}'.format(idx + 1)] = np.diag(np.full(dim, np.sum(item)))
    if config == 'gen':
        return
    # expand other configurations
    for k, value in map2pandas(get_eqtlbma_configurations(dim, True if config == "all" else False),
                               "dv").items():
        for idx, item in enumerate(data['gridS']):
            res["cfg.{}.{}".format(k, idx + 1)] = np.full((dim, dim), item[1])
            np.fill_diagonal(res["cfg.{}.{}".format(k, idx + 1)], np.sum(item))
            res["cfg.{}.{}".format(k, idx + 1)][np.where(np.outer(value, value) == 0)] = 0
    return res
