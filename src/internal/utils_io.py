#! /usr/bin/env python3
# utils_io.py
# Gao Wang (c) 2015
import numpy as np
import deepdish as dd
import pandas as pd
import re
from .utils import env, rename_stamp, flatten_dict, is_null
import tables as tb
from .pyeqtlbma import get_eqtlbma_configurations, dict_x2_vectors


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
    for filename in filenames if type(filenames) is list else [filenames]:
        if verbose:
            env.log(filename, flush = True)
        with tb.open_file(filename) as f:
            names.update([node._v_name for node in \
                          (f.root if group_name is None else getattr(f.root, '{}'.format(group_name)))])
        if verbose:
            env.log('%s unique groups identified from %s files\n' \
                    % (len(names), len(filenames) if type(filenames) is list else 1), flush = True)
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

def str2list(value):
    if value is None:
        return []
    else:
        return [x.strip() for x in re.split(" |\+|,", value) if x.strip()]

class ConfigReader(dict):
    def __init__(self, params):
        self.update({'association_data': None,
                     'input_sumstats_data': None,
                     'mixture_weights_data': None,
                     'posterior_data': None,
                     'extract_average_bf_per_class': 0,
                     'permsep': 0,
                     'trick': 0,
                     'nperm': 100,
                     'tricut': 10,
                     'maxbf': 0,
                     'pbf': 'all',
                     'gcoord': None,
                     'snp': None,
                     'geno': None,
                     'sbgrp': None,
                     'covar': None,
                     'output_sumstats_data': None,
                     'exp': None,
                     'prior_data': None,
                     'scoord': None,
                     'bf_config': 'sin',
                     'nb_groups': None,
                     'verbose': 1,
                     'seed': 10086,
                     'thread': 4,
                     'fiterr': 0.5,
                     'anchor': 'TSS',
                     'error': None,
                     'qnorm': 1,
                     'maf': None,
                     'analys': None,
                     'lik': 'normal',
                     'bfs': 'sin',
                     'cis': 1000,
                     'wrtsize': 10})
        self.update(flatten_dict(params))
        self.check_input()
        self.check_output()
        self.check_analysis()
        self.check_permutation()
        self.check_runtime()

    def check_input(self):
        self["sbgrp"] = str2list(self["sbgrp"])

    def check_output(self):
        for item in ['output_sumstats_data', 'association_data', 'posterior_data']:
            try:
                self[item] = rename_stamp(re.match(r'stamp\((.*)\)', self[item]).group(1))
            except AttributeError:
                pass

    def check_analysis(self):
        self["bfs"] = str2list(self["bfs"])

    def check_permutation(self):
        pass

    def check_runtime(self):
        if not 'optimizer_control' in self:
            self['optimizer_control'] = {}
        if "thread" in self:
            self['optimizer_control']['iparam.num_threads'] = self['thread']
