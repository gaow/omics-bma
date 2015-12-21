#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015

def convert_data(data, to_obj):
    if type(data) in (dict, list, tuple):
        return Py2Swig(to_obj).convert(data)
    elif type(data) in (str, int, float):
        return data
    else:
        return Swig2Py(to_obj).convert(data)

class Swig2Py(object):
    '''
    Convert data from SWIG objects to Python objects
    '''
    def __init__(self, data_type):
        if data_type == "ddm":
            self.convert = self.get_dict_dict_matrix
        else:
            raise ValueError("Unknown data type {}".format(data_type))

    def get_dict_dict_matrix(self, value):
        res = {}
        for k1, val1 in list(dict(value).items()):
            res[k1] = {}
            for k2, val2 in list(dict(val1).items()):
                res[k1][k2] = list(val2)
        return res

class Py2Swig(object):
    '''
    Convert data from Python objects to SWIG interface compatible objects
    '''
    def __init__(self, data_type):
        if data_type == "dict":
            self.convert = self.get_dict

    def get_dict(self, value):
        params = {}
        for item in ["string", "float", "int", "vectors", "vectorf", "vectori", "dict"]:
            params[item] = {}
        for k, val in list(value.items()):
            if type(val) == str:
                params["string"][k] = val
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
            else:
                pass
        return params
