#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015

def convert_data(data, to_obj = None):
    if to_obj is not None:
        return Map2Dict(to_obj).convert(data)
    elif type(data) == dict:
        return Dict2Map().convert(data)
    else:
        return data

class Map2Dict(object):
    '''
    Convert data from SWIG map objects to Python dict objects
    '''
    def __init__(self, data_type):
        if data_type == "ddm":
            self.convert = self.get_dict_dict_matrix
        elif data_type == "dm":
            self.convert = self.get_dict_matrix
        else:
            raise ValueError("Unknown data type {}".format(data_type))

    def get_dict_dict_matrix(self, value):
        res = {}
        for k1, val1 in list(dict(value).items()):
            res[k1] = {}
            for k2, val2 in list(dict(val1).items()):
                res[k1][k2] = list(val2)
        return res

    def get_dict_matrix(self, value):
        res = {}
        for k1, val1 in list(dict(value).items()):
            res[k1] = list(val1)
        return res

class Dict2Map(object):
    '''
    Convert data from Python heterogeneous dict to SWIG homogeneous map objects
    '''
    def __init__(self):
        pass

    def convert(self, value):
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
