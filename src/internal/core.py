#! /usr/bin/env python3
# core.py
# Gao Wang (c) 2015
from OmicsBMA.utils_io import *
from OmicsBMA.pyeqtlbma import *
import deepdish as dd

def run_eqtlbma_bf(params, to_file):
    params = convert_data_(params)
    exe = eQtlBma()
    exe.eqtlbma_bf(params["string"], params["int"], params["float"], params["vectors"],
                   {'':{'':[[]]}} if "sumstats" in params["None"]
                   else load_sumstats_(params["string"]["sumstats"]))
    res = {"SumStats": convert_data_(exe.GetSstats(), "ddm", rownames = exe.GetSstatsRownames()),
           "Abfs": convert_data_(exe.GetAbfs(), "dm", rownames = exe.GetAbfsNames()),
           "SepPermPvals": convert_data_(exe.GetSepPermPvals(), "dm",
                                         rownames = exe.GetSepPermPvalsRownames()),
           "JoinPermPvals": convert_data_(exe.GetJoinPermPvals(), "m",
                                          rownames = exe.GetJoinPermPvalsRownames())
           }
    dd.io.save(to_file, dict((k, v) for k, v in res.items() if not is_empty_(v)),
               compression=("zlib", 9))
