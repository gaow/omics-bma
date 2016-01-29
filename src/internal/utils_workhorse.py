#! /usr/bin/env python3
# utils_workhorse.py
# Gao Wang (c) 2015
import pandas as pd
import deepdish as dd
from .utils_io import get_tb_groups, load_priors
from .utils_math import SnpPosteriorCalculator, BlockPosteriorCalculator
import os

class PosteriorController(object):
    def __init__(self, sumstats_path, cov_path, bf_path, prior_path, weights_path, params_path):
        '''each path variable is a tuple of (HDF5 file name, table name)
        the controller should get data, make sure they are properly matched,
        feed into calculators and collect results'''
        # FIXME: check path, should not be too deep!
        self.params = dd.io.load(params_path[0], params_path[1])
        priors = load_priors(prior_path, self.params['nb_groups'], self.params['bf_config'])
        weights = dd.io.load(weights_path[0], weights_path[1]).to_dict()
        self.Uk_names = sorted(self.weights.keys())
        self.SnpCalculator = SnpPosteriorCalculator([priors[k] for k in self.Uk_names],
                                                    [weights[k] for k in self.Uk_names])
        self.BlockCalculator = BlockPosteriorCalculator(weights['null'])
        self.sumstats = sumstats_path
        self.cov = cov_path
        self.bf = bf_path
        self.tmp_snps = None

    def ScanBlocks(self, block_list):
        '''
        Input
        -----
        A list of names of blocks, e.g., genes, to compute

        Output
        ------
        Various posterior quantities for the block of interest
        '''
        # FIXME: single core implementation for now. will move to multiple cores down the line
        res = {}
        for item in block_list:
            res[item] = self.__CalcBlock(item)
        return res

    def __CalcBlock(self, block_name):
        '''
        For each block extract BFs, betahat and vhat and match them by SNP_ID
        '''
        # bfs is a matrix with SNP ID as row names
        bfs = dd.io.load(self.bf[0], os.path.join(self.bf[1], block_name))
        # vhat is a dictionary with SNP ID as keys
        vhat = dd.io.load(self.cov[0], os.path.join(self.cov[1], block_name))
        betahat = {}
        # for sbgrp, ss in dd.io.load(self.sumstats[0], os.path.join(self.sumstats[1], block_name)).items():
            # for
