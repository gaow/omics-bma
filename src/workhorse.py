#! /usr/bin/env python3
# utils_workhorse.py
# Gao Wang (c) 2015
import numpy as np
import deepdish as dd
import os
from .io import get_tb_groups, load_priors
from .core import SnpPosteriorCalculator, BlockPosteriorCalculator

class PosteriorController(object):
    def __init__(self, sumstats_path, bf_path, prior_path, weights_path, output_path, config):
        '''each path variable is a tuple of (HDF5 file name, table name)
        the controller should get data, make sure they are properly matched,
        feed into calculators and collect results'''
        # FIXME: check path, should not be too deep!
        priors = load_priors(prior_path, config['nb_groups'], config['bf_config'])
        weights = dd.io.load(weights_path[0], weights_path[1]).to_dict()
        self.Uk_names, self.Uk_classes, self.Uk_partitions = self.__SortUkNames(list(weights.keys()))
        self.SnpCalculator = SnpPosteriorCalculator(np.array([priors[k] for k in self.Uk_names]),
                                                    np.array([weights[k] for k in self.Uk_names]),
                                                    self.Uk_classes, self.Uk_partitions)
        self.BlockCalculator = BlockPosteriorCalculator(weights['null'])
        self.sumstats = sumstats_path
        self.bf = bf_path
        self.output = output_path
        self.blocks = get_tb_groups(sumstats_path[0], sumstats_path[1])

    def ScanBlocks(self):
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
        for item in self.blocks:
            res[item] = self.__CalcBlock(item)
        dd.io.save(self.output[0], res)

    def __CalcBlock(self, block_name):
        '''
        For each block extract BFs, betahat and vhat and match them by SNP_ID
        '''
        res = {}
        # bfs is a matrix with SNP ID as row names and Uk as columns, rearranged here column-wise
        bfs = dd.io.load(self.bf[0], os.path.join(self.bf[1], block_name)).\
          rename(columns = {'nb_groups' : 'null'})
        bfs = bfs[self.Uk_names]
        # FIXME: need BF input in original scale not log scale
        bfs['null'] = 0
        # sumstats is a dictionary with SNP ID as keys
        sumstats = dd.io.load(self.sumstats[0], os.path.join(self.sumstats[1], block_name))
        for snp in list(bfs.index):
            # FIXME: need BF input in original scale not log scale
            bfs_snp = np.power(10, bfs.loc[snp]).as_matrix()
            betahat = np.matrix(sumstats[snp][0:1])
            vhat = np.matrix(sumstats[snp][1:])
            self.SnpCalculator.Calclikelihood(betahat, vhat)
            self.SnpCalculator.CalcPosterior(betahat, vhat, bfs_snp)
            res[snp] = self.SnpCalculator.GetResults()
        # FIXME: block level posterior not saved yet
        return res

    def __SortUkNames(self, names):
        res = {}
        res_names = []
        res_counts = []
        for x in names:
            key = os.path.splitext(x)[0]
            if key not in res:
                res[key] = []
            res[key].append(x)
        #
        classes = sorted(res.keys())
        classes.insert(0, classes.pop(classes.index('null')))
        for key in classes:
            res_names.extend(sorted(res[key]))
            res_counts.append(len(res[key]))
        return res_names, classes, np.cumsum(res_counts[:-1])
