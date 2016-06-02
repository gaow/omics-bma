#! /usr/bin/env python3
__author__ = "Gao Wang"
__copyright__ = "Copyright 2016, Stephens lab"
__email__ = "gaow@uchicago.edu"
__license__ = "MIT"
__version__ = "0.1.0"
import os
import numpy as np
import pandas as pd
from numpy.linalg import inv
from scipy.stats import multivariate_normal

class SnpPosteriorCalculator(object):
    '''Calculate SNP specific posterior quantities'''
    def __init__(self, Uks, w_Uk, Uk_classes, partition_by):
        '''
        Input
        -----
        Uks, {U_kl}: a (K * L) array of R x R matrices
        w_Uk, {\hat{w}_kl}: a (K * L) array of weight estimates corresponding to each Uk
        partition_by, {p_k}: a K array of integers each stands for the size of one partition.
          \sum{p_k} = (K * L)
        '''
        self.Uks, self.w_Uk, self.Uk_classes, self.partition_by = Uks, w_Uk, Uk_classes, partition_by
        self.average_bayes_factor = self.posterior_beta = self.posterior_cov = None
        self.posterior_weights = self.posterior_betas = self.posterior_covs = None
        self.likelihood = self.Uk_weights = None

    def __CalcAvgBF(self, bfs):
        '''
        Input
        -----
        bfs, {BF_kl}: a (K * L) array of Bayes factors (scalars)

        Output
        ------
        self.average_bayes_factor, BF_avg: a scalar of averaged Bayes Factor
        '''
        self.average_bayes_factor = np.average(bfs, weights = self.w_Uk)

    def __CalcPosteriorWeights(self, bfs):
        '''
        Input
        -----
        bfs, {BF_kl}: a (K * L) array of Bayes factors (scalars)
        self.average_bayes_factor, BF_avg: a scalar of averaged Bayes Factor

        Output
        ------
        self.posterior_weights, {\tilda{w}_kl}: a (K * L) array of weights (scalars)
        '''
        self.posterior_weights = np.multiply(bfs, self.w_Uk) / self.average_bayes_factor

    def __CalcPosteriorEffectSizes(self, betahat, vhat):
        '''
        Input
        -----
        betahat, \hat{beta}: a 1 X R (betahat.T is R x 1) vector of effect size estimates
        vhat, \hat{V}: an R x R variance-covariance matrix estimate

        Output
        ------
        self.posterior_covs, {\tilda{V}_kl}: a (K * L) array of R x R posterior variance-covariance matrices
        self.posterior_betas, {\tilda{beta}_kl}: a (K * L) array of R x 1 posterior vectors
        '''
        def compute_cov(vhat_inv):
            return [np.dot(U, inv(np.dot(vhat_inv, U) + np.identity(U.shape[0]))) for U in self.Uks]

        def compute_mean(betahat, vhat_inv, Us):
            return [np.dot(np.dot(U, vhat_inv), betahat.T) for U in Us]

        vhat_inv = inv(vhat)
        self.posterior_covs = compute_cov(vhat_inv)
        self.posterior_betas = compute_mean(betahat, vhat_inv, self.posterior_covs)

    def __CalcUkWeights(self, bfs):
        '''
        Input
        -----
        bfs, see above
        self.average_bayes_factor = np.average(bfs, weights = self.w_Uk)

        Output
        ------
        self.Uk_weights, {w_k}: a K array of weights corresponding to each U_k matrices
        '''
        weighted_bfs = [np.dot(bf, w) for bf, w in
                        zip(np.split(bfs, self.partition_by), np.split(self.w_Uk, self.partition_by))]
        self.Uk_weights = pd.Series(np.array(weighted_bfs) / self.average_bayes_factor,
                                    index = self.Uk_classes)

    def Calclikelihood(self, betahat, vhat):
        '''
        Input
        -----
        See above

        Output
        ------
        L: likelihood
        '''
        likelihoods = [multivariate_normal.pdf(betahat, np.zeros(vhat.shape[0]), p + vhat)
                       for p in self.Uks]
        self.likelihood = np.average(likelihoods, weights = self.w_Uk)
        # Note: likelihoods / self.likelihood should equal self.posterior_weights

    def CalcPosterior(self, betahat, vhat, bfs):
        '''
        Input
        -----
        See above

        Output
        ------
        self.posterior_beta: \tilda{beta}: an R x 1 vector of posterior effect size
        self.posterior_cov: \tilda{V}: an R x R matrix of posterior variance-covariance matrix
        '''
        self.__CalcAvgBF(bfs)
        self.__CalcUkWeights(bfs)
        self.__CalcPosteriorEffectSizes(betahat, vhat)
        self.__CalcPosteriorWeights(bfs)
        # FIXME: double-check the math
        self.posterior_beta = np.average(self.posterior_betas, axis = 0, weights = self.posterior_weights)
        # law of total variance
        tmp_covs = [U + np.dot(b - self.posterior_beta, (b - self.posterior_beta).T)
                    for U, b in zip(self.posterior_covs, self.posterior_betas)]
        # FIXME: double-check the math
        self.posterior_cov = np.average(tmp_covs, axis = 0, weights = self.posterior_weights)

    def GetResults(self):
        return {'post_beta': self.posterior_beta,
                'post_cov': self.posterior_cov,
                'bf_avg': self.average_bayes_factor,
                'likelihood': self.likelihood,
                'Uk_weights': self.Uk_weights}


class BlockPosteriorCalculator(object):
    '''Calculate Block specific posterior quantities'''
    def __init__(self, pi0):
        # FIXME: documentation
        self.posterior_has_signal = None
        self.posterior_is_signal = None
        self.posterior_signal = None
        self.posterior_block_bayes_factor = None
        self.pi0 = pi0

    def Calculate(self, bfs):
        m = len(bfs)
        self.posterior_block_bayes_factor = np.sum(bfs) / m
        tmp_prop = (1 - self.pi0) * self.posterior_block_bayes_factor
        self.posterior_has_signal = tmp_prop / (self.pi0 + tmp_prop)
        self.posterior_is_signal = bfs / (m * self.posterior_block_bayes_factor)
        self.posterior_signal = bfs / (bfs + m - 1)

    def GetResults(self):
        return {'bf': self.posterior_block_bayes_factor,
                'has_signal': self.posterior_has_signal,
                'is_signal': self.posterior_is_signal,
                'signal': self.posterior_signal}

from .io import get_tb_groups, load_priors
from .deepdishio import load, save

class PosteriorController(object):
    def __init__(self, sumstats_path, bf_path, prior_path, weights_path, output_path, config):
        '''each path variable is a tuple of (HDF5 file name, table name)
        the controller should get data, make sure they are properly matched,
        feed into calculators and collect results'''
        # FIXME: check path, should not be too deep!
        priors = load_priors(prior_path, config['nb_groups'], config['bf_config'])
        weights = load(weights_path[0], weights_path[1]).to_dict()
        self.Uk_names, self.Uk_classes, self.Uk_partitions = self.__SortUkNames(list(weights.keys()))
        self.SnpCalculator = SnpPosteriorCalculator(np.array([priors[k] for k in self.Uk_names]),
                                                    np.array([weights[k] for k in self.Uk_names]),
                                                    self.Uk_classes, self.Uk_partitions)
        self.BlockCalculator = BlockPosteriorCalculator(weights['null'])
        self.sumstats = sumstats_path
        self.bf = bf_path
        self.output = output_path
        self.blocks = get_tb_groups(sumstats_path[0], sumstats_path[1])

    def __call__(self, nb_threads = 2):
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
        save(self.output[0], res)

    def __CalcBlock(self, block_name):
        '''
        For each block extract BFs, betahat and vhat and match them by SNP_ID
        '''
        res = {}
        # bfs is a matrix with SNP ID as row names and Uk as columns, rearranged here column-wise
        bfs = load(self.bf[0], os.path.join(self.bf[1], block_name)).rename(columns = {'nb_groups' : 'null'})
        bfs = bfs[self.Uk_names]
        # FIXME: need BF input in original scale not log scale
        bfs['null'] = 0
        # sumstats is a dictionary with SNP ID as keys
        sumstats = load(self.sumstats[0], os.path.join(self.sumstats[1], block_name))
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
