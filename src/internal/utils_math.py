#! /usr/bin/env python3
# utils_math.py
# Gao Wang (c) 2015
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


def test_almost_equal_recursive(x, y, level = 6, level_cutoff = 3):
    try:
      np.testing.assert_almost_equal(x, y, level)
    except AssertionError:
      if level <= level_cutoff:
          raise
      else:
          level = level - 1
          test_almost_equal_recursive(x, y, level)
    return level
