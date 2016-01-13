#!/usr/bin/env python3
# Gao Wang (c) 2015
import numpy as np
import mosek as mk
from scipy.sparse import coo_matrix as make_aij_sparse
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

def mixIP(f_data, f_prior):
    # Load data
    matrix_lik, prior = load_data(f_data, f_prior)
    n, k = matrix_lik.shape
    # Add in observations corresponding to prior
    matrix_lik = np.vstack((np.diag(np.ones(len(prior))), matrix_lik))
    w = np.concatenate((prior - 1, np.ones(n)))
    # Remove zero weight entries
    matrix_lik = matrix_lik[w != 0,:]
    w = w[w != 0]
    # Optimize
    res = kw_dual(matrix_lik, np.ones(k), normalize(w))


def load_data(lik, prior):
    return np.loadtxt(lik), np.loadtxt(prior)

def normalize(v):
    return v / np.sum(v)

def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def kw_dual(A, d, w, rtol = 1E-6, verb = 0, **kwargs):
    '''
    Dual Kiefer-Wolfowitz MLE for Mixture Problems

       This version implements a class of density estimators solving:

       min_x  {F(x) := sum -log (x_i)}  s.t. A' x <= d, 0 <= x,

	where e.g.  A = phi(outer(Y,g,"fun")), with Y data and g a grid on the support of Y,
	and "fun" is some function representing the dependence of the base distribution.

    This Python implementation is adapted from R function REBayes::KWDual()

    Parameters
    ----------
    A : matrix
      Linear constraint matrix
    d : vector
      constraint vector
    w : vector
      weights for `x`, should sum to one
    rtol :
      the relative tolerance for dual gap convergence criterion
    verb : Integer from 0 to 5
      to control verbosity desired from mosek
    control : dictionary
      various mosek control parameters, see mosek Python API documentation for details
      http://docs.mosek.com/7.1/pythonapi/Parameters.html
      http://docs.mosek.com/7.1/pythonapi/Conventions_employed_in_the_API.html

    Returns
    -------
    data : anything
    '''
    # Since the value of infinity is ignored, we define it solely
    # for symbolic purposes
    inf = 0.0
    # number of variables and number of constraints
    numvar, numcon = A.shape
    # Make sparse matrix to satisfy Mosek requirement
    A_sparse = make_aij_sparse(A)
    # Bound key for constraints, <= d
    bkc = [ mk.boundkey.up ] * numcon
    # Lower / upper bounds for constraints
    blc = [ 0.0 ] * numcon
    buc = d
    # Bound key for variables, >= 0
    bkx = [ mk.boundkey.lo ] * numvar
    blx = [ 0.0 ] * numvar
    bux = [ inf ] * numvar
    # separable terms of the objective
    opro = [ mk.scopr.log ] * numvar
    oprjo = np.arange(numvar)
    oprfo = -1 * w
    oprgo = np.ones(numvar)
    oprho = np.zeros(numvar)
    #
    res = [ 0.0 ] * numcon
    # set mosek params
    # FIXME: all kwargs should be set here
    mk.dparam.intpnt_nl_tol_rel_gap = rtol

    with mk.Env() as mk_env:
        mk_env.set_Stream (mk.streamtype.log, streamprinter)
        with mk_env.Task(0,0) as task:
            task.set_Stream (mk.streamtype.log, streamprinter)
            task.appendvars(numvar)
            task.appendcons(numcon)
            # logic follows from REBayes::KWDual hereafter
            task.putobjsense(mk.objsense.minimize)
            ## objective coefficients
            task.putclist(np.arange(numvar), np.zeros(numvar))
            ## constraint matrix coefficients
            task.putaijlist(A_sparse.col, A_sparse.row, A_sparse.data)
            ## bounds
            task.putconboundslice(0, numcon, bkc, blc, buc)
            task.putvarboundslice(0, numvar, bkx, blx, bux)
            ## solver setup
            task.putSCeval(opro = opro, oprjo = oprjo, oprfo = oprfo,
                           oprgo = oprgo, oprho = oprho)
            ## run
            task.optimize()
            ## get results
            task.getsolutionslice(mk.soltype.itr, mk.solitem.xc, 0, numcon, res)
            converged = task.getsolsta(mk.soltype.itr) in [mk.solsta.optimal, mk.solsta.near_optimal]
    print(oprfo[0:10], normalize(res), converged)
    return normalize(res)

if __name__ == '__main__':
    mixIP("data.txt", "prior.txt")
