#!/usr/bin/env python3
__author__ = "Gao Wang"
__copyright__ = "Copyright 2016, Stephens lab"
__email__ = "gaow@uchicago.edu"
__license__ = "MIT"
__version__ = "0.1.0"
import sys
import numpy as np
import mosek as mk
from mosek import iparam, dparam, sparam
from scipy.sparse import csc_matrix as make_sparse
import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

def mixIP(matrix_lik, prior = None, rtol = None, control = {}):
    n, k = matrix_lik.shape
    w = np.ones(n)
    if prior is not None:
        # Add in observations corresponding to prior
        matrix_lik = np.vstack((np.diag(np.ones(len(prior))), matrix_lik))
        w = np.concatenate((np.array(prior) - 1, w))
        # Remove zero weight entries
        matrix_lik = matrix_lik[w != 0,:]
        w = w[w != 0]
    # Optimize
    res = kw_dual(matrix_lik, np.ones(k), normalize(w), rtol = rtol, control = control)
    return res

def normalize(v):
    return v / np.sum(v)

def streamprinter(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def kw_dual(A, d, w, rtol = None, control = {}):
    '''
    Dual Kiefer-Wolfowitz MLE for Mixture Problems

       This version implements a class of density estimators solving:

       min_x  {F(x) := sum -log (x_i)}  s.t. A' x <= d, 0 <= x,

	where e.g.  A = phi(outer(Y,g,"fun")), with Y data and g a grid on the support of Y,
	and "fun" is some function representing the dependence of the base distribution.

    This Python implementation is adapted from R function REBayes::KWDual()

    Parameters
    ----------
    A : double[][]
      Linear constraint matrix
    d : double[]
      constraint vector
    w : double[]
      weights for `x`, should sum to one
    rtol : double
      the relative tolerance for dual gap convergence criterion
      default to min(1E-7, 0.1 / numvars)
    control  : dictionary
      various mosek control parameters, e.g. control = {'iparam.log' : 10, 'iparam.num_threads' : 8}
      see mosek Python API documentation for details:
      http://docs.mosek.com/7.1/pythonapi/Parameters.html
      http://docs.mosek.com/7.1/pythonapi/Conventions_employed_in_the_API.html

    Returns
    -------
    c : double[]
      solution for the constraints
    converged : bool
      whether or not the optimizer converged
    '''
    # Since the value of infinity is ignored, we define it solely
    # for symbolic purposes
    inf = 0.0
    # number of variables and number of constraints
    numvar, numcon = A.shape
    # Make sparse matrix to satisfy Mosek requirement
    ## aptrb = AT_sparse.indptr[:-1]
    ## aptre = AT_sparse.indptr[1:]
    ## asub = AT_sparse.indices
    ## aval = AT_sparse.data
    AT_sparse = make_sparse(A.T)
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
    if rtol is None:
        rtol = min(1E-7, 0.1 / float(numvar))
    res = [ 0.0 ] * numcon

    with mk.Env() as mk_env:
        mk_env.set_Stream (mk.streamtype.log, streamprinter)
        with mk_env.Task(0,0) as task:
            task.set_Stream (mk.streamtype.log, streamprinter)
            # set mosek params
            task.putdouparam(dparam.intpnt_nl_tol_rel_gap, rtol)
            for (param, value) in control.items():
                if param[:6] == "iparam":
                    task.putintparam(eval(param), value)
                elif param[:6] == "dparam":
                    task.putdouparam(eval(param), value)
                elif str(param)[:6] == "sparam":
                    task.putstrparam(eval(param), value)
                else:
                    raise ValueError("invalid MOSEK parameter: {}".format(param))
            task.appendvars(numvar)
            task.appendcons(numcon)
            # logic follows from REBayes::KWDual hereafter
            task.putobjsense(mk.objsense.minimize)
            ## objective coefficients
            task.putclist(np.arange(numvar), np.zeros(numvar))
            ## constraint matrix coefficients
            task.putacolslice(0, numvar,
                              AT_sparse.indptr[:-1], AT_sparse.indptr[1:],
                              AT_sparse.indices, AT_sparse.data)
            ## bounds
            task.putconboundslice(0, numcon, bkc, blc, buc)
            task.putvarboundslice(0, numvar, bkx, blx, bux)
            ## solver setup
            task.putSCeval(opro = opro, oprjo = oprjo, oprfo = oprfo,
                           oprgo = oprgo, oprho = oprho)
            ## run
            task.optimize()
            ## get results
            task.getsolutionslice(mk.soltype.itr, mk.solitem.suc, 0, numcon, res)
            converged = task.getsolsta(mk.soltype.itr) in [mk.solsta.optimal, mk.solsta.near_optimal]
    return normalize(res), converged

if __name__ == '__main__':
    res = mixIP(np.loadtxt("data.txt"), np.loadtxt("prior.txt"))
    print('\t'.join(list(map(str, res[0]))))
    print(res[1])
