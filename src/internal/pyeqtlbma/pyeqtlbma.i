// pyeqtlbma.i
// Gao Wang (c) 2015

%module pyeqtlbma

%{
#include "pyeqtlbma.hpp"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_cdf.h"
%}

// gsl functions initialization
%inline %{
void my_error_handler(const char * reason, const char * file,
                       int line, int gsl_errno)
{
    fprintf(stderr, "GSL Error %d:\t%s", gsl_errno, reason);
}
int gsl_initialize()
{
    gsl_set_error_handler(&my_error_handler);
    return 1;
}
%}

%init
%{
    pyeqtlbma::initialize();
    gsl_initialize();
%}

%include exception.i

%exception
{
    try
    {
        $function
    }
    catch(pyeqtlbma::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(pyeqtlbma::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(pyeqtlbma::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(pyeqtlbma::RuntimeError e)
    {
        SWIG_exception(SWIG_RuntimeError, e.message());
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


%newobject *::clone;

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

namespace std
{
  %template(vectors)    vector<string>;
  %template(vectorf)    vector<double>;
  %template(vectori)    vector<int>;
  %template(matrixi)    vector<vector<int> >;
  %template(matrixf)    vector<vector<double> >;
  %template(dict_vectors) map<string, vector<string> >;
  %template(dict_vectorf) map<string, vector<double> >;
  %template(dict_vectori) map<string, vector<int> >;
  %template(dict_matrixi) map<string, vector<vector<int> > >;
  %template(dict_matrixf) map<string, vector<vector<double> > >;
  %template(dict_dict_matrixf) map<string, map<string, vector<vector<double> > > >;
  %template(dict_dict_vectors) map<string, map<string, vector<string> > >;
  %template(dict_string) map<string, string>;
  %template(dict_int) map<string, int>;
  %template(dict_float) map<string, double>;
}

%ignore pyeqtlbma::PyAction::PyAction(const PyAction & rhs);
%ignore pyeqtlbma::PyFunc;

%include "pyeqtlbma.hpp"

// gsl functions
extern double gsl_cdf_gaussian_P(double x, double sigma);
extern double gsl_cdf_gaussian_Q(double x, double sigma);
extern double gsl_cdf_gaussian_Pinv(double P, double sigma);
extern double gsl_cdf_gaussian_Qinv(double Q, double sigma);
extern double gsl_cdf_ugaussian_P(double x);
extern double gsl_cdf_ugaussian_Q(double x);
extern double gsl_cdf_ugaussian_Pinv(double P);
extern double gsl_cdf_ugaussian_Qinv(double Q);
