// pyeqtlbma.i
// Gao Wang (c) 2015

%module pyeqtlbma

%{
#include "pyeqtlbma.h"
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
    catch(vtools::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(vtools::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(vtools::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(vtools::RuntimeError e)
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
    %template(vectora)    vector<vtools::BaseAction * >;
}

%ignore vtools::PyAction::PyAction(const PyAction & rhs);
%ignore vtools::PyFunc;

%include "pyeqtlbma.h"

// gsl functions
extern double gsl_cdf_gaussian_P(double x, double sigma);
extern double gsl_cdf_gaussian_Q(double x, double sigma);
extern double gsl_cdf_gaussian_Pinv(double P, double sigma);
extern double gsl_cdf_gaussian_Qinv(double Q, double sigma);
extern double gsl_cdf_ugaussian_P(double x);
extern double gsl_cdf_ugaussian_Q(double x);
extern double gsl_cdf_ugaussian_Pinv(double P);
extern double gsl_cdf_ugaussian_Qinv(double Q);
