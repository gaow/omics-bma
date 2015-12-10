// pyeqtlbma.hpp
// Gao Wang (c) 2015

#ifndef __PYEQTLBMA_H_
#define __PYBEQTLBMA_H_

#include "libeqtlbma.hpp"

namespace pyeqtlbma {

typedef std::vector<int> vectori;
typedef std::vector<std::vector<int> > matrixi;
typedef std::vector<std::string> vectors;
typedef std::vector<double> vectorf;
typedef std::vector<std::vector<double> > matrixf;
typedef std::map<std::string, std::vector<double> > dict_vectorf;
typedef std::map<std::string, std::vector<std::vector<double> > > dict_matrixf;
typedef std::map<std::string, std::vector<int> > dict_vectori;
typedef std::map<std::string, std::vector<std::vector<int> > > dict_matrixi;
typedef std::map<std::string, std::vector<std::string> > dict_vectors;
typedef std::map<std::string, std::string> dict_string;
typedef std::map<std::string, int> dict_int;
typedef std::map<std::string, double> dict_float;
typedef std::map<std::string,
                 std::map<std::string,
                          std::vector<std::vector<double> > > >
    dict_dict_matrixf;

int eqtlbma_bf(
	const dict_string & param_s,
	const dict_int & param_i,
	const dict_float & param_f,
	const dict_dict_matrixf & sstats,
	const vectors & subgroups_tokeep);

}
#endif
