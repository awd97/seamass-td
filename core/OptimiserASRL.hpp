//
// Original author: Andrew Dowsey <andrew.dowsey <a.t> liverpool.ac.uk>
//
// Copyright (C) 2015  biospi Laboratory, EEE, University of Liverpool, UK
//
// This file is part of seaMass-TD.
//
// seaMass-TD is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// seaMass is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with seaMass-TD.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef _SEAMASS_CORE_OPTIMISERASRL_HPP_
#define _SEAMASS_CORE_OPTIMISERASRL_HPP_


#include "Basis.hpp"


class OptimiserASRL
{
public:
    const std::vector<Basis*>& bases;
    const Matrix& g;
 
	std::vector<Matrix> cs;
	std::vector<Matrix> l1s;
	std::vector<Matrix> l2s;
    
    ii accell;
	std::vector<Matrix> c0s;
	std::vector<Matrix> u0s;
	std::vector<Matrix> q0s;
    
    ii iteration;
	Basis::ErrorInfo info;

	double synthesis_wtime;
	double error_wtime;
	double analysis_wtime;
	double shrinkage_wtime;

	Matrix f;

public:    
	OptimiserASRL(const std::vector<Basis*>& bases, const Matrix& g, ii accell = 2);
	~OptimiserASRL();
    
    double step(ii iteration, double strength);

    void threshold(double threshold);

	std::vector<Matrix>& get_cs() { return cs; }
	const std::vector<Matrix>& get_l1s() const { return l1s; }
	const std::vector<Matrix>& get_l2s() const { return l2s; }
	const Basis::ErrorInfo& get_info() const { return info; }
};


#endif

