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
protected:
    const std::vector<Basis*>& bases;
    std::vector<fp>& gs;
 
    std::vector< std::vector<fp> > cs;
    std::vector< std::vector<fp> > wcs;
    std::vector< std::vector<fp> > l2;
    
    ii accell;
    std::vector< std::vector<fp> > c0s;
    std::vector< std::vector<fp> > u0s;
    std::vector< std::vector<fp> > q0s;
    
    ii iteration;
	Info info;

public:    
	OptimiserASRL(const std::vector<Basis*>& _bases,
                  std::vector<fp>& gs,
                  ii accell = 2);
	~OptimiserASRL();
    
    double step(ii iteration, double strength);

    void threshold(double threshold);

    std::vector< std::vector<fp> >& get_cs() { return cs; }
	const std::vector< std::vector<fp> >& get_l2() const { return l2; }
	const Info& get_info() const { return info; }
};


#endif

