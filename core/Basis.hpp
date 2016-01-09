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


#ifndef _SEAMASS_CORE_BASIS_HPP_
#define _SEAMASS_CORE_BASIS_HPP_


#include "core.hpp"
#include <iostream>


/*struct CoeffsMetadata
{
	CoeffsMetadata(ii d);
	~CoeffsMetadata();

	li size() const;
	void print(std::ostream& out) const;
	void operator=(const CoeffsMetadata& cm);

	ii d;        // dimension of the output coefficients
	std::vector<ii> l; // coefficient dyadic level for each dimension
	std::vector<ii> o; // coefficient offset
	std::vector<ii> n; // number of coefficients per set
};*/


class Basis
{
protected:
    ii index;       // index of this basis in the serialised tree
    Basis* parent;  // parent node
    ii child_count; // how many children synthesise to this node
    bool transient; // if transient, coefficients not part of fitting
    //CoeffsMetadata cm;
  
    double volume; // statistic of the last error() call
    double discrep; // statistic of the last error() call
    double erro; // statistic of the last error() call
    double maxerr; // statistic of the last error() call

public:
	Basis(std::vector<Basis*>& bases,
          ii dimensions,
          Basis* parent = NULL, bool
          transient = false);
    virtual ~Basis() {}
    
    virtual void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true) = 0;
    virtual void analysis(std::vector<fp>& es, const std::vector<fp>& fs) = 0;
    virtual void l2norm(std::vector<fp>& es, const std::vector<fp>& fs) = 0;
	virtual void error(std::vector<fp>& fs, const std::vector<fp>& gs);
	virtual void shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage);
	virtual li get_nc() = 0;

    ii get_index() { return index; }
    Basis* get_parent() { return parent; }
    bool is_transient() { return transient; }
    
    // some statistics of the last error() call
    double get_volume() { return volume; }
    double get_discrep() { return discrep; }
    double get_error() { return erro; }
    double get_maxerror() { return maxerr; }
};


#endif

