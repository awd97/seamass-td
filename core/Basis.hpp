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


struct Info
{
	double volume;
	double discrep;
	double error;
	double max_error;
};


class Basis
{
private:
    ii index;       // index of this basis in the serialised tree
    Basis* parent;  // parent node
    ii child_count; // how many children synthesise to this node
    bool transient; // if transient, coefficients not part of fitting

public:
	Basis(std::vector<Basis*>& bases,
          Basis* parent = NULL, bool
          transient = false);
    virtual ~Basis() {}
    
    virtual void synthesis(std::vector<fp>& fs, const std::vector<fp>& cs, bool accum = true) const = 0;
	virtual void analysis(std::vector<fp>& es, const std::vector<fp>& fs) const = 0;
	virtual void l2norm(std::vector<fp>& es, const std::vector<fp>& fs) const = 0;
	virtual Info error(std::vector<fp>& fs, const std::vector<fp>& gs) const;
	virtual void shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const;
	virtual li get_nc() const = 0;

    ii get_index() const { return index; }
    Basis* get_parent() const { return parent; }
    bool is_transient() const { return transient; }
};


#endif

