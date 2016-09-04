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


#include "Matrix.hpp"


class Basis
{
public:
	struct ErrorInfo
	{
		double volume;
		double discrep;
		double error;
		double max_error;
	};

private:
	ii index;       // index of this basis in the serialised tree
	ii parent_index;  // parent node
	ii child_count; // how many children synthesise to this node
	bool transient; // if transient, coefficients not part of fitting

public:
	Basis(std::vector<Basis*>& bases, bool transient = false, ii parent_index = -1);
	virtual ~Basis() {}

	virtual void synthesis(Matrix& f, const Matrix& c, bool accum = true) const = 0;
	virtual void analysis(Matrix& c_err, const Matrix& f_err, bool a_sqrd = false) const = 0;

	virtual ErrorInfo error(Matrix& f_err, const Matrix& f, const Matrix& g) const;
	virtual void shrink(Matrix& c_err, const Matrix& c, const Matrix& l1, const Matrix& l2, fp shrinkage) const;

	virtual li get_nc() const = 0;

	ii get_index() { return index; }
	ii get_parent_index() { return parent_index; }
	bool is_transient() { return transient; }
};


#endif

