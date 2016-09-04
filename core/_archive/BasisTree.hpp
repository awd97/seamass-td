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


#ifndef _SEAMASS_CORE_BASISTREE_HPP_
#define _SEAMASS_CORE_BASISTREE_HPP_


#include "Basis.hpp"


class BasisTree : public Basis
{
private:
	std::vector<Basis*> node;    // nodes in the tree
	std::vector<ii> parent;  // parent node
	std::vector<ii> child_count; // how many children synthesise to this node
	std::vector<bool> transient; // if transient, coefficients not part of fitting

public:
	BasisTree();
	virtual ~BasisTree() {}

	virtual void synthesis(SparseMatrix& fs, const SparseMatrix& cs, bool accum = true) const = 0;
	virtual void analysis(SparseMatrix& es, const SparseMatrix& fs) const = 0;
	virtual void l2norm(SparseMatrix& es, const SparseMatrix& fs) const = 0;
	virtual ErrorInfo error(SparseMatrix& fs, const SparseMatrix& gs) const;
	virtual void shrink(SparseMatrix& es, const SparseMatrix& cs, const SparseMatrix& l2, const SparseMatrix& wcs, double shrinkage) const;
	virtual li get_nc() const = 0;
};


#endif

