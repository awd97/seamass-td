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
	ii child_count; // how many children synthesise to this node
	Basis* parent;  // parent node
	bool transient; // if transient, coefficients not part of fitting

public:
	Basis(std::vector<Basis*>& bases, Basis* parent = NULL, bool transient = false);
	virtual ~Basis() {}

	virtual void synthesis(SparseMatrix& fs, const SparseMatrix& cs, bool accum = true) const = 0;
	virtual void analysis(SparseMatrix& es, const SparseMatrix& fs) const = 0;
	virtual void l2norm(SparseMatrix& es, const SparseMatrix& fs) const = 0;
	virtual ErrorInfo error(SparseMatrix& fs, const SparseMatrix& gs) const;
	virtual void shrink(SparseMatrix& es, const SparseMatrix& cs, const SparseMatrix& l2, const SparseMatrix& wcs, double shrinkage) const;
	virtual li get_nc() const = 0;

	ii get_index() const { return index; }
	Basis* get_parent() const { return parent; }
	bool is_transient() const { return transient; }
};


class BasisBSpline : public Basis
{
public:
	struct MeshInfo
	{
		MeshInfo(ii d);
		~MeshInfo();

		li size() const;
		void print(std::ostream& out) const;
		void operator=(const MeshInfo& bo);

		ii d;				 // dimension of the b-spline coefficients mesh
		std::vector<ii> ls;  // dyadic level for each dimension
		std::vector<ii> os;  // coefficient offset for each dimension
		std::vector<ii> ns;  // number of coefficients for each dimension
		ii m;                // number of meshes
	};

protected:
	MeshInfo mi;

public:
	BasisBSpline(std::vector<Basis*>& bases, ii dims, Basis* parent = NULL, bool transient = false);
	virtual ~BasisBSpline() {}

	li get_nc() { return mi.size(); }
};


#endif

