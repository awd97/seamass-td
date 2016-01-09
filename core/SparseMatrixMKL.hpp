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


#ifndef _SEAMASS_CORE_SPARSEMATRIX_HPP_
#define _SEAMASS_CORE_SPARSEMATRIX_HPP_


#include <mkl.h>
#include <iostream>
#include <vector>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected indexing integer size (32 or 64 bit)
typedef long long li; // 64bit integer


class SparseMatrix
{
private:
	ii m;
	ii n;

	fp* a;
	ii* ia;
	ii* ja;

public:
	SparseMatrix();
	~SparseMatrix();

	void init(ii m, ii n, std::vector<fp>& acoo, std::vector<ii>& rowind, std::vector<ii>& colind);

	void mult(fp* bs, const fp* xs, bool transpose = false, bool accum = false) const;
	void sqr_mult(fp* bs, const fp* xs, bool transpose = false, bool accum = false) const;

	ii get_m() const { return m; }
	ii get_n() const { return n; }

	void print(std::ostream& out) const;
};


#endif

