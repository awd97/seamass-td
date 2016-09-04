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


#ifndef _SEAMASS_MATRIX_SPARSEMKL_HPP_
#define _SEAMASS_MATRIX_SPARSEMKL_HPP_


#include <mkl.h>
#include <mkl_spblas.h>
#include <iostream>
#include <vector>


typedef float fp; // fp is the selected floating point precision (float or double)
typedef MKL_INT ii; // ii is the selected indexing integer size (32 or 64 bit)
typedef long long li; // 64bit integer


class MatrixSparseMKL
{
private:
	sparse_matrix_t mat;

public:
	MatrixSparseMKL();
	~MatrixSparseMKL();

	void init(const std::vector<fp>& xs);
	void init(ii n, fp v = 0.0);
	void init(ii m, ii n, const std::vector<fp>& acoo, const std::vector<ii>& rowind, const std::vector<ii>& colind);
	
	bool operator!() const;
	void copy(const MatrixSparseMKL& m);
	void free();

	ii m() const;
	ii n() const;

	void mult(const MatrixSparseMKL& a, const MatrixSparseMKL& x, bool accum, bool a_sqrd);
	void elem_div(const MatrixSparseMKL& a, const MatrixSparseMKL& b);
	void sqrt();
	void shrink(const MatrixSparseMKL& c, const MatrixSparseMKL& l1, const MatrixSparseMKL& l2, fp shrinkage);

	double sum_sqrd() const;
	double sum_sqrd_diffs(const MatrixSparseMKL& a) const;

	friend std::ostream& operator<<(std::ostream& os, const MatrixSparseMKL& mat);

	li mem() const;
	void save(const std::string filename) const;
};

std::ostream& operator<<(std::ostream& os, const MatrixSparseMKL& mat);


#endif

