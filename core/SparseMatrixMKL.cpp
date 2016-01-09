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


#include "SparseMatrixMKL.hpp"
#include <iomanip>
using namespace std;


SparseMatrix::
SparseMatrix() :
	a(0)
{
}


SparseMatrix::
~SparseMatrix()
{
	if (a)
	{
		delete[] a;
		delete[] ia;
		delete[] ja;
	}
}


void
SparseMatrix::
init(ii _m, ii _n, std::vector<fp>& acoo, std::vector<ii>& rowind, std::vector<ii>& colind)
{
	m = _m;
	n = _n;

	if (a)
	{
		delete[] a;
		delete[] ia;
		delete[] ja;
	}

	ii nnz = acoo.size();
	a = new fp[nnz];
	ia = new ii[m + 1];
	ja = new ii[nnz];

	ii job[] = { 2, 0, 0, 0, nnz, 0 };
	ii info;
	mkl_scsrcoo(job, &m, a, ja, ia, &nnz, acoo.data(), rowind.data(), colind.data(), &info);
}


void
SparseMatrix::
mult(fp* bs, const fp* xs, bool transpose, bool accum) const
{
	fp alpha = 1.0;
	fp beta = accum ? 1.0 : 0.0;
	if (transpose)
	{
		mkl_scsrmv("T", &m, &n, &alpha, "G**C", a, ja, ia, &(ia[1]), const_cast<fp*>(xs), &beta, bs);
	}
	else
	{
		mkl_scsrmv("N", &m, &n, &alpha, "G**C", a, ja, ia, &(ia[1]), const_cast<fp*>(xs), &beta, bs);
	}
}


void
SparseMatrix::
sqr_mult(fp* bs, const fp* xs, bool transpose, bool accum) const
{
	fp* sqr_a = new fp[ia[m]];
	vsSqr(ia[m], a, sqr_a);

	fp alpha = 1.0;
	fp beta = accum ? 1.0 : 0.0;
	if (transpose)
	{
		mkl_scsrmv("T", &m, &n, &alpha, "G**C", sqr_a, ja, ia, &(ia[1]), const_cast<fp*>(xs), &beta, bs);
	}
	else
	{
		mkl_scsrmv("N", &m, &n, &alpha, "G**C", sqr_a, ja, ia, &(ia[1]), const_cast<fp*>(xs), &beta, bs);
	}

	delete[] sqr_a;
}


void
SparseMatrix::
print(ostream& out) const
{
	li size = sizeof(this);
	size += sizeof(fp) * ia[m];
	size += sizeof(ii) * (m + 1);
	size += sizeof(ii) * ia[m];

	cout << "[" << m << "," << n << "]:" << ia[m] << " mem=" << setprecision(2) << fixed << size / (1024.0*1024.0) << "Mb";
}