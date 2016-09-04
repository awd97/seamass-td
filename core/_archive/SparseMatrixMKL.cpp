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
	free();
}


void
SparseMatrix::
init(const std::vector<fp>& xs)
{
	free();

	ii nnz = 0;
	for (ii j = 0; j < (ii)xs.size(); j++)
	{
		if (xs[j] != 0.0) nnz++;
	}

	std::vector<fp> acoo(nnz);
	std::vector<ii> rowind(nnz, 0);
	std::vector<ii> colind(nnz);

	ii nz = 0;
	for (ii j = 0; j < (ii)xs.size(); j++)
	{
		if (xs[j] != 0.0)
		{
			acoo[nz] = xs[j];
			colind[nz] = j;
			nz++;
		}
	}

	sparse_matrix_t t;
	mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, 1, (ii)xs.size(), (ii)acoo.size(), (ii*)rowind.data(), (ii*)colind.data(), (fp*)acoo.data());
	mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &a);
	mkl_sparse_destroy(t);
}


void
SparseMatrix::
init(ii n, fp v)
{
	free();

	std::vector<fp> xs(n, v);
	std::vector<ii> rowind(n, 0);
	std::vector<ii> colind(n); for (ii i = 0; i < n; i++) colind[i] = i;

	sparse_matrix_t t;
	mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, 1, n, (ii)xs.size(), (ii*)rowind.data(), (ii*)colind.data(), (fp*)xs.data());
	mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &a);
	mkl_sparse_destroy(t);
}


void
SparseMatrix::
init(ii m, ii n, const std::vector<fp>& acoo, const std::vector<ii>& rowind, const std::vector<ii>& colind)
{
	free();

	sparse_matrix_t t;
	mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m, n, (ii)acoo.size(), (ii*)rowind.data(), (ii*)colind.data(), (fp*)acoo.data());
	mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &a);
	mkl_sparse_destroy(t);
}


void
SparseMatrix::
operator=(const SparseMatrix& m)
{
	free();
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	mkl_sparse_copy(m.a, descr, &a);
}


bool
SparseMatrix::
operator!() const
{
	return a == 0;
}


void
SparseMatrix::
free()
{
	if (a)
	{
		mkl_sparse_destroy(a);
		a = 0;
	}
}


void
SparseMatrix::
mult(SparseMatrix& b, const SparseMatrix& x, bool accum, bool a_sqrd) const
{
	if (!!b && accum)
	{
		SparseMatrix t;
		mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, x.a, a, &t.a);
		mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, b.a, 1.0f, t.a, &b.a);
	}
	else
	{
		mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, x.a, a, &b.a);
	}
}


/*fp* sqr_acoo = new fp[acoo.size()];
vsSqr(acoo.size(), acoo.data(), sqr_acoo);
mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m, n, acoo.size(), (ii*)rowind.data(), (ii*)colind.data(), sqr_acoo);
mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &sqr_a);
mkl_sparse_destroy(t);
delete sqr_acoo;*/


ostream&
operator<<(ostream& os, const SparseMatrix& mat)
{
	if (mat.a)
	{
		ii m, n;
		ii* i0s, *i1s, *js;
		fp* vs;
		sparse_index_base_t indexing;
		mkl_sparse_s_export_csr(mat.a, &indexing, &m, &n, &i0s, &i1s, &js, &vs);
		ii nnz = i1s[m - 1] - i0s[0];
		os << "[" << m << "," << n << "]:" << nnz << " mem=" << (2 * nnz + m + 1) / 1024.0 / 1024.0 << "Mb";
	}
	else
	{
		os << "[]";
	}

	return  os;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
SparseMatrixOld::
SparseMatrixOld() :
a(0)
{
}


SparseMatrixOld::
~SparseMatrixOld()
{
	if (a)
	{
		delete[] a;
		delete[] ia;
		delete[] ja;
	}
}


void
SparseMatrixOld::
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
SparseMatrixOld::
mult(fp* bs, const fp* xs, bool transpose, bool accum) const
{
	fp alpha = 1.0;
	fp beta = (fp) accum ? 1.0f : 0.0f;
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
SparseMatrixOld::
sqr_mult(fp* bs, const fp* xs, bool transpose, bool accum) const
{
	fp* sqr_a = new fp[ia[m]];
	vsSqr(ia[m], a, sqr_a);

	fp alpha = 1.0;
	fp beta = (fp) accum ? 1.0f : 0.0f;
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
SparseMatrixOld::
print(ostream& out) const
{
	li size = sizeof(this);
	size += sizeof(fp) * ia[m];
	size += sizeof(ii) * (m + 1);
	size += sizeof(ii) * ia[m];

	cout << "[" << m << "," << n << "]:" << ia[m] << " mem=" << setprecision(2) << fixed << size / (1024.0*1024.0) << "Mb";
}
*/