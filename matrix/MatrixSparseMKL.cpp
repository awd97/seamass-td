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


#include "MatrixSparseMKL.hpp"
#include <iomanip>
#include <fstream>
using namespace std;


MatrixSparseMKL::
MatrixSparseMKL() :
	mat(0)
{
}


MatrixSparseMKL::
~MatrixSparseMKL()
{
	free();
}


void
MatrixSparseMKL::
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
	mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat);
	mkl_sparse_destroy(t);
}


void
MatrixSparseMKL::
init(ii n, fp v)
{
	free();

	std::vector<fp> xs(n, v);
	std::vector<ii> rowind(n, 0);
	std::vector<ii> colind(n); for (ii i = 0; i < n; i++) colind[i] = i;

	sparse_matrix_t t;
	mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, 1, n, (ii)xs.size(), (ii*)rowind.data(), (ii*)colind.data(), (fp*)xs.data());
	mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat);
	mkl_sparse_destroy(t);
}


void
MatrixSparseMKL::
init(ii m, ii n, const std::vector<fp>& acoo, const std::vector<ii>& rowind, const std::vector<ii>& colind)
{
	free();

	sparse_matrix_t t;
	mkl_sparse_s_create_coo(&t, SPARSE_INDEX_BASE_ZERO, m, n, (ii)acoo.size(), (ii*)rowind.data(), (ii*)colind.data(), (fp*)acoo.data());
	mkl_sparse_convert_csr(t, SPARSE_OPERATION_NON_TRANSPOSE, &mat);
	mkl_sparse_destroy(t);
}


void
MatrixSparseMKL::
copy(const MatrixSparseMKL& a)
{
	free();
	matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	mkl_sparse_copy(a.mat, descr, &mat);
}


bool
MatrixSparseMKL::
operator!() const
{
	return mat == 0;
}


void
MatrixSparseMKL::
free()
{
	if (mat)
	{
		mkl_sparse_destroy(mat);
		mat = 0;
	}
}


void
MatrixSparseMKL::
mult(const MatrixSparseMKL& a, const MatrixSparseMKL& x, bool accum, bool a_sqrd)
{
	if (a_sqrd)
	{
		if (mat && accum)
		{
			cerr << "error: not implemented!" << endl;
		}
		else
		{
			MatrixSparseMKL a2; a2.copy(a);
			ii m, n, *i0s, *i1s, *js; fp* a2s; sparse_index_base_t indexing;
			mkl_sparse_s_export_csr(a2.mat, &indexing, &m, &n, &i0s, &i1s, &js, &a2s);
			ii nnz = i1s[m - 1] - i0s[0];

			vsSqr(nnz, a2s, a2s);
			mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, x.mat, a2.mat, &mat);
		}
	}
	else
	{
		if (mat && accum)
		{
			MatrixSparseMKL t;
			mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, x.mat, a.mat, &t.mat);
			mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat, 1.0f, t.mat, &mat);
		}
		else
		{
			mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, x.mat, a.mat, &mat);
		}
	}
}


void
MatrixSparseMKL::
elem_div(const MatrixSparseMKL& a, const MatrixSparseMKL& b)
{
	ii a_m, a_n, *a_i0s, *a_i1s, *a_js; fp* as; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(a.mat, &indexing, &a_m, &a_n, &a_i0s, &a_i1s, &a_js, &as);
	ii a_nnz = a_i1s[a_m - 1] - a_i0s[0];

	ii b_m, b_n, *b_i0s, *b_i1s, *b_js; fp* bs;
	mkl_sparse_s_export_csr(b.mat, &indexing, &b_m, &b_n, &b_i0s, &b_i1s, &b_js, &bs);
	ii b_nnz = b_i1s[b_m - 1] - b_i0s[0];

	vector<fp> xs(a_n);
	ii a_j = 0;
	ii b_j = 0;
	for (ii j = 0; j < a_n; j++)
	{
		fp a_x, b_x;

		if (a_j < a_nnz && j == a_js[a_j])
		{
			a_x = as[a_j];
			a_j++;
		}
		else
		{
			a_x = 0.0;
		}
		//cout << b_j << endl;
		if (b_j < b_nnz && j == b_js[b_j])
		{
			b_x = bs[b_j];
			b_j++;
		}
		else
		{
			b_x = 0.0;
		}
		cout << "A" << a << " B" << b << endl;
		cout << setprecision(10) << a_x << "/" << b_x << endl;
		xs[j] = b_x > 0.0f ? a_x / b_x : 1.0f;
	}

	init(xs);

	// nb: only first row processed at present
}


void
MatrixSparseMKL::
sqrt()
{
	ii m, n, *i0s, *i1s, *js; fp* as; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(mat, &indexing, &m, &n, &i0s, &i1s, &js, &as);
	ii nnz = i1s[m - 1] - i0s[0];

	vsSqrt(nnz, as, as);
}

void
MatrixSparseMKL::
shrink(const MatrixSparseMKL& c, const MatrixSparseMKL& l1, const MatrixSparseMKL& l2, fp shrinkage)
{
	// w should probably be precomputed...
	MatrixSparseMKL w;
	mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, l2.mat, shrinkage, l1.mat, &w.mat);

	ii w_m, w_n, *w_i0s, *w_i1s, *w_js; fp* ws; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(w.mat, &indexing, &w_m, &w_n, &w_i0s, &w_i1s, &w_js, &ws);
	ii w_nnz = w_i1s[w_m - 1] - w_i0s[0];

	ii c_m, c_n, *c_i0s, *c_i1s, *c_js; fp* cs;
	mkl_sparse_s_export_csr(c.mat, &indexing, &c_m, &c_n, &c_i0s, &c_i1s, &c_js, &cs);
	ii c_nnz = c_i1s[c_m - 1] - c_i0s[0];

	ii m, n, *i0s, *i1s, *js; fp* xs;
	mkl_sparse_s_export_csr(mat, &indexing, &m, &n, &i0s, &i1s, &js, &xs);
	ii nnz = i1s[m - 1] - i0s[0];


	/*ii l1_m, l1_n, *l1_i0s, *l1_i1s, *l1_js; fp* l1s;
	mkl_sparse_s_export_csr(l1.mat, &indexing, &l1_m, &l1_n, &l1_i0s, &l1_i1s, &l1_js, &l1s);
	ii l1_nnz = l1_i1s[w_m - 1] - l1_i0s[0];

	ii l2_m, l2_n, *l2_i0s, *l2_i1s, *l2_js; fp* l2s;
	mkl_sparse_s_export_csr(l2.mat, &indexing, &l2_m, &l2_n, &l2_i0s, &l2_i1s, &l2_js, &l2s);
	ii l2_nnz = l2_i1s[l2_m - 1] - l2_i0s[0];*/


	vector<fp> ys(c_n, 0.0);
	for (ii i = 0; i < c_nnz; i++)
	{
		//cout << endl << setprecision(10) << l1s[i] << "," << l2s[i] << "=" << ws[i] << ":" << shrinkage << endl;

		ys[c_js[i]] = xs[c_js[i]] * cs[i] / ws[c_js[i]];
		//cout << ys[c_js[i]] << "=" << xs[c_js[i]] << "*" << cs[i] << "/" << ws[c_js[i]] << endl;
		//if (i == 1) exit(0);
	}
	init(ys);
}


double
MatrixSparseMKL::
sum_sqrd() const
{
	MatrixSparseMKL a2; a2.copy(*this);
	ii m, n, *i0s, *i1s, *js; fp* a2s; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(a2.mat, &indexing, &m, &n, &i0s, &i1s, &js, &a2s);
	ii nnz = i1s[m - 1] - i0s[0];

	vsSqr(nnz, a2s, a2s);
	double sum = 0.0;
	for (ii i = 0; i < nnz; i++) sum += a2s[i];

	return sum;
}

double
MatrixSparseMKL::
sum_sqrd_diffs(const MatrixSparseMKL& a) const
{
	MatrixSparseMKL a2;
	mkl_sparse_s_add(SPARSE_OPERATION_NON_TRANSPOSE, mat, -1.0, a.mat, &a2.mat);

	ii m, n, *i0s, *i1s, *js; fp* a2s; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(a2.mat, &indexing, &m, &n, &i0s, &i1s, &js, &a2s);
	ii nnz = i1s[m - 1] - i0s[0];

	vsSqr(nnz, a2s, a2s);
	double sum = 0.0;
	for (ii i = 0; i < nnz; i++) sum += a2s[i];

	return sum;
}


ii
MatrixSparseMKL::
n() const
{
	ii m, n, *i0s, *i1s, *js; fp* vs; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(mat, &indexing, &m, &n, &i0s, &i1s, &js, &vs);

	return n;
}


ii
MatrixSparseMKL::
m() const
{
	ii m, n, *i0s, *i1s, *js; fp* vs; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(mat, &indexing, &m, &n, &i0s, &i1s, &js, &vs);

	return m;
}


li
MatrixSparseMKL::
mem() const
{
	if (mat)
	{
		ii m, n, *i0s, *i1s, *js; fp* vs; sparse_index_base_t indexing;
		mkl_sparse_s_export_csr(mat, &indexing, &m, &n, &i0s, &i1s, &js, &vs);
		ii nnz = i1s[m - 1] - i0s[0];

		return (ii) (sizeof(fp) * nnz + sizeof(ii) * (nnz + m + 1));
	}
	else
	{
		return 0;
	}
}


void
MatrixSparseMKL::
save(const std::string filename) const
{
	ofstream ofs(filename);

	ii m, n, *i0s, *i1s, *js; fp* as; sparse_index_base_t indexing;
	mkl_sparse_s_export_csr(mat, &indexing, &m, &n, &i0s, &i1s, &js, &as);
	ii nnz = i1s[m - 1] - i0s[0];

	for (ii i = 0; i < nnz; i++)
	{
		ofs << js[i] << "," << as[i] << endl;
	}
}



ostream&
operator<<(ostream& os, const MatrixSparseMKL& a)
{
	if (a.mat)
	{
		ii m, n, *i0s, *i1s, *js; fp* vs; sparse_index_base_t indexing;
		mkl_sparse_s_export_csr(a.mat, &indexing, &m, &n, &i0s, &i1s, &js, &vs);
		ii nnz = i1s[m - 1] - i0s[0];

		os << "[" << m << "," << n << "]:" << nnz;
	}
	else
	{
		os << "[]";
	}

	return  os;
}

