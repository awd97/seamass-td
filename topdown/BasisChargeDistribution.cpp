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


#include "BasisChargeDistribution.hpp"
#include <limits>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


double proton_mass = 1.007276466879; // in Daltons
double isotope_interval = 1.00286084990559; // average difference between monoisotope and second isotope peaks


BasisUniformChargeDistribution::
BasisUniformChargeDistribution(vector<Basis*>& bases,
                        const vector< vector<double> >& mzs,
                        const vector<fp>& gs, const vector<li>& _is, const vector<ii>& js,
                        ii mass_res, ii max_z, double peak_fwhm,
                        bool transient) :
    Basis(bases, 2, 0, transient),
    is(_is),
	rc(isotope_interval / pow(2.0, (double) mass_res))
{
	BSpline bspline(3, 65535);

    ///////////////////////////////////////////////////////////////////////
    // create A as a temporary COO matrix
    
    // calculate indicies of non-empty spectra
    cm.l[1] = numeric_limits<ii>::min();
    cm.o[1] = numeric_limits<ii>::min();
    cm.n[1] = js.size();

    // init arrays
    a.resize(cm.n[1]);
    ia.resize(cm.n[1]);
    ja.resize(cm.n[1]);
    m.resize(cm.n[1]);
    nnz.resize(cm.n[1], 0);

    // find min and max m/z across spectra
    double mz0 = DBL_MAX;
    double mz1 = 0.0;
    for (ii j = 0; j < cm.n[1]; j++)
    {
        m[j] = (ii) (is[j+1] - is[j]);
        mz0 = mzs[js[j]].front() < mz0 ? mzs[js[j]].front() : mz0;
        mz1 = mzs[js[j]].back() > mz1 ? mzs[js[j]].back() : mz1;
    }
	mz0 -= 0.5*peak_fwhm;
	mz1 += 0.5*peak_fwhm;

	// calculate output mass range
	double mass0 = mz0 - proton_mass;
	double mass1 = max_z * (mz1 - proton_mass);
    cm.l[0] = mass_res;
	cm.o[0] = ceil(mass0 * rc);
	cm.n[0] = floor(mass1 * rc) - cm.o[0] + 1;

    // figure out nnz (number of non-zeros)
    #pragma omp parallel for
	for (ii j = 0; j < cm.n[1]; ++j)
	for (ii i = 0; i < m[j]; i++)
	if (gs[is[j] + i] >= 0.0)
	for (ii z = 1; z <= max_z; z++)
	{
		double cf0 = z * (mzs[js[j]][i] - 0.5*peak_fwhm - proton_mass) * rc;
		double cf1 = z * (mzs[js[j]][i + 1] + 0.5*peak_fwhm - proton_mass) * rc;

		ii ci0 = (ii) ceil(cf0);
		ii ci1 = (ii) floor(cf1);

		nnz[j] += ci1 - ci0 + 1;
	}

	// populate coo matrix
    ii done = 0;
    #pragma omp parallel for
    for(ii j = 0; j < cm.n[1]; ++j)
    {
        vector<fp> acoo(nnz[j]);
        vector<ii> rowind(nnz[j]);
        vector<ii> colind(nnz[j]);
        
        ii k = 0;
		for (ii i = 0; i < m[j]; i++)
		{
			if (gs[is[j] + i] >= 0.0)
			for (ii z = 1; z <= max_z; z++)
			{
				double cf0 = z * (mzs[js[j]][i] - 0.5*peak_fwhm - proton_mass) * rc;
				double cf1 = z * (mzs[js[j]][i + 1] + 0.5*peak_fwhm - proton_mass) * rc;

				ii ci0 = (ii)ceil(cf0);
				ii ci1 = (ii)floor(cf1);

				// work out basis coefficients
				for (ii ci = ci0; ci <= ci1; ci++)
				{
					double bin0 = z * (mzs[js[j]][i] - proton_mass) * rc;
					double bin1 = z * (mzs[js[j]][i + 1] - proton_mass) * rc;

					double basis0 = ci - z*rc*0.5*peak_fwhm;
					double basis1 = ci + z*rc*0.5*peak_fwhm;

					// intersection of bin and basis
					double b0 = bin0 > basis0 ? (bin0 - basis0) / (basis1 - basis0) : 0.0;
					double b1 = bin1 < basis1 ? (bin1 - basis0) / (basis1 - basis0) : 1.0;

					// basis coefficient b is _integral_ of area under b-spline basis
					double b = bspline.ibasis(b1) - bspline.ibasis(b0);
					if (b <= FLT_MIN) b = FLT_MIN;

					acoo[k] = b;
					rowind[k] = i;
					colind[k] = ci - cm.o[0];

					k++;
				}
			}

			// display progress update
			#pragma omp critical
			{
				done++;
				if (done % 10000 == 0)
				{
					for (int i = 0; i < 256; ++i) cout << '\b';
					cout << index << " BasisUniformChargeDistribution " << setw(1 + (int)(log10((float)m[j] * cm.n[1]))) << done << "/" << m[j]*cm.n[1] << " " << flush;
				}
			}
		}
        
        // create A and free coo
        a[j].resize(nnz[j]);
        ia[j].resize(m[j]+1);
        ja[j].resize(nnz[j]);
        
        ii job[] = {2, 0, 0, 0, nnz[j], 0}; ii info;
        mkl_scsrcoo(job, &m[j], a[j].data(), ja[j].data(), ia[j].data(), &(nnz[j]), acoo.data(), rowind.data(), colind.data(), &info);
    }
    for (int i = 0; i < 256; ++i) cout << '\b';
    
    cout << index << " BasisUniformChargeDistribution ";
    cm.print(cout);
    li size = is.back();
    for (ii j = 0; j < a.size(); j++) size += 2*nnz[j]+1;
    cout << " mem=" << setprecision(2) << fixed << (sizeof(this)+size*sizeof(fp))/(1024.0*1024.0) << "Mb";
    if (transient) cout << " (t)";
    cout << endl;
}


BasisUniformChargeDistribution::~BasisUniformChargeDistribution()
{
}


void
BasisUniformChargeDistribution::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum)
{
    static fp alpha = 1.0;
    fp beta = accum ? 1.0 : 0.0;
    # pragma omp parallel for
    for (li j = 0; j < cm.n[1]; j++)
    {
        fp* c = const_cast<fp*>(&(cs.data()[j*cm.n[0]]));
        mkl_scsrmv("N", &(m[j]), &(cm.n[0]), &alpha, "G**C", a[j].data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), c, &beta, &(fs.data()[is[j]]));
    }
}


void
BasisUniformChargeDistribution::
analysis(vector<fp>& es, const vector<fp>& fs)
{
    static fp alpha = 1.0, beta = 0.0;
    //# pragma omp parallel for
    for (li j = 0; j < cm.n[1]; j++)
    {
        fp* f = const_cast<fp*>(&(fs.data()[is[j]]));
         mkl_scsrmv("T", &(m[j]), &(cm.n[0]), &alpha, "G**C", a[j].data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
    }
}


void
BasisUniformChargeDistribution::
l2norm(vector<fp>& es, const vector<fp>& fs)
{
    static fp alpha = 1.0, beta = 0.0;
    //# pragma omp parallel for
    for (li j = 0; j < cm.n[1]; j++)
    {
        vector<fp> a2(nnz[j]);
        vsSqr(nnz[j], a[j].data(), a2.data());
        fp* f = const_cast<fp*>(&(fs.data()[is[j]]));
        mkl_scsrmv("T", &(m[j]), &(cm.n[0]), &alpha, "G**C", a2.data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
    }
}


void
BasisUniformChargeDistribution::
write_cs(const std::vector<fp>& cs)
{
	ofstream ofs("cs.csv");

	ofs << "index,mass,intensity" << endl;

	for (ii i = 0; i < cs.size(); ++i)
	{
		if (cs[i] > 0.0 || i > 0 && cs[i - 1] > 0.0 || i < cs.size()-1 && cs[i + 1] > 0.0)
			ofs << setprecision(10) << i << "," << (cm.o[0] + i) / rc << "," << cs[i] << endl;
	}

	ofs << endl;
}


BasisFreeformChargeDistribution::
BasisFreeformChargeDistribution(vector<Basis*>& bases,
                                const vector< vector<double> >& mzs,
								const vector<fp>& gs, const vector<li>& _is, const vector<ii>& js,
								ii out_res, ii _max_z, double peak_fwhm,
								bool transient) :
	Basis(bases, 2, 0, transient),
	is(_is),
	out_interval(isotope_interval / pow(2.0, (double)out_res)),
	max_z(_max_z)
{
	BSpline bspline(3, 65535);

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	// calculate indicies of non-empty spectra
	cm.l[1] = numeric_limits<ii>::min();
	cm.o[1] = numeric_limits<ii>::min();
	cm.n[1] = js.size();

	// init arrays
	a.resize(cm.n[1]);
	ia.resize(cm.n[1]);
	ja.resize(cm.n[1]);
	m.resize(cm.n[1]);
	nnz.resize(cm.n[1], 0);

	// find min and max m/z across spectra
	double mz0 = DBL_MAX;
	double mz1 = 0.0;
	for (ii j = 0; j < cm.n[1]; j++)
	{
		m[j] = (ii)(is[j + 1] - is[j]);
		mz0 = mzs[js[j]].front() < mz0 ? mzs[js[j]].front() : mz0;
		mz1 = mzs[js[j]].back() > mz1 ? mzs[js[j]].back() : mz1;
	}
	mz0 -= 0.5*peak_fwhm;
	mz1 += 0.5*peak_fwhm;

	// calculate output mass range
	double mass0 = mz0 - proton_mass;
	double mass1 = max_z * (mz1 - proton_mass);
	cm.l[0] = out_res;
	cm.o[0] = ceil(mass0 / out_interval);
	n = floor(mass1 / out_interval) - cm.o[0] + 1;
	cm.n[0] = n * max_z;

	// figure out nnz (number of non-zeros)
	#pragma omp parallel for
	for (ii j = 0; j < cm.n[1]; ++j)
	{
		for (ii i = 0; i < m[j]; i++)
		if (gs[is[j] + i] >= 0.0)
		for (ii z = 0; z < max_z; z++)
		{
			double cf0 = (z + 1) * (mzs[js[j]][i] - 0.5*peak_fwhm - proton_mass) / out_interval;
			double cf1 = (z + 1) * (mzs[js[j]][i + 1] + 0.5*peak_fwhm - proton_mass) / out_interval;

			ii ci0 = (ii)ceil(cf0);
			ii ci1 = (ii)floor(cf1);

			nnz[j] += ci1 - ci0 + 1;
		}
	}

	// populate coo matrix
	ii done = 0;
	#pragma omp parallel for
	for (ii j = 0; j < cm.n[1]; ++j)
	{
		vector<fp> acoo(nnz[j]);
		vector<ii> rowind(nnz[j]);
		vector<ii> colind(nnz[j]);

		ii k = 0;
		for (ii i = 0; i < m[j]; i++)
		{
			if (gs[is[j] + i] >= 0.0)
			for (ii z = 0; z < max_z; z++)
			{
				// compute the range of peak coefficient indicies [ci0, ci1] that overlap with the raw data bin [mzs[js[j]][i], mzs[js[j]][i+1]]   
				double cf0 = (z + 1) * (mzs[js[j]][i] - 0.5*peak_fwhm - proton_mass) / out_interval;
				double cf1 = (z + 1) * (mzs[js[j]][i + 1] + 0.5*peak_fwhm - proton_mass) / out_interval;
				ii ci0 = (ii)ceil(cf0);
				ii ci1 = (ii)floor(cf1);

				// for each coefficient index, determine what protortion of that peak overlaps with the raw data bin
				// (Gaussian peak shape is approximated by a cubic b-spline basis function) 
				for (ii ci = ci0; ci <= ci1; ci++)
				{
					double bin0 = (z + 1) * (mzs[js[j]][i] - proton_mass) / out_interval;
					double bin1 = (z + 1) * (mzs[js[j]][i + 1] - proton_mass) / out_interval;

					double basis0 = ci - (z + 1)*0.5*peak_fwhm / out_interval;
					double basis1 = ci + (z + 1)*0.5*peak_fwhm / out_interval;

					// intersection of bin and basis
					double b0 = bin0 > basis0 ? (bin0 - basis0) / (basis1 - basis0) : 0.0;
					double b1 = bin1 < basis1 ? (bin1 - basis0) / (basis1 - basis0) : 1.0;

					// basis coefficient b is _integral_ of area under b-spline basis
					double b = bspline.ibasis(b1) - bspline.ibasis(b0);
					if (b <= FLT_MIN) b = FLT_MIN;

					acoo[k] = b;
					rowind[k] = i;
					colind[k] = max_z * (ci - cm.o[0]) + z;

					//ofs << setprecision(10) << b << "," << 0.5*(mzs[js[j]][i] + mzs[js[j]][i + 1]) << "," << get_mass(ci - cm.o[0]) << endl;

					k++;
				}
			}

			// display progress update
			#pragma omp critical
			{
				done++;
				if (done % 10000 == 0)
				{
					for (int i = 0; i < 256; ++i) cout << '\b';
					cout << index << " BasisFreeformChargeDistribution " << setw(1 + (int)(log10((float)m[j] * cm.n[1]))) << done << "/" << m[j] * cm.n[1] << " " << flush;
				}
			}
		}

		// create A and free coo
		a[j].resize(nnz[j]);
		ia[j].resize(m[j] + 1);
		ja[j].resize(nnz[j]);

		ii job[] = { 2, 0, 0, 0, nnz[j], 0 }; ii info;
		mkl_scsrcoo(job, &m[j], a[j].data(), ja[j].data(), ia[j].data(), &(nnz[j]), acoo.data(), rowind.data(), colind.data(), &info);
	}
	for (int i = 0; i < 256; ++i) cout << '\b';

	cout << index << " BasisFreeformChargeDistribution ";
	cm.print(cout);
	li size = is.back();
	for (ii j = 0; j < a.size(); j++) size += 2 * nnz[j] + 1;
	cout << " mem=" << setprecision(2) << fixed << (sizeof(this) + size*sizeof(fp)) / (1024.0*1024.0) << "Mb";
	if (transient) cout << " (t)";
	cout << endl;
}


BasisFreeformChargeDistribution::~BasisFreeformChargeDistribution()
{
}


void
BasisFreeformChargeDistribution::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum)
{
	static fp alpha = 1.0;
	fp beta = accum ? 1.0 : 0.0;
	# pragma omp parallel for
	for (li j = 0; j < cm.n[1]; j++)
	{
		fp* c = const_cast<fp*>(&(cs.data()[j*cm.n[0]]));
		mkl_scsrmv("N", &(m[j]), &(cm.n[0]), &alpha, "G**C", a[j].data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), c, &beta, &(fs.data()[is[j]]));
	}
}


void
BasisFreeformChargeDistribution::
analysis(vector<fp>& es, const vector<fp>& fs)
{
	static fp alpha = 1.0, beta = 0.0;
	//# pragma omp parallel for
	for (li j = 0; j < cm.n[1]; j++)
	{
		fp* f = const_cast<fp*>(&(fs.data()[is[j]]));
		mkl_scsrmv("T", &(m[j]), &(cm.n[0]), &alpha, "G**C", a[j].data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
	}
}


void
BasisFreeformChargeDistribution::
l2norm(vector<fp>& es, const vector<fp>& fs)
{
	static fp alpha = 1.0, beta = 0.0;
	//# pragma omp parallel for
	for (li j = 0; j < cm.n[1]; j++)
	{
		vector<fp> a2(nnz[j]);
		vsSqr(nnz[j], a[j].data(), a2.data());
		fp* f = const_cast<fp*>(&(fs.data()[is[j]]));
		mkl_scsrmv("T", &(m[j]), &(cm.n[0]), &alpha, "G**C", a2.data(), ja[j].data(), ia[j].data(), &(ia[j].data()[1]), f, &beta, &(es.data()[j*cm.n[0]]));
	}
}


void
BasisFreeformChargeDistribution::
shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage)
{
	// GROUP-WISE SHRINKAGE! (note - intuitive implemention, not mathematically verified yet)
	#pragma omp parallel for
	for (ii i = 0; i < n; ++i)
	{
		// total of all the coefficients in this group
		double sum_c = 0.0;
		for (ii z = 0; z < max_z; z++)
		{
			if (es[i*max_z + z] >= FLT_MIN && cs[i*max_z + z] >= FLT_MIN) sum_c += cs[i*max_z + z];
		}

		// scale the shrinkage to be protortional to the contribution of this coefficient to the group total
		for (ii z = 0; z < max_z; z++)
		{
			double scale = cs[i*max_z + z] / sum_c;
			if (es[i*max_z + z] >= FLT_MIN && cs[i*max_z + z] >= FLT_MIN)
			{
				es[i*max_z + z] *= cs[i*max_z + z] / (scale * shrinkage * l2[i*max_z + z] + wcs[i*max_z + z]);
			}
			else
			{
				es[i*max_z + z] = 0.0;
			}
		}
	}
}


void
BasisFreeformChargeDistribution::
write_cs(const std::vector<fp>& cs)
{
	ofstream ofs("cs.csv");

	vector<fp> sum(n, 0.0);
	for (ii i = 0; i < n; ++i)
	{
		for (ii z = 0; z < max_z; z++) sum[i] += cs[i*max_z + z];
	}

	ofs << "index,mass";
	for (ii z = 0; z < max_z; z++) ofs << ",z=" << z+1;
	ofs << setprecision(10) << endl;

	for (ii i = 0; i < n; ++i)
	{
		if (sum[i] > 0.0 || i > 0 && sum[i - 1] > 0.0 || i < n - 1 && sum[i + 1] > 0.0)
		{
			ofs << i << "," << (cm.o[0] + i) * out_interval;
			for (ii z = 0; z < max_z; z++) ofs << "," << cs[i*max_z + z];
			ofs << endl;
		}
	}

	ofs << endl;
}
