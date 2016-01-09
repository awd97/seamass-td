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


#include "BasisIsotopeDistribution.hpp"
#include <limits>
#include <map>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cfloat>
using namespace std;


BasisIsotopeDistribution::
BasisIsotopeDistribution(vector<Basis*>& bases, BasisChargeDistribution* parent,
	                     ii out_res, ii factor_res, ii max_z, bool transient) :
	Basis(bases, parent, transient),
	ns(parent->get_ns()),
	nc(parent->get_nc())
{
	cout << get_index() << " BasisIsotopeDistribution " << flush;

	ifstream ifs("seamass-td.factors");
	//ofstream ofs("grr2");
	multimap< ii, vector<fp> > factors;
	ii ci0 = parent->get_ci0s().front() / (1 << (out_res - factor_res));
	ii ci1 = parent->get_ci1s().back() / (1 << (out_res - factor_res));

	while (ifs.good())
	{
		char comma;
		ii ci, z, id;
		ifs >> ci;
		if (ci < ci0) { ifs.ignore(65535, '\n'); continue; }
		if (ci > ci1) break;

		ifs >> comma >> z >> comma >> id;
		map< ii, vector<fp> >::iterator elem = factors.insert(pair< ii, vector<fp> >(ci, vector<fp>()));
		for (ii j = 0; j < 21; j++)
		{
			fp isotope;
			ifs >> comma >> isotope;
			if (elem->second.size() > 0 && isotope <= elem->second.back() && isotope < 0.00001)
			{
				ifs.ignore(65535, '\n');
				break;
			}
			else
			{
				elem->second.push_back(isotope);
			}
		}
		//ofs << ci << " " << z << " " << id << " " << elem->second.size() << endl;

		if ((ci-ci0) % 1000 == 0)
		{
			for (int i = 0; i < 256; ++i) cout << '\b';
			cout << get_index() << " BasisIsotopeDistribution " << setw(1 + (int)(log10((float)ci1))) << ci-ci0 << "/" << ci1-ci0 << " " << flush;
		}
	}
	for (int i = 0; i < 256; ++i) cout << '\b';

	///////////////////////////////////////////////////////////////////////
	// create A as a temporary COO matrix

	ii nnz = 0;
	for (ii z = 0; z < parent->get_cos().size() - 1; z++)
	for (ii i = parent->get_cos()[z]; i < parent->get_cos()[z + 1]; i++)
	{
		ii ci = i - parent->get_cos()[z] + parent->get_ci0s()[z];
		ii fi = ci / (1 << (out_res - factor_res));
		
		bool yes = true;
		for (pair<multimap< ii, vector<fp> >::iterator, multimap< ii, vector<fp> >::iterator> fs = factors.equal_range(fi); fs.first != fs.second; ++fs.first)
		{
			//if (yes) ofs << "z=" << z << " ci=" << ci << " fi=" << fi << endl;
			//yes = false;
			for (ii j = 0; j < fs.first->second.size(); j++)
			{
				ii oi = i + j * (1 << out_res);
				if (oi < parent->get_cos()[z + 1])
				{
					nnz++;
					//ofs << " " << fs.first->second[i];
				}
			}
			break;
			//ofs << endl;
		}
		//cs[j*cos.back() + i];
	}

	cout << endl << "nnz=" << nnz << endl << endl;;
	vector<fp> acoo(nnz);
	vector<ii> rowind(nnz);
	vector<ii> colind(nnz);

	ii k = 0;
	for (ii z = 0; z < parent->get_cos().size() - 1; z++)
	for (ii i = parent->get_cos()[z]; i < parent->get_cos()[z + 1]; i++)
	{
		ii ci = i - parent->get_cos()[z] + parent->get_ci0s()[z];
		ii fi = ci / (1 << (out_res - factor_res));

		for (pair<multimap< ii, vector<fp> >::iterator, multimap< ii, vector<fp> >::iterator> fs = factors.equal_range(fi); fs.first != fs.second; ++fs.first)
		{
			for (ii j = 0; j < fs.first->second.size(); j++)
			{
				ii oi = i + j * (1 << out_res);
				if (oi < parent->get_cos()[z + 1])
				{
					acoo[k] = fs.first->second[j];
					rowind[k] = oi;
					colind[k] = i;
					k++;
				}
			}
			break;
		}
	}

	// create A
	a.init(parent->get_cos().back(), parent->get_cos().back(), acoo, rowind, colind);

	cout << "HELLO" << endl;

	cout << get_index() << " BasisIsotopeDistribution ";
	cout << " A=";
	a.print(cout);	
	if (transient) cout << " (t)" << endl; else cout << " nc=" << nc << endl;
}


BasisIsotopeDistribution::~BasisIsotopeDistribution()
{
}


void
BasisIsotopeDistribution::
synthesis(vector<fp>& fs, const vector<fp>& cs, bool accum) const
{
	for (li j = 0; j < ns; j++)
	{
		a.mult(&(fs.data()[j*a.get_m()]), &(cs.data()[j*a.get_n()]), false, accum);
	}
}


void
BasisIsotopeDistribution::
analysis(vector<fp>& es, const vector<fp>& fs) const
{
	for (li j = 0; j < ns; j++)
	{
		a.mult(&(es.data()[j*a.get_n()]), &(fs.data()[j*a.get_m()]), true);
	}
}


void
BasisIsotopeDistribution::
l2norm(vector<fp>& es, const vector<fp>& fs) const
{
	for (li j = 0; j < ns; j++)
	{
		a.sqr_mult(&(es.data()[j*a.get_n()]), &(fs.data()[j*a.get_m()]), true);
	}
}


void
BasisIsotopeDistribution::
shrink(std::vector<fp>& es, const std::vector<fp>& cs, const std::vector<fp>& l2, const std::vector<fp>& wcs, double shrinkage) const
{
	get_parent()->shrink(es, cs, l2, wcs, shrinkage);
}
