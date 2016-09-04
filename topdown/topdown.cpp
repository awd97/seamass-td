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


#include "topdown.hpp"
#include "../core/DataMS.hpp"
#include "../core/BasisBSplineMZ.hpp"
#include "../core/BasisBSplineScale.hpp"
#include "../core/OptimiserASRL.hpp"
#include <fstream>
#include <sstream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <omp.h>
using namespace std;


namespace seamass
{


	void
	topdown_notice()
	{
		cout << endl;
		cout << "seaMass-TD : Copyright (C) 2015 : biospi Laboratory, EEE, University of Liverpool, UK" << endl;
		cout << "This program comes with ABSOLUTELY NO WARRANTY." << endl;
		cout << "This is free software, and you are welcome to redistribute it under certain conditions." << endl;
		cout << endl;
	}


	void
	topdown(const std::string& id, const std::string& config_id, int instrument_type,
		    vector<double>& sts, vector< vector<double> >& mzs, vector< vector<double> >& intensities,
		    double out_mass0, double out_mass1, int out_res, int max_z, double max_peak_width, int shrink, int tol,
		    int threads, int debug)
	{
        cout << "seaMass matrix_capacity=" << 8 * sizeof(ii) << "bit" << endl;
        
		double start = omp_get_wtime();
		int _threads = omp_get_num_threads();
		omp_set_num_threads(threads);

		////////////////////////////////////////////////////////////////////////////////////
		// INIT RAW DATA AND BASIS FUNCTIONS

		DataMS input(id, config_id, instrument_type, sts, mzs, intensities);
		for (ii j = 0; j < (ii)intensities.size(); j++) std::vector<double>().swap(intensities[j]); // free intensities

		vector<Basis*> bases;
		BasisBSplineMZ bRoot(bases, input, mzs, out_res);
		for (ii j = 0; j < (ii)mzs.size(); j++) std::vector<double>().swap(mzs[j]); // free mzs

		BasisBSplineScale bScale1(bases, bRoot.get_index(), 0);
		BasisBSplineScale bScale2(bases, bScale1.get_index(), 0);
		BasisBSplineScale bScale3(bases, bScale2.get_index(), 0);

		////////////////////////////////////////////////////////////////////////////////////
		// OPTIMISATION
		double thres = 1.0; // L0 threshold

		OptimiserASRL optimiser(bases, input.get_g(), 2);

		double shrinkage = pow(2.0, (double)shrink);
		double tolerance = pow(2.0, (double)tol);
		double grad = DBL_MAX;

		// l1
		//li nc = 0;
		//for (ii j = 0; j < (ii)bases.size(); j++)
		//if (!bases[j]->is_transient())
		//	nc += bases[j]->get_nc();
		//cout << " L1 nc=" << nc << " shrinkage=" << shrink << ":" << fixed << setprecision(2) << shrinkage << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;

		//bIsotopeDistribution.restrict_range(optimiser.get_cs()[bIsotopeDistribution.get_index()], out_mass0, out_mass1);

		for (ii i = 0; grad > tolerance; i++)
		{
			grad = optimiser.step(i, shrinkage);
            
			cout << "  f: " << fixed << setw(5) << i;
			cout << "  grad: " << setiosflags(ios::fixed) << setprecision(6) << setw(8) << grad << endl;

			/*if (isnan(grad)) grad = 1.0;

			li nnz = 0;
			for (ii j = 0; j < (ii)bases.size(); j++)
			if (!bases[j]->is_transient())
			for (ii i = 0; i < (ii)bases[j]->get_nc(); i++)
			//if (optimiser.get_cs()[j][i] > thres)
			{
				nnz++;
			}

			cout << "  f: " << fixed << setw(5) << i;
			cout << "  err: " << setw(8) << setprecision(5) << optimiser.get_info().error;
			cout << "  max: " << setw(8) << setprecision(1) << optimiser.get_info().max_error;
			cout << "  dis: " << setw(8) << setprecision(5) << optimiser.get_info().discrep;
			cout << "  vol: " << setw(8) << setprecision(5) << optimiser.get_info().volume;
			cout << "  nnz: " << setw(8) << setprecision(5) << nnz;
			cout << "  grad: " << setiosflags(ios::fixed) << setprecision(6) << setw(8) << grad << endl;*/
		}

		input.get_g().save("g.csv");
		optimiser.f.save("f.csv");

		// l0
		/*cout << " L0 threshold=" << fixed << setprecision(2) << thres << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;

		optimiser.threshold(thres);
		grad = DBL_MAX;
		for (ii i = 0; grad > tolerance; i++)
		{
			grad = optimiser.step(i, 0.0);
            if (isnan(grad)) grad = 1.0;

			li nnz = 0;
			for (ii j = 0; j < (ii)bases.size(); j++)
			if (!bases[j]->is_transient())
			for (ii i = 0; i < (ii)bases[j]->get_nc(); i++)
			//if (optimiser.get_cs()[j][i] > thres)
			{
				nnz++;
			}

			cout << "  f: " << fixed << setw(5) << i;
			cout << "  err: " << setw(8) << setprecision(5) << optimiser.get_info().error;
			cout << "  max: " << setw(8) << setprecision(1) << optimiser.get_info().max_error;
			cout << "  dis: " << setw(8) << setprecision(5) << optimiser.get_info().discrep;
			cout << "  vol: " << setw(8) << setprecision(5) << optimiser.get_info().volume;
			cout << "  nnz: " << setw(8) << setprecision(5) << nnz;
			cout << "  grad: " << setiosflags(ios::fixed) << setprecision(6) << setw(8) << grad << endl;
		}*/
        //cout << "a" << endl;
		//////////////////////////////////////////////////////////////////////////////////
		// OUTPUT
		//bIsotopeDistribution.write_cs(optimiser.get_cs()[bIsotopeDistribution.get_index()]);
		//cout << "b" << endl;
		//vector<fp> ts(bChargeDistribution.get_nc());
		//bIsotopeDistribution.synthesis(ts, optimiser.get_cs()[bIsotopeDistribution.get_index()]);
		//cout << "c" << endl;
		//bChargeDistribution.write_cs(optimiser.get_cs()[bChargeDistribution.get_index()]);
        //cout << "d" << endl;
		//vector<fp> fs(gs.size());
		//bChargeDistribution.synthesis(fs, ts);
		//ofstream ofs("fs.csv");
		//for (ii i = 0; i < fs.size(); i++) ofs << fs[i] << endl;
		//ofstream ofs2("gs.csv");
		//for (ii i = 0; i < gs.size(); i++) ofs2 << gs[i] << endl;
        //cout << "e" << endl;
		////////////////////////////////////////////////////////////////////////////////////
		omp_set_num_threads(_threads);
		cout << "Duration: " << (omp_get_wtime() - start) << "seconds" << endl;
	}


}

