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
#include "BasisChargeDistribution.hpp"
#include "BasisIsotopeDistribution.hpp"
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
		    vector<double>& rts, vector< vector<double> >& mzs, vector< vector<double> >& intensities,
		    int out_res, int max_z, double max_peak_width, int shrink, int tol,
		    int threads, int debug)
	{
        cout << "SIZE: " << sizeof(ii) << endl;
        
		double start = omp_get_wtime();
		int _threads = omp_get_num_threads();
		omp_set_num_threads(threads);

		cout << "Input id=" << id << " config_id=" << config_id << " instrument_type=" << instrument_type << endl;
		cout << endl;

		////////////////////////////////////////////////////////////////////////////////////
		// INIT RAW DATA AND MZ BASIS

		// Ensure the raw data is in binned format and compute exposures
		vector<fp> exposures;
		utils::bin_mzs_intensities(mzs, intensities, instrument_type, exposures);

		// Convert intensities into format used in algorithm for gs
		vector<fp> gs; vector<li> is; vector<ii> js;
		utils::create_gs(gs, is, js, intensities);
		for (ii j = 0; j < (ii)intensities.size(); j++) vector<double>().swap(intensities[j]);

		// Create our tree of bases
		vector<Basis*> bases;
		ii order = 3; // B-spline order

		// Construct BasisChargeDistribution, which gets added as the root node of bases
		BasisChargeDistribution bChargeDistribution(bases, mzs, gs, is, js, 1.00286084990559, out_res, 0, max_z, max_peak_width, true);
		for (ii j = 0; j < (ii)mzs.size(); j++) vector<double>().swap(mzs[j]);

		// Construct BasisIsotopeDistribution, which gets added to bases with bChargeDistribution as parent
		BasisIsotopeDistribution bIsotopeDistribution(bases, &bChargeDistribution, out_res, 0, max_z);

		////////////////////////////////////////////////////////////////////////////////////
		// OPTIMISATION
		double thres = 1.0; // L0 threshold

		OptimiserASRL optimiser(bases, gs, 2);

		double shrinkage = pow(2.0, (double)shrink);
		double tolerance = pow(2.0, (double)tol);
		double grad = DBL_MAX;

		// l1
		li nc = 0;
		for (ii j = 0; j < (ii)bases.size(); j++)
		if (!bases[j]->is_transient())
			nc += bases[j]->get_nc();
		cout << " L1 nc=" << nc << " shrinkage=" << shrink << ":" << fixed << setprecision(2) << shrinkage << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;

		for (ii i = 0; grad > tolerance; i++)
		{
			grad = optimiser.step(i, shrinkage);
            if (isnan(grad)) grad = 1.0;

			li nnz = 0;
			for (ii j = 0; j < (ii)bases.size(); j++)
			if (!bases[j]->is_transient())
			for (ii i = 0; i < (ii)bases[j]->get_nc(); i++)
			if (optimiser.get_cs()[j][i] > thres)
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
		}

		// l0
		cout << " L0 threshold=" << fixed << setprecision(2) << thres << " tolerance=" << tol << ":" << setprecision(6) << tolerance << endl;

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
			if (optimiser.get_cs()[j][i] > thres)
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
		}

		//////////////////////////////////////////////////////////////////////////////////
		// OUTPUT
		bIsotopeDistribution.write_cs(optimiser.get_cs()[bIsotopeDistribution.get_index()]);
		
		vector<fp> ts(bChargeDistribution.get_nc());
		bIsotopeDistribution.synthesis(ts, optimiser.get_cs()[bIsotopeDistribution.get_index()]);
		bChargeDistribution.write_cs(ts);

		vector<fp> fs(gs.size());
		bChargeDistribution.synthesis(fs, ts);
		ofstream ofs("fs.csv");
		for (ii i = 0; i < fs.size(); i++) ofs << fs[i] << endl;
		ofstream ofs2("gs.csv");
		for (ii i = 0; i < gs.size(); i++) ofs2 << gs[i] << endl;

		////////////////////////////////////////////////////////////////////////////////////
		omp_set_num_threads(_threads);
		cout << "Duration: " << (omp_get_wtime() - start) / 60.0 << "mins" << endl;
	}


}

