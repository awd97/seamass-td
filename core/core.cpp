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


#include "core.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
using namespace std;


namespace utils
{

	void
	bin_mzs_intensities(vector< vector<double> >& mzs,
	                    vector< vector<double> >& intensities,
						ii instrument_type, vector<fp>& exposures)
	{
		// This modifies the raw data for some limitations of the mzML spec and makes
		// sure the intensities are treated as binned between m/z datapoints.
		//
		// For ToF data, it interpolates the mzs to represent the edges of the bins rather
		//   than the centres
		// For FT data, it converts the sampled data to binned counts by treating the
		//   mz values as the bin edges, and using trapezoid rule to integrate intensity values
		//
		// This is all a bit rough at the moment, should be fitting splines to the data

		// initialise exposures to default of 1 (unit exposure)
		exposures.assign(intensities.size(), 1);

		// if more than one spectrum, ignore last as we do not know its scan end time
		if (mzs.size() > 1) mzs.resize(mzs.size() - 1);
		if (intensities.size() > 1) intensities.resize(intensities.size() - 1);

		if (instrument_type == 1) // ToF
		{
			#pragma omp parallel for
			for (ii j = 0; j < mzs.size(); j++)
			if (mzs[j].size() >= 2)
			{
				// dividing by minimum to get back to ion counts for SWATH data which appears to be scaled to correct for dynamic range restrictions (hack!)
				double minimum = std::numeric_limits<double>::max();

				// we drop the first and last m/z datapoint as we don't know both their bin edges
				for (ii i = 1; i < mzs[j].size(); i++)
				{
					if (intensities[j][i] > 0) minimum = minimum < intensities[j][i] ? minimum : intensities[j][i];

					// linear interpolation of mz extent (probably should be cubic)
					mzs[j][i - 1] = 0.5 * (mzs[j][i - 1] + mzs[j][i]);
					intensities[j][i - 1] = intensities[j][i];
				}
				mzs[j].resize(mzs[j].size() - 1);
				intensities[j].resize(intensities[j].size() - 2);
				exposures[j] = 1.0 / minimum;
				for (ii i = 0; i < intensities[j].size(); i++)
				{
					intensities[j][i] *= exposures[j];
				}
			}
			else
			{
				mzs[j].resize(0);
				intensities[j].resize(0);
			}
		}
		/*else if (instrument_type == 2) // Orbitrap - but disabled for now as doesn't play nicely with merge_bins function
		{
			#pragma omp parallel for
			for (ii j = 0; j < mzs.size(); j++)
			if (mzs[j].size() >= 2)
			{
				// for Orbitrap only, mark zeros as missing data
				for (ii i = 1; i < mzs[j].size(); i++)
				if (intensities[j][i-1] <= 0.0 || intensities[j][i] <= 0.0)
				{
					intensities[j][i-1] = -1.0;
				}
				else
				{
					intensities[j][i-1] = (mzs[j][i] - mzs[j][i-1]) * 0.5 * (intensities[j][i] + intensities[j][i-1]);
				}
				intensities[j].resize(intensities[j].size() - 1);
			}
			else
			{
				mzs[j].resize(0);
				intensities[j].resize(0);
			}
		}*/
		else // FT-ICR (Orbitrap is a type of FT-ICR)
		{
			#pragma omp parallel for
			for (ii j = 0; j < mzs.size(); j++)
			if (mzs[j].size() >= 2)
			{
				for (ii i = 1; i < mzs[j].size(); i++)
				{
					intensities[j][i - 1] = (mzs[j][i] - mzs[j][i - 1]) * 0.5 * (intensities[j][i] + intensities[j][i - 1]);
				}
				intensities[j].resize(intensities[j].size() - 1);
			}
			else
			{
				mzs[j].resize(0);
				intensities[j].resize(0);
			}
		}
	}


	void
	create_gs(vector<fp>& gs, vector<li>& is, vector<ii>& js,
	          const vector< vector<double> >& intensities)
	{
		ii nj = 0;
		for (ii j = 0; j < intensities.size(); j++)
		if (intensities[j].size() > 0) nj++;

		js.resize(nj);
		is.resize(nj + 1);
		is.front() = 0;
		for (ii i = 0, j = 0; j < intensities.size(); j++)
		if (intensities[j].size() > 0)
		{
			js[i] = j;
			is[i + 1] = is[i] + intensities[j].size();
			i++;
		}

		gs.resize(is.back());
		#pragma omp parallel for
		for (ii j = 0; j < nj; j++)
		for (ii i = 0; i < intensities[js[j]].size(); i++)
		{
			gs[is[j] + i] = intensities[js[j]][i];
		}

		cout << "Raw gs=[" << js.size() << "/" << intensities.size() << "]:" << gs.size() << " mem=" << fixed << setprecision(2) << (is.size()*sizeof(li) + js.size()*sizeof(ii) + gs.size()*sizeof(fp)) / (1024.0*1024.0) << "Mb";
		cout << endl;
	}

}
