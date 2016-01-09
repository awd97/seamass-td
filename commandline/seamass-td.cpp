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


#include "../topdown/topdown.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <H5Cpp.h>
using namespace std;


struct spectrum
{
	size_t index;
	unsigned short preset_config;
	double precursor_mz;
	double scan_start_time;
	size_t scan_start_time_index;
};


bool scan_start_time_order(const spectrum& lhs, const spectrum& rhs)
{
	return lhs.scan_start_time <= rhs.scan_start_time;
}


bool seamass_order(const spectrum& lhs, const spectrum& rhs)
{
	if (lhs.preset_config == rhs.preset_config)
	{
		if (lhs.precursor_mz == rhs.precursor_mz)
		{
			return lhs.scan_start_time <= rhs.scan_start_time;
		}
		else
		{
			return lhs.precursor_mz < rhs.precursor_mz;
		}
	}
	else
	{
		return lhs.preset_config < rhs.preset_config;
	}
}


int main(int argc, char *argv[])
{
	cout << argv[0] << endl;


	H5::Exception::dontPrint();
	seamass::topdown_notice();
	if (argc != 9)
	{
		cout << "seaMass-TD (Top Down) Usage" << endl;
		cout << "-----" << endl;
		cout << "seamass-td <in_file> <mz_res> <shrink> <tol> <threads> <out_type>" << endl;
		cout << endl;
		cout << "<in_file>:   Raw input file in seaMass Input format (smj)" << endl;
		cout << "             guidelines: Use pwiz-seamass to convert from mzML or vendor format" << endl;
		cout << "<out_res>:   Output resolution interval between output data points is ~1.00286 / 2^out_res" << endl;
		cout << "             guidelines: around 5" << endl;
		cout << "<max_z>:     Maximum charge state to consider" << endl;
		cout << "             guidelines: as low as possible but must be higher than what you expect is visible" << endl;
		cout << "<peak_fwhm>: The mass resolution (fwhm @ 400 m/z)" << endl;
		cout << "             note: only Orbitrap peak model supported currently" << endl;
		cout << "<shrink>:    Amount of denoising given as: \"L1 shrinkage = 2^shrinkage\"" << endl;
		cout << "             guidelines: around -4" << endl;
		cout << "<tol>:       Convergence tolerance, given as: \"gradient <= 2^tol\"" << endl;
		cout << "             guidelines: around -9" << endl;
		cout << "<threads>:   Number of OpenMP threads to use" << endl;
		cout << "             guidelines: set to amount of CPU cores or 4, whichever is smaller" << endl;
		cout << "<out_type>:  Type of output desired" << endl;
		cout << "             guidelines: 0 = just viz_client input; 1 = also smo; 2 = also debug" << endl;
		return 0;
	}
	string in_file = argv[1];
	int mass_res = atoi(argv[2]);
	int max_z = atoi(argv[3]);
	double max_peak_width = atof(argv[4]);
	int shrink = atoi(argv[5]);
	int tol = atoi(argv[6]);
	int threads = atoi(argv[7]);
	int out_type = atoi(argv[8]);

	int lastdot = in_file.find_last_of("."); 
	string id = (lastdot == string::npos) ? in_file : in_file.substr(0, lastdot); 

	// open H5 file
	H5::H5File file(in_file, H5F_ACC_RDONLY);

	// query necessary metadata
    hsize_t ns;

    H5::DataSet start_time_ds = file.openDataSet("StartTime");
    start_time_ds.getSpace().getSimpleExtentDims(&ns);
    cout << "Querying metadata from " << ns << " spectra..." << endl;
    vector<double> start_times(ns);
    start_time_ds.read(start_times.data(), H5::DataType(H5::PredType::NATIVE_DOUBLE));

    H5::DataSet preset_config_ds = file.openDataSet("PresetConfig");
    preset_config_ds.getSpace().getSimpleExtentDims(&ns);
    vector<unsigned long> config_indices(ns);
    preset_config_ds.read(config_indices.data(), H5::DataType(H5::PredType::NATIVE_ULONG));
    H5::DataSet precursor_mz_ds = file.openDataSet("PrecursorMZ");
    precursor_mz_ds.getSpace().getSimpleExtentDims(&ns);
    vector<double> precursor_mzs(ns);
    precursor_mz_ds.read(precursor_mzs.data(), H5::DataType(H5::PredType::NATIVE_DOUBLE));
    vector<spectrum> spectra(ns);
    for (size_t i = 0; i < ns; i++)
    {
        spectra[i].index = i;
        spectra[i].preset_config = config_indices[i];
        spectra[i].precursor_mz = precursor_mzs[i];
        spectra[i].scan_start_time = start_times[i];
    }
	// determine start_scan_time order of spectra
	sort(spectra.begin(), spectra.end(), scan_start_time_order);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		spectra[i].scan_start_time_index = i;
	}
	// save scan_start_times and sort into seamass processing order:
	// preset_config -> ms_level -> scan_start_time
	vector<double> scan_start_times(ns);
	for (size_t i = 0; i < spectra.size(); i++)
	{
		scan_start_times[i] = spectra[i].scan_start_time;
	}
	sort(spectra.begin(), spectra.end(), seamass_order);
	// load spectra and process
    H5::DataSet mzs_ds = file.openDataSet("SpectrumMZ");
    H5::DataSet intensities_ds = file.openDataSet("SpectrumIntensity");
    // read instrument type
    unsigned long instrument_type = 1;
    H5::Attribute att = intensities_ds.openAttribute("instrumentType");
    att.read(H5::IntType(H5::PredType::NATIVE_USHORT), &instrument_type);
    hsize_t ns1;
    H5::DataSet index_ds = file.openDataSet("SpectrumIndex");
    index_ds.getSpace().getSimpleExtentDims(&ns1);
    vector<double> indices(ns1);
    index_ds.read(indices.data(), H5::DataType(H5::PredType::NATIVE_DOUBLE));
	vector< std::vector<double> > mzs(ns);
	vector< std::vector<double> > intensities(ns);
	int loaded = 0;
    bool precursor_mz_is_constant = true;
	for (size_t i = 0; i < spectra.size(); i++)
	{
        if (loaded > 1 && spectra[i].precursor_mz != spectra[i-1].precursor_mz)
            precursor_mz_is_constant = false;
        
		// load spectra into mzs and intensities vectors if precursor_mz
        hsize_t offset = indices[spectra[i].index];
        hsize_t count = indices[spectra[i].index+1] - offset;
        
        H5::DataSpace memspace(1, &count);

        H5::DataSpace mzs_dsp = mzs_ds.getSpace();
        mzs_dsp.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        mzs[spectra[i].scan_start_time_index].resize(count);
        mzs_ds.read(mzs[spectra[i].scan_start_time_index].data(), H5::DataType(H5::PredType::NATIVE_DOUBLE), memspace, mzs_dsp);

        H5::DataSpace intensities_dsp = intensities_ds.getSpace();
        intensities_dsp.selectHyperslab(H5S_SELECT_SET, &count, &offset);
        intensities[spectra[i].scan_start_time_index].resize(count);
        intensities_ds.read(intensities[spectra[i].scan_start_time_index].data(), H5::DataType(H5::PredType::NATIVE_DOUBLE), memspace, intensities_dsp);
       
        loaded++;
		if ((i == spectra.size()-1 ||
			spectra[i].preset_config != spectra[i+1].preset_config ||
            (spectra[i].precursor_mz != spectra[i+1].precursor_mz && spectra[i].precursor_mz == 0.0)))
		{
            if (precursor_mz_is_constant)
			{
				ostringstream oss; oss << spectra[i].preset_config << "_" << spectra[i].precursor_mz;
				seamass::topdown(id,
					oss.str().c_str(),
					instrument_type,
					scan_start_times, mzs, intensities,
					mass_res, max_z, max_peak_width, shrink, tol,
					threads, out_type);
			}
            loaded = 0;
            precursor_mz_is_constant = true;
        }
	}

	return 0;
}
