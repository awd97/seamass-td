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


#include "OptimiserASRL.hpp"
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <omp.h>
using namespace std;


OptimiserASRL::
OptimiserASRL(const std::vector<Basis*>& _bases, const Matrix& _g, ii _accell) :
	g(_g),
    bases(_bases),
    accell(_accell),
    iteration(0),
	cs(bases.size()),
	l1s(bases.size()),
	l2s(bases.size()),
	c0s(accell >= 1 ? bases.size() : 0),
	u0s(accell >= 1 ? bases.size() : 0),
	q0s(accell >= 2 ? bases.size() : 0),
	synthesis_wtime(0),
	error_wtime(0),
	analysis_wtime(0),
	shrinkage_wtime(0)
{
	cout << "Accelerated Sparse Richardson-Lucy accell=" << accell << endl;

	// compute basis function L1 norms 
	l1s.resize(bases.size());
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		if (i == 0)
		{
			Matrix t; t.init(g.n(), 1.0);
			bases[i]->analysis(l1s.front(), t);
		}
		else
		{
			bases[i]->analysis(l1s[i], l1s[bases[i]->get_parent_index()]);
		}
	}
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		cout << l1s[i] << endl;
		if (bases[i]->is_transient()) l1s[i].free();
	}

	// compute basis function L2 norms
	l2s.resize(bases.size());
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		if (i == 0)
		{
			Matrix t; t.init(g.n(), 1.0);
			bases[i]->analysis(l2s.front(), t, true);
		}
		else
		{
			bases[i]->analysis(l2s[i], l2s[bases[i]->get_parent_index()], true);
		}
	}
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		if (bases[i]->is_transient())
		{
			l2s[i].free();
		}
		else
		{
			l2s[i].sqrt();
		}
	}

	vector<Matrix> ts(bases.size() + 1);
	for (ii i = (ii)bases.size() - 1; i >= 0; i--)
	{
		ii parent_i = bases[i]->get_parent_index();

		// if parent Matrix not initialised and not transient then we need to copy in its cs before synthesis
		if (!ts[parent_i + 1] && parent_i >= 0 && !bases[parent_i]->is_transient())
		{
			ts[parent_i + 1].init(cs[parent_i].n(), 1.0);
		}

		// perform the synthesis
		if (!ts[i + 1])
		{
			ts[i + 1].init(cs[i].n, 1.0);
		}
		bases[i]->synthesis(ts[parent_i + 1], ts[i + 1]);
		ts[i + 1].free();
	}
	ts[0].save("t.csv");
	exit(0);

	// initialising C from G
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		if (i == 0)
		{
			bases[i]->analysis(cs[0], g);
		}
		else
		{
			bases[i]->analysis(cs[i], cs[bases[i]->get_parent_index()]);
		}
	}
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		if (bases[i]->is_transient())
		{
			cs[i].free();
		}
		else
		{
			cs[i].elem_div(cs[i], l1s[i]);
		}
	}



	cout << "Initialised C from G" << endl;
}


OptimiserASRL::
~OptimiserASRL()
{
}


void
OptimiserASRL::
threshold(double thresh)
{
    /*for (ii j = 0; j < (ii) bases.size(); j++)
    if (!bases[j]->is_transient())
    #pragma omp parallel for
    for (li i = 0; i < (li) bases[j]->get_nc(); i++)
    {
        if (cs[j][i] < thresh) cs[j][i] = 0.0;
    }*/
    
    //for (ii j = 0; j < (ii) bases.size(); j++)
    //if (!bases[j]->is_transient())
    //#pragma omp parallel for
    //for (li i = 0; i < (li) bases[j]->get_nc(); i++)
    //{
    //    cs[j][i] /= info.volume;
    //}
}


double
OptimiserASRL::
step(ii iteration, double shrinkage)
{
	//////////////////////////////////////////////////////////////////////////////////////
	// SYNTHESIS

	double synthesis_stime = omp_get_wtime();

	vector<Matrix> ts(bases.size() + 1);
	for (ii i = (ii)bases.size() - 1; i >= 0; i--)
	{
		ii parent_i = bases[i]->get_parent_index();

		// if parent Matrix not initialised and not transient then we need to copy in its cs before synthesis
		if (!ts[parent_i + 1] && parent_i >= 0 && !bases[parent_i]->is_transient())
		{
			ts[parent_i + 1].copy(cs[parent_i]);
		}

		// perform the synthesis
		if (!ts[i + 1])
		{
			bases[i]->synthesis(ts[parent_i + 1], cs[i]);
		}
		else
		{
			bases[i]->synthesis(ts[parent_i + 1], ts[i + 1]);
		}
		ts[i + 1].free();
	}

	synthesis_wtime += omp_get_wtime() - synthesis_stime;

	f.copy(ts[0]);

	//////////////////////////////////////////////////////////////////////////////////////
	// POISSON LIKELIHOOD ERROR

	double error_stime = omp_get_wtime();

    info = bases.front()->error(ts[0], ts[0], g);

	error_wtime += omp_get_wtime() - error_stime;

	//////////////////////////////////////////////////////////////////////////////////////
	// ANALYSIS

    double analysis_stime = omp_get_wtime();

	// initialising cs from gs
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		bases[i]->analysis(ts[i + 1], ts[bases[i]->get_parent_index() + 1]);
	}
	ts.front().free();
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		if (bases[i]->is_transient())
		{
			ts[i + 1].free();
		}
	}

	analysis_wtime += omp_get_wtime() - analysis_stime;
	
	//////////////////////////////////////////////////////////////////////////////////////
	// SHRINKAGE

    double shrinkage_stime = omp_get_wtime();

    // shrinkage
    for (ii i = 0; i < (ii) bases.size(); i++)
    {
        if (!bases[i]->is_transient())
        {
			bases[i]->shrink(ts[i + 1], cs[i], l1s[i], l2s[i], shrinkage);
        }
    }

	shrinkage_wtime += omp_get_wtime() - shrinkage_stime;

	double sum = 0.0;
	double sumd = 0.0;
	for (ii i = 0; i < (ii)bases.size(); i++)
	{
		sum += ts[i + 1].sum_sqrd();
		sumd += cs[i].sum_sqrd_diffs(ts[i + 1]);

		cs[i].copy(ts[i + 1]);
	}

	return sqrt(sumd) / sqrt(sum);

/*		if (!bases[j]->is_transient())
#pragma omp parallel for reduction(+:sum,sumd)
			for (li i = 0; i < (li)bases[j]->get_nc(); i++)
			{
		sum += cs[j][i] * cs[j][i];
		sumd += (es[j][i] - cs[j][i])*(es[j][i] - cs[j][i]);

		cs[j][i] = es[j][i];
			}*/

	/*
	//////////////////////////////////////////////////////////////////////////////////////
	// ACCELERATION

    double accel_start = omp_get_wtime();
    // unaccelerated update
    double sum = 0.0;
    double sumd = 0.0;
    fp a = 0.0;
    if (accell == 0)
    {
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        #pragma omp parallel for reduction(+:sum,sumd)
        for (li i = 0; i < (li) bases[j]->get_nc(); i++)
        {
            sum += cs[j][i] * cs[j][i];
            sumd += (es[j][i] - cs[j][i])*(es[j][i] - cs[j][i]);

            cs[j][i] = es[j][i];
        }        
    }
    else // accelerated update    
    {
        // init/update u0s and compute acceleration factor a
        if (iteration == 0)
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            #pragma omp parallel for
            for (li i = 0; i < (li) bases[j]->get_nc(); i++)
            if (cs[j][i] > 0.0)
            {
                u0s[j][i] = es[j][i]/cs[j][i];
            }
            a = 0.0;
        }
        else
        {
            double sum_u0u = 0.0;
            double sum_u0u0 = 0.0;
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            #pragma omp parallel for reduction(+:sum_u0u,sum_u0u0)
            for (li i = 0; i < (li) bases[j]->get_nc(); i++)
            if (cs[j][i] > 0.0)
            {
                double old_u0 = u0s[j][i];
                u0s[j][i] = es[j][i]/cs[j][i];
                
                old_u0 = old_u0 > 0.0 ? c0s[j][i]*log(old_u0) : 0.0;
                sum_u0u += old_u0 * (u0s[j][i] > 0.0 ? es[j][i]*log(u0s[j][i]) : 0.0);
                sum_u0u0 += old_u0 * old_u0;
            }
            a = sqrt(sum_u0u/sum_u0u0);
            a = a > 0.0f ? a : 0.0f;
            a = a < 1.0f ? a : 1.0f;
        }    
        
        if (iteration == 0) // unaccelerated this time, but save es as c0s
        {            
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (li i = 0; i < (li) bases[j]->get_nc(); i++)
                if (cs[j][i] > 0.0)
                {
                    sum += cs[j][i] * cs[j][i];
                    sumd += (es[j][i] - cs[j][i])*(es[j][i] - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = es[j][i];
                    // for next itteration
                    c0s[j][i] = es[j][i];
                }
                vector<fp>().swap(es[j]);
            }
        }
        else if (accell == 1) // linear vector extrapolation
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (li i = 0; i < (li) bases[j]->get_nc(); i++)
                if (cs[j][i] > 0.0)
                {
                    fp c1 = es[j][i] * powf(es[j][i]/c0s[j][i], a);
                    
                    sum += cs[j][i] * cs[j][i];
                    sumd += (c1 - cs[j][i])*(c1 - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = c1;
                    
                    // for next itteration
                    c0s[j][i] = es[j][i];
                }
                vector<fp>().swap(es[j]);
            }
        }
        else if (iteration == 1) // linear vector extrapolation this time, but save the qs
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (li i = 0; i < (li) bases[j]->get_nc(); i++)
                if (cs[j][i] > 0.0)
                {
                    fp q = es[j][i]/c0s[j][i];
                    fp c1 = es[j][i] * powf(es[j][i]/c0s[j][i], a);
                    
                    sum += cs[j][i] * cs[j][i];
                    sumd += (c1 - cs[j][i])*(c1 - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = c1;
                    
                    // for next itteration
                    c0s[j][i] = es[j][i];
                    q0s[j][i] = q;
                }
                vector<fp>().swap(es[j]);
            }
        }
        else // quadratic vector extrapolation
        {
            for (ii j = 0; j < (ii) bases.size(); j++)
            if (!bases[j]->is_transient())
            {
                #pragma omp parallel for reduction(+:sum,sumd)
                for (li i = 0; i < (li) bases[j]->get_nc(); i++)
                if (cs[j][i] > 0.0)
                {
                    fp q = es[j][i]/c0s[j][i];
                    fp c1 = es[j][i] * powf(q, a) * powf(q/q0s[j][i], 0.5f*a*a);
                    
                    sum += cs[j][i] * cs[j][i];
                    sumd += (c1 - cs[j][i])*(c1 - cs[j][i]);
                    
                    // for this itteration
                    cs[j][i] = c1;
                    
                    // for next itteration
                    c0s[j][i] = es[j][i];
                    q0s[j][i] = q;
                }
                vector<fp>().swap(es[j]);
            }
        }
    }
    double accel_duration = omp_get_wtime() - accel_start;
    static double culm_accel_duration = 0.0;
    culm_accel_duration += accel_duration;
    
    //cout << "Iteration   Durations: synthesis=" << setprecision(2) << synthesis_duration << " error=" << error_duration << " analysis=" << analysis_duration << " shrinkage=" << shrinkage_duration << " accel=" << accel_duration << " all=" << synthesis_duration+error_duration+analysis_duration+shrinkage_duration+accel_duration << endl;
    //cout << "Culminative Durations: synthesis=" << culm_synthesis_duration << " error=" << culm_error_duration << " analysis=" << culm_analysis_duration << " shrinkage=" << culm_shrinkage_duration << " accel=" << culm_accel_duration<< " all=" << culm_synthesis_duration+culm_error_duration+culm_analysis_duration+culm_shrinkage_duration +culm_accel_duration << endl;
    
    return sqrt(sumd)/sqrt(sum);*/
}
