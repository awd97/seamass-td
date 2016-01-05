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
OptimiserASRL(const vector<Basis*>& _bases,
              vector<fp>& _gs,
              ii _accell) :
    bases(_bases),
    gs(_gs),
    accell(_accell),
    iteration(0)
{
    // pre-compute weights
    wcs.resize(bases.size());
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        //cout << j << endl;
        wcs[j].resize(bases[j]->get_cm().size());
        if (j == 0)
        {
            vector<fp> ts(gs.size(), 1.0);
            bases[j]->analysis(wcs[0], ts);
        }
        else
        {
            bases[j]->analysis(wcs[j], wcs[bases[j]->get_parent()->get_index()]);
        }
    }
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient()) vector<fp>().swap(wcs[j]);
    }
    
    // pre-compute l2norm
    l2.resize(bases.size());
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        l2[j].resize(bases[j]->get_cm().size());
        if (j == 0)
        {
            vector<fp> ts(gs.size(), 1.0);
            bases[j]->l2norm(l2[0], ts);
        }
        else
        {
            bases[j]->l2norm(l2[j], l2[bases[j]->get_parent()->get_index()]);
        }
    }
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient())
        {
            vector<fp>().swap(l2[j]);
        }
        else
        {   // sqrt
            // not 64 bit compliant atm // vsSqrt(bases[j]->get_cm().size(), &(l2[j][0]), &(l2[j][0]));
            for (li i = 0; i < (li) bases[j]->get_cm().size(); i++) l2[j][i] = sqrt(l2[j][i]);
        }
    }
    
    // starting cs from ng
    cs.resize(bases.size());
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        cs[j].resize(bases[j]->get_cm().size());
        if (j == 0)
        {
            bases[j]->analysis(cs[0], gs);
        }
        else
        {
            bases[j]->analysis(cs[j], cs[bases[j]->get_parent()->get_index()]);
        }
    }
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient()) vector<fp>().swap(cs[j]);
    }
    // any bases that are too small must be deleted to avoid numerical problems
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (!bases[j]->is_transient())
        for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
        if (wcs[j][i] < 0.001) cs[j][i] = 0.0;
    }
    
    // temporaries required for acceleration
    if (accell >= 1)
    {
        c0s.resize(bases.size());
        u0s.resize(bases.size());
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        {
            c0s[j].resize(bases[j]->get_cm().size());
            u0s[j].resize(bases[j]->get_cm().size());
        }
    }
    if (accell >= 2)
    {
        q0s.resize(bases.size());
        for (ii j = 0; j < (ii) bases.size(); j++)
        if (!bases[j]->is_transient())
        {
            q0s[j].resize(bases[j]->get_cm().size());
        }
    }
    
    // how much memory are we using?
    li size = 0;
    for (ii j = 0; j < cs.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < wcs.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < l2.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < c0s.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < u0s.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    for (ii j = 0; j < q0s.size(); j++) if (!bases[j]->is_transient()) size += bases[j]->get_cm().size();
    cout << endl << "Vector Extrapolated Sparse Richardson-Lucy mem=" << fixed << setprecision(2) << (size*sizeof(fp) + sizeof(this))/(1024.0*1024.0) << "Mb" << endl;
}


OptimiserASRL::
~OptimiserASRL()
{
}


void
OptimiserASRL::
threshold(double thresh)
{
    for (ii j = 0; j < (ii) bases.size(); j++)
    if (!bases[j]->is_transient())
    #pragma omp parallel for
    for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
    {
        if (cs[j][i] < thresh) cs[j][i] = 0.0;
    }
    
    double volume = bases.front()->get_volume();
    
    for (ii j = 0; j < (ii) bases.size(); j++)
    if (!bases[j]->is_transient())
    #pragma omp parallel for
    for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
    {
        cs[j][i] /= volume;
    }
}


double
OptimiserASRL::
step(ii iteration, double shrinkage)
{
    double synthesis_start = omp_get_wtime();
    // synthesis except root
    vector< vector<fp> > ts(bases.size());
    for (ii j = (ii) bases.size() - 1; j > 0; j--)
    {
        ii pj = bases[j]->get_parent()->get_index();
        if (ts[pj].size() == 0)
        if (bases[pj]->is_transient())
        {
            ts[pj].resize(bases[pj]->get_cm().size(), 0.0);
        }
        else
        {
            ts[pj].resize(bases[pj]->get_cm().size());
            for (li i = 0; i < (li) bases[pj]->get_cm().size(); i++) ts[pj][i] = cs[pj][i];
        }

        if (ts[j].size() == 0)
        {
            bases[j]->synthesis(ts[pj], cs[j]);
        }
        else
        {
            bases[j]->synthesis(ts[pj], ts[j]);
            
        }
        vector<fp>().swap(ts[j]);
    }
    //synthesis root
    vector<fp> fs(gs.size());
    if (ts.front().size() == 0)
    {
        bases.front()->synthesis(fs, cs.front());
    }
    else
    {
        bases.front()->synthesis(fs, ts.front());
    }
    vector<fp>().swap(ts.front());
    // timing
    double synthesis_duration = omp_get_wtime() - synthesis_start;
    static double culm_synthesis_duration = 0.0;
    culm_synthesis_duration += synthesis_duration;
    
    double error_start = omp_get_wtime();
    // Poisson likelihood error
    bases.front()->error(fs, gs);
    // timing
    double error_duration = omp_get_wtime() - error_start;
    static double culm_error_duration = 0.0;
    culm_error_duration += error_duration;
    
    double analysis_start = omp_get_wtime();
    // analysis root
    vector< vector<fp> > es(bases.size());
    es.front().resize(bases.front()->get_cm().size());
    bases.front()->analysis(es.front(), fs);
    vector<fp>().swap(fs);
    // analysis except root
    for (ii j = 1; j < (ii) bases.size(); j++)
    {
        es[j].resize(bases[j]->get_cm().size());
        bases[j]->analysis(es[j], es[bases[j]->get_parent()->get_index()]);
        // use child count here to delete transient es sooner
    }
    // timing
    double analysis_duration = omp_get_wtime() - analysis_start;
    static double culm_analysis_duration = 0.0;
    culm_analysis_duration += analysis_duration;

    double shrinkage_start = omp_get_wtime();
    // shrinkage
    for (ii j = 0; j < (ii) bases.size(); j++)
    {
        if (bases[j]->is_transient())
        {
            vector<fp>().swap(es[j]);
        }
        else
        {
			bases[j]->shrink(es[j], cs[j], l2[j], wcs[j], shrinkage);
        }
    }
    double shrinkage_duration = omp_get_wtime() - shrinkage_start;
    static double culm_shrinkage_duration = 0.0;
    culm_shrinkage_duration += shrinkage_duration;
    
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
        for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
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
            for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
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
            for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
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
                for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
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
                for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
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
                for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
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
                for (li i = 0; i < (li) bases[j]->get_cm().size(); i++)
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
    
    return sqrt(sumd)/sqrt(sum);
}
