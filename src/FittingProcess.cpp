#include "FittingProcess.h"

#include <gsl/gsl_randist.h>
#include <cmath>

namespace EERAModel {
namespace FittingProcess {

void parameter_select_first_step(std::vector<double> &selected_param, std::vector<double> flag1,
	std::vector<double> flag2, gsl_rng * r,int nPar) {
 	//pick randomly parameters' value from priors with gamma shape
    for (int xx = 0; xx < nPar; ++xx) {
        double tmpval = 0.0;
		if(xx == 2 ){
			tmpval = (double)gsl_ran_poisson(r, flag1[xx]);
		} else if( xx == (nPar-1)){
			tmpval = gsl_ran_flat(r, flag1[xx], flag2[xx]);
		} else if(xx > 4 && xx < (nPar-1) ){
			tmpval = gsl_ran_gamma(r, flag1[xx], flag2[xx]);
		} else {
			tmpval = gsl_ran_beta(r, flag1[xx], flag2[xx]);
		}
    	//return a vector with all selection measures
    	selected_param[xx] = tmpval;
	}
}

void parameter_select_nsteps(std::vector<double> &selected_param, int nPar, gsl_rng * r,
	std::vector<particle> particleList, int pick_val, double vlimitKernel[], double vect_min[],
	double vect_Max[]) {

	particle perturbList;

    //define the elements of the chosen particle
    perturbList = particleList[pick_val];

	//perturb all elements of the chosen particle with a kernel following a uniform distribution +/- kernelFactor

    for (int xx = 0; xx < nPar; ++xx) {
        double tmpval = perturbList.parameter_set[xx]+vlimitKernel[xx]*gsl_ran_flat(r, -1, 1);
        while(std::isnan(tmpval) || tmpval<=vect_min[xx] || tmpval>=vect_Max[xx]){
        	tmpval = perturbList.parameter_set[xx]+vlimitKernel[xx]*gsl_ran_flat(r, -1, 1);
        }


    	//return a vector with all selection measures
    	selected_param[xx] = tmpval;
	}
}

void weight_calc(int smc,int pastNpart, std::vector<EERAModel::particle> pastPart,
	EERAModel::particle &currentPart, double vlimitKernel[], int nPar) {

	//calculate the weight of the accepted particle here
	if(smc == 0){
		currentPart.weight = 1.0;
	} else {
		double denom = 0;
		for (int jj = 0; jj < pastNpart; ++jj) {
			int rval_dist=0;
			double m = 0.0;
			//for each parameter of each particle
			for (int ii = 0; ii < nPar; ++ii) {
				//compute the distance between the jjth previous particle and the chosen particle for the iith parameters
				m= std::fabs(currentPart.parameter_set[ii] - pastPart[jj].parameter_set[ii]);
				//add 1 to rval_dist if iith parameter is within the range of the corresponding parameters of the jjth particle
				if(m <= vlimitKernel[ii]){
					++rval_dist;
				}
			}
			//if all values of parameter of the particle are within the range of perturbation of all parameters of a given source particle
			if(rval_dist==nPar){
				//add the weight of this particle to denom
				denom += pastPart[jj].weight;
			}
		}

		if(denom == 0) denom = 1.0 ;
		//the weight of each particle is 1 over the sum of the weights from the particles at s-1 that may generate this particle
		currentPart.weight = 1.0 / denom;
	}
}

} // namespace FittingProcess
} // namespace EERAModel