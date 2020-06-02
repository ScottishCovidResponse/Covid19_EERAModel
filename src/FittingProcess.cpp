#include "FittingProcess.h"

#include <gsl/gsl_randist.h>
#include <cmath>

namespace EERAModel {
namespace FittingProcess {

void weight_calc(int smc,int pastNpart, std::vector<EERAModel::particle> pastPart,
	EERAModel::particle &currentPart, const std::vector<double>& vlimitKernel, int nPar) {

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