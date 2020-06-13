#include "Model.h"
#include "IO.h"
#include "ModelTypes.h"
#include "FittingProcess.h"
#include "Observations.h"
#include "InferenceParameters.h"

#include <algorithm>
#include <iostream>
#include <ostream>
#include <numeric>
#include <vector>
#include <gsl/gsl_sf_gamma.h>

#include <functional>

namespace EERAModel {
namespace Model {

std::vector<double> BuildPopulationSeed(const std::vector<int>& age_nums,  ModelStructureId structure)
{
    unsigned int end_age = age_nums.size() - 1;
    unsigned int start_age = 1;
    if (ModelStructureId::IRISH == structure)
    {
        start_age = 0;
    }

    std::vector<double> _temp;
    for (unsigned int age = start_age; age < end_age; ++age)
    {
        _temp.push_back(static_cast<double>(age_nums[age]));
    }

    return _temp;
}

std::vector<Compartments> BuildPopulationArray(Random::RNGInterface::Sptr rng,
    const std::vector<int>& age_nums, const seed& seedlist, ModelStructureId structure)
{
    unsigned int distribution_size = 6;
    if (ModelStructureId::IRISH == structure)
    {
        distribution_size = 7;
    }

    std::vector<Compartments> _temp(age_nums.size(), Compartments());
    for (int age = 0; age < age_nums.size(); ++age) 
    {
        //set up the starting population as fully susceptible
         _temp[age].S = age_nums[age];
    }

    if(seedlist.seedmethod != "background")
    {        
        std::vector<int> startdist(distribution_size, 0);
        startdist[static_cast<int>(rng->Flat(0, distribution_size))] = seedlist.nseed;

        if (ModelStructureId::ORIGINAL == structure)
        {
            for (unsigned int age = 1; age < age_nums.size() - 1; ++age)
            {
                _temp[age].S -=  startdist[age - 1];	// take diseased out of S
                _temp[age].I_p +=  startdist[age - 1];	// put diseased in I	
            }
        }
        else if (ModelStructureId::IRISH == structure)
        {
            for (unsigned int age = 0; age < age_nums.size() - 1; ++age)
            {
                _temp[age].S -=  startdist[age];	// take diseased out of S
                _temp[age].I_p +=  startdist[age];	// put diseased in I	
            }
        }
    }

    return _temp;
}

void GenerateDiseasedPopulation(Random::RNGInterface::Sptr rng,
    std::vector<Compartments>& poparray, std::vector<double>& seedarray,
    const double& bkg_lambda, ModelStructureId structure)
{
    unsigned int start_age = 1;
    if (ModelStructureId::IRISH == structure)
    {
        start_age = 0;
    }
	
	int n_susc = 0;
    for (int age = start_age; age < poparray.size() - 1; ++age) 
    {
        seedarray[age - start_age] = static_cast<double>(poparray[age].S);
        n_susc += poparray[age].S;	
    }

    // how many diseased is introduced in each given day before lockdown
    // as a proportion of number of Susceptible available for background infection
    // (not total population, only 20-70 individuals)
    int startdz = rng->Poisson(static_cast<double>(n_susc) * bkg_lambda);

    //unsigned int startdist[seedarray.size()];
    std::vector<unsigned int> startdist(seedarray.size(), 0);
    rng->Multinomial(seedarray.size(), startdz, &seedarray[0], &startdist[0]); //distribute the diseased across the older age categories
    
    for (int age = start_age; age < poparray.size() - 1; ++age) {
    	int nseed = startdist[age - start_age];
    	poparray[age].S -=  nseed;	// take diseased out of S
    	poparray[age].E +=  nseed;	// put diseased in E
    }

}

Status RunModel(std::vector<double> parameter_set, std::vector<::EERAModel::params> fixed_parameters,
				AgeGroupData per_age_data, seed seedlist, int day_shut, std::vector<int> agenums, 
				int n_sim_steps, ModelStructureId structure, Random::RNGInterface::Sptr rng) {


///	std::cout<< "top0..\n";

	//initialize saving of the detection for t=0 in simulations, deaths, hospital deaths
	/**@todo status.deaths is unused anywhere (although it is populated in RunModel) 
	* TP: yes, I had issues in memory to output it and write it in the output_abc-smc_simu file. 
	* an output of this should be created similarly to status.deaths_hospital.
	*/
	Status status = {{0}, {0}, {0}, {}};

	const int n_agegroup = per_age_data.waifw_norm.size();

	// Start without lockdown
	bool inLockdown = false;

	// Assumes that the number of age groups matches the size of the 'agenums' vector
	std::vector<double> seed_pop = BuildPopulationSeed(agenums, structure);

	std::vector<std::vector<double>> parameter_fit(per_age_data.waifw_norm.size(), parameter_set);	
	parameter_fit[0][5] = fixed_parameters[0].juvp_s;

//	std::cout<< "top1..\n";
	std::vector<Compartments> poparray = BuildPopulationArray(rng, agenums, seedlist, structure);
	

//	std::cout<< "top2..\n";

//	std::cout<< "top3..\n";
	
	//run the simulation
	for (int tt{1}; tt < n_sim_steps; ++tt) {
		//initialize return value
		InfectionState infection_state;	
		
//		std::cout<< "top4..\n";
		//identify the lock down
		if(tt > day_shut){ inLockdown = true; }
		
	
		//introduce disease from background infection until lockdown
		if(!inLockdown && seedlist.seedmethod == "background")
		{
			const double bkg_lambda = parameter_set[parameter_set.size()-1];

			//compute the total number of susceptible and the number of susceptible per age class			
			GenerateDiseasedPopulation(rng, poparray, seed_pop, bkg_lambda, structure);
		}

        //compute the forces of infection
        std::vector<double> lambda = GenerateForcesOfInfection(infection_state.hospitalised, parameter_set, fixed_parameters[0].inf_asym, per_age_data,
            poparray, inLockdown);	

        // step each agegroup through infections
        for ( int age{0}; age < n_agegroup; ++age) {	
            InfectionState new_spread;
            if (ModelStructureId::ORIGINAL == structure) {
                new_spread = 
                    GenerateInfectionSpreadOriginal(rng, poparray[age], infection_state.hospitalised,
                        fixed_parameters[age],parameter_fit[age],
                        per_age_data.cfr_byage[age], per_age_data.pf_byage[age],lambda[age]);
            } else {
                new_spread = 
                    GenerateInfectionSpreadIrish(rng, poparray[age], infection_state.hospitalised,
                        fixed_parameters[age],parameter_fit[age],
                        per_age_data.cfr_byage[age], per_age_data.pf_byage[age],lambda[age]);
            }

            infection_state.deaths += new_spread.deaths;
            infection_state.hospital_deaths += new_spread.hospital_deaths;
            infection_state.detected += new_spread.detected;
        }

        status.simulation.push_back(infection_state.detected); 
        status.deaths.push_back(infection_state.deaths); 
        status.hospital_deaths.push_back(infection_state.hospital_deaths);
    }

	//save the population in each epi compt for the last day
	for ( int age{0}; age < n_agegroup; ++age) {
		status.ends.push_back(poparray[age]);
	}

	return status;
}

InfectionState GenerateInfectionSpreadOriginal(Random::RNGInterface::Sptr rng, Compartments& pop,
    const int& n_hospitalised, params fixed_parameters, std::vector<double> parameter_set,
    std::vector<double> cfr_tab, double pf_val, double lambda)
{
    Compartments newpop(pop);

    // Fixed parameters
    const double T_lat= fixed_parameters.T_lat;
    const double T_inf= fixed_parameters.T_inf;
    const double T_rec= fixed_parameters.T_rec;
    const double T_sym= fixed_parameters.T_sym;
    const double T_hos= fixed_parameters.T_hos;
    const double K=fixed_parameters.K;

    double capacity = n_hospitalised / K;
    capacity = std::min(1.0, capacity);

    const double p_h= cfr_tab[0];
    const double p_d= cfr_tab[2];	
    const double p_s= parameter_set[5];
    const double rrd= parameter_set[6];

    // hospitalized  - non-frail
    const int newdeathsH= Flow(rng, pop.H, newpop.H, p_d * (1.0 / T_hos));
    newpop.H -= newdeathsH;
    newpop.D += newdeathsH;

    const int recoverH = Flow(rng, pop.H, newpop.H, (1.0 - p_d) * (1.0 / T_hos));
    newpop.H -= recoverH;
    newpop.R += recoverH;

    // symptomatic - non-frail
    const int hospitalize = Flow(rng, pop.I_s4, newpop.I_s4, p_h  * (1.0 - capacity) * (4.0 / T_sym));
    newpop.I_s4 -= hospitalize;
    newpop.H += hospitalize;

    const int newdeathsI_s = Flow(rng, pop.I_s4, newpop.I_s4, p_h  * p_d * rrd * capacity * ( 4.0 / T_sym));
    newpop.I_s4 -= newdeathsI_s;
    newpop.D += newdeathsI_s;

    const int recoverI_s = Flow(rng, pop.I_s4, newpop.I_s4, ( (1.0 - p_h) + p_h  * (1 - p_d * rrd) * capacity) * ( 4.0 / T_sym));
    newpop.I_s4 -= recoverI_s;
    newpop.R += recoverI_s;

    const int Is_from3_to_4 = Flow(rng, pop.I_s3, newpop.I_s3, (4.0 / T_sym));
    newpop.I_s3 -= Is_from3_to_4;
    newpop.I_s4 += Is_from3_to_4;

    const int Is_from2_to_3 = Flow(rng, pop.I_s2, newpop.I_s2, (4.0 / T_sym));
    newpop.I_s2 -= Is_from2_to_3;
    newpop.I_s3 += Is_from2_to_3;
        
    const int Is_from1_to_2 = Flow(rng, pop.I_s1, newpop.I_s1, (4.0 / T_sym));
    newpop.I_s1 -= Is_from1_to_2;
    newpop.I_s2 += Is_from1_to_2;

    // asymptomatic
    const int recoverI = Flow(rng, pop.I4, newpop.I4, ( 4.0 / T_rec ));
    newpop.I4 -= recoverI;
    newpop.R += recoverI;

    const int I_from3_to_4 = Flow(rng, pop.I3, newpop.I3, ( 4.0 / T_rec ));
    newpop.I3 -= I_from3_to_4;
    newpop.I4 += I_from3_to_4;

    const int I_from2_to_3=Flow(rng, pop.I2, newpop.I2, ( 4.0 / T_rec ));
    newpop.I2 -= I_from2_to_3;
    newpop.I3 += I_from2_to_3;

    const int I_from1_to_2 = Flow(rng, pop.I1, newpop.I1, ( 4.0 / T_rec ));
    newpop.I1 -= I_from1_to_2;
    newpop.I2 += I_from1_to_2;

    // infectious - pre-clinical
    int newasymptomatic = Flow(rng, pop.I_p, newpop.I_p, (1.0 - p_s) * ( 1.0 / T_inf ));
    newpop.I_p -= newasymptomatic;
    newpop.I1 += newasymptomatic;

    const int newsymptomatic = Flow(rng, pop.I_p, newpop.I_p, p_s * ( 1.0 / T_inf ));
    newpop.I_p -= newsymptomatic;
    newpop.I_s1 += newsymptomatic;

    // latent
    const int infectious = Flow(rng, pop.E, newpop.E, ( 1.0 / T_lat ));
    newpop.E -= infectious;
    newpop.I_p += infectious;

    const int infectious_t = Flow(rng, pop.E_t, newpop.E_t, ( 1.0 / T_lat ));
    newpop.E_t -= infectious_t;
    newpop.I_t += infectious_t;

    // susceptible
    const int newinfection = Flow(rng, pop.S, newpop.S, lambda);
    newpop.S -= newinfection;
    newpop.E += newinfection;

    pop = newpop;

    // Define infection state for this spread
    InfectionState infection_state;
    infection_state.deaths              += (newdeathsH + newdeathsI_s);
    infection_state.hospital_deaths     +=  newdeathsH;
    infection_state.detected            += (hospitalize + infectious_t);

    return infection_state;
}

InfectionState GenerateInfectionSpreadIrish(Random::RNGInterface::Sptr rng, Compartments& pop,
    const int& n_hospitalised, params fixed_parameters, std::vector<double> parameter_set,
    std::vector<double> cfr_tab, double pf_val, double lambda)
{
    Compartments newpop(pop);

    // Fixed parameters
    const double T_lat= fixed_parameters.T_lat;
    const double T_inf= fixed_parameters.T_inf;
    const double T_rec= fixed_parameters.T_rec;
    const double T_sym= fixed_parameters.T_sym;
    const double T_hos= fixed_parameters.T_hos;
    const double K=fixed_parameters.K;

    double capacity = n_hospitalised / K;
    capacity = std::min(1.0, capacity);

    const double p_h= cfr_tab[0];
    const double p_d= cfr_tab[2];	
    const double p_s= parameter_set[5];
    const double rrd= parameter_set[6];

    // hospitalized  - non-frail
    const int newdeathsH= Flow(rng, pop.H, newpop.H, p_d * (1.0 / T_hos));
    newpop.H -= newdeathsH;
    newpop.D += newdeathsH;

    const int recoverH = Flow(rng, pop.H, newpop.H, (1.0 - p_d) * (1.0 / T_hos));
    newpop.H -= recoverH;
    newpop.R += recoverH;

    // symptomatic - non-frail
    const int hospitalize = Flow(rng, pop.I_s4, newpop.I_s4, p_h  * (1.0 - capacity) * (4.0 / T_sym));
    newpop.I_s4 -= hospitalize;
    newpop.H += hospitalize;

    const int newdeathsI_s = Flow(rng, pop.I_s4, newpop.I_s4, p_h  * p_d * rrd * capacity * ( 4.0 / T_sym));
    newpop.I_s4 -= newdeathsI_s;
    newpop.D += newdeathsI_s;

    const int recoverI_s = Flow(rng, pop.I_s4, newpop.I_s4, ( (1.0 - p_h) + p_h  * (1 - p_d * rrd) * capacity) * ( 4.0 / T_sym));
    newpop.I_s4 -= recoverI_s;
    newpop.R += recoverI_s;

    const int Is_from3_to_4 = Flow(rng, pop.I_s3, newpop.I_s3, (4.0 / T_sym));
    newpop.I_s3 -= Is_from3_to_4;
    newpop.I_s4 += Is_from3_to_4;

    const int Is_from2_to_3 = Flow(rng, pop.I_s2, newpop.I_s2, (4.0 / T_sym));
    newpop.I_s2 -= Is_from2_to_3;
    newpop.I_s3 += Is_from2_to_3;
        
    const int Is_from1_to_2 = Flow(rng, pop.I_s1, newpop.I_s1, (4.0 / T_sym));
    newpop.I_s1 -= Is_from1_to_2;
    newpop.I_s2 += Is_from1_to_2;

    //pre-clinical
    const int showsymptoms = Flow(rng, pop.I_p, newpop.I_p, (1.0 / T_inf));
    newpop.I_p -= showsymptoms;
    newpop.I_s1 += showsymptoms;

    // asymptomatic
    const int recoverI = Flow(rng, pop.I4, newpop.I4, ( 4.0 / T_rec ));
    newpop.I4 -= recoverI;
    newpop.R += recoverI;

    const int I_from3_to_4 = Flow(rng, pop.I3, newpop.I3, ( 4.0 / T_rec ));
    newpop.I3 -= I_from3_to_4;
    newpop.I4 += I_from3_to_4;

    const int I_from2_to_3=Flow(rng, pop.I2, newpop.I2, ( 4.0 / T_rec ));
    newpop.I2 -= I_from2_to_3;
    newpop.I3 += I_from2_to_3;

    const int I_from1_to_2 = Flow(rng, pop.I1, newpop.I1, ( 4.0 / T_rec ));
    newpop.I1 -= I_from1_to_2;
    newpop.I2 += I_from1_to_2;

    // latent
    const int newsymptomatic = Flow(rng, pop.E, newpop.E, p_s * ( 1.0 / T_lat ));
    newpop.E -= newsymptomatic;
    newpop.I_p += newsymptomatic;

    const int newasymptomatic = Flow(rng, pop.E, newpop.E, (1.0 - p_s) * ( 1.0 / T_lat ));
    newpop.E -= newasymptomatic;
    newpop.I1 += newasymptomatic;

    const int infectious_t = Flow(rng, pop.E_t, newpop.E_t, ( 1.0 / T_lat ));
    newpop.E_t -= infectious_t;
    newpop.I_t += infectious_t;

    // susceptible
    const int newinfection = Flow(rng, pop.S, newpop.S, lambda);
    newpop.S -= newinfection;
    newpop.E += newinfection;

    pop = newpop;

    // Define infection state for this spread
    InfectionState infection_state;
    infection_state.deaths              += (newdeathsH + newdeathsI_s);
    infection_state.hospital_deaths     +=  newdeathsH;
    infection_state.detected            += (hospitalize + infectious_t);

    return infection_state;
}

std::vector<double> GenerateForcesOfInfection(int& inf_hosp, const std::vector<double>& parameter_set, double u_val, 
			const AgeGroupData& age_data, const std::vector<Compartments>& pops, bool shut) 
{

	double p_i = parameter_set[0];	
	double p_hcw = parameter_set[1];	
	double c_hcw = parameter_set[2];	
	double d_val = parameter_set[3];					
	double q_val = parameter_set[4];
	
	//lambda=rep(NA,nrow(pops))
	int n_agegroup = age_data.waifw_norm.size();

	std::vector<double> lambda(n_agegroup);

	double quarantined = 0.0;	// mean proportion reduction in contacts due to quarantine

	std::vector<double> I_mat(n_agegroup); 

	//contact matrix
	std::vector<std::vector<double>> waifw = (shut) ? std::vector<std::vector<double>>(n_agegroup, std::vector<double> (n_agegroup)) : age_data.waifw_norm; 

	//age-dependent transmission rate
	std::vector<std::vector<double>> beta(n_agegroup, std::vector<double> (n_agegroup));

  if(shut)
  {
    for ( int from{0}; from < n_agegroup; ++from) 
    {
      for ( int to{0}; to < n_agegroup; ++to) 
      {
        //contact network during shutdown period, assuming a proportion d will properly do it
        waifw[from][to] = (1.0-d_val) * age_data.waifw_sdist[from][to] + (d_val) * age_data.waifw_home[from][to];
      }
    }
  }
  
  for ( int from{0}; from < n_agegroup; ++from) {
    for ( int to{0}; to < n_agegroup; ++to) 
    {
      if(waifw[from][to]>0)
      {
        quarantined += age_data.waifw_home[from][to]/waifw[from][to];
      } 
			
      else
      {
         quarantined += 0;
      }
      beta[from][to] = waifw[from][to] * p_i ;
	}
  } 
  
  //mean proportion reduction in contacts due to quarantine accounting for compliance
  quarantined = 1 - (quarantined / (double)(n_agegroup * n_agegroup) ) * (1 - q_val);
  
  //compute the pressure from infectious groups, normalizes by group size
  //compute the pressure from hospitalised
  for ( int from{0}; from < n_agegroup; ++from) {
	if(from < (n_agegroup-1)) { //n-agegroup-1 to not account for hcw (as different contact)
	 	// are considered those that are infectious pre-clinical, asymp, tested or symptomatic only
		// asym are infectious pre-clinical and asymptomatic
		// sym are infectious plus those that are tested positive
		int tot_asym = pops[from].I_p + u_val * (pops[from].I1 + pops[from].I2 + pops[from].I3 + pops[from].I4 );
		int tot_sym = pops[from].I_t + pops[from].I_s1 + pops[from].I_s2 + pops[from].I_s3 + pops[from].I_s4;
		
		I_mat[from] = (double)tot_asym + (1-quarantined) * (double)tot_sym;
		//normalisation shouldnt account for dead individuals	
		const int tot_pop = accumulate_compartments(pops[from]);
		I_mat[from] = (double)I_mat[from] / (double)(tot_pop); 
				
	}
	//sum up of number of hospitalised cases (frail and non-frail)
	inf_hosp+= pops[from].H ;
  }

  //sum up infectious pressure from each age group
  for ( int to{0}; to < (n_agegroup); ++to) {
	  for(int from{0}; from < (n_agegroup-1); ++from){ //assume perfect self isolation and regular testing of infected HCW if infected
	  	lambda[to] += ( beta[to][from]*I_mat[from] );
	  }
  }

  const int tot_pop = accumulate_compartments(pops[n_agegroup-1]);
  //int rem_pop = pops[n_agegroup-1][4]+pops[n_agegroup-1][5]+pops[n_agegroup-1][6]+pops[n_agegroup-1][8];
  lambda[n_agegroup-1] = p_hcw * c_hcw * ((double)inf_hosp/(double)(tot_pop) );

  return lambda;
}

/*void Flow(Random::RNGInterface::Sptr rng, unsigned int& pop_from, unsigned int& pop_to, double rate, unsigned int& outs){
    outs = rng->Poisson(rate * (double)pop_from); //symptomatic become hospitalized
    outs = std::min(pop_from, outs);
	pop_to += outs;
	pop_from -= outs;
}*/

int Flow(Random::RNGInterface::Sptr rng, int pops_from_val, int pops_new_from_val, double rate)
{
	int outs = rng->Poisson(rate * static_cast<double>(pops_from_val)); //symptomatic become hospitalized
	outs = std::min(pops_new_from_val, outs);
	return outs;
}

std::vector<int> ComputeAgeNums(int shb_id, int Npop, int N_hcw, const InputObservations& obs) {
	std::vector<int> agenums;
	
	// define age structure of the shb of interest. the -1 is to account for difference in number of
	// rows between two datasets (age_pop does not have a row with column name)
	const auto& agedist = obs.age_pop[shb_id - 1];
	
	for (const auto& var : agedist) {
		// modulate the population of non-hcw now to proportion in each age group recorded in 2011	
		agenums.push_back(round(var * (Npop - N_hcw))); 

	}		
	agenums.push_back(N_hcw);

	return agenums;
}

int GetPopulationOfRegion(const InputObservations& obs, int region_id)
{
	return obs.cases[region_id][0];
}

int ComputeNumberOfHCWInRegion(int regionalPopulation, int totalHCW, const InputObservations& observations)
{
    int scotlandPopulation = 0;
	for (unsigned int region = 0; region < observations.cases.size() - 1; ++region) {
		scotlandPopulation += observations.cases[region][0];
	}
	double regionalProportion = static_cast<double>(regionalPopulation) / scotlandPopulation;
	
    return static_cast<int>(round(totalHCW * regionalProportion)); 
}

std::vector<params> BuildFixedParameters(unsigned int size, params parameters)
{
    return std::vector<params>(size, parameters);
}

int accumulate_compartments(const Compartments& comp)
{
	int _total = 0;
	_total += comp.S + comp.E + comp.E_t + comp.I_p;
	_total += comp.I_t + comp.I1 + comp.I2 + comp.I3;
	_total += comp.I4 + comp.I_s1 + comp.I_s2 + comp.I_s3;
	_total += comp.I_s4 + comp.H + comp.R + comp.D;

	return _total;
}

std::vector<std::vector<int>> compartments_to_vector(const std::vector<Compartments>& cmps_vec)
{
	std::vector<std::vector<int>> _temp;

	for(auto cmps : cmps_vec)
	{
		_temp.push_back({cmps.S, cmps.E, cmps.E_t, cmps.I_p,
						cmps.I_t, cmps.I1, cmps.I2, cmps.I3,
						cmps.I4, cmps.I_s1, cmps.I_s2, cmps.I_s3,
						cmps.I_s4, cmps.H, cmps.R, cmps.D});
	}

	return _temp;
}

} // namespace Model
} // namespace EERAModel
