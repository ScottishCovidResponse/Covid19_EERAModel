#include "OriginalModel.h"

namespace EERAModel {
namespace Model {

OriginalModel::OriginalModel(const CommonModelInputParameters& commonParameters,
    HealthBoardData& observations, Random::RNGInterface::Sptr rng, Utilities::logging_stream::Sptr log) 
    : rng_(rng) {
    
    fixedParameters_ = BuildFixedParameters(
        observations.waifw_norm.size(), commonParameters.paramlist
    );
        
    ageGroupData_ = AgeGroupData{
            observations.waifw_norm,
            observations.waifw_home,
            observations.waifw_sdist,
            observations.cfr_byage,
            observations.pf_byage
    };
    
    int regionalPopulation = GetPopulationOfRegion(
        observations
    );

    int healthCareWorkers = ComputeNumberOfHCWInRegion(
        regionalPopulation, commonParameters.totN_hcw, observations
    );
    
    ageNums_ = ComputeAgeNums(
        regionalPopulation, healthCareWorkers, observations
    );

    (*log) << "[Model settings]" << std::endl;
    (*log) << "    Model Structure: Original" << std::endl;
    (*log) << "    Population size: " << regionalPopulation << std::endl;
	(*log) << "    Number of HCW: " << healthCareWorkers << std::endl;
}

std::vector<double> OriginalModel::BuildPopulationSeed(const std::vector<int>& age_nums)
{
    unsigned int end_age = age_nums.size() - 1;
    unsigned int start_age = 1;

    std::vector<double> _temp;
    for (unsigned int age = start_age; age < end_age; ++age)
    {
        _temp.push_back(static_cast<double>(age_nums[age]));
    }

    return _temp;
}

std::vector<Compartments> OriginalModel::BuildPopulationArray(const std::vector<int>& age_nums,
    const seed& seedlist)
{
    unsigned int distribution_size = 6;

    std::vector<Compartments> _temp(age_nums.size(), Compartments());
    for (unsigned int age = 0; age < age_nums.size(); ++age) 
    {
        //set up the starting population as fully susceptible
         _temp[age].S = age_nums[age];
    }

    if(seedlist.seedmethod != "background")
    {        
        std::vector<int> startdist(distribution_size, 0);
        startdist[static_cast<int>(rng_->Flat(0, distribution_size))] = seedlist.nseed;

        for (unsigned int age = 1; age < age_nums.size() - 1; ++age)
        {
            _temp[age].S -=  startdist[age - 1];	// take diseased out of S
            _temp[age].I_p +=  startdist[age - 1];	// put diseased in I	
        }
    }

    return _temp;
}

void OriginalModel::GenerateDiseasedPopulation(std::vector<Compartments>& poparray,
    std::vector<double>& seedarray, const double& bkg_lambda)
{
    size_t start_age = 1;
	int n_susc = 0;
    for (size_t age = start_age; age < poparray.size() - 1; ++age) 
    {
        seedarray[age - start_age] = static_cast<double>(poparray[age].S);
        n_susc += poparray[age].S;	
    }

    // how many diseased is introduced in each given day before lockdown
    // as a proportion of number of Susceptible available for background infection
    // (not total population, only 20-70 individuals)
    int startdz = static_cast<int>(rng_->Poisson(static_cast<double>(n_susc) * bkg_lambda));

    //unsigned int startdist[seedarray.size()];
    std::vector<unsigned int> startdist(seedarray.size(), 0);
    rng_->Multinomial(seedarray.size(), startdz, &seedarray[0], &startdist[0]); //distribute the diseased across the older age categories
    
    for (size_t age = start_age; age < poparray.size() - 1; ++age) {
    	int nseed = startdist[age - start_age];
    	poparray[age].S -=  nseed;	// take diseased out of S
    	poparray[age].E +=  nseed;	// put diseased in E
    }

}

Status OriginalModel::Run(std::vector<double> parameter_set, const seed& seedlist, int day_shut, int n_sim_steps) {
	Status status = {{0}, {0}, {0}, {}, {}};

	const int n_agegroup = ageGroupData_.waifw_norm.size();

	// Start without lockdown
	bool inLockdown = false;

	// Assumes that the number of age groups matches the size of the 'agenums' vector
	std::vector<double> seed_pop = BuildPopulationSeed(ageNums_);

	std::vector<std::vector<double>> parameter_fit(ageGroupData_.waifw_norm.size(), parameter_set);	
	parameter_fit[0][5] = fixedParameters_[0].juvp_s;

	std::vector<Compartments> poparray = BuildPopulationArray(ageNums_, seedlist);

	for (int tt{1}; tt < n_sim_steps; ++tt) {
		//initialize return value
		InfectionState infection_state;	
		
		//identify the lock down
		if(tt > day_shut){ inLockdown = true; }
		
	
		//introduce disease from background infection until lockdown
		if(!inLockdown && seedlist.seedmethod == "background")
		{
			//compute the total number of susceptible and the number of susceptible per age class			
			GenerateDiseasedPopulation(poparray, seed_pop, parameter_set[ModelParameters::LAMBDA]);
		}

        //compute the forces of infection
        std::vector<double> lambda = GenerateForcesOfInfection(infection_state.hospitalised, parameter_set, fixedParameters_[0].inf_asym, ageGroupData_,
            poparray, inLockdown);	

        // step each agegroup through infections
        for ( int age{0}; age < n_agegroup; ++age) {	
            InfectionState new_spread = 
                GenerateInfectionSpread(poparray[age], infection_state.hospitalised,
                    fixedParameters_[age],parameter_fit[age],
                    ageGroupData_.cfr_byage[age], lambda[age]);
            std::cout << age << "\t" << tt << "\t" << n_sim_steps << std::endl;
            infection_state.deaths += new_spread.deaths;
            infection_state.hospital_deaths += new_spread.hospital_deaths;
            infection_state.detected += new_spread.detected;
        }

        status.pop_array.push_back(poparray);
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

void OriginalModel::SetFixedParameters(const params& paramlist)
{
    fixedParameters_ = BuildFixedParameters(fixedParameters_.size(), paramlist);
}

InfectionState OriginalModel::GenerateInfectionSpread(Compartments& pop,
    const int& n_hospitalised, params fixed_parameters, std::vector<double> parameter_set,
    std::vector<double> cfr_tab, double lambda)
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

    const double p_h = cfr_tab[0];
    const double p_d = cfr_tab[2];	
    const double p_s = parameter_set[ModelParameters::PS];
    const double rrd = parameter_set[ModelParameters::RRD];

    // hospitalized  - non-frail
	const int outpatient = Flow(rng_, pop.H, newpop.H, (1.0 / T_hos));
    const int newdeathsH = static_cast<int>(rng_->Binomial(p_d,outpatient));
    const int recoverH = outpatient - newdeathsH;	
	
    newpop.H -= outpatient;
    newpop.D += newdeathsH;
    newpop.R += recoverH;
	

    // symptomatic - non-frail
	
	const int outClinical = Flow(rng_, pop.I_s4, newpop.I_s4, (4.0 / T_sym));
	const int severe = static_cast<int>(rng_->Binomial(p_h,outClinical));
	const int mild = outClinical - severe;
	
	const int hospitalize = static_cast<int>(rng_->Binomial((1.0 - capacity),severe));
	const int nothospitalize = severe - hospitalize;
	
	const int newdeathsI_s = static_cast<int>(rng_->Binomial(p_d * rrd,nothospitalize));
	const int recoverI_s = nothospitalize - newdeathsI_s;
	
	newpop.I_s4 -= outClinical;
	newpop.H += hospitalize;
	newpop.D += newdeathsI_s;
	newpop.R += mild+recoverI_s;	
	
	

    const int Is_from3_to_4 = Flow(rng_, pop.I_s3, newpop.I_s3, (4.0 / T_sym));
    newpop.I_s3 -= Is_from3_to_4;
    newpop.I_s4 += Is_from3_to_4;

    const int Is_from2_to_3 = Flow(rng_, pop.I_s2, newpop.I_s2, (4.0 / T_sym));
    newpop.I_s2 -= Is_from2_to_3;
    newpop.I_s3 += Is_from2_to_3;
        
    const int Is_from1_to_2 = Flow(rng_, pop.I_s1, newpop.I_s1, (4.0 / T_sym));
    newpop.I_s1 -= Is_from1_to_2;
    newpop.I_s2 += Is_from1_to_2;

    // asymptomatic
    const int recoverI = Flow(rng_, pop.I4, newpop.I4, ( 4.0 / T_rec ));
    newpop.I4 -= recoverI;
    newpop.R += recoverI;

    const int I_from3_to_4 = Flow(rng_, pop.I3, newpop.I3, ( 4.0 / T_rec ));
    newpop.I3 -= I_from3_to_4;
    newpop.I4 += I_from3_to_4;

    const int I_from2_to_3=Flow(rng_, pop.I2, newpop.I2, ( 4.0 / T_rec ));
    newpop.I2 -= I_from2_to_3;
    newpop.I3 += I_from2_to_3;

    const int I_from1_to_2 = Flow(rng_, pop.I1, newpop.I1, ( 4.0 / T_rec ));
    newpop.I1 -= I_from1_to_2;
    newpop.I2 += I_from1_to_2;

    // infectious - pre-clinical
	
	const int outpreclin = Flow(rng_, pop.I_p, newpop.I_p,  ( 1.0 / T_inf ));	
	const int newsymptomatic = static_cast<int>(rng_->Binomial(p_s,outpreclin));
	int newasymptomatic = 	outpreclin - newsymptomatic;	
			
    newpop.I_p -= outpreclin;
    newpop.I_s1 += newsymptomatic;
    newpop.I1 += newasymptomatic;

    // latent
    const int infectious = Flow(rng_, pop.E, newpop.E, ( 1.0 / T_lat ));
    newpop.E -= infectious;
    newpop.I_p += infectious;

    const int infectious_t = Flow(rng_, pop.E_t, newpop.E_t, ( 1.0 / T_lat ));
    newpop.E_t -= infectious_t;
    newpop.I_t += infectious_t;

    // susceptible
    const int newinfection = Flow(rng_, pop.S, newpop.S, lambda);
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

std::vector<double> OriginalModel::GenerateForcesOfInfection(int& inf_hosp, const std::vector<double>& parameter_set, double u_val, 
			const AgeGroupData& age_data, const std::vector<Compartments>& pops, bool shut) 
{
	double p_i = parameter_set[ModelParameters::PINF];	
	double p_hcw = parameter_set[ModelParameters::PHCW];	
	double c_hcw = parameter_set[ModelParameters::CHCW];	
	double d_val = parameter_set[ModelParameters::D];					
	double q_val = parameter_set[ModelParameters::Q];
	
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
    quarantined = 1 - (quarantined / static_cast<double>(n_agegroup * n_agegroup) ) * (1 - q_val);

    //compute the pressure from infectious groups, normalizes by group size
    //compute the pressure from hospitalised
    for ( int from{0}; from < n_agegroup; ++from) {
        if(from < (n_agegroup-1)) { //n-agegroup-1 to not account for hcw (as different contact)
            // are considered those that are infectious pre-clinical, asymp, tested or symptomatic only
            // asym are infectious pre-clinical and asymptomatic
            // sym are infectious plus those that are tested positive
            int tot_asym = pops[from].I_p + static_cast<int>(u_val * (pops[from].I1 + pops[from].I2 + pops[from].I3 + pops[from].I4 ));
            int tot_sym = pops[from].I_t + pops[from].I_s1 + pops[from].I_s2 + pops[from].I_s3 + pops[from].I_s4;
            
            I_mat[from] = static_cast<double>(tot_asym) + (1-quarantined) * static_cast<double>(tot_sym);
            //normalisation shouldnt account for dead individuals	
            const int tot_pop = accumulate_compartments(pops[from]);
            I_mat[from] = static_cast<double>(I_mat[from]) / static_cast<double>(tot_pop);            
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
    lambda[n_agegroup-1] = p_hcw * c_hcw * (static_cast<double>(inf_hosp)/static_cast<double>(tot_pop) );

    return lambda;
}

} // namespace Model
} // namespace EERAModel
