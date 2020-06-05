#include "IO.h"
#include "IniFile.h"
#include "CSV.h"

#include <fstream>
#include <sstream>

namespace EERAModel {
namespace IO {

EERAModel::ModelInputParameters ReadParametersFromFile(const std::string& filePath, const Utilities::logging_stream::Sptr& log)
{
	ModelInputParameters modelInputParameters;

	CIniFile parameters;

	//Settings
	modelInputParameters.herd_id = atoi(parameters.GetValue("shb_id", "Settings", filePath).c_str());

	//time step for the modelling process
	modelInputParameters.tau = atof(parameters.GetValue("tau", "Settings", filePath).c_str());

	//number of threads used for computation
	modelInputParameters.num_threads = atoi(parameters.GetValue("num_threads", "Settings", filePath).c_str());

	//Seed settings
	modelInputParameters.seedlist.seedmethod = parameters.GetValue("seedmethod", "Seed settings", filePath).c_str();
	if(modelInputParameters.seedlist.seedmethod == "random"){
		modelInputParameters.seedlist.nseed = atoi(parameters.GetValue("nseed", "Seed settings", filePath).c_str());
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		modelInputParameters.seedlist.hrp = atoi(parameters.GetValue("hrp", "Seed settings", filePath).c_str());
	} else {
		(*log) << "Warning!!! Unknown method - using random seed method instead." << endl;
		modelInputParameters.seedlist.nseed = atoi(parameters.GetValue("nseed", "Seed settings", filePath).c_str());
	}
	modelInputParameters.seedlist.use_fixed_seed = static_cast<bool>(
		atoi(parameters.GetValue("use_fixed_seed", "Seed settings", filePath).c_str())
	);
	modelInputParameters.seedlist.seed_value = strtoul(
		parameters.GetValue("seed_value", "Seed settings", filePath).c_str(), NULL, 0
	);
	
	//Fit settings
	modelInputParameters.nsteps = atoi(parameters.GetValue("nsteps", "Fit settings", filePath).c_str());
	modelInputParameters.nParticalLimit = atoi(parameters.GetValue("nParticLimit", "Fit settings", filePath).c_str());
	modelInputParameters.nSim = atoi(parameters.GetValue("nSim", "Fit settings", filePath).c_str());
	modelInputParameters.kernelFactor = atof(parameters.GetValue("kernelFactor", "Fit settings", filePath).c_str());

	//Tolerance settings
	for (int ii = 1; ii <= modelInputParameters.nsteps; ++ii) {
		modelInputParameters.toleranceLimit.push_back(0.0);
	}

	for (int ii = 1; ii <= modelInputParameters.nsteps; ++ii) {
		if(ii==1) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key1", "Tolerance settings", filePath).c_str());
		if(ii==2) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key2", "Tolerance settings", filePath).c_str());
		if(ii==3) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key3", "Tolerance settings", filePath).c_str());
		if(ii==4) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key4", "Tolerance settings", filePath).c_str());
		if(ii==5) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key5", "Tolerance settings", filePath).c_str());
		if(ii==6) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key6", "Tolerance settings", filePath).c_str());
		if(ii==7) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key7", "Tolerance settings", filePath).c_str());
		if(ii==8) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key8", "Tolerance settings", filePath).c_str());
		if(ii==9) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key9", "Tolerance settings", filePath).c_str());
		if(ii==10) modelInputParameters.toleranceLimit[ii-1] = atof(parameters.GetValue("Key10", "Tolerance settings", filePath).c_str());
	}

	//Fixed parameters
	modelInputParameters.paramlist.T_lat = atof(parameters.GetValue("T_lat", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.juvp_s = atof(parameters.GetValue("juvp_s", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_inf = atof(parameters.GetValue("T_inf", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_rec = atof(parameters.GetValue("T_rec", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_sym = atof(parameters.GetValue("T_sym", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.T_hos = atof(parameters.GetValue("T_hos", "Fixed parameters", filePath).c_str());
	modelInputParameters.day_shut = atoi(parameters.GetValue("day_shut", "Fixed parameters", filePath).c_str());
	modelInputParameters.totN_hcw = atoi(parameters.GetValue("totN_hcw", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.K = atoi(parameters.GetValue("K", "Fixed parameters", filePath).c_str());
	modelInputParameters.paramlist.inf_asym = atof(parameters.GetValue("inf_asym", "Fixed parameters", filePath).c_str());

	//priors settings
	modelInputParameters.nPar = atoi(parameters.GetValue("nPar", "Priors settings", filePath).c_str());
	modelInputParameters.prior_pinf_shape1 = atof(parameters.GetValue("prior_pinf_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_pinf_shape2 = atof(parameters.GetValue("prior_pinf_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_phcw_shape1 = atof(parameters.GetValue("prior_phcw_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_phcw_shape2 = atof(parameters.GetValue("prior_phcw_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_chcw_mean = atof(parameters.GetValue("prior_chcw_mean", "Priors settings", filePath).c_str());
	modelInputParameters.prior_d_shape1 = atof(parameters.GetValue("prior_d_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_d_shape2 = atof(parameters.GetValue("prior_d_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_q_shape1 = atof(parameters.GetValue("prior_q_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_q_shape2 = atof(parameters.GetValue("prior_q_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_lambda_shape1 = atof(parameters.GetValue("prior_lambda_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_lambda_shape2 = atof(parameters.GetValue("prior_lambda_shape2", "Priors settings", filePath).c_str());
	
	modelInputParameters.prior_ps_shape1 = atof(parameters.GetValue("prior_ps_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_ps_shape2 = atof(parameters.GetValue("prior_ps_shape2", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrd_shape1 = atof(parameters.GetValue("prior_rrd_shape1", "Priors settings", filePath).c_str());
	modelInputParameters.prior_rrd_shape2 = atof(parameters.GetValue("prior_rrd_shape2", "Priors settings", filePath).c_str());

	// modelInputParameters.run_type = parameters.GetValue("run_type", "Run type", filePath).c_str();
	modelInputParameters.run_type = "Inference";
	// modelInputParameters.run_type = "Prediction";

	return modelInputParameters;
}

EERAModel::PriorParticleParameters ReadPriorParametersFromFile(const std::string& filePath, const Utilities::logging_stream::Sptr& log)
{
	EERAModel::PriorParticleParameters priorParticleParameters;

	std::ifstream infile(filePath.c_str());
<<<<<<< HEAD
	std::string line;
	
	char delimiter = '\n';
	while (std::getline(infile, line, delimiter))
	{
		priorParticleParameters.prior_param_list.push_back(atof(line.c_str()));
=======
	// if (infile.fail()) { std::cout << "Prior Parameters File" << std::endl;}
	std::string line;
	int linepos = 0;
	// int inquotes = false;
	char c;
	int linemax = line.length();
	std::string curstring;
	char delimiter = '\n';
	while (std::getline(infile, line, delimiter))
	{
		while (line[linepos] != 0 && linepos < linemax)
		{
			if (c == delimiter)
			{
				priorParticleParameters.prior_param_list.push_back(atof(curstring.c_str()));
			}
			else
			{
				curstring.push_back(c);
			}
		}
>>>>>>> 825517cf9a21a1eb7c4d5e86292efbcdfc899c20
	}
	
	return priorParticleParameters;
}

EERAModel::InputObservations ReadObservationsFromFiles(const Utilities::logging_stream::Sptr& log)
{
	EERAModel::InputObservations observations;
	(*log) << "[Observations Files]:" << std::endl;

	const std::string scot_data_file = std::string(ROOT_DIR)+"/data/scot_data.csv";
	const std::string scot_deaths_file = std::string(ROOT_DIR)+"/data/scot_deaths.csv";
	const std::string scot_ages_file = std::string(ROOT_DIR)+"/data/scot_age.csv";
	const std::string waifw_norm_file = std::string(ROOT_DIR)+"/data/waifw_norm.csv";
	const std::string waifw_home_file = std::string(ROOT_DIR)+"/data/waifw_home.csv";
	const std::string waifw_sdist_file = std::string(ROOT_DIR)+"/data/waifw_sdist.csv";
	const std::string cfr_byage_file = std::string(ROOT_DIR)+"/data/cfr_byage.csv";
	const std::string scot_frail_file = std::string(ROOT_DIR)+"/data/scot_frail.csv";
	
	//Uploading observed disease data
	//Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
	//rows from 1 are indivudual health board
	//last row is for all of scotland
	
	(*log) << "\t- " << scot_data_file << std::endl;
	EERAModel::CSV::read_csv_int(observations.cases,scot_data_file,',');
	
	//Uploading observed death data
	//Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
	//rows from 1 are indivudual health board
	//last row is for all of scotland
	
	(*log) << "\t- " << scot_deaths_file << std::endl;
	EERAModel::CSV::read_csv_int(observations.deaths,scot_deaths_file,',');
	
	//Uploading population per age group
	//columns are for each individual Health Borad
	//last column is for Scotland
	//rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	(*log) << "\t- " << scot_ages_file << std::endl;
	EERAModel::CSV::read_csv_double(observations.age_pop,scot_ages_file,',');	
	
	//mean number of daily contacts per age group (overall)	
	(*log) << "\t- " << waifw_norm_file << std::endl;
	EERAModel::CSV::read_csv_double(observations.waifw_norm,waifw_norm_file,',');

	//mean number of daily contacts per age group (home only)		
	(*log) << "\t- " << waifw_home_file << std::endl;
	EERAModel::CSV::read_csv_double(observations.waifw_home,waifw_home_file,',');
	
	//mean number of daily contacts per age group (not school, not work)			
	(*log) << "\t- " << waifw_sdist_file << std::endl;
	EERAModel::CSV::read_csv_double(observations.waifw_sdist,waifw_sdist_file,',');	
	
	//Upload cfr by age group
	//col0: p_h: probability of hospitalisation
	//col1: cfr: case fatality ratio
	//col2: p_d: probability of death, given hospitalisation
	//rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	(*log) << "\t- " << cfr_byage_file << std::endl;
	EERAModel::CSV::read_csv_double(observations.cfr_byage,cfr_byage_file,',');	
		
	//Upload frailty probability p_f by age group
	//columns are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	//rows are for each individual Health Borad
	//last row is for Scotland
	(*log) << "\t- " << scot_frail_file << std::endl;
	EERAModel::CSV::read_csv_double(observations.pf_pop,scot_frail_file,',');	

	return observations;
}

void WriteOutputsToFiles(int smc, int herd_id, int Nparticle, int nPar, 
	const std::vector<EERAModel::particle>& particleList, const std::string& outDirPath, const Utilities::logging_stream::Sptr& log)
{
	std::stringstream namefile, namefile_simu, namefile_ends;
	namefile << (outDirPath + "/output_abc-smc_particles_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";
	namefile_simu << (outDirPath + "/output_abc-smc_simu_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";
	namefile_ends << (outDirPath + "/output_abc-smc_ends_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";		

	std::ofstream output_step (namefile.str().c_str());
	std::ofstream output_simu (namefile_simu.str().c_str());
	std::ofstream output_ends (namefile_ends.str().c_str());
	
	//add the column names for each output list of particles
	output_step << "iterID,nsse_cases,nsse_deaths,p_inf,p_hcw,c_hcw,d,q,p_s,rrd,intro,weight\n";

	//add the column names for each output list of chosen simulations
	output_simu << "iterID" << "," << "day" << "," << "inc_case" << "," << "inc_death\n";
	
	//add the column names for each output list of the compartment values of the last day of the chosen simulations
	output_ends << "iterID" << "," << "age_group" << "," << "comparts" << "," << "value\n";		

	// outputs the list of particles (and corresponding predictions) that were accepted at each steps of ABC-smc
	for (int kk = 0; kk < Nparticle; ++kk) {
		output_step << particleList[kk].iter << ", " << particleList[kk].nsse_cases <<  ", " << particleList[kk].nsse_deaths <<  ", " ;

		for (int var = 0; var < nPar; ++var) {
			output_step << particleList[kk].parameter_set[var] << ", ";
		}
		output_step	<< particleList[kk].weight<< '\n';
		
		for (unsigned int var = 0; var < particleList[kk].simu_outs.size(); ++var) {
//				cout  << particleList[kk].iter << ", " << var << ", " <<  particleList[kk].simu_outs[var] << ", " <<  particleList[kk].death_outs[var] << '\n';
			output_simu << particleList[kk].iter << ", " << var << ", " <<  particleList[kk].simu_outs[var] << ", " <<  particleList[kk].death_outs[var] << '\n';
		}
		
		for (unsigned int age = 0; age < particleList[kk].end_comps.size(); ++age) {
			for (unsigned int var = 0; var < particleList[kk].end_comps[0].size(); ++var) {
				output_ends << particleList[kk].iter << ", " << age << ", " << var << ", " <<  particleList[kk].end_comps[age][var] << '\n';
			}
		}
	}
	
	output_step.close();
	output_simu.close();
	output_ends.close();
}

} // namespace IO
} // namespace EERAModel
