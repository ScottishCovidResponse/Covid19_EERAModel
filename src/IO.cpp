#include "IO.h"
#include "ModelCommon.h"
#include "Utilities.h"

#include <valarray>
#include <fstream>
#include <sstream>
#include <iostream>

namespace EERAModel {
namespace IO {

ModelInputParameters ReadParametersFromFile(const DataSourcing::DataFiles& data_files, const Utilities::logging_stream::Sptr& log)
{
	ModelInputParameters modelInputParameters;

	CIniFile parameters;

	//Settings
	//modelInputParameters.herd_id = atoi(parameters.GetValue_modif(&SettingName, &SettingCategory, &filePath).c_str());
<<<<<<< HEAD
	modelInputParameters.herd_id = ReadNumberFromFile<int>("shb_id", "Settings", filePath);

	//time step for the modelling process
	// modelInputParameters.tau = atof(parameters.GetValue("tau", "Settings", filePath).c_str());
	modelInputParameters.tau = ReadNumberFromFile<double>("tau", "Settings", filePath);

	//number of threads used for computation
	modelInputParameters.num_threads = ReadNumberFromFile<int>("num_threads", "Settings", filePath);
=======
	modelInputParameters.herd_id = ReadNumberFromFile<int>("shb_id", "Settings", data_files.parameters);

	//time step for the modelling process
	// modelInputParameters.tau = atof(parameters.GetValue("tau", "Settings", filePath).c_str());
	modelInputParameters.tau = ReadNumberFromFile<double>("tau", "Settings", data_files.parameters);

	//number of threads used for computation
	modelInputParameters.num_threads = ReadNumberFromFile<int>("num_threads", "Settings", data_files.parameters);
>>>>>>> kzscisoft/SCRC-440-abstraction-layer

    // Model structure
	std::string SettingName = "model";
	std::string SettingCategory = "Settings";
<<<<<<< HEAD
    std::string model_structure(CIniFile::GetValue(SettingName, SettingCategory, filePath));
=======
    std::string model_structure(CIniFile::GetValue(SettingName, SettingCategory, data_files.parameters));
>>>>>>> kzscisoft/SCRC-440-abstraction-layer
    if ("irish" == model_structure) {
        modelInputParameters.model_structure = ModelStructureId::IRISH;
    } else if("original" == model_structure) {
        modelInputParameters.model_structure = ModelStructureId::ORIGINAL;
    } else {
    	modelInputParameters.model_structure = ModelStructureId::TEMP;
    }

	//Seed settings
	SettingName = "seedmethod";
	SettingCategory = "Seed settings";
<<<<<<< HEAD
	modelInputParameters.seedlist.seedmethod = CIniFile::GetValue(SettingName, SettingCategory, filePath);
	if(modelInputParameters.seedlist.seedmethod == "random"){
		modelInputParameters.seedlist.nseed = ReadNumberFromFile<int>("nseed", "Seed settings", filePath);
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		modelInputParameters.seedlist.hrp = ReadNumberFromFile<int>("hrp", "Seed settings", filePath);
	} else {
		(*log) << "Warning!!! Unknown method - using random seed method instead." << endl;
		modelInputParameters.seedlist.nseed = ReadNumberFromFile<int>("nseed", "Seed settings", filePath);
	}
	modelInputParameters.seedlist.use_fixed_seed = static_cast<bool>(
		ReadNumberFromFile<int>("use_fixed_seed", "Seed settings", filePath)
	);
	modelInputParameters.seedlist.seed_value = ReadNumberFromFile<int>(
		"seed_value", "Seed settings", filePath
	);
	
	//Fit settings
	modelInputParameters.nsteps = ReadNumberFromFile<int>("nsteps", "Fit settings", filePath);
	modelInputParameters.nParticalLimit = ReadNumberFromFile<int>("nParticLimit", "Fit settings", filePath);
	modelInputParameters.nSim = ReadNumberFromFile<int>("nSim", "Fit settings", filePath);
	modelInputParameters.kernelFactor = ReadNumberFromFile<double>("kernelFactor", "Fit settings", filePath);
=======
	modelInputParameters.seedlist.seedmethod = CIniFile::GetValue(SettingName, SettingCategory, data_files.parameters);
	if(modelInputParameters.seedlist.seedmethod == "random"){
		modelInputParameters.seedlist.nseed = ReadNumberFromFile<int>("nseed", "Seed settings", data_files.parameters);
	} else if(modelInputParameters.seedlist.seedmethod == "background"){
		modelInputParameters.seedlist.hrp = ReadNumberFromFile<int>("hrp", "Seed settings", data_files.parameters);
	} else {
		(*log) << "Warning!!! Unknown method - using random seed method instead." << endl;
		modelInputParameters.seedlist.nseed = ReadNumberFromFile<int>("nseed", "Seed settings", data_files.parameters);
	}
	modelInputParameters.seedlist.use_fixed_seed = static_cast<bool>(
		ReadNumberFromFile<int>("use_fixed_seed", "Seed settings", data_files.parameters)
	);
	modelInputParameters.seedlist.seed_value = ReadNumberFromFile<int>(
		"seed_value", "Seed settings", data_files.parameters
	);
	
	//Fit settings
	modelInputParameters.nsteps = ReadNumberFromFile<int>("nsteps", "Fit settings", data_files.parameters);
	modelInputParameters.nParticalLimit = ReadNumberFromFile<int>("nParticLimit", "Fit settings", data_files.parameters);
	modelInputParameters.nSim = ReadNumberFromFile<int>("nSim", "Fit settings", data_files.parameters);
	modelInputParameters.kernelFactor = ReadNumberFromFile<double>("kernelFactor", "Fit settings", data_files.parameters);
>>>>>>> kzscisoft/SCRC-440-abstraction-layer

	//Tolerance settings
	for (int ii = 1; ii <= modelInputParameters.nsteps; ++ii) {
		modelInputParameters.toleranceLimit.push_back(0.0);
	}

	for (int ii = 0; ii < modelInputParameters.nsteps; ++ii) {
		std::stringstream KeyName;
		KeyName << "Key" << (ii + 1);
<<<<<<< HEAD
		modelInputParameters.toleranceLimit[ii] = ReadNumberFromFile<double>(KeyName.str(), "Tolerance settings", filePath);
	}

	//Fixed parameters
	modelInputParameters.paramlist.T_lat = ReadNumberFromFile<double>("T_lat", "Fixed parameters", filePath);
	modelInputParameters.paramlist.juvp_s = ReadNumberFromFile<double>("juvp_s", "Fixed parameters", filePath);
	modelInputParameters.paramlist.T_inf = ReadNumberFromFile<double>("T_inf", "Fixed parameters", filePath);
	modelInputParameters.paramlist.T_rec = ReadNumberFromFile<double>("T_rec", "Fixed parameters", filePath);
	modelInputParameters.paramlist.T_sym = ReadNumberFromFile<double>("T_sym", "Fixed parameters", filePath);
	modelInputParameters.paramlist.T_hos = ReadNumberFromFile<double>("T_hos", "Fixed parameters", filePath);
	modelInputParameters.day_shut = ReadNumberFromFile<int>("day_shut", "Fixed parameters", filePath);
	modelInputParameters.totN_hcw = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", filePath);
	modelInputParameters.paramlist.K = ReadNumberFromFile<int>("K", "Fixed parameters", filePath);
	modelInputParameters.paramlist.inf_asym = ReadNumberFromFile<double>("inf_asym", "Fixed parameters", filePath);

	//priors settings
	modelInputParameters.nPar = ReadNumberFromFile<int>("nPar", "Priors settings", filePath);
	modelInputParameters.prior_pinf_shape1 = ReadNumberFromFile<double>("prior_pinf_shape1", "Priors settings", filePath);
	modelInputParameters.prior_pinf_shape2 = ReadNumberFromFile<double>("prior_pinf_shape2", "Priors settings", filePath);
	modelInputParameters.prior_phcw_shape1 = ReadNumberFromFile<double>("prior_phcw_shape1", "Priors settings", filePath);
	modelInputParameters.prior_phcw_shape2 = ReadNumberFromFile<double>("prior_phcw_shape2", "Priors settings", filePath);
	modelInputParameters.prior_chcw_mean = ReadNumberFromFile<double>("prior_chcw_mean", "Priors settings", filePath);
	modelInputParameters.prior_d_shape1 = ReadNumberFromFile<double>("prior_d_shape1", "Priors settings", filePath);
	modelInputParameters.prior_d_shape2 = ReadNumberFromFile<double>("prior_d_shape2", "Priors settings", filePath);
	modelInputParameters.prior_q_shape1 = ReadNumberFromFile<double>("prior_q_shape1", "Priors settings", filePath);
	modelInputParameters.prior_q_shape2 = ReadNumberFromFile<double>("prior_q_shape2", "Priors settings", filePath);
	modelInputParameters.prior_lambda_shape1 = ReadNumberFromFile<double>("prior_lambda_shape1", "Priors settings", filePath);
	modelInputParameters.prior_lambda_shape2 = ReadNumberFromFile<double>("prior_lambda_shape2", "Priors settings", filePath);
	
	modelInputParameters.prior_ps_shape1 = ReadNumberFromFile<double>("prior_ps_shape1", "Priors settings", filePath);
	modelInputParameters.prior_ps_shape2 = ReadNumberFromFile<double>("prior_ps_shape2", "Priors settings", filePath);
	modelInputParameters.prior_rrd_shape1 = ReadNumberFromFile<double>("prior_rrd_shape1", "Priors settings", filePath);
	modelInputParameters.prior_rrd_shape2 = ReadNumberFromFile<double>("prior_rrd_shape2", "Priors settings", filePath);

	// modelInputParameters.run_type = parameters.GetValue("run_type", "Run type", filePath).c_str();
	modelInputParameters.run_type = ModelModeId::INFERENCE;
	// modelInputParameters.run_type = ModelModeId::PREDICTION;

	modelInputParameters.posterior_parameter_select = ReadNumberFromFile<int>("posterior_parameter_select", "Posterior Parameters Select", filePath);
=======
		modelInputParameters.toleranceLimit[ii] = ReadNumberFromFile<double>(KeyName.str(), "Tolerance settings", data_files.parameters);
	}

	//Fixed parameters
	modelInputParameters.paramlist.T_lat = ReadNumberFromFile<double>("T_lat", "Fixed parameters", data_files.parameters);
	modelInputParameters.paramlist.juvp_s = ReadNumberFromFile<double>("juvp_s", "Fixed parameters", data_files.parameters);
	modelInputParameters.paramlist.T_inf = ReadNumberFromFile<double>("T_inf", "Fixed parameters", data_files.parameters);
	modelInputParameters.paramlist.T_rec = ReadNumberFromFile<double>("T_rec", "Fixed parameters", data_files.parameters);
	modelInputParameters.paramlist.T_sym = ReadNumberFromFile<double>("T_sym", "Fixed parameters", data_files.parameters);
	modelInputParameters.paramlist.T_hos = ReadNumberFromFile<double>("T_hos", "Fixed parameters", data_files.parameters);
	modelInputParameters.day_shut = ReadNumberFromFile<int>("day_shut", "Fixed parameters", data_files.parameters);
	modelInputParameters.totN_hcw = ReadNumberFromFile<int>("totN_hcw", "Fixed parameters", data_files.parameters);
	modelInputParameters.paramlist.K = ReadNumberFromFile<int>("K", "Fixed parameters", data_files.parameters);
	modelInputParameters.paramlist.inf_asym = ReadNumberFromFile<double>("inf_asym", "Fixed parameters", data_files.parameters);

	//priors settings
	modelInputParameters.nPar = ReadNumberFromFile<int>("nPar", "Priors settings", data_files.parameters);
	modelInputParameters.prior_pinf_shape1 = ReadNumberFromFile<double>("prior_pinf_shape1", "Priors settings", data_files.parameters);
	modelInputParameters.prior_pinf_shape2 = ReadNumberFromFile<double>("prior_pinf_shape2", "Priors settings", data_files.parameters);
	modelInputParameters.prior_phcw_shape1 = ReadNumberFromFile<double>("prior_phcw_shape1", "Priors settings", data_files.parameters);
	modelInputParameters.prior_phcw_shape2 = ReadNumberFromFile<double>("prior_phcw_shape2", "Priors settings", data_files.parameters);
	modelInputParameters.prior_chcw_mean = ReadNumberFromFile<double>("prior_chcw_mean", "Priors settings", data_files.parameters);
	modelInputParameters.prior_d_shape1 = ReadNumberFromFile<double>("prior_d_shape1", "Priors settings", data_files.parameters);
	modelInputParameters.prior_d_shape2 = ReadNumberFromFile<double>("prior_d_shape2", "Priors settings", data_files.parameters);
	modelInputParameters.prior_q_shape1 = ReadNumberFromFile<double>("prior_q_shape1", "Priors settings", data_files.parameters);
	modelInputParameters.prior_q_shape2 = ReadNumberFromFile<double>("prior_q_shape2", "Priors settings", data_files.parameters);
	modelInputParameters.prior_lambda_shape1 = ReadNumberFromFile<double>("prior_lambda_shape1", "Priors settings", data_files.parameters);
	modelInputParameters.prior_lambda_shape2 = ReadNumberFromFile<double>("prior_lambda_shape2", "Priors settings", data_files.parameters);
	
	modelInputParameters.prior_ps_shape1 = ReadNumberFromFile<double>("prior_ps_shape1", "Priors settings", data_files.parameters);
	modelInputParameters.prior_ps_shape2 = ReadNumberFromFile<double>("prior_ps_shape2", "Priors settings", data_files.parameters);
	modelInputParameters.prior_rrd_shape1 = ReadNumberFromFile<double>("prior_rrd_shape1", "Priors settings", data_files.parameters);
	modelInputParameters.prior_rrd_shape2 = ReadNumberFromFile<double>("prior_rrd_shape2", "Priors settings", data_files.parameters);

	// modelInputParameters.run_type = parameters.GetValue("run_type", "Run type", data_files.parameters).c_str();
	modelInputParameters.run_type = ModelModeId::INFERENCE;
	// modelInputParameters.run_type = ModelModeId::PREDICTION;

	modelInputParameters.posterior_parameter_select = ReadNumberFromFile<int>("posterior_parameter_select", "Posterior Parameters Select", data_files.parameters);
>>>>>>> kzscisoft/SCRC-440-abstraction-layer

	return modelInputParameters;
}

<<<<<<< HEAD
std::vector<double> ReadPosteriorParametersFromFile(const std::string& filePath, int set_selection)
=======
std::vector<double> ReadPosteriorParametersFromFile(const DataSourcing::DataFiles& data_files, int set_selection)
>>>>>>> kzscisoft/SCRC-440-abstraction-layer
{
	// Temporary matrix to hold data from input file
	std::vector<std::vector<double>> lines;
	char delimiter = ',';
	
<<<<<<< HEAD
	lines = Utilities::read_csv<double>(filePath, delimiter);
=======
	lines = Utilities::read_csv<double>(data_files.parameters, delimiter);
>>>>>>> kzscisoft/SCRC-440-abstraction-layer

	// Select line from input file and store result in another temporary vector
	if (set_selection >= static_cast<int>(lines.size())){
		std::stringstream SetSelectError;
		SetSelectError << "Parameter set selection out of bounds! Please select between 0-" << (lines.size() - 1) << "..." << std::endl;
		throw std::overflow_error(SetSelectError.str());
	}
	std::vector<double> line_select = lines[set_selection];

	/** Extract posterior parameters from selected line
	 * the slicing is relevant to the format of the output file
	 * *output_abc-smc_particles_step* coming from the inference framework
	 */
	auto first = line_select.cbegin() + 3;
	auto last = line_select.cend() - 1;
	int nPar = static_cast<int>(line_select.size()) - 4;
	if ((last - line_select.begin()) - (first - line_select.begin()) < nPar) {
		std::stringstream PosteriorFileFormatError;
		PosteriorFileFormatError << "Please check formatting of posterios parameter input file, 8 parameter values are needed..." << std::endl;
		throw std::runtime_error(PosteriorFileFormatError.str());
	}

	std::vector<double> parameter_sets(first, last);

	return parameter_sets;
}

<<<<<<< HEAD
InputObservations ReadObservationsFromFiles(const Utilities::logging_stream::Sptr& log)
=======
EERAModel::InputObservations ReadObservationsFromFiles(const DataSourcing::DataFiles& data_files, const Utilities::logging_stream::Sptr& log)
>>>>>>> kzscisoft/SCRC-440-abstraction-layer
{
	InputObservations observations;
	(*log) << "[Observations Files]:" << std::endl;
	
	//Uploading observed disease data
	//Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
	//rows from 1 are indivudual health board
	//last row is for all of scotland
	
	(*log) << "\t- " << data_files.data << std::endl;
	observations.cases = Utilities::read_csv<int>(data_files.data,',');
	
	//Uploading observed death data
	//Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
	//rows from 1 are indivudual health board
	//last row is for all of scotland
	
	(*log) << "\t- " << data_files.deaths << std::endl;
	observations.deaths = Utilities::read_csv<int>(data_files.deaths,',');
	
	//Uploading population per age group
	//columns are for each individual Health Borad
	//last column is for Scotland
	//rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	(*log) << "\t- " << data_files.ages << std::endl;
	observations.age_pop = Utilities::read_csv<double>(data_files.ages,',');	
	
	//mean number of daily contacts per age group (overall)	
	(*log) << "\t- " << data_files.waifw.norm << std::endl;
	observations.waifw_norm = Utilities::read_csv<double>(data_files.waifw.norm,',');

	//mean number of daily contacts per age group (home only)		
	(*log) << "\t- " << data_files.waifw.home << std::endl;
	observations.waifw_home = Utilities::read_csv<double>(data_files.waifw.home,',');
	
	//mean number of daily contacts per age group (not school, not work)			
	(*log) << "\t- " << data_files.waifw.sdist << std::endl;
	observations.waifw_sdist = Utilities::read_csv<double>(data_files.waifw.sdist,',');	
	
	//Upload cfr by age group
	//col0: p_h: probability of hospitalisation
	//col1: cfr: case fatality ratio
	//col2: p_d: probability of death, given hospitalisation
	//rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	(*log) << "\t- " << data_files.cfr_byage << std::endl;
	observations.cfr_byage = Utilities::read_csv<double>(data_files.cfr_byage,',');	
		
	//Upload frailty probability p_f by age group
	//columns are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	//rows are for each individual Health Borad
	//last row is for Scotland
	(*log) << "\t- " << data_files.frail << std::endl;
	observations.pf_pop = Utilities::read_csv<double>(data_files.frail,',');	

	return observations;
}

void WriteOutputsToFiles(int smc, int herd_id, int Nparticle, int nPar, 
	const std::vector<particle>& particleList, const std::string& outDirPath, const Utilities::logging_stream::Sptr& log)
{
	std::stringstream namefile, namefile_simu, namefile_ends;
	namefile << (outDirPath + "/output_abc-smc_particles_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";
	namefile_simu << (outDirPath + "/output_abc-smc_simu_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";
	namefile_ends << (outDirPath + "/output_abc-smc_ends_step") << smc << "_shb"<< herd_id << "_" << log->getLoggerTime() << ".txt";		

	std::ofstream output_step (namefile.str().c_str());
	std::ofstream output_simu (namefile_simu.str().c_str());
	std::ofstream output_ends (namefile_ends.str().c_str());
	
	//add the column names for each output list of particles
	output_step << "iterID,nsse_cases,nsse_deaths,p_inf,p_hcw,c_hcw,d,q,p_s,rrd,intro,weight" << std::endl;

	//add the column names for each output list of chosen simulations
	output_simu << "iterID" << "," << "day" << "," << "inc_case" << "," << "inc_death_hospital" << "," << "inc_death" << std::endl;
	
	//add the column names for each output list of the compartment values of the last day of the chosen simulations
	output_ends << "iterID" << "," << "age_group" << "," << "comparts" << "," << "value" << std::endl;		

	// outputs the list of particles (and corresponding predictions) that were accepted at each steps of ABC-smc
	for (int kk = 0; kk < Nparticle; ++kk) {
		output_step << particleList[kk].iter << ", " << particleList[kk].nsse_cases <<  ", " << particleList[kk].nsse_deaths <<  ", " ;

		for (int var = 0; var < nPar; ++var) {
			output_step << particleList[kk].parameter_set[var] << ", ";
		}
		output_step	<< particleList[kk].weight<< '\n';
		
		for (unsigned int var = 0; var < particleList[kk].simu_outs.size(); ++var) {
			output_simu << particleList[kk].iter << ", " << var << ", " <<  particleList[kk].simu_outs[var] << ", " \
						<<  particleList[kk].hospital_death_outs[var] << ", " << particleList[kk].death_outs[var] << '\n';
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

void WritePredictionsToFiles(Status status, std::vector<std::vector<int>>& end_comps, const std::string& outDirPath, const Utilities::logging_stream::Sptr& log)
{
	std::stringstream namefile_simu, namefile_ends, namefile_full;
	namefile_simu << (outDirPath + "/output_prediction_simu") << "_" << log->getLoggerTime() << ".txt";
	namefile_ends << (outDirPath + "/output_prediction_ends") << "_" << log->getLoggerTime() << ".txt";
	namefile_full << (outDirPath + "/output_prediction_full") << "_" << log->getLoggerTime() << ".txt";

	std::ofstream output_simu (namefile_simu.str().c_str());
	std::ofstream output_ends (namefile_ends.str().c_str());
	std::ofstream output_full (namefile_full.str().c_str());

	output_simu << "day" << "," << "inc_case" << "," << "inc_death_hospital" << "," << "inc_death" << std::endl;

	output_ends << "age_group" << "," << "comparts" << "," << "value" << std::endl;

	for (unsigned int var = 0; var < status.simulation.size(); ++var) {
		output_simu << var << ", " << status.simulation[var] << ", "\
					<< status.hospital_deaths[var] << ", " << status.deaths[var] << '\n';
	}
	
	for (unsigned int age = 0; age < end_comps.size(); ++age) {
		for (unsigned int var = 0; var < end_comps[0].size(); ++var) {
			output_ends << age << ", " << var << ", " << end_comps[age][var] << '\n';
		}
	}

	/** The first two lines are just markers for age group and compartments
	 * similar to the first and second columns of output_ends. These can be 
	 * removed if it'll make parsing easier.
	 */
	for (unsigned int age = 0; age < end_comps.size(); ++age) {
		for (unsigned int comp = 0; comp < end_comps[age].size(); ++comp) {
			output_full << age;
			if (comp < end_comps[0].size() - 1) { output_full << ", "; }
		}
		output_full << "\t";
	}
	output_full << std::endl;

	for (unsigned int age = 0; age < end_comps.size(); ++age) {
		for (unsigned int comp = 0; comp < end_comps[age].size(); ++comp) {
			output_full << comp;
			if (comp < end_comps[0].size() - 1) { output_full << ", "; }
		}
		output_full << "\t";
	}

	output_full << std::endl;
	

	for (unsigned int var = 0; var < status.pop_array.size(); ++var) {
		std::vector<std::vector<int>> pop_array_compartment_to_vector = Model::compartments_to_vector(status.pop_array[var]);
		for (auto & age: pop_array_compartment_to_vector) {
			for (unsigned int comp = 0; comp < age.size(); ++comp) {
				output_full << age[comp];
				if (comp < age.size() - 1) { output_full << ", "; }
			}
			output_full << "\t";
		}
		if (var < status.pop_array.size() - 1) { output_full << '\n'; }
	}

	output_simu.close();
	output_ends.close();
	output_full.close();
}

} // namespace IO
} // namespace EERAModel
