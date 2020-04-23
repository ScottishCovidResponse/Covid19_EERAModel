/* Created on: 17 04 2020
 * Authors: Thibaud Porphyre
 *
 *
 *version 0.3.2
 *
 * version report: 
 *		  - fix the transmission flow between Is,H, D and R relative to %frail
 * 
 *  
 * fitting procedure: ABS-SMC
 * model to fit: spread, SEIIsHRD
 * number of herd: single
 * model type: stochastic, age-structured population, tau-leap
 *
 * time-step = 1 day
 *
 * selection measures: normalise sum squared error 
 *
 * fitted parameters: p_i, p_hcw, c_hcw, q and d
 *
 * main.cpp
 *
 *
 */


#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <omp.h>
#include <random>
#include <array>
#include "IniFile.h"

using namespace std;



//create structure for particles generated during inference
struct particle {
	double nsse_cases;
	double nsse_deaths;
	vector<double> parameter_set;
	int iter;
	double weight;
	vector<int> simu_outs;
	vector<int> death_outs;
	vector<vector<int>> end_comps;
};


//create structure for fixed parameters.
struct params {
	double T_lat;
	double p_s;
	double juvp_s;
	double T_inf;
	double T_rec;
	double T_sym;
	double T_hos;
};

//create structure for seed settings.
struct seed {
	string seedmethod;
	int nseed;
	int hrp;
	int day_intro;
	double lambda;
};



//declare vectors of inputs as global variables
vector<int > obsHosp;
vector<int > obsDeaths;

//declare functions
	//tools
void read_csv_int(vector<vector<int> > &data, const string &inputfile, char delimiter);
void read_csv_double(vector<vector<double> > &data, const string &inputfile, char delimiter);

void read_parameters(int& herd_id, double& tau, int&num_threads,int& nsteps, int& nParticLimit, int& nSim, double& kernelFactor, 
		vector<double>& toleranceLimit, params& paramlist, seed& seedlist, int& day_shut, int& totN_hcw,
		int& nPar, double& prior_pinf_shape1,double& prior_pinf_shape2,double& prior_phcw_shape1,double& prior_phcw_shape2,
		double& prior_chcw_mean,double& prior_d_shape1,double& prior_d_shape2,double& prior_q_shape1, double& prior_q_shape2,
		double& prior_rrdh_shape1,double& prior_rrdh_shape2,double& prior_rrdc_shape1,double& prior_rrdc_shape2,
		double& prior_rrh_shape1,double& prior_rrh_shape2,double& prior_lambda_shape1, double& prior_lambda_shape2);
		
double mean_calc_double( vector<double> v);
void cumsum_calc(vector<double> v,vector<double>& r_val);
double max_calc(double a,double b);
double sum_calc(vector<double> v);
double colsum_calc(vector<vector<double>> m, int colid);

double mean_calc_int( vector<int> v);
void cumsum_calc_int(vector<int> v,vector<int>& r_val);
int max_calc_int(int a,int b);
int sum_calc_int(vector<int> v);
int colsum_calc_int(vector<vector<int>> m, int colid);

void compute_incidence(vector<int> v,vector<int>& r_val);
void correct_incidence(vector<int>& v,vector<int> cumv);

void select_obs(int& Npop, int& t_index, int& duration, int& day_intro, int& day_shut, 
	vector<int>& obsHosp_tmp, vector<int>& obsDeaths_tmp, 
	vector<vector<int> > data_tmp, vector<vector<int> > death_tmp, int herd_id, 
	int time_back);

	//fitting process functions
void model_select(int smc,particle &outvec, vector<params> fixed_parameters,vector<vector<double>> cfr_byage, 
	vector<double> pf_byage,
	vector<vector<double>> waifw_norm, vector<vector<double>> waifw_sdist, vector<vector<double>> waifw_home, 
	vector <int> agenums, double tau, int duration, seed seedlist, int day_shut, int Npop, gsl_rng * r);
void parameter_select_first_step(vector<double> &selected_param,vector<double> flag1, vector<double> flag2, gsl_rng * r,int nPar);
void parameter_select_nsteps(vector<double> &selected_param, int nPar, gsl_rng * r, vector<particle> particleList, int pick_val,
		double vlimitKernel[], double vect_min[], double vect_Max[]);
void weight_calc(int smc,int pastNpart, vector<particle > pastPart, particle & currentPart, double vlimitKernel[],int nPar);


	//distance computation functions
double sse_calc(int Npop, vector<double> simHosp);
double sse_calc_int(vector<int> simval, vector<int> obsval);
double nsse_calc_int(vector<int> simHosp);
double sst_calc(int Npop, vector<double> simHosp);
double sst_calc_int(int Npop, vector<int> simHosp);

	//model functions
void my_model(vector<double> parameter_set, vector<params> fixed_parameters, vector<vector<double>> cfr_byage, 
				vector<double> pf_byage,
				vector<vector<double>> waifw_norm, vector<vector<double>> waifw_sdist, vector<vector<double>> waifw_home,
				int duration, seed seedlist, int day_shut, int Npop, vector<int> agenums, 
				double tau, gsl_rng * r, vector<int> &sim_status,vector<vector<int>> &ends,
				vector<int> &death_status,vector<int> &deathH_status);

void infspread(gsl_rng * r, vector<int>& pop, int& deaths, int& deathsH, int& detected, 
				params fixed_parameters, vector<double> parameter_set, vector<double> cfr_tab, double pf_val, double lambda);

void Lambda(vector<double> &lambda, vector<double> parameter_set,vector<vector<double>> waifw_norm,
				vector<vector<double>> waifw_sdist,vector<vector<double>> waifw_home, 
				vector<vector<int>> pops, int shut);



int main(int argc, char **argv) {

	/*---------------------------------------
	 * Model parameters and fitting settings
	 *---------------------------------------*/
	//get the start time
	clock_t startTime = clock();
	clock_t time_taken1=0;


	//variable declaration
	int smc,kk;
	double time_taken;
	double valMax, valmin;
	int herd_id,nsteps,nParticLimit,nSim,num_threads;
	double tau,kernelFactor;
	vector<double> toleranceLimit;
	
	params paramlist;
	seed seedlist;
	double prior_pinf_shape1,prior_pinf_shape2, prior_phcw_shape1,prior_phcw_shape2;
	double prior_chcw_mean,prior_d_shape1,prior_d_shape2,prior_q_shape1, prior_q_shape2;
	
	double prior_rrdh_shape1,prior_rrdh_shape2,prior_rrdc_shape1,prior_rrdc_shape2;
	double prior_rrh_shape1,prior_rrh_shape2,prior_lambda_shape1,prior_lambda_shape2;
	int nPar=0, day_shut=0, totN_hcw=0;
	
	//read settings information for fitting procedure
	read_parameters(herd_id, tau, num_threads, nsteps, nParticLimit,nSim,kernelFactor,toleranceLimit,
			paramlist,seedlist,day_shut, totN_hcw,
			nPar,prior_pinf_shape1,prior_pinf_shape2, prior_phcw_shape1,prior_phcw_shape2,
			prior_chcw_mean,prior_d_shape1,prior_d_shape2,prior_q_shape1, prior_q_shape2,
			prior_rrdh_shape1,prior_rrdh_shape2,prior_rrdc_shape1,prior_rrdc_shape2,
			prior_rrh_shape1,prior_rrh_shape2,prior_lambda_shape1,prior_lambda_shape2);

    cout << "[Settings]:\n";
	cout<< "number of parameters tested: "<<nPar<<endl;
    cout<< "seeding method: "<<seedlist.seedmethod<<endl;
	if(seedlist.seedmethod == "random"){
		cout<< "number of seed: " << seedlist.nseed <<endl;
	} else if(seedlist.seedmethod == "background"){
		cout<< "duration of the high risk period: " << seedlist.hrp <<endl;
	}
	
	/*---------------------------------------
	 * Observed data
	 *---------------------------------------*/

	//Uploading observed disease data
	//Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
	//rows from 1 are indivudual health board
	//last row is for all of scotland
	vector<vector<int> > data_tmp;
	read_csv_int(data_tmp,"./DATA/scot_data.csv",',');
	
	//Uploading observed death data
	//Note: first vector is the vector of time. value of -1 indicate number of pigs in the herd
	//rows from 1 are indivudual health board
	//last row is for all of scotland
	vector<vector<int> > death_tmp;
	read_csv_int(death_tmp,"./DATA/scot_deaths.csv",',');
	
	//Uploading population per age group
	//columns are for each individual Health Borad
	//last column is for Scotland
	//rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	vector<vector<double> > age_pop;
	read_csv_double(age_pop,"./DATA/scot_age.csv",',');	
	//Uploading contact matrices
	vector<vector<double> > waifw_norm, waifw_home, waifw_sdist;

	//mean number of daily contacts per age group (overall)	
	read_csv_double(waifw_norm,"./DATA/waifw_norm.csv",',');

	//mean number of daily contacts per age group (home only)		
	read_csv_double(waifw_home,"./DATA/waifw_home.csv",',');
	
	//mean number of daily contacts per age group (not school, not work)			
	read_csv_double(waifw_sdist,"./DATA/waifw_sdist.csv",',');	
	
	
	//Upload cfr by age group
	//col0: p_h: probability of hospitalisation
	//col1: cfr: case fatality ratio
	//col2: p_d: probability of death, given hospitalisation
	//rows are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	vector<vector<double> > cfr_byage;
	read_csv_double(cfr_byage,"./DATA/cfr_byage.csv",',');	
		
	//Upload frailty probability p_f by age group
	//columns are for each age group: [0] Under20,[1] 20-29,[2] 30-39,[3] 40-49,[4] 50-59,[5] 60-69,[6] Over70,[7] HCW
	//rows are for each individual Health Borad
	//last row is for Scotland
	vector<vector<double> > pf_pop;
	read_csv_double(pf_pop,"./DATA/scot_frail.csv",',');	

	//keep information for the health board if interest
	vector<double> pf_byage = pf_pop[herd_id-1];//define frailty structure of the shb of interest.

	//create vector of fixed parameters		
	vector<params> fixed_parameters(waifw_norm.size());	
	for (unsigned int var = 0; var < fixed_parameters.size(); ++var) {
		fixed_parameters[var] = paramlist;
	}
	fixed_parameters[0].p_s = fixed_parameters[0].juvp_s;

	//Separate case information for each herd_id
	vector<int> obsHosp_tmp, obsDeaths_tmp;
	int Npop=0;
	int t_index=-1;
	seedlist.day_intro=0;
	
	int duration = 0;
	
	if(seedlist.seedmethod == "background"){
		select_obs(	Npop, t_index, duration, seedlist.day_intro, day_shut, 
			obsHosp_tmp, obsDeaths_tmp, data_tmp, death_tmp, herd_id, 
			seedlist.hrp);
	} else{
		select_obs(	Npop, t_index, duration, seedlist.day_intro, day_shut, 
			obsHosp_tmp, obsDeaths_tmp, data_tmp, death_tmp, herd_id, 
			paramlist.T_inf + paramlist.T_sym);
		if(seedlist.seedmethod!= "random"){
			cout << "Warning!! Unknown seeding method - applying _random_ seed method\n";
		}
	}

	obsHosp = obsHosp_tmp;
	obsDeaths = obsDeaths_tmp;
	obsHosp_tmp.clear();	
	obsDeaths_tmp.clear();	
	
	//define age structure and number of hcw of the population at risk
	//compute the number of hcw in the shb
	int N_scot = 0;
	for (unsigned int nn = 0; nn < data_tmp.size()-1; ++nn) {
		N_scot += data_tmp[nn][0];  //compute number of scots in scotland
	}
	double prop_scot = (double)Npop / (double)N_scot;  //proportion of Scots in each shb
	int N_hcw = round(totN_hcw * prop_scot); // modulate total number of hcw in Scotland to population in shb

	//adjust population structure to the current population size
	vector<double> agedist = age_pop[herd_id-1];//define age structure of the shb of interest. the -1 is to account for difference in number of rows between two datasets (age_pop does not have a row with column name)
	vector<int> agenums;
	for (unsigned int var = 0; var < agedist.size(); ++var) {
		// modulate the population of non-hcw now to proportion in each age group recorded in 2011	
		agenums.push_back(round(agedist[var] * (Npop - N_hcw))); 

	}		
	//push back the number of hcw in the shb of interest
	agenums.push_back(N_hcw);
	
    cout << "[Health Board settings]:\n";
	cout << "    SHB id: " << herd_id <<'\n';
	cout << "    Population size: " << Npop << '\n';
	cout << "    Number of HCW: " << N_hcw << '\n';
	cout << "    Simulation period: " << duration << "days\n";
	cout << "    time step: " << tau << "days\n";


	/*---------------------------------------
	 * Model settings
	 *---------------------------------------*/
	//set output file name
	stringstream namefile, namefile_simu, namefile_ends;

	//initialise the gsl random number generator with a seed depending of time of the run
	gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937); //“Mersenne Twister” random number generator
	gsl_rng_set(r, time(NULL));

	//initialise the random number generator for importance sampling
    //default_random_engine gen;
	mt19937 gen (time(NULL));

	//declare vectors for priors
	vector<double> flag1, flag2;
	
	flag1 = {prior_pinf_shape1, prior_phcw_shape1, prior_chcw_mean, prior_d_shape1, prior_q_shape1, 
		prior_rrdh_shape1,prior_rrdc_shape1,prior_rrh_shape1,prior_lambda_shape1};	
	flag2 = { prior_pinf_shape2, prior_phcw_shape2, prior_chcw_mean, prior_d_shape2, prior_q_shape2,
		prior_rrdh_shape2,prior_rrdc_shape2,prior_rrh_shape2,prior_lambda_shape2};			

	//declare vectors of outputs/inputs for ABC process
	vector<particle > particleList,particleList1;
	vector<double> weight_val;

	//declare intermediate vectors for ABC smc process
	int pick_val;
	double vlimitKernel[nPar];
	double vect_Max[nPar];
	double vect_min[nPar];

	//initialise the number of accepted particles
	int Nparticle = 0;
	int nsim_val = 0;

	/*--------------------------------------
	 * abc-smc loop
	 *-------------------------------------*/
	 cout << "[Simulations]:\n";
	for (smc = 0; smc < nsteps; ++smc) {//todo:

		//set the output files for each steps
		namefile.str(std::string());
		namefile_simu.str(std::string());
		namefile_ends.str(std::string());
		namefile << "./OUTPUTS/output_abc-smc_particles_step" << smc << "_shb"<< herd_id << ".txt";
		namefile_simu << "./OUTPUTS/output_abc-smc_simu_step" << smc << "_shb"<< herd_id << ".txt";
		namefile_ends << "./OUTPUTS/output_abc-smc_ends_step" << smc << "_shb"<< herd_id << ".txt";		
		ofstream output_step (namefile.str().c_str());
		ofstream output_simu (namefile_simu.str().c_str());
		ofstream output_ends (namefile_ends.str().c_str());
		//add the column names for each output list of particles
		output_step << "iterID,nsse_cases,nsse_deaths,p_inf,p_hcw,c_hcw,d,q,rrdh,rrdc,rrh,intro,weight\n";

		//add the column names for each output list of chosen simulations
		output_simu << "iterID" << "," << "day" << "," << "inc_case" << "," << "inc_death\n";
		
		//add the column names for each output list of the compartment values of the last day of the chosen simulations
		output_ends << "iterID" << "," << "age_group" << "," << "comparts" << "," << "value\n";		

		//the abort statement for keeping the number of particles less than 1000
		bool aborting = false;

		//initialise the counter for the number of simulation prior maximum particle is reached
		int nsim_count = 0;

		//initialise the number of accepted particles
		int counter = 0;

		//initialise the weight and particle lists
		if (smc>0) {
			//update the vectors
			particleList = particleList1;
			particleList1.clear();

			//compute the kernel window
			for (int i = 0; i < nPar; ++i) {
				particle valMax1 = *max_element(particleList.begin(),particleList.end(), [&i](particle a , particle b) { return a.parameter_set[i] < b.parameter_set[i]; } ); //find the element of a vector that is the biggest and return its value
				particle valmin1 = *min_element(particleList.begin(),particleList.end(), [&i](particle a , particle b) { return a.parameter_set[i] < b.parameter_set[i]; } ); //find the element of a vector  that is the smallest and return its value
				valMax = vect_Max[i] = valMax1.parameter_set[i];
				valmin = vect_min[i] = valmin1.parameter_set[i];
				vlimitKernel[i] = kernelFactor*fabs((valMax)-(valmin));
			}
		}


		//create the discrete distribution of the weights for the "importance sampling" process
		for (unsigned int i = 0; i < particleList.size(); ++i) {
			weight_val.push_back(particleList[i].weight);
		}
		//weight_val = weight_val1;
		discrete_distribution<int> weight_distr(weight_val.begin(), weight_val.end());
		weight_val.clear();

	/*---------------------------------------
	 * simulate the infection data set
	 *---------------------------------------*/
		//omp_set_num_threads(num_threads); // Use maximum num_threads threads for all consecutive parallel regions
		//#pragma omp parallel for default(shared), private(pick_val), firstprivate(weight_distr)
		for (kk = 0; kk < nSim; ++kk) {//todo
			//abort statement if number of accepted particles reached nParticLimit particles
			//#pragma omp flush (aborting)
			if (!aborting) {

				//Update progress
				if (counter >= (nParticLimit)) {
					aborting = true;
					//#pragma omp flush (aborting)
				}

				//declare and initialise output variables
				particle outs_vec;
				outs_vec.iter = kk;
				outs_vec.nsse_cases = 0.0;
				outs_vec.nsse_deaths = 0.0;
//				outs_vec.sum_sq = 1.e06;
				for (int i = 0; i < nPar; ++i) {
					outs_vec.parameter_set.push_back(0.0);
				}
				//pick the values of each particles
				if (smc==0) {
					//pick randomly and uniformly parameters' value from priors
					parameter_select_first_step(outs_vec.parameter_set,flag1, flag2, r,nPar);
				} else {
					//sample 1 particle from the previously accepted particles and given their weight (also named "importance sampling")
				    pick_val = weight_distr(gen);
				    parameter_select_nsteps(outs_vec.parameter_set, nPar, r, particleList, pick_val, vlimitKernel,vect_min,vect_Max);
				}
				//run the model and compute the different measures for each potential parameters value
				model_select(smc,outs_vec, fixed_parameters, cfr_byage,pf_byage,waifw_norm,waifw_sdist,waifw_home,
							agenums,tau, duration, seedlist, day_shut, Npop, r);
				//count the number of simulations that were used to reach the maximum number of accepted particles
				//#pragma omp critical
					{
					if (counter < nParticLimit) ++nsim_count;
					}
				//if the particle agrees with the different criteria defined for each ABC-smc step
				//if(counter < nParticLimit && outs_vec.nsse_cases <= toleranceLimit[smc]){
				//if(counter < nParticLimit && outs_vec.nsse_deaths <= toleranceLimit[smc]){
				if(counter < nParticLimit && outs_vec.nsse_cases <= toleranceLimit[smc] && outs_vec.nsse_deaths <= toleranceLimit[smc]*1.5){
					
					//#pragma omp critical
					{
						weight_calc(smc,Nparticle, particleList, outs_vec, vlimitKernel,nPar);
						particleList1.push_back(outs_vec);
						++counter;		//count the number of accepted particles
						if(counter % 10 == 0) cout << "|" << flush;
						//cout << counter << " " ;
					}
				}
			}
		}


		Nparticle = counter;
		nsim_val = nsim_count;

		//time taken per step
		if(smc == 0){
			time_taken = double( clock() - startTime ) / (double)CLOCKS_PER_SEC;
			time_taken1 = clock();
		} else {
			time_taken = double( clock() - time_taken1 ) / (double)CLOCKS_PER_SEC;
			time_taken1 = clock();
		}

		/*---------------------------------------
		 * Outputs
		 *---------------------------------------*/
		//output on screen of the number of accepted particles, the number of simulations and the computation time at each step
		cout << "\nStep:" << smc << ", <number of accepted particles> " << Nparticle << "; <number of simulations> " << nsim_val << "; <computation time> " <<  time_taken << " seconds.\n";

		//break the ABC-smc at the step where no particles were accepted
		if(Nparticle==0){
			break;
		}

		// outputs the list of particles (and corresponding predictions) that were accepted at each steps of ABC-smc
		for (kk = 0; kk < Nparticle; ++kk) {
			output_step << particleList1[kk].iter << ", " << particleList1[kk].nsse_cases <<  ", " << particleList1[kk].nsse_deaths <<  ", " ;

			for (int var = 0; var < nPar; ++var) {
				output_step << particleList1[kk].parameter_set[var] << ", ";
			}
			output_step	<< particleList1[kk].weight<< '\n';
			
			for (unsigned int var = 0; var < particleList1[kk].simu_outs.size(); ++var) {
//				cout  << particleList1[kk].iter << ", " << var << ", " <<  particleList1[kk].simu_outs[var] << ", " <<  particleList1[kk].death_outs[var] << '\n';
				output_simu << particleList1[kk].iter << ", " << var << ", " <<  particleList1[kk].simu_outs[var] << ", " <<  particleList1[kk].death_outs[var] << '\n';
			}
			
			for (unsigned int age = 0; age < particleList1[kk].end_comps.size(); ++age) {
				for (unsigned int var = 0; var < particleList1[kk].end_comps[0].size(); ++var) {
					output_ends << particleList1[kk].iter << ", " << age << ", " << var << ", " <<  particleList1[kk].end_comps[age][var] << '\n';
				}
			}
		}
		output_step.close();
		output_simu.close();
		output_ends.close();
	}

	//output on screen the overall computation time
	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
}
///##########################################################################################################//todo:


void read_csv_int(vector<vector<int> > &data, const string &inputfile, char delimiter)
{
	ifstream infile(inputfile.c_str());
	if (infile.fail())  { cout << "Input file not found" << endl; return; }
	string line;
	vector<int> record;
	
	while (getline(infile, line,'\n') )
	{
		int linepos=0;
		int inquotes=false;
		char c;
		int linemax=line.length();
		string curstring;
		record.clear();
	while(line[linepos]!=0 && linepos < linemax)
		{
			c = line[linepos];
			
			if (!inquotes && curstring.length()==0 && c=='"'){ inquotes=true;}	//beginquotechar
			else if (inquotes && c=='"')
			{																	//quotechar
				if ( (linepos+1 <linemax) && (line[linepos+1]=='"') )
				{
					curstring.push_back(c);										//encountered 2 double quotes in a row (resolves to 1 double quote)
					linepos++;
				}
				else {inquotes=false;}											//endquotechar
			}
			else if (!inquotes && c==delimiter)
			{																//end of field
				record.push_back(atoi(curstring.c_str()) );
				curstring="";
			}
/*			else if (!inquotes && (c=='\r' || c=='\n') )
			{
				record.push_back( atoi(curstring.c_str()) );
				break;
			}
*/			else
			{
				curstring.push_back(c);
			}
			linepos++;
		}
		record.push_back( atoi(curstring.c_str()) );
		data.push_back(record);
	}
}

void read_csv_double(vector<vector<double> > &data, const string &inputfile, char delimiter)
{
	ifstream infile(inputfile.c_str());
	if (infile.fail())  { cout << "Input file not found" << endl; return; }
	string line;
	vector<double> record;
	
	while (getline(infile, line,'\n') )
	{
		int linepos=0;
		int inquotes=false;
		char c;
		int linemax=line.length();
		string curstring;
		record.clear();
	while(line[linepos]!=0 && linepos < linemax)
		{
			c = line[linepos];
			
			if (!inquotes && curstring.length()==0 && c=='"'){ inquotes=true;}	//beginquotechar
			else if (inquotes && c=='"')
			{																	//quotechar
				if ( (linepos+1 <linemax) && (line[linepos+1]=='"') )
				{
					curstring.push_back(c);										//encountered 2 double quotes in a row (resolves to 1 double quote)
					linepos++;
				}
				else {inquotes=false;}											//endquotechar
			}
			else if (!inquotes && c==delimiter)
			{																//end of field
				record.push_back(atof(curstring.c_str()) );
				curstring="";
			}
/*			else if (!inquotes && (c=='\r' || c=='\n') )
			{
				record.push_back( atoi(curstring.c_str()) );
				break;
			}
*/			else
			{
				curstring.push_back(c);
			}
			linepos++;
		}
		record.push_back( atof(curstring.c_str()) );
		data.push_back(record);
	}
}


double sse_calc(int Npop, vector<double> simHosp){
	double sum_sq=0.0;
	//verify that the 2 vectors have the same size
	if(simHosp.size() != obsHosp.size()){
		cout<< "Warning!!! Simulated and observed vectors of different sizes." << endl;
	}
	//compute the sum of squared errors
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = simHosp[xx]/Npop - obsHosp[xx]/Npop;
		sum_sq += pow(error,2);
	}
	return sum_sq;
}

double sse_calc_int(vector<int> simval, vector<int> obsval){
	
	//verify that the 2 vectors have the same size
	if(simval.size() != obsval.size()){
		cout<< "Warning!!! Simulated and observed vectors of different sizes." << endl;
	}

	//compute the sum of squared errors
	double sum_sq=0.0;
	for (unsigned int xx = 0; xx < obsval.size(); ++xx) {
		double error = 0.0;
		error = (double)( obsval[xx] - simval[xx] ) ;
		sum_sq += pow(error,2);
	}

	return sum_sq;
	
}
double nsse_calc_int(vector<int> simHosp){
	
	//verify that the 2 vectors have the same size
	if(simHosp.size() != obsHosp.size()){
		cout<< "Warning!!! Simulated and observed vectors of different sizes." << endl;
	}
	
	//maximum numbers of cases
	int max_obs = accumulate(obsHosp.begin(), obsHosp.end(), 0);
	cout << "max_obs: " << max_obs <<'\n';
	//compute the sum of squared errors
	double sum_sq=0.0;
	//compute the sum of squared errors
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = (double)simHosp[xx] - (double)obsHosp[xx];
		sum_sq += pow(error,2);
	}
	
	//comute the normalised sum of squared error by the total number of cases
	double nsse = sum_sq / (double)max_obs;

	return nsse;
	
}


double sst_calc(int Npop, vector<double> simHosp){
	double sum_tot=0.0;
	double mean_obs=0.0;
	//verify that the 2 vectors have the same size

	if(simHosp.size() != obsHosp.size()){
		cout << obsHosp.size() << "vs. " << simHosp.size() << endl;
		cout<< "Warning!!! Simulated and observed vectors of different sizes." << endl;
	}

	//compute the mean of the observed value
	mean_obs=mean_calc_int(obsHosp);
	//compute the sum of squared total
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = simHosp[xx]/Npop - mean_obs/Npop;
		sum_tot += pow(error,2);
	}
	return sum_tot;
}

double sst_calc_int(int Npop, vector<int> simHosp){
	double sum_tot=0.0;
	double mean_obs=0.0;
	//verify that the 2 vectors have the same size

	if(simHosp.size() != obsHosp.size()){
		cout << obsHosp.size() << "vs. " << simHosp.size() << endl;
		cout<< "Warning!!! Simulated and observed vectors of different sizes." << endl;
	}

	//compute the mean of the observed value
	mean_obs=mean_calc_int(obsHosp);
	//compute the sum of squared total
	for (unsigned int xx = 0; xx < obsHosp.size(); ++xx) {
		double error = 0.0;
		error = ((double)simHosp[xx]/Npop) - (mean_obs/Npop);
		sum_tot += pow(error,2);
	}
	return sum_tot;
}


void model_select(int smc,particle &outvec, vector<params> fixed_parameters,vector<vector<double>> cfr_byage, 
				vector<double> pf_byage,
				vector<vector<double>> waifw_norm, vector<vector<double>> waifw_sdist, vector<vector<double>> waifw_home, 
				vector <int> agenums, double tau, int duration, seed seedlist, int day_shut, int Npop, gsl_rng * r){

	//---------------------------------------
	// the root model
	//---------------------------------------
	//declare the vector of outputs		int t, t_int, i, j, yr,yr_int,increment;
	vector<int> sim_status;
	vector<int> death_status;
	vector<int> deathH_status;
	vector<vector<int>> ends;
	
	//main loop
	my_model(outvec.parameter_set, fixed_parameters, cfr_byage, pf_byage,waifw_norm, waifw_sdist, waifw_home,
					duration, seedlist, day_shut, Npop,agenums, tau, r, sim_status, ends,
					death_status,deathH_status);

	//---------------------------------------
	// compute the  sum of squared errors
	//---------------------------------------
	double sum_sq_cases = 1000000.0,sum_sq_deaths = 1000000.0,nsse_cases=1000000.0,nsse_deaths=1000000.0;
	vector<int> obs_case = obsHosp;
	vector<int> obs_death = obsDeaths;
	sum_sq_cases= sse_calc_int(sim_status,obs_case);
		
	sum_sq_deaths= sse_calc_int(deathH_status,obs_death);
	//---------------------------------------
	// compute the deviation of the fit expressed as % total number of cases
	//---------------------------------------
	//total numbers of cases and deaths (at hospital)
	int sum_obs_cases = accumulate(obsHosp.begin(), obsHosp.end(), 0);
	int sum_obs_deaths = accumulate(obsDeaths.begin(), obsDeaths.end(), 0);
	
	nsse_cases = sqrt(sum_sq_cases)/ (double)sum_obs_cases;
	nsse_deaths = sqrt(sum_sq_deaths)/ (double)sum_obs_deaths;
 	
	//cout << "nsse_cases: " << nsse_cases << " , nsse_deaths: " << nsse_deaths <<'\n';

	//---------------------------------------
	// Return a vector with all selection measures
	//---------------------------------------
	//	outvec.sum_sq = sum_sq;
	outvec.nsse_cases = nsse_cases;
	outvec.nsse_deaths = nsse_deaths;
	for (unsigned int ii = 0; ii < sim_status.size(); ++ii) {
		outvec.simu_outs.push_back(sim_status[ii]);
		outvec.death_outs.push_back(deathH_status[ii]);
	}
	for (unsigned int ii = 0; ii < ends.size(); ++ii) {
		outvec.end_comps.push_back(ends[ii]);
	}		

}

void parameter_select_first_step(vector<double> &selected_param, vector<double> flag1, vector<double> flag2, gsl_rng * r,int nPar){
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

void parameter_select_nsteps(vector<double> &selected_param, int nPar, gsl_rng * r, vector<particle> particleList, int pick_val,
		double vlimitKernel[], double vect_min[], double vect_Max[]){

	particle perturbList;

    //define the elements of the chosen particle
    perturbList = particleList[pick_val];

	//perturb all elements of the chosen particle with a kernel following a uniform distribution +/- kernelFactor

    for (int xx = 0; xx < nPar; ++xx) {
        double tmpval = perturbList.parameter_set[xx]+vlimitKernel[xx]*gsl_ran_flat(r, -1, 1);
        while(isnan(tmpval) || tmpval<=vect_min[xx] || tmpval>=vect_Max[xx]){
        	tmpval = perturbList.parameter_set[xx]+vlimitKernel[xx]*gsl_ran_flat(r, -1, 1);
        }


    	//return a vector with all selection measures
    	selected_param[xx] = tmpval;
	}
}

void weight_calc(int smc,int pastNpart, vector<particle > pastPart, particle & currentPart, double vlimitKernel[],int nPar){

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
				m=fabs(currentPart.parameter_set[ii] - pastPart[jj].parameter_set[ii]);
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

//transform the timeseries of cummulative cases into incidence
void compute_incidence(vector<int> v,vector<int>& r_val){
	for (unsigned int xx = 0; xx < v.size(); ++xx) {
		int tmp_inc = 0;
		if(xx>0) tmp_inc = v[xx] - v[xx-1];
		r_val.push_back(tmp_inc);
	}	
}


//correct the timeseries to avoid negative incidence records
void correct_incidence(vector<int>& v,vector<int> cumv){
	if ( any_of(v.begin(), v.end(), [](int i){return i<0;}) ){
		for (unsigned int ii = 1; ii < v.size(); ++ii) {
			if(v[ii]<0){
				int case_bef = cumv[ii-1];//look next day number of cases
				int case_aft = cumv[ii+1];//look at the previous number of cases	
				
      		  	//check if the next tot cases is bigger than tot cases before
    			//if(bigger) compute difference and remove it from day before
		        if(case_bef<case_aft){
		          if(v[ii-1] >= abs(v[ii])){
		            cumv[ii-1] += v[ii];
		          } else {
		          	cumv[ii]=0;
		          }
		        } else{
		          //if smaller
		          //check when the day prior today was smaller than the next day and remove extra cases to this day
		          int counter_check=1;
		          while( cumv[ii-counter_check]>case_aft){
		            ++counter_check ;
		          }
		          if(v[ii-counter_check+1] >= abs(v[ii])){
				  	for ( int jj = (ii-counter_check+1); jj < (ii); ++jj) {
						cumv[jj] += v[ii];
					}			
		          } else {
		            cumv[ii]=0;
		          }
		        }
			}

		}
		//recompute the incident cases
		v.clear();
		compute_incidence(cumv,v);
	}
}
  

void select_obs(int& Npop, int& t_index, int& duration, int& day_intro, int& day_shut, 
	vector<int>& obsHosp_tmp, vector<int>& obsDeaths_tmp, 
	vector<vector<int> > data_tmp, vector<vector<int> > death_tmp, int herd_id, 
	int time_back){
		
		
	int maxTime=0;
	vector<int> seqTime;
		
	//define population size
	Npop = data_tmp[herd_id][0];
	
	//create the vector of cases (cummulative) and define time of first detection (index)
	for (unsigned int iv = 1; iv < data_tmp[0].size(); ++iv) {
		if(data_tmp[herd_id][iv]>=0){
			seqTime.push_back(data_tmp[0][iv]);
			maxTime = max_calc_int(maxTime,(int)data_tmp[0][iv]);
			obsHosp_tmp.push_back(data_tmp[herd_id][iv]);
			obsDeaths_tmp.push_back(death_tmp[herd_id][iv]);
			if(data_tmp[herd_id][iv]>0 && t_index<0){
				t_index = (int)data_tmp[0][iv];//identify when the first case is detected/hospitalised 
			}		
		}
	}

	//identify the first day of infectiousness for the index case, which will be the start of our simulation, and add days in the observation, and define duration of the disease process
	int intro = (t_index+1)-(time_back);

	//define the duration of the study period
	if(intro < 0) duration = maxTime - intro + 1;
	else duration = maxTime + 1;

	//define the day of the incursion
	if(intro < 0 ) day_intro = 0;
	else day_intro=intro;
		
	//add the extra information on the observations
	if(duration > obsHosp_tmp.size()){
		int extra_time = duration - obsHosp_tmp.size();
		vector<int> extra_cases, extra_deaths;
		for (int extra = 0; extra < extra_time; ++extra) {
			extra_cases.push_back(0);
			extra_deaths.push_back(0);
		}		
		for (unsigned int iv = 0; iv < obsHosp_tmp.size(); ++iv) {
			extra_cases.push_back(obsHosp_tmp[iv]);
			extra_deaths.push_back(obsDeaths_tmp[iv]);
		}		
		day_shut = day_shut+ extra_time;
		obsHosp_tmp = extra_cases;
		obsDeaths_tmp = extra_deaths;
		extra_cases.clear();	
		extra_deaths.clear();	
	}
	
	cout << "Number of days of obs cases: "<< obsHosp_tmp.size()<<'\n';
	cout << "Number of days of obs deaths: "<< obsDeaths_tmp.size()<<'\n';
	//transform  cumulative numbers into incident cases
	vector<int> obsHosp_tmp2;
	compute_incidence(obsHosp_tmp,obsHosp_tmp2);
	correct_incidence(obsHosp_tmp2,obsHosp_tmp);	
	obsHosp_tmp = obsHosp_tmp2;
	obsHosp_tmp2.clear();	
	
	//transform  cumulative numbers into incident deaths
	vector<int> obsDeaths_tmp2;
	compute_incidence(obsDeaths_tmp,obsDeaths_tmp2);
	correct_incidence(obsDeaths_tmp2,obsDeaths_tmp);	
	obsDeaths_tmp = obsDeaths_tmp2;
	obsDeaths_tmp2.clear();	
}
	

void Lambda(vector<double> &lambda, vector<double> parameter_set,vector<vector<double>> waifw_norm,
			vector<vector<double>> waifw_sdist,vector<vector<double>> waifw_home, 
			vector<vector<int>> pops, int shut){

  double p_i = parameter_set[0];	
  double p_hcw = parameter_set[1];	
  double c_hcw = parameter_set[2];	
  double d_val = parameter_set[3];					
  double q_val = parameter_set[4];					
	

  //lambda=rep(NA,nrow(pops))
  int n_agegroup = waifw_norm.size();
  double quarantined = 0.0;//mean proportion reduction in contacts due to quarantine
  int inf_hosp=0;
  vector<double> I_mat(n_agegroup); 
  //contact matrix
  vector<vector<double>> waifw(n_agegroup, vector<double> (n_agegroup)); 
  //age-dependent transmission rate
  vector<vector<double>> beta(n_agegroup, vector<double> (n_agegroup)); 
  
  for ( int from = 0; from < n_agegroup; ++from) {
	  for ( int to = 0; to < n_agegroup; ++to) {
		  if(shut==0){
		  	//contact network during normal period
		    waifw[from][to]  = waifw_norm[from][to] ;
		  } else {
			//contact network during shutdown period, assuming a proportion d will properly do it
		  	waifw[from][to] = (1.0-d_val) * waifw_sdist[from][to] + (d_val) * waifw_home[from][to];
		  }
		  beta[from][to] = waifw[from][to] * p_i ;
		  if(waifw[from][to]>0) {
			  quarantined += waifw_home[from][to]/waifw[from][to];
		  } else {
			  quarantined += 0;
		  }
	  }
  } 
  
  //mean proportion reduction in contacts due to quarantine accounting for compliance
  quarantined = 1 - (quarantined / (double)(n_agegroup * n_agegroup) ) * (1 - q_val);
  
  //compute the pressure from infectious groups, normalizes by group size
  //compute the pressure from hospitalised, normalized by group size
  for ( int from = 0; from < n_agegroup; ++from) {
	if(from < (n_agegroup-1)) { //n-agegroup-1 to not account for hcw (as different contact)
	 	// are considered those that are infectious asymp, tested or symptomatic only
		I_mat[from] = (double)pops[from][3] + (1-quarantined) * ( (double)pops[from][4] + (double)pops[from][5] );
		//normalisation shouldnt account for dead individuals	
		int tot_pop = accumulate(pops[from].begin(), pops[from].end(), 0);
		I_mat[from] = (double)I_mat[from] / (double)(tot_pop- pops[from][pops[from].size()-1] ); 
				
	}
	//sum up of number of hospitalised cases 
	inf_hosp+= pops[from][6];
  }

  //sum up infectious pressure from each age group
  for ( int from = 0; from < (n_agegroup); ++from) {
	  for(int ii = 0; ii < (n_agegroup-1); ++ii){
	  	lambda[from] += ( beta[from][ii]*I_mat[ii] );
	  }
  }

  int tot_pop = accumulate(pops[n_agegroup-1].begin(), pops[n_agegroup-1].end(), 0);
  int rem_pop = pops[n_agegroup-1][4]+pops[n_agegroup-1][5]+pops[n_agegroup-1][6]+pops[n_agegroup-1][8];
  lambda[n_agegroup-1] = p_hcw * c_hcw * ((double)inf_hosp/(double)(tot_pop - rem_pop) );
 
}


void infspread(gsl_rng * r, vector<int>& pop, int& deaths, int& deathsH, int& detected, 
	params fixed_parameters, vector<double> parameter_set, vector<double> cfr_tab, double pf_val, double lambda){
		
    unsigned int S=pop[0], E=pop[1], E_t=pop[2], I=pop[3], I_t=pop[4], I_s=pop[5], H=pop[6], R=pop[7], D=pop[8] ;
	
	double T_lat= fixed_parameters.T_lat;
	double p_s= fixed_parameters.p_s;
	double T_inf= fixed_parameters.T_inf;
	double T_rec= fixed_parameters.T_rec;
	double T_sym= fixed_parameters.T_sym;
	double T_hos= fixed_parameters.T_hos;
	
	double p_h= cfr_tab[0];
	double cfr= cfr_tab[1];	
	double p_d= cfr_tab[2];	
	//double p_dc= cfr_tab[3];
		
	double rrdh = parameter_set[5];
	double rrdc = parameter_set[6];
	double rrh = parameter_set[7];	
	
	double A_val = (pf_val * rrdh + (1 - pf_val ) ) * p_d;
	double B_val = pf_val * (1 - rrdh * p_d) + (1 - pf_val ) * (1 - p_d);
	double C_val = pf_val  * rrh * p_h ;
	double D_val = (1 - pf_val) * p_h;
	double E_val = (( pf_val - C_val ) * rrdc + ((1 - pf_val ) - D_val ) ) * cfr;
	double F_val = ( pf_val - C_val ) * (1 - rrdc * cfr ) + ((1 - pf_val ) - D_val ) * ( 1 - cfr ) ;
	
    // hospitalized
    unsigned int newdeathsH = gsl_ran_poisson(r, A_val * (1 / T_hos) * (double)H);
	newdeathsH = min(H,newdeathsH);
    D += newdeathsH;
	H -= newdeathsH;
    unsigned int recoverH = gsl_ran_poisson(r,  B_val * (1 / T_hos) * (double)H);//hospitalized recover
	recoverH=min(H,recoverH) ;
    R += recoverH;
	H -= recoverH;
	
    // symptomatic
    unsigned int hospitalize = gsl_ran_poisson(r, (C_val + D_val) * (1 / T_sym) * (double)I_s); //symptomatic become hospitalized
    hospitalize = min(I_s, hospitalize);
	H += hospitalize;
	I_s -= hospitalize;
	
    unsigned int newdeathsC = gsl_ran_poisson(r, E_val * ( 2 / T_hos) * (double)I_s); //symptomatic die at home / in communities	
    newdeathsC = min(I_s, newdeathsC);
    D += newdeathsC;
	I_s -= newdeathsC;
	
    unsigned int recoverI_s = gsl_ran_poisson(r, F_val * ( 1 / T_rec) * (double)I_s); //symptomatic recover at home / in communities	
    recoverI_s = min(I_s, recoverI_s);
    R += recoverI_s;
	I_s -= recoverI_s;
	
    // infectious
    unsigned int symptomsI = gsl_ran_poisson(r, p_s * ( 1/T_inf ) * (double)I ) ; //infectious become symptomatic
    symptomsI = min(I, symptomsI);
	unsigned int symptomsI_t=gsl_ran_poisson(r, p_s * ( 1/T_inf ) * (double)I_t); //tested infectious become symptomatic
    symptomsI_t=min(I_t,symptomsI_t);
	I_s += (symptomsI+symptomsI_t);
	I -= symptomsI;
	I_t -= symptomsI_t;
	
    unsigned int recoverI = gsl_ran_poisson(r, (1 - p_s) * ( 1/T_inf ) * (double)I); //infectious recover
    recoverI = min(I, recoverI);
	unsigned int recoverI_t = gsl_ran_poisson(r, (1 - p_s) * ( 1/T_inf ) * (double)I_t); //tested infectious recover
	recoverI_t = min(I_t, recoverI_t);
    R += (recoverI+recoverI_t) ;
	I -= recoverI;
	I_t -= recoverI_t;
	
    // latent
    unsigned int infectious = gsl_ran_poisson(r, (double)E/T_lat); //latent become infectious
	infectious = min(E, infectious);
    I += infectious;
	E -= infectious;
    unsigned int infectious_t = gsl_ran_poisson(r, (double)E_t/T_lat); //tested latent become infectious
    infectious_t = min(E_t, infectious_t);
	I_t += infectious_t;
	E_t -= infectious_t;
	
    // susceptible
    unsigned int newinfection= gsl_ran_poisson(r, (double)S*lambda); //susceptible become infected
    newinfection= min(S, newinfection);
	S -= newinfection;
	E += newinfection;
	
    // recover population
	vector<int> newpop= {(int)S,(int)E,(int)E_t,(int)I,(int)I_t,(int)I_s,(int)H,(int)R,(int)D};
	pop = newpop;
	deaths += (newdeathsC + newdeathsH);
	deathsH += newdeathsH;
	detected += (H + I_t);	
	
}


void my_model(vector<double> parameter_set, vector<params> fixed_parameters, vector<vector<double>> cfr_byage, 
				vector<double> pf_byage,
				vector<vector<double>> waifw_norm, vector<vector<double>> waifw_sdist, vector<vector<double>> waifw_home,
				int duration, seed seedlist, int day_shut, int Npop,
				vector<int> agenums, double tau, gsl_rng * r, vector<int> &sim_status, vector<vector<int>> &ends,
				vector<int> &death_status,vector<int> &deathH_status){

	int n_agegroup = waifw_norm.size();
	int inLockdown = 0;
	double seed_pop[6];
	

	//sets up an array for the population at each timestep in each age and disease category	
	//also set up the age distribution of old ages as target for disease introduction
	vector<vector<int>> poparray(n_agegroup, vector<int>(9));	  
	for ( int age = 0; age < (n_agegroup); ++age) {
		for ( int st = 0; st < 9; ++st) {
			if(st ==0) {
				poparray[age][st]=agenums[age]; //set up the starting population as fully susceptible
			} else {
				poparray[age][st] = 0;
			}
		}
		if(age >0 & age < (n_agegroup-1) ){ //seed only occurs in >20yo and not in HCW
			seed_pop[age-1] = (double)agenums[age];
		}
	}


	//introduce disease at t=0. if seedmethod != "background"
	//WARNING!! MINOR BUG HERE: gsl multinomial sometimes freeze for unknown reasons. replaced by random flat, picking a value where to seed infection. will only works because seed=1.
	if(seedlist.seedmethod != "background"){
	//	size_t k = sizeof(seed_pop);
	//	unsigned int startdist[k];
		vector<int> startdist = {0,0,0,0,0,0};
		startdist[(int)gsl_ran_flat(r, 0, 6)] = seedlist.nseed;
	//	gsl_ran_multinomial(r, k, startdz, seed_pop,startdist); //distribute the diseased across the older age categories

		for ( int age =1; age < (n_agegroup-1); ++age) {
			poparray[age][0] -=  startdist[age-1];//take diseased out of S
			poparray[age][3] +=  startdist[age-1];//put diseased in I	
		}	
	}


  	//initialize saving of the detection for t=0
	sim_status.push_back(0); 
	death_status.push_back(0); 
	deathH_status.push_back(0); 
	//run the simulation
	for (int tt = 1; tt < (int)ceil(duration/tau); ++tt) {
		//initialize return value
		int deaths=0;
		int deathsH=0;
		int detected=0;
		vector<int> popcomp;
		
		//identify the lock down
		if(tt >= day_shut)  inLockdown = 1;
		
		
		//introduce disease from background infection until lockdown
	    if(inLockdown<1){
			if(seedlist.seedmethod == "background"){
				double bkg_lambda = parameter_set[parameter_set.size()-1];
				unsigned int startdz = gsl_ran_poisson(r, (double)Npop * bkg_lambda);	//how many diseased is introduced in each given day before lockdown
		
			//	size_t k = sizeof(seed_pop);
			//	unsigned int startdist[k];
				vector<int> startdist = {0,0,0,0,0,0};
				int pickcomp = (int)gsl_ran_flat(r, 0, 6);
				while(poparray[pickcomp][0]<1){
					pickcomp = (int)gsl_ran_flat(r, 0, 6);
				} 
				startdist[pickcomp] = startdz;
			//	gsl_ran_multinomial(r, k, startdz, seed_pop,startdist); //distribute the diseased across the older age categories

				for ( int age =1; age < (n_agegroup-1); ++age) {
					int nseed = min(poparray[age][0], startdist[age-1]);
					poparray[age][0] -=  nseed;//take diseased out of S
					poparray[age][1] +=  nseed;//put diseased in E	
				}	
			}
	    }


		//compute the forces of infection
		vector<double> lambda(n_agegroup);
		Lambda(lambda, parameter_set,waifw_norm, waifw_sdist,waifw_home, poparray, inLockdown);	
	
		//step each agegroup through infections
		for ( int age = 0; age < (n_agegroup); ++age) {	
			infspread(r, poparray[age], deaths, deathsH, detected,fixed_parameters[age],parameter_set,
				cfr_byage[age],pf_byage[age],lambda[age]);
		}

		//cout << tt << " , " <<  detected << " , " << obsHosp[tt]<< " , "<< deaths << '\n';

		sim_status.push_back(detected); 
		death_status.push_back(deaths); 
		deathH_status.push_back(deathsH);
    }	
	//save the population in each epi compt for the last day
	for ( int age = 0; age < n_agegroup; ++age) {
		ends.push_back(poparray[age]);
	}
}




void read_parameters(int& herd_id, double& tau, int&num_threads,int& nsteps, int& nParticLimit, int& nSim, double& kernelFactor, 
		vector<double>& toleranceLimit, params& paramlist, seed& seedlist, int& day_shut, int& totN_hcw,
		int& nPar, double& prior_pinf_shape1,double& prior_pinf_shape2,double& prior_phcw_shape1,double& prior_phcw_shape2,
		double& prior_chcw_mean,double& prior_d_shape1,double& prior_d_shape2,double& prior_q_shape1, double& prior_q_shape2,
		double& prior_rrdh_shape1,double& prior_rrdh_shape2,double& prior_rrdc_shape1,double& prior_rrdc_shape2,
		double& prior_rrh_shape1,double& prior_rrh_shape2,double& prior_lambda_shape1, double& prior_lambda_shape2){
	string parameterfile = "parameters.ini";
	CIniFile parameters; //Create a new variable of type CIniFile (see Inifile.h for type definition)

	//Settings
	herd_id = atoi(parameters.GetValue("shb_id", "Settings", parameterfile).c_str());

	//time step for the modelling process
	tau = atof(parameters.GetValue("tau", "Settings", parameterfile).c_str());

	//number of threads used for computation
	num_threads = atoi(parameters.GetValue("num_threads", "Settings", parameterfile).c_str());

	//Seed settings
	seedlist.seedmethod = parameters.GetValue("seedmethod", "Seed settings", parameterfile).c_str();
	if(seedlist.seedmethod == "random"){
		seedlist.nseed = atoi(parameters.GetValue("nseed", "Seed settings", parameterfile).c_str());
	} else if(seedlist.seedmethod == "background"){
		seedlist.hrp = atoi(parameters.GetValue("hrp", "Seed settings", parameterfile).c_str());
	} else {
		cout<< "Warning!!! Unknown method - using random seed method instead." << endl;
		seedlist.nseed = atoi(parameters.GetValue("nseed", "Seed settings", parameterfile).c_str());
	}
	
	//Fit settings
	nsteps = atoi(parameters.GetValue("nsteps", "Fit settings", parameterfile).c_str());
	nParticLimit = atoi(parameters.GetValue("nParticLimit", "Fit settings", parameterfile).c_str());
	nSim = atoi(parameters.GetValue("nSim", "Fit settings", parameterfile).c_str());
	kernelFactor = atof(parameters.GetValue("kernelFactor", "Fit settings", parameterfile).c_str());

	//Tolerance settings
	for (int ii = 1; ii <= nsteps; ++ii) {
		toleranceLimit.push_back(0.0);
	}

	for (int ii = 1; ii <= nsteps; ++ii) {
		if(ii==1) toleranceLimit[ii-1] = atof(parameters.GetValue("Key1", "Tolerance settings", parameterfile).c_str());
		if(ii==2) toleranceLimit[ii-1] = atof(parameters.GetValue("Key2", "Tolerance settings", parameterfile).c_str());
		if(ii==3) toleranceLimit[ii-1] = atof(parameters.GetValue("Key3", "Tolerance settings", parameterfile).c_str());
		if(ii==4) toleranceLimit[ii-1] = atof(parameters.GetValue("Key4", "Tolerance settings", parameterfile).c_str());
		if(ii==5) toleranceLimit[ii-1] = atof(parameters.GetValue("Key5", "Tolerance settings", parameterfile).c_str());
		if(ii==6) toleranceLimit[ii-1] = atof(parameters.GetValue("Key6", "Tolerance settings", parameterfile).c_str());
		if(ii==7) toleranceLimit[ii-1] = atof(parameters.GetValue("Key7", "Tolerance settings", parameterfile).c_str());
		if(ii==8) toleranceLimit[ii-1] = atof(parameters.GetValue("Key8", "Tolerance settings", parameterfile).c_str());
		if(ii==9) toleranceLimit[ii-1] = atof(parameters.GetValue("Key9", "Tolerance settings", parameterfile).c_str());
		if(ii==10) toleranceLimit[ii-1] = atof(parameters.GetValue("Key10", "Tolerance settings", parameterfile).c_str());
	}


	//Fixed parameters
	
	paramlist.T_lat = atof(parameters.GetValue("T_lat", "Fixed parameters", parameterfile).c_str());
	paramlist.p_s = atof(parameters.GetValue("p_s", "Fixed parameters", parameterfile).c_str());
	paramlist.juvp_s = atof(parameters.GetValue("juvp_s", "Fixed parameters", parameterfile).c_str());
	paramlist.T_inf = atof(parameters.GetValue("T_inf", "Fixed parameters", parameterfile).c_str());
	paramlist.T_rec = atof(parameters.GetValue("T_rec", "Fixed parameters", parameterfile).c_str());
	paramlist.T_sym = atof(parameters.GetValue("T_sym", "Fixed parameters", parameterfile).c_str());
	paramlist.T_hos = atof(parameters.GetValue("T_hos", "Fixed parameters", parameterfile).c_str());
	day_shut = atoi(parameters.GetValue("day_shut", "Fixed parameters", parameterfile).c_str());
	totN_hcw = atoi(parameters.GetValue("totN_hcw", "Fixed parameters", parameterfile).c_str());

	//priors settings

	nPar = atoi(parameters.GetValue("nPar", "Priors settings", parameterfile).c_str());
	prior_pinf_shape1 = atof(parameters.GetValue("prior_pinf_shape1", "Priors settings", parameterfile).c_str());
	prior_pinf_shape2 = atof(parameters.GetValue("prior_pinf_shape2", "Priors settings", parameterfile).c_str());
	prior_phcw_shape1 = atof(parameters.GetValue("prior_phcw_shape1", "Priors settings", parameterfile).c_str());
	prior_phcw_shape2 = atof(parameters.GetValue("prior_phcw_shape2", "Priors settings", parameterfile).c_str());
	prior_chcw_mean = atof(parameters.GetValue("prior_chcw_mean", "Priors settings", parameterfile).c_str());
	prior_d_shape1 = atof(parameters.GetValue("prior_d_shape1", "Priors settings", parameterfile).c_str());
	prior_d_shape2 = atof(parameters.GetValue("prior_d_shape2", "Priors settings", parameterfile).c_str());
	prior_q_shape1 = atof(parameters.GetValue("prior_q_shape1", "Priors settings", parameterfile).c_str());
	prior_q_shape2 = atof(parameters.GetValue("prior_q_shape2", "Priors settings", parameterfile).c_str());
	prior_rrdh_shape1= atof(parameters.GetValue("prior_rrdh_shape1", "Priors settings", parameterfile).c_str());
	prior_rrdh_shape2=atof(parameters.GetValue("prior_rrdh_shape2", "Priors settings", parameterfile).c_str());
	prior_rrdc_shape1=atof(parameters.GetValue("prior_rrdc_shape1", "Priors settings", parameterfile).c_str());
	prior_rrdc_shape2=atof(parameters.GetValue("prior_rrdc_shape2", "Priors settings", parameterfile).c_str());
	prior_rrh_shape1=atof(parameters.GetValue("prior_rrh_shape1", "Priors settings", parameterfile).c_str());
	prior_rrh_shape2=atof(parameters.GetValue("prior_rrh_shape2", "Priors settings", parameterfile).c_str());	
	
	prior_lambda_shape1 = atof(parameters.GetValue("prior_lambda_shape1", "Priors settings", parameterfile).c_str());
	prior_lambda_shape2 = atof(parameters.GetValue("prior_lambda_shape2", "Priors settings", parameterfile).c_str());
}

double mean_calc_double( vector<double> v){
	double return_value = 0.0;
	int n = v.size();
	for(vector<double>::iterator j=v.begin();j!=v.end();++j)  return_value += *j;
    return ( return_value / n);
}

double mean_calc_int( vector<int> v){
	int return_value = 0;
	int n = v.size();
	for(vector<int>::iterator j=v.begin();j!=v.end();++j)  return_value += *j;
    return ( (double)return_value / n);
}

void cumsum_calc(vector<double> v,vector<double>& r_val){
	int size_vec = v.size();
	double tmp_sum=0.0;
	for (int xx = 0; xx < size_vec; ++xx) {
		tmp_sum += v[xx];
		r_val.push_back(tmp_sum);
	}
}

void cumsum_calc_int(vector<int> v,vector<int>& r_val){
	int size_vec = v.size();
	int tmp_sum=0;
	for (int xx = 0; xx < size_vec; ++xx) {
		tmp_sum += v[xx];
		r_val.push_back(tmp_sum);
	}
}

double max_calc(double a,double b){
	double rval;
	if(a>=b) rval = a;
	else rval = b;
	return rval;
}

int max_calc_int(int a,int b){
	int rval;
	if(a>=b) rval = a;
	else rval = b;
	return rval;
}

double sum_calc(vector<double> v){
    double rval = 0.0;
	for(vector<double>::iterator j=v.begin();j!=v.end();++j)  rval += *j;
	return rval;
}

int sum_calc_int(vector<int> v){
    int rval = 0;
	for(vector<int>::iterator j=v.begin();j!=v.end();++j)  rval += *j;
	return rval;
}

int colsum_calc_int(vector<vector<int>> m, int colid){
    int rval = 0;
	for (unsigned int j = 0; j < m.size(); ++j) {
		rval += m[j][colid];
	}
	return rval;
}

double colsum_calc(vector<vector<double>> m, int colid){
    double rval = 0.0;
	for (unsigned int j = 0; j < m.size(); ++j) {
		rval += m[j][colid];
	}
	return rval;
}


