//PowerLawDiskCO.cpp
//Basic disk modeling example, described in section 2.2.1
#include "magrathea/magrathea.h"
#include <random>

double goodnessOfFit(image& model, const image& obs){
	double sum=0;
	for(int i=0;i<obs.hpix;i++){
		for(int j=0;j<obs.vpix;j++){
			double realPoint=obs.data[0][i][j];
			double modelPoint=model.data[0][i][j];
			sum+=(realPoint-modelPoint)*(realPoint-modelPoint);
		}
	}
	//and just for good measure, the residuals
	/*image resid=obsFixed;
	for(int i=0;i<resid.hpix;i++){
		for(int j=0;j<resid.vpix;j++){
			resid.data[0][i][j]-=final.data[0][i][j];
		}
	}
	resid.printToFits("resid.fits");*/
	return -sum;
}

int main(int argc, char* argv[]){
	const size_t nThreadsMax=std::thread::hardware_concurrency();
	ThreadPool pool(nThreadsMax);
	
	ParameterSet params;
	params.addParameter("rc"); //units of au
	params.setParameterLowerLimit("rc",40); params.setParameterUpperLimit("rc",500);
	params.addParameter("P"); //unitless
	params.setParameterLowerLimit("P",-1); params.setParameterUpperLimit("P",1.5);
	
	std::mt19937 rng(137);
	auto randInRange=[&rng](const double min, const double max){
		double range = max-min;
		double randNum=rng()/(double)rng.max();
		return min+randNum*range;
	};
	std::vector<std::vector<double>> initialEnsemble; //the starting positions
	unsigned int numWalkers=10;
	//initialize the ensemble
	for(int i=0;i<numWalkers;i++){
		std::vector<double> ensembleMember;
		ensembleMember.push_back(randInRange(50,55));//rc
		ensembleMember.push_back(randInRange(1.1,1.3));//P
		initialEnsemble.push_back(ensembleMember);
	}
	
	image data=fitsExtract("data/fitData.fits",{2.30538e11},3.18e20);
	
	astroParams diskData(246.6,-24.72,3.18e18*100,pi/4,pi/4);
	std::shared_ptr<powerLawDisk> dptr=std::make_shared<powerLawDisk>(1,100*AU,5*AU,1.0,1.25);
	std::ifstream gridfile("data/amr_grid.inp");
	std::ifstream datafile("data/dust_temperature_phi0.ascii");
	std::shared_ptr<fileTemp> tptr=std::make_shared<fileTemp>(gridfile,datafile);
	std::shared_ptr<dustOpacity> doptr=std::make_shared<dustOpacity>("data/dustopac.txt");
	
	grid g(0.1*AU,150*AU,75*pi/180,105*pi/180, 2.18*mSun, dptr, tptr, doptr);
	image model(500, 500,250*AU, 250*AU, {2.30538e11}, diskData);
	vect offset(0,0,0);
	
	auto LLH=[&data,&params,&dptr,&g,&model,&offset,&diskData,&pool](const std::vector<double>& c){
		double rc=params.extractParameter("rc",c);
		double P=params.extractParameter("P",c);
		
		dptr->rc=rc*AU;
		dptr->P=P;
		
		model.propagate(g, grid::continuum, offset, pool, false, 50, 2);
		double fit=-goodnessOfFit(model,data)*10000000;
		
		std::cout << rc << "\t" << P << "\t" << fit << std::endl;
		
		return fit;
	};
	
	int nSamples=1000;
	auto samples=markovSample(LLH,stretchMove(),params,nSamples,0,0,rng,initialEnsemble);
	
	std::cout << "done" << std::endl;
	
	//print the final set
	std::ofstream outfile("samples.txt");
	outfile.precision(10);
	for(int i=0;i<samples.size();i++){
		for(int j=0;j<params.numberOfParameters();j++){
			outfile << samples[i].coordinates[j] << "\t";
		}
		outfile <<  samples[i].value << std::endl;
	}
	outfile.close();

	functionEvaluation best;
	best.value=1e9;//start very large
	for(auto sample : samples){
		if(sample.value < best.value)
			best = sample;
	}
	std::cout << "best fit is: ";
	for(int i=0; i< best.coordinates.size(); i++){
		std::cout << best.coordinates[i] << "\t";
	}
	std::cout << best.value << std::endl;
	
	return 0;
}