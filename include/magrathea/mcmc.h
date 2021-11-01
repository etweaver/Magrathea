//mcmc.h
//support structs and functions for sampling, and for non-distributed fitting

#include <iostream>
#include <random>
#include <vector>

#include "ParameterSet.h"

struct functionEvaluation{
	std::vector<double> coordinates;
	double value;
};

struct moveProposal{
	std::vector<double> coordinates;
	double logAcceptanceFactor;
	
	moveProposal()=default;
	moveProposal(std::size_t n):coordinates(n){}
	moveProposal(moveProposal&& mp):
	coordinates(std::move(mp.coordinates)),logAcceptanceFactor(mp.logAcceptanceFactor){}
};

template<typename EnsembleMember>
struct coordinateExtractor{
	//we only need a version of this for the special case of std::vector<double>
};

template<>
struct coordinateExtractor<std::vector<double> >{
	const std::vector<double>& operator()(const std::vector<double>& ensembleMember)const{
		return ensembleMember;
	}
};

template<typename EnsembleMember = std::vector<double> >
struct stretchMove{

	//square root distribution
	struct sqrtDist{
		double range; //extends from 1/r to r
		double norm;
		mutable std::uniform_real_distribution<double> uniform; //mutable is needed to override a const later because the
		//uniform distribution has no const call operator. Ugly but necessary.

		sqrtDist(double range): range(range),norm(1/((sqrt(range)-sqrt(1/range)))),uniform(norm*sqrt(1/range),norm*sqrt(range)){
			if(range<=1)
				throw std::domain_error("square_root_distribution requires range>1");
		}

		template<typename RNG>
		double operator()(RNG& rng) const{
			double v=uniform(rng)/(norm);
			return(v*v);
		}
	};

	sqrtDist jumpDist;
	stretchMove():jumpDist(2){}

	template<typename RNG>
	moveProposal operator()(const std::vector<double>& coords, const std::vector<EnsembleMember>& ensemble, RNG& rng) const {
		coordinateExtractor<EnsembleMember> coordExtractor;
		assert(!coords.empty());
		assert(ensemble.size() > 1);
		moveProposal mp(coords.size());

		//choose a random member of the ensemble, but not the one we are already using
		std::uniform_int_distribution<int> idxDist(0,ensemble.size()-1);
		int idx=idx=idxDist(rng);;
		unsigned int maxTrials=100; //if we need to keep trying
		//std::cout << coords.size() << "\t" << ensemble.size() << "\t" << idx << std::endl;
		while(std::equal(coords.begin(),coords.end(),coordExtractor(ensemble[idx]).begin())){
			idx=idxDist(rng);
			if(!--maxTrials)
				throw std::runtime_error("StretchMove failed too many times to find a distinct ensemble member. "
					"Ensmeble may have become degenerate.");
		}

		//jump distance
		double z=jumpDist(rng);

		//direction
		for(std::size_t i=0; i<coords.size(); i++)
			mp.coordinates[i]=coordExtractor(ensemble[idx])[i]+z*(coords[i]-coordExtractor(ensemble[idx])[i]);
		//and the penalty associated with going there
		mp.logAcceptanceFactor=(coords.size()-1)*log(z);

		return mp;
	}
};

template<typename Dist, typename Jumper, typename RNG>
std::vector<functionEvaluation>
markovSample(const Dist& distribution, const Jumper& jumper, 
			 const ParameterSet& parameters,
			 int nSamples, int burnIn, int skip, RNG& rng,
			 const std::vector<std::vector<double>>& initialEnsemble){

	int ensembleSize=initialEnsemble.size();
	std::vector<functionEvaluation> results;
	results.reserve(nSamples);
	int skipCounter=skip;
	unsigned int totalSteps=0;
	unsigned int acceptedSteps=0;
	
	std::vector<std::vector<double>> coordinates(ensembleSize);
	for(size_t i=0; i<ensembleSize; i++){
		//this is sort of a trick to hide parameters that you don't want to fit for right now
		//that way the minimizer only sees what it needs.
		coordinates[i].resize(parameters.numberOfFreeParameters());
		parameters.extractFreeParameters(initialEnsemble[i],coordinates[i]);
	}
	std::vector<std::vector<double>> oldCoordinates(ensembleSize);
	std::vector<double> logLikelihood(ensembleSize), lastLogLikelihood(ensembleSize);
	
	std::uniform_real_distribution<double> acceptDist(0,1);
	std::vector<double> baseCoordinates=parameters.getParameterValues();
	std::vector<double> proposedCoordinates(parameters.numberOfParameters());
	
	//initialize the logLikelihood vector
	for(unsigned int i=0; i<ensembleSize; i++){
		std::copy(baseCoordinates.begin(),baseCoordinates.end(),proposedCoordinates.begin());
		parameters.insertFreeParameters(coordinates[i],proposedCoordinates);
		logLikelihood[i]=distribution(proposedCoordinates);
	}
	while(results.size()<nSamples){
		oldCoordinates=coordinates;
		lastLogLikelihood=logLikelihood;
		
		//update the whole ensemble
		//std::cout << "updating ensemble" << std::endl;
		for(unsigned int i=0; i<ensembleSize; i++){
			moveProposal mp=jumper(coordinates[i],oldCoordinates,rng);
			
			std::copy(baseCoordinates.begin(),baseCoordinates.end(),proposedCoordinates.begin());
			parameters.insertFreeParameters(mp.coordinates,proposedCoordinates);
			if(!parameters.inBounds(proposedCoordinates)){
				//std::cout << "Step rejected. Out of bounds." << std::endl;
				/*for(auto coord : proposedCoordinates){
					std::cout << coord/AU << "\t";
				}
				std::cout << std::endl;*/
				continue; //check bounds
			}
			double newLogLikelihood=distribution(proposedCoordinates);
			double log_ratio=(lastLogLikelihood[i]-newLogLikelihood)+mp.logAcceptanceFactor;
			//std::cout << mp.coordinates[0] << "\t" << mp.coordinates[1] << "\t" << mp.coordinates[2] << "\t" << log_ratio << std::endl;
			
			double randomAcceptThreshold=log(acceptDist(rng));
			bool accept=(log_ratio>=0) || (log_ratio>randomAcceptThreshold);
			/*for(auto coord : mp.coordinates){
				std::cout << coord << "\t";
			}
			std::cout << std::endl;*/
			totalSteps++;
			//std::cout << "log ratio: " << log_ratio << "\trandom accept threshold: " << randomAcceptThreshold << std::endl;
			if(accept){
				std::cout << "accepted" << std::endl;
				acceptedSteps++;
				std::copy(mp.coordinates.begin(),mp.coordinates.end(),coordinates[i].begin());
				logLikelihood[i]=newLogLikelihood;
			}
		}
		//std::cout << std::endl;
		//decide whether to use the current positions
		if(burnIn){
			//std::cout << "burn in" << std::endl;
			burnIn--;
			continue;
		}
		if(skipCounter){
			//std::cout << "skip" << std::endl;
			skipCounter--;
			continue;
		}
		else
			skipCounter=skip;
		
		for(unsigned int i=0; i<ensembleSize && results.size()<nSamples; i++){
			//augment internal coordinates with all fixed parameters before
			//returning to the user
			std::copy(baseCoordinates.begin(),baseCoordinates.end(),proposedCoordinates.begin());
			parameters.insertFreeParameters(coordinates[i],proposedCoordinates);
			//std::cout << "result: " << proposedCoordinates[0] << "\t" << proposedCoordinates[1] << "\t" << proposedCoordinates[2] << "\t" << logLikelihood[i] << std::endl;
			results.push_back(functionEvaluation{proposedCoordinates,logLikelihood[i]});
		}
	}
	
	std::cout << "Sampling complete. Acceptance ratio: " << (double)acceptedSteps/(double)totalSteps << std::endl;
	return results;			 
}