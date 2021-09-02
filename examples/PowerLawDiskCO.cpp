//PowerLawDiskCO.cpp
//Basic disk modeling example, described in section 2.2.1
#include "magrathea/magrathea.h"

int main(int argc, char* argv[]){
	const size_t nThreadsMax=std::thread::hardware_concurrency();
	ThreadPool pool(nThreadsMax);
	astroParams diskData(246.6,-24.72,3.18e18*100,pi/4,pi/4);
	std::shared_ptr<powerLawDisk> dptr=std::make_shared<powerLawDisk>(0.01,100*AU,1*AU,0.5,0.5);
	std::shared_ptr<powerLawDisk> gptr=std::make_shared<powerLawDisk>(14,100*AU,5*AU,1.0,1.25);
	std::ifstream gridfile("data/amr_grid.inp");
	std::ifstream datafile("data/dust_temperature_phi0.ascii");
	std::shared_ptr<fileTemp> tptr=std::make_shared<fileTemp>(gridfile,datafile);
	std::shared_ptr<dustOpacity> doptr=std::make_shared<dustOpacity>("data/dustopac.txt");
	std::shared_ptr<COopac> goptr=std::make_shared<COopac>(1,COopac::iso_12CO);
	
	grid g(0.1*AU,150*AU,75*pi/180,105*pi/180, 2.18*mSun, gptr, dptr, tptr, doptr, goptr, false, 0);
	//set up the frequency
	//center around the 12CO 2-1 line, with 50KHz spacing
	unsigned int nfreqs=20;
	double freqRange=2e7; //1MHz
	double freqStep=freqRange/nfreqs;
	std::vector<double> frequencies;
	std::cout.precision(10);
	for(int i=0;i<nfreqs;i++){
		double freq=2.30538e11+freqStep*i-freqRange/2.0;
		frequencies.push_back(freq);
	}
	image img(2000, 2000,250*AU, 250*AU, frequencies, diskData);
	vect offset(0,0,0);
	
	img.propagate(g, grid::normal, offset, pool, true, 100, 5);
	img.printToFits("PowerLawDiskCO.fits");
	
	return 0;
}