//PowerLawDisk.cpp
//Basic disk modeling example, described in section 2.2.1
#include "magrathea/magrathea.h"

struct advancedTemp:public temperature_base{
	interpSurface tempSurf;
	advancedTemp():tempSurf(){}
	advancedTemp(const advancedTemp& other):tempSurf(other.tempSurf){}
	
	advancedTemp(std::ifstream& gridFile, std::ifstream& dataFile):tempSurf(gridFile, dataFile) {}
	
	double operator()(double r, double theta, double phi) const{
		return(tempSurf.interp2d(r,theta));
	}
	
};

int main(int argc, char* argv[]){
	const size_t nThreadsMax=std::thread::hardware_concurrency();
	ThreadPool pool(nThreadsMax);
	astroParams diskData(246.6,-24.72,3.18e18*100,pi/4,pi/4);
	std::shared_ptr<powerLawDisk> dptr=std::make_shared<powerLawDisk>(1,100*AU,5*AU,1.0,1.25);
	std::ifstream gridfile("data/amr_grid.inp");
	std::ifstream datafile("data/dust_temperature_phi0.ascii");
	std::shared_ptr<advancedTemp> tptr=std::make_shared<advancedTemp>(gridfile,datafile);
	std::shared_ptr<dustOpacity> doptr=std::make_shared<dustOpacity>("data/dustopac.txt");
	
	grid g(0.1*AU,150*AU,75*pi/180,105*pi/180, 2.18*mSun, dptr, tptr, doptr);
	image img(500, 500,250*AU, 250*AU, {2.30538e11}, diskData);
	vect offset(0,0,0);
	
	img.propagate(g, grid::continuum, offset, pool, false, 100, 5);
	img.printToFits("PowerLawDisk.fits");
	
	return 0;
}