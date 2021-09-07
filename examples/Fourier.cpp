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
	unsigned int nfreqs=240;
	double freqRange=2e7; //1MHz
	double freqStep=freqRange/nfreqs;
	std::vector<double> frequenciesBig;	
	std::vector<double> frequenciesSmall;
	std::cout.precision(10);
	for(int i=0;i<nfreqs;i++){
		double freq=2.30538e11+freqStep*i-freqRange/2.0;
		frequenciesBig.push_back(freq);
		if(!((i-10)%20)) frequenciesSmall.push_back(freq);
		
	}
	image img(500, 500,250*AU, 250*AU, frequenciesBig, diskData);
	vect offset(0,0,0);
	
	img.propagate(g, grid::normal, offset, pool, true, 100, 5);
	//img.printToFits("FourierDisk.fits");
	
	std::cout << "averaging" << std::endl;
	image imgAvg(500, 500,250*AU, 250*AU, frequenciesSmall, diskData);
	for(int f=0;f<frequenciesSmall.size();f++){
		for(int i=0;i<img.vpix;i++){
			for(int j=0;j<img.hpix;j++){
				//first we average the values in the full image
				double value=0;
				for(int subF=0;subF<20;subF++){
					value+=img.data[f*20+subF][i][j];
				}
				value/=20;
				//then add to the final image
				imgAvg.data[f][i][j]=value;
			}
		}
	}
	imgAvg.printToFits("FourierDisk.fits");
	//calculate the FFT
	fourierImage FFTs=FFT(imgAvg);
	FFTs.printToFits("test.fits");
	
	
	//now we do the beam smoothing
	beam bm(imgAvg, 0.030,0.020,pi/4);
	for(int f=0;f<frequenciesSmall.size();f++){
		for(int i=0;i<img.vpix;i++){
			for(int j=0;j<img.hpix;j++){
				FFTs.realPart[f][i][j]*=bm.realPart[i][j];
				FFTs.imaginaryPart[f][i][j]*=bm.realPart[i][j];
			}
		}
	}
	FFTs.printToFits("test2.fits");
	//now transform back
	image final=backTransform(FFTs);
	
	final.printToFits("FourierDisk2.fits");
	
	return 0;
}