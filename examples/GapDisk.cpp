//PowerLawDisk.cpp
//Basic disk modeling example, described in section 2.2.1
#include "magrathea/magrathea.h"
struct gapDisk: public density_base {
	double Sigma0;
	double rc;
	double h0;
	double P;//density index
	double S;//scale height index
	double p1, p2; //ring positions (au)
	double d1, d2; //ring depths (0 to 1)
	double w1, w2; //ring widths (au)

	gapDisk(): Sigma0(0), rc(0), h0(0), P(0), S(0), p1(0), p2(0), d1(0), d2(0), w1(0), w2(0){ }
	gapDisk(const gapDisk& other): Sigma0(other.Sigma0), rc(other.rc), h0(other.h0),
	P(other.P), S(other.S), p1(other.p1), p2(other.p2), d1(other.d1), d2(other.d2), 
	w1(other.w1), w2(other.w2) {  }
	gapDisk(double Sigma0, double rc, double h0, double P, double S, double p1, double p2,
	double d1, double d2, double w1, double w2):
	Sigma0(Sigma0), rc(rc), h0(h0), P(P), S(S), p1(p1), p2(p2), d1(d1), d2(d2),
	w1(w1), w2(w2){	}

	gapDisk& operator= (const gapDisk& other){
		Sigma0=other.Sigma0; rc=other.rc; h0=other.h0; P=other.P; S=other.S;
		p1=other.p1; p2=other.p2; d1=other.d1; d2=other.d2; w1=other.w1; w2=other.w2;
		return *this;
	}

	double surfaceMassDensity(const double r) const{
		return (Sigma0*pow(r/rc, -P)) * exp(-(pow(r/rc, 2-P)));
	}

	double scaleHeight(const double r) const{
		return (h0*pow(r/rc,S));
	}

	double operator()(double r, double theta, double phi) const{
		double r_cyl=r*sin(theta);
		double z=r*cos(theta);
		double h=scaleHeight(r_cyl);
		double gap1=(1-d1*(gaussianNotNorm(r_cyl,p1,w1)));
		double gap2=(1-d2*(gaussianNotNorm(r_cyl,p2,w2)));

		return(((surfaceMassDensity(r_cyl))/(sqrt(2*pi)*h)) * exp(-z*z/(2*h*h))*gap1*gap2);
	}
};

int main(int argc, char* argv[]){
	const size_t nThreadsMax=std::thread::hardware_concurrency();
	ThreadPool pool(nThreadsMax);
	astroParams diskData(246.6,-24.72,3.18e18*100,pi/4,pi/4);
	std::shared_ptr<gapDisk> dptr=std::make_shared<gapDisk>(1,100*AU,5*AU,0.5,1.25,25*AU,50*AU,0.8,0.9,5*AU,10*AU);
	std::ifstream gridfile("data/amr_grid.inp");
	std::ifstream datafile("data/dust_temperature_phi0.ascii");
	std::shared_ptr<fileTemp> tptr=std::make_shared<fileTemp>(gridfile,datafile);
	std::shared_ptr<dustOpacity> doptr=std::make_shared<dustOpacity>("data/dustopac.txt");
	
	grid g(0.1*AU,150*AU,75*pi/180,105*pi/180, 2.18*mSun, dptr, tptr, doptr);
	image img(500, 500,250*AU, 250*AU, {2.30538e11}, diskData);
	vect offset(0,0,0);
	
	img.propagate(g, grid::continuum, offset, pool, false, 100, 5);
	img.printToFits("GapDisk.fits");
	
	return 0;
}