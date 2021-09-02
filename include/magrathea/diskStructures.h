#ifndef DISKSTRUCTURES_H
#define DISKSTRUCTURES_H

#include "diskPhysics.h"

//the interpSurface is an r/theta/phi grid used for temperatures or desnities
//that come in a table format
//The main use is temperatures that come from RADMC output
//Currently, I've written this to assume axial symmetry, so the setup and interp2d
//functions don't consider phi
//For a very complicated 3d surface, this can be expanded
class interpSurface {
private:
	struct interpCell{
		double r;
		double theta;
		double phi;
		double value; //this will be either the temperature or the density
		interpCell(double r, double theta, double phi, double value): r(r), theta(theta), phi(phi), value(value) {	}
		interpCell(const interpCell& other): r(other.r), theta(other.theta), phi(other.phi), value(other.value)	{	}
	};
	
	std::vector<interpCell> data;
	std::vector<double> rCenters,tCenters;
	void initCenterLists(){
		for(size_t i=0; i<nr; i++){
			rCenters.push_back(data[i].r);
		}
		for(size_t i=0; i<data.size(); i+=nr){
			tCenters.push_back(data[i].theta);
		}
	}
	double rmin, rmax, tmin, tmax, pmin, pmax;
	int ndim;
	int nr; //used for 2d interpolation. Stride length.
	
public:
	interpSurface():data(std::vector<interpCell>()), rmin(0), rmax(0), tmin(0), tmax(0), pmin(0), pmax(0), ndim(0), nr(0)	{
		initCenterLists();
	}
	interpSurface(std::vector<interpCell> data, double rmin, double rmax, double tmin, double tmax, double pmin, double pmax , int ndim, int nr):
	data(data), rmin(rmin), rmax(rmax), tmin(tmin), tmax(tmax), pmin(pmin), pmax(pmax), ndim(ndim), nr(nr)	{
		initCenterLists();
	}
	interpSurface(const interpSurface& other):
	data(other.data), rmin(other.rmin), rmax(other.rmax), tmin(other.tmin), tmax(other.tmax), pmin(other.pmin), pmax(other.pmax), ndim(other.ndim), nr(other.nr)	{
		rCenters.clear();
		tCenters.clear();
		initCenterLists();
	}
	
	interpSurface& operator= (const interpSurface& other){
		data=other.data;
		rmin=other.rmin; rmax=other.rmax;
		tmin=other.tmin; tmax=other.tmax;
		pmin=other.pmin; pmax=other.pmax;
		ndim=other.ndim; nr=other.nr;
		rCenters.clear();
		tCenters.clear();
		initCenterLists();
		return *this;
	}
	
	//special constructor for reading in 2d radmc output files.
	interpSurface(std::istream& gridFile, std::istream& dataFile){
		unsigned int rtotal, ttotal, ptotal;
		std::string line;
		if(!gridFile){
			std::cout << "Could not open grid file " << std::endl;
			exit(1);
		}
		if(!dataFile){
			std::cout << "Could not open data file " << std::endl;
			exit(1);
		}
		std::stringstream stream;
		int linenum=0;
	
		for (linenum=0; linenum<5; linenum++) {
			getline(gridFile, line);
		}
	
		getline(gridFile, line);
		stream << line;
		stream >> rtotal >> ttotal >> ptotal;
		nr=rtotal; ndim=2;
		//std::cout << rtotal << "\t" << ttotal << "\t" << ptotal << std::endl;
		std::vector<double> thetaAndr(rtotal+ttotal+ptotal+3);
		double dummy;
		for(int i=0;i<rtotal+ttotal+ptotal+3;i++){
			getline(gridFile, line);
			//std::cout << line << std::endl;
			thetaAndr[i] = strtod(line.c_str(), NULL);
			//std::cout << thetaAndr[i] << std::endl;
		}	
		rmin=thetaAndr[0];
		rmax=thetaAndr[rtotal];
		tmin=thetaAndr[rtotal+1];
		tmax=thetaAndr[rtotal+ttotal+1];
		pmin=thetaAndr[rtotal+ttotal+2];
		pmax=thetaAndr[rtotal+ttotal+ptotal+2];
		//std::cout << rmin << "\t" << rmax << std::endl;
		//std::cout << tmin << "\t" << tmax << std::endl;
		//std::cout << pmin << "\t" << pmax << std::endl;
	
		double rLowerEdge, rUpperEdge, tLowerEdge, tUpperEdge;
		double rCenter, tCenter;
			
		tLowerEdge=thetaAndr[rtotal+1];
		for(int i=1;i<=ttotal;i++){
			tUpperEdge=thetaAndr[rtotal+1+i];
			tCenter=(tLowerEdge+tUpperEdge)/2;
			rLowerEdge=thetaAndr[0];
			for(int j=1;j<=rtotal;j++){
				rUpperEdge=thetaAndr[j];
				rCenter=(rLowerEdge+rUpperEdge)/2;
				data.push_back(interpCell(rCenter,tCenter,0,0));
				rLowerEdge=rUpperEdge;
			}
			tLowerEdge=tUpperEdge;
		}

		//getline(infile, line);
		//getline(infile, line);
		//unsigned int numcells = strtod(line.c_str(), NULL);
		//std::cout << data.size() << std::endl;
		//assert(numcells == data.size());
		//getline(infile, line);
		double temp;
		for(int i=0; i<data.size();i++){
			std::stringstream newstream;
			getline(dataFile, line);
			if(line != ""){
				newstream << line;
				//std::cout << line << std::endl;
				newstream >> dummy >> dummy >> dummy >> temp;
				data[i].value=temp;
				//std::cout << temp << std::endl;
			} else {
				i--;
			}
			//std::cout << i << " ";
		}
		initCenterLists();
	}
	
	double interp2d(double r, double theta) const{
	//note that theta is measured from the pole, not the midplane
		if(data.size()==0)
			return 0;
		if(theta>pi/2){
			theta = pi - theta;	//the file only goes to pi/2 because it is midplane symmetric
		}
		double rmin, rmax, tmin, tmax;
		//double tol=1e-6;
		rmin=data.front().r;
		tmin=data.front().theta;
		rmax=data.back().r;
		tmax=data.back().theta;
		/*if(r<surf.rmin-tol || r>surf.rmax+tol || theta<surf.tmin-tol || theta>surf.tmin*2+tol){
			std::cout << "Error: point (" << r << "," << theta << ") not inside region" << std::endl;
			std::cout << "rmin: " << surf.rmin << ", rmax: " << surf.rmax << ", tmin: " << surf.tmin << ", tmax " << surf.tmin*2 << std::endl;
			//std::cout << "(" << x << "," << y << "," << z << ")" << std::endl;
			exit(1);
		}*/
		//the points are ordered by theta then by r.
		//we need the nearest 4 points
		int indexT=0; int indexR=0;
		indexT=(std::lower_bound(tCenters.begin(),tCenters.end(),theta)-tCenters.begin())*nr;
		indexR=std::lower_bound(rCenters.begin(),rCenters.end(),r)-rCenters.begin();
		//because of the way we assign values to the centers of bins, our interpolation
		//surface is actually smaller than the grid, so we need a limited ability to 
		//extrapolate, provided the region is still inside the grid.
		if(indexT==0)
			indexT+=nr;
		if(indexT>=data.size())
			indexT-=nr;
		if(indexR==0)
			indexR++;
		if(indexR==nr)
			indexR--;
		//std::cout << indexT << '\t' << indexR << std::endl;
		//std::cout << indexR << "\t" << indexT << std::endl;
		double r0,r1,t0,t1;
		interpCell c1=data[indexT-nr+indexR-1];
		interpCell c2=data[indexT-nr+indexR];
		interpCell c3=data[indexT+indexR-1];
		interpCell c4=data[indexT+indexR];
		//std::cout << c1.value << "\t" << c2.value << "\t" << c3.value << "\t" << c4.value << std::endl;
	
		r0=c1.r;
		r1=c2.r;
		t0=c1.theta;
		t1=c3.theta;
		//std::cout << r0 << '\t' << r1 << '\t' << t0 << '\t' << t1 << std::endl;
		double n1= (r1-r)*(theta-t0)/((r1-r0)*(t1-t0));
		double n2= (r-r0)*(theta-t0)/((r1-r0)*(t1-t0));
		double n3= (r1-r)*(t1-theta)/((r1-r0)*(t1-t0));
		double n4= (r-r0)*(t1-theta)/((r1-r0)*(t1-t0));
		return c3.value*n1+c4.value*n2+c1.value*n3+c2.value*n4;	
	}
};

//the density_base provides the abstract template for user defined density classes
//The primary requirement is dervided classes provide an overloaded call operator
//which takes an r, theta, and phi coordinate and returns the density as a double.
//There is also a virtual columnDensityAbove function, which is available for 
//calculation of photodisassociation
class density_base{
public:
	virtual double operator()(double r, double theta, double phi)const{
		return 0;
	}
	//the main use of the scaleHeight function is checking if we should contract the stepsize.
	//If the user makes their own density that doesn't provide this, we just get this default
	//version, whichn returns 1cm. We can't return 0, since we divide by this value, but at 1cm,
	//it won't do the contraction.
	virtual double scaleHeight(double r)const{
		return 1;
	}
	//provides the column density above a given point. This is primarily useful for 
	//calculating photodisassociation.
	virtual double colDensAbove(double r, double theta, double phi)const{
		return 0;
	}
};

struct powerLawDisk : public density_base {
private:
	double Sigma0;
	double rc;
	double h0;
	double P;//density index
	double S;//scale height index
	
public:
	powerLawDisk(): Sigma0(0), rc(0), h0(0), P(0), S(0){ }
	powerLawDisk(const powerLawDisk& other): Sigma0(other.Sigma0), rc(other.rc), h0(other.h0),
	P(other.P), S(other.S){  }
	powerLawDisk(double Sigma0, double rc, double h0, double P, double S):
	Sigma0(Sigma0), rc(rc), h0(h0), P(P), S(S){	}

	powerLawDisk& operator= (const powerLawDisk& other){
		Sigma0=other.Sigma0; rc=other.rc; h0=other.h0; P=other.P; S=other.S;
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

		return(((surfaceMassDensity(r_cyl))/(sqrt(2*pi)*h)) * exp(-z*z/(2*h*h)));
	}
	
	//for a given point in the disk, return the column density ABOVE that position
	//needed to calculate the CO photodisociation
	//Here, the z that matters is the absolute value, since we want this symmetric
	//about the midplane;
	double colDensAbove(double r, double theta, double phi) const{
		double r_cyl=r*sin(theta);
		double z=abs(r*cos(theta));
		double result=surfaceMassDensity(r_cyl)/2;
		double arg=z/(sqrt(2)*scaleHeight(r_cyl));
		result *= erfc(arg);
		return result;
	}
};

//Because far fewer assumptions can be made of the temperature structure, there can't 
//necessarily be analogues of the scale height, etc, and the interface is much simpler.
class temperature_base{
public:
	virtual double operator()(double r, double theta, double phi)const{
		return 0;
	}
};

//the basic temp is primarily for testing.
//It's an r^-0.5 radial profile, and vertically isothermal.
//Too simple for a real project, but a close enough approximation for simple tests
struct basicTemp: public temperature_base{
	double r0, t0, index;
	basicTemp():r0(0),t0(0),index(0){}
	basicTemp(const basicTemp& other):r0(other.r0),t0(other.t0),index(other.index){}
	basicTemp(double r0, double t0, double index):r0(r0),t0(t0),index(index){}
	
	double operator()(double r, double theta, double phi) const{

		return(t0*pow(r/r0,index));
	}
};

//The fileTemp is designed to work with a temperature from an external table, 
//in this case one from RADMC. You need to supply both the RADMC grid file 
//as well as the temperature file.
struct fileTemp:public temperature_base{
	interpSurface tempSurf;
	fileTemp():tempSurf(){}
	fileTemp(const fileTemp& other):tempSurf(other.tempSurf){}
	
	fileTemp(std::ifstream& gridFile, std::ifstream& dataFile):tempSurf(gridFile, dataFile) {}
	
	double operator()(double r, double theta, double phi) const{
		return(tempSurf.interp2d(r,theta));
	}
	
};

#endif