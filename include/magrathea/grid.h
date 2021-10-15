/*
 *  grid.h
 *  Created by Erik on 11/9/14.
 *  
 *  This file defines the procedures for ray tracing through a grid,
 *  As well as basic grid structure
 *  
 *  
 *  Note that for the spherical coordinates, r is radial component,
 *  theta is the polar angle, and phi is the azimuthal angle.
 *  
 */
#ifndef GRID_H
#define GRID_H

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <sstream>

#include "diskPhysics.h"
#include "diskStructures.h"
#include "vectorclass.h"
#include "vectormath_trig.h"
#include "aligned_alloc.h"


class grid{
public:
	double rmin, rmax, tmin, tmax, pmin, pmax;	//starting/ending values for each dimension
	double starMass;
	std::shared_ptr<density_base> dens;
	std::shared_ptr<density_base> dustDens;//used when the dust and gas structures don't match
	std::shared_ptr<temperature_base> diskTemp;
	std::shared_ptr<opacity_base> dustOpac;
	std::shared_ptr<opacity_base> lineOpac;
	bool freezeout;
	double turbulence;
	
	enum prop_type {normal, continuum, attenuated_continuum, subtracted_bad, subtracted_real};
	
	//empty
	grid(){
		rmin=0; rmax=0; tmin=0; tmax=0; pmin=0; pmax=0;
		dens=NULL; dustDens=NULL; diskTemp=NULL;
		freezeout=false;
		turbulence=0;
	}
	
	//standard, including a dust and gas structure
	grid(double rmin,double rmax,double tmin,double tmax, double starMass, 
		std::shared_ptr<density_base> dens, std::shared_ptr<density_base> dustDens, 
		std::shared_ptr<temperature_base> diskTemp, std::shared_ptr<opacity_base> dustOpac,
		std::shared_ptr<opacity_base> lineOpac, bool freezeout, double turb):
		rmin(rmin), rmax(rmax), tmin(tmin), tmax(tmax), pmin(0),pmax(2*pi),
		starMass(starMass), dens(dens), dustDens(dustDens), diskTemp(diskTemp), 
		dustOpac(dustOpac), lineOpac(lineOpac), freezeout(freezeout), turbulence(turb){	}
	
	//simple continuum-only constructor
	grid(double rmin, double rmax, double tmin, double tmax, double starMass,
		std::shared_ptr<density_base> dustDens, std::shared_ptr<temperature_base> diskTemp,
		std::shared_ptr<opacity_base> dustOpac):
		rmin(rmin), rmax(rmax), tmin(tmin), tmax(tmax), pmin(0),pmax(2*pi),
		starMass(starMass), dens(std::make_shared<density_base>()), dustDens(dustDens), diskTemp(diskTemp), 
		dustOpac(dustOpac), lineOpac(std::make_shared<opacity_base>()), freezeout(false), turbulence(0){	}


	//copy
	//Currently real broken
	grid(const grid& other):rmin(other.rmin), rmax(other.rmax),tmin(other.tmin),tmax(other.tmax),pmin(other.pmin),pmax(other.pmax), 
	dens(other.dens), dustDens(other.dustDens),	dustOpac(other.dustOpac), lineOpac(other.lineOpac), 
	freezeout(other.freezeout), turbulence(other.turbulence){	}
	
	//assignment
	grid& operator= (const grid& other){
		if(this!=&other){
			rmin=other.rmin; rmax=other.rmax; tmin=other.tmin;
			tmax=other.tmax; pmin=other.pmin; pmax=other.pmax;
			dens=other.dens, dustDens=other.dustDens, 
			diskTemp=other.diskTemp, 
			dustOpac=other.dustOpac; lineOpac=other.lineOpac;
			freezeout=other.freezeout, turbulence=other.turbulence;
		}
		return *this;
	}
	
	//the widthInfo struct is used for calculating the width of the emission region, among other things
	//a widthInfo is stored for each step through the disk. Calculations are done on the set at the end.
	struct widthInfo{
		vect pos;
		double intensity;
		double tau;
		widthInfo():pos(0,0,0),intensity(0),tau(0){	}
		widthInfo(vect pos, double intensity, double tau): pos(pos), intensity(intensity), tau(tau) { }
		widthInfo(const widthInfo& other):pos(other.pos),intensity(other.intensity),tau(other.tau) {	}
	};
	
	std::vector<std::pair<intersection,intersection> > entriesAndExits(const line& l) const{
		const vect origin(0,0,0);
		std::vector<std::pair<intersection,intersection> > diskPairs;
		
		//make the geometric objects for the boundaries
		sphere s1(origin,rmin);
		sphere s2(origin,rmax);
		cone c1(origin,vect(0,0,1),tmin);
		cone c2(origin,vect(0,0,-1),pi-tmax);
		shortIntersectionList<8> intersects; //store all intersections here, there should be 2 or 4, but I'm being paranoid.
		//std::cout << "rmin: " << rmin << ", rmax: " << rmax << ",tmin: " << tmin << ",tmax: " << tmax << std::endl;
		
		//first get the sphere intersection points, and check to make
		//sure that they are within the correct theta and phi boundaries
		shortIntersectionList<2> temp; //use this for each stage of intersection checking
		temp=lineSphereIntersect(l, s1);	//inner sphere
		for(const auto &isect: temp){
			if(!((isect.location.theta()<tmin&&isect.location.theta()<tmax) ||
				 (isect.location.theta()>tmin&&isect.location.theta()>tmax))){
				vect v(isect.location-l.position);
				if(v*l.direction >=0){//make sure the point is past the start of the line
					intersects.push_back(isect);
					intersects.back().side=0;
					//std::cout << "Intersection with inner sphere at " << isect.location.r() << ", "
					//<< isect.location.theta() << ", " << isect.location.phi() << std::endl;
				}
			}
		}
		temp=lineSphereIntersect(l, s2);
		for(int i=0;i<temp.size();i++){
			if(!((temp[i].location.theta()<tmin&&temp[i].location.theta()<tmax)||
				 (temp[i].location.theta()>tmin&&temp[i].location.theta()>tmax))){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					intersects.back().side=1;
					//std::cout << "Intersection with outer sphere at " << temp[i].location.r() << ", "
					//<< temp[i].location.theta() << ", " <<temp[i].location.phi() << std::endl;
				}
			}
		}
		
		//next find the intersections between the line and the cones,
		temp=lineConeIntersect(l, c1);
		for(int i=0;i<temp.size();i++){
			if(!( (temp[i].location.r()<rmin && temp[i].location.r()<rmax) || (temp[i].location.r()>rmin && temp[i].location.r()>rmax))){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					intersects.back().side=0;
					//std::cout << "Intersection with upper cone at " << temp[i].location.r() << ", "
					//<< temp[i].location.theta() << ", " <<temp[i].location.phi() << std::endl;
				}
			}
		}
		temp=lineConeIntersect(l, c2);
		for(int i=0;i<temp.size();i++){
			if(!( (temp[i].location.r()<rmin && temp[i].location.r()<rmax) || (temp[i].location.r()>rmin && temp[i].location.r()>rmax))){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					intersects.back().side=1;
					//std::cout << "Intersection with lower cone at " << temp[i].location.r() << ", "
					//<< temp[i].location.theta() << ", " <<temp[i].location.phi() << std::endl;
				}
			}
		}
		
		//sort the intersections in order of distance to the start of the ray
		std::sort(intersects.begin(),intersects.end(),
				  [&](const intersection& i1, const intersection& i2)->bool{
					  return(((i1.location-l.position)*l.direction)>((i2.location-l.position)*l.direction));
				  });
		
		//if there are two intersections, there is just one pair
		if(intersects.size()==2){
			std::pair<intersection,intersection> testpair (intersects[0],intersects[1]);
			diskPairs.push_back(testpair);
		}
		//if there are three intersections, it means that the line just grazes the inner sphere
		//or one of the cones, and we can ignore the middle intersection
		if(intersects.size()==3){
			std::pair<intersection,intersection> testpair (intersects[0],intersects[2]);
			diskPairs.push_back(testpair);
		}
		//if there are four intersections, the line enters the cell, leaves again, reenters, and
		//then leaves. We only care about the path between the first and second, and third and
		//fourth intersections
		if(intersects.size()==4){
			std::pair<intersection,intersection> testpair1 (intersects[0],intersects[1]);
			std::pair<intersection,intersection> testpair2 (intersects[2],intersects[3]);
			diskPairs.push_back(testpair1);
			diskPairs.push_back(testpair2);
		}
		//Six intersections is rare but possible if the ray comes in almost parallel to a cone, but
		//doesn't go though the center. The cone then has some curvature, so the ray can enter through
		//the outer sphere, exit through the inner sphere, enter through the inner sphere, exit through
		//the cone, enter through the cone, and exit through the outer sphere again. This isn't a floating
		//point roundoff issue, it is just geometric, though very rare.
		if(intersects.size()==6){
			//std::cout << "6 intersections" << std::endl;
			std::pair<intersection,intersection> testpair1 (intersects[0],intersects[1]);
			std::pair<intersection,intersection> testpair2 (intersects[2],intersects[3]);
			std::pair<intersection,intersection> testpair3 (intersects[4],intersects[5]);
			diskPairs.push_back(testpair1);
			diskPairs.push_back(testpair2);
			diskPairs.push_back(testpair3);
		}
		//I don't think there's any way to have more than 6 intersections. If there are 5, it is likely a
		//degenerate version of the 6 case, but where is only grazes the cone or the inner sphere.
		//If that is true it can be handled, but only once I find a ray with 5 intersections so I can confirm this.
		if(intersects.size() == 5 || intersects.size() > 6){
			std::cout << "Something has gone horribly wrong in propagate() finding the disk intersections" << std::endl;
			std::cout << (int)intersects.size() << " intersections?!?" << std::endl;
			std::cout << "Ray start: " << l.position << ", ray direction: " << l.direction << std::endl;
			std::cout << "Intersection\t\ttype" << std::endl;
			for(auto isect : intersects){
				std::cout << isect.location << "\t\t" << isect.type << std::endl;
			}
			exit(1);
		}
		
		return(diskPairs);
	}
	
	//std::vector<std::vector<widthInfo> > propagateRay(const line l, const std::vector<double> frequencies, const vect cameraPosition, prop_type type) const{
	std::vector<double> propagateRay(const line l, const std::vector<double> frequencies, const vect cameraPosition, prop_type type, 
	unsigned int nsteps, double contractionFactor) const{
		std::vector<std::pair<intersection,intersection> > diskPairs=entriesAndExits(l);

		vect dir=l.direction/-l.direction.mag(); //ensure unit vector
		//std::cout << "pos:" << l.position << "\tdir: " << l.direction << std::endl;
		//std::cout << "dir: " << dir << std::endl;
		std::vector<double> value (frequencies.size(),0);	//final results will go here
		std::vector<double> opticalDepth (frequencies.size(),0); //for diagnostics
		std::vector<double> temperatures; //temps at each step through the disk
	
		for(auto pair : diskPairs){
			int stepcount=0;
			vect pos=pair.first.location; //starting point
			double remainingDist=(pair.second.location-pair.first.location).mag(); //the total distance along which to integrate
			double stepsize = remainingDist/nsteps;	//currently nonadaptive stepsize
			std::vector<double> valNew(frequencies.size(),0);
			pos+=dir*stepsize/2;
			bool isShrunk=false;
			
			//double numScaleHeights; 
			double scaleHeightLimit=2.5; //how many scale heights above/below the midplane we increase the stepsize for
			double numScaleHeights=(pos).z/dustDens->scaleHeight((pos).r()*sin((pos).theta()));
			while (remainingDist >= stepsize/2) {
				double numScaleHeightsNext=(pos+dir*stepsize).z/dustDens->scaleHeight((pos+dir*stepsize).r()*sin((pos+dir*stepsize).theta()));
				if(!isShrunk && ((numScaleHeightsNext > -scaleHeightLimit && numScaleHeights < -scaleHeightLimit) || (numScaleHeightsNext < scaleHeightLimit && numScaleHeights > scaleHeightLimit))){
					//std::cout << "shrinking step size" << std::endl;
					isShrunk=true;
					stepsize /= contractionFactor;
					numScaleHeightsNext=(pos+dir*stepsize).z/dustDens->scaleHeight((pos+dir*stepsize).r()*sin((pos+dir*stepsize).theta()));
				}
				if(isShrunk && ((numScaleHeightsNext > scaleHeightLimit && numScaleHeights < scaleHeightLimit) || (numScaleHeightsNext < -scaleHeightLimit && numScaleHeights > -scaleHeightLimit))){
					//std::cout << "expanding step size" << std::endl;
					isShrunk=false;
					stepsize *= contractionFactor;
					numScaleHeightsNext=(pos+dir*stepsize).z/dustDens->scaleHeight((pos+dir*stepsize).r()*sin((pos+dir*stepsize).theta()));
				}
				numScaleHeights=numScaleHeightsNext;
				//std::cout << stepcount << "\t" << pos.z/AU << "\t" << numScaleHeights << "\t" << numScaleHeightsNext << std::endl;
				//densities of gas and dust.
				//For CO, which we're using here, we add a factor to account for photodissociation 
				double d=(*dens)(pos.r(),pos.theta(),pos.phi());	//gas density
				double ddust=(*dustDens)(pos.r(),pos.theta(),pos.phi());	//seperate dust density
				
				//Previously, I just turned off the CO if the column (number) density of H2 
				//was less than 10^21. This is from C. Qi et al 2011, and Rosenfeld 2013
				//Now, it gets smoothly decreased over an order of magnitude around 10^21
				if(type!=continuum){
					double colNumDens=dens->colDensAbove(pos.r(),pos.theta(),pos.phi())/amu/2.0;
					//I'm further complicating this by smoothly phasing out the CO
					double factor=lineOpac->disassociationFactor(colNumDens);
					d*=factor;
				}
				double temp=(*diskTemp)(pos.r(),pos.theta(),pos.phi()); //starting temp
				temperatures.push_back(temp);
		
				double r_cyl=pos.r()*sin(pos.theta());	//cylindrical radius
				double azdif=(cameraPosition.phi()-pos.phi());//difference between the camera orientation and the point in question
				
				double velocity=sqrt(gravConst*starMass/r_cyl/r_cyl/r_cyl)*r_cyl*sin(cameraPosition.theta())*cos(-pi/2+azdif);//velocity with respect to LOS
				double soundspeed=sqrt(temp*kboltzmann/(3.819239518e-24));
			
				for(int i=0;i<frequencies.size();i++){
					double newfreq=doppler(frequencies[i],velocity);
					double dustOp=(*dustOpac)(newfreq);
					double BB=blackBody(temp,newfreq);
					double opacity;
					if(type==continuum){
						opacity=ddust*dustOp;	
						valNew[i]=(stepsize*opacity*BB+value[i]*(1-stepsize*opacity/2))/(1+stepsize*opacity/2);
						//std::cout << "dens: " << ddust << "\tTemp: " << temp << "\tOpac: " << dustOp <<	std::endl;	
					}else if(type==normal || type==subtracted_bad || type==subtracted_real){
						opacity=d*(*lineOpac)(temp,newfreq,turbulence*soundspeed,freezeout)+ddust*dustOp;
						valNew[i]=(stepsize*opacity*BB+value[i]*(1-stepsize*opacity/2))/(1+stepsize*opacity/2);
						//std::cout << stepcount << "\t" << pos.x/AU << "\t" << pos.y/AU << "\t" << pos.z/AU << "\tdens: " << ddust << "\tTemp: " << temp << "\tOpac: " << opacity <<	std::endl;
					}else if(type==attenuated_continuum){
						double opacityDust=ddust*dustOp;
						double opacityBoth=d*(*lineOpac)(temp,newfreq,turbulence*soundspeed,freezeout)+ddust*dustOp;
						opacity=opacityDust; //for calculating the optical depth
						valNew[i]=(stepsize*opacityDust*BB+value[i]*(1-stepsize*opacityBoth/2))/(1+stepsize*opacityBoth/2);
					}
					if(type==subtracted_bad || type==subtracted_real){
						opacity=d*(*lineOpac)(temp,newfreq,turbulence*soundspeed,freezeout);//in the subtracted cases, the optical depth is just that of the line
					}
					//std::cout << value[i] << "\t" << valueSlow[i] << "\t" << opticalDepth[i] << "\t" << stepsize*opacity*BB << "\t" << exp(-opticalDepth[i]) << std::endl;
					//valueSlow[i]+=stepsize*opacity*BB * (exp(-opticalDepth[i]));
					opticalDepth[i]+=(stepsize*opacity);
					
					//std::cout << pos/AU << "\t" << opacity << "\t" << opticalDepth[i] << std::endl;
				}
				value=valNew;
				pos+=dir*stepsize;
				remainingDist-=stepsize;
				stepcount++;
			}
		}
	
		if(type==subtracted_bad || type==subtracted_real){
			std::vector<double> valueSpecial(frequencies.size(),0);	//final results will go here
			for(auto pair : diskPairs){
				int stepcount=0;
				vect pos=pair.first.location; //starting point
				double remainingDist=(pair.second.location-pair.first.location).mag(); //the total distance along which to integrate
				double stepsize = remainingDist/nsteps;	//currently nonadaptive stepsize
				std::vector<double> valNew(frequencies.size(),0);
				pos+=dir*stepsize/2;
				while (remainingDist >= stepsize/2) {
					double d=(*dens)(pos.r(),pos.theta(),pos.phi());	//gas density
					double ddust=(*dustDens)(pos.r(),pos.theta(),pos.phi());	//seperate dust density
					double temp=(*diskTemp)(pos.r(),pos.theta(),pos.phi());	//starting temp
					double r_cyl=pos.r()*sin(pos.theta());	//cylindrical radius
					double azdif=(cameraPosition.phi()-pos.phi());//difference between the camera orientation and the point in question
					double velocity=sqrt(gravConst*starMass/r_cyl/r_cyl/r_cyl)*r_cyl*sin(cameraPosition.theta())*cos(-pi/2+azdif);//velocity with respect to LOS
					double soundspeed=sqrt(temp*kboltzmann/(3.819239518e-24));
					for(int i=0;i<frequencies.size();i++){
						double newfreq=doppler(frequencies[i],velocity);
						double dustOp=(*dustOpac)(newfreq);
						double BB=blackBody(temp,newfreq);
						if(type==subtracted_bad){
							double opacity= ddust*dustOp;					
							valNew[i]=(stepsize*opacity*BB+valueSpecial[i]*(1-stepsize*opacity/2))/(1+stepsize*opacity/2);
						}else if(type==subtracted_real){
							double opacityDust=ddust*dustOp;
							double opacityBoth=d*(*lineOpac)(temp,newfreq,turbulence*soundspeed,freezeout)+ddust*dustOp;
							valNew[i]=(stepsize*opacityDust*BB+valueSpecial[i]*(1-stepsize*opacityBoth/2))/(1+stepsize*opacityBoth/2);
						}
					}
					valueSpecial=valNew;
					pos+=dir*stepsize;
					remainingDist-=stepsize;
					stepcount++;
				}
			}
			for(int i=0;i<value.size();i++){
				value[i]-=valueSpecial[i];
			}
		}
		
		//optional conversion to brightness temperature
		/*assert(value.size()==frequencies.size());
		for(int f=0;f<frequencies.size();f++){
			value[f]=tempFromSurfaceBrightness(value[f], frequencies[f]);
		}*/
		
		return(value);
	}
	
	//same as the basic progagateRay, but with frequencies packed into sets of 4, so we can do AVX operations
	//for this, I'm taking out the special versions of prop_type, since I probably will never need to do
	//the incorrect continuum subtraction again.
	std::vector<double> propagateRayAVX(const line l, const std::vector<Vec4d,aligned_allocator<Vec4d> > freqsAVX, const vect cameraPosition, prop_type type, 
	unsigned int nsteps, double contractionFactor) const{
		if(type==attenuated_continuum || type==subtracted_bad || type==subtracted_real){
			std::cout << "Error: AVX ray tracing does not support propagation of type " << type << std::endl;
			exit(1);
		}
		
		std::vector<std::pair<intersection,intersection> > diskPairs=entriesAndExits(l);

		vect dir=l.direction/-l.direction.mag(); //ensure unit vector
		//std::cout << "pos:" << l.position << "\tdir: " << l.direction << std::endl;
		//std::cout << "dir: " << dir << std::endl;
		std::vector<Vec4d,aligned_allocator<Vec4d> > value (freqsAVX.size(),0,aligned_allocator<Vec4d>(5));	//final results will go here
		std::vector<Vec4d,aligned_allocator<Vec4d> > opticalDepth (freqsAVX.size(),0,aligned_allocator<Vec4d>(5)); //for diagnostics
		std::vector<double> temperatures; //temps at each step through the disk
		
		for(auto pair : diskPairs){
			int stepcount=0;
			vect pos=pair.first.location; //starting point
			double remainingDist=(pair.second.location-pair.first.location).mag(); //the total distance along which to integrate
			double stepsize = remainingDist/nsteps;	//currently nonadaptive stepsize
			std::vector<Vec4d,aligned_allocator<Vec4d> > valNew(freqsAVX.size(),0,aligned_allocator<Vec4d>(5));
			pos+=dir*stepsize/2;
			bool isShrunk=false;
			
			//double numScaleHeights; 
			double scaleHeightLimit=2.5; //how many scale heights above/below the midplane we increase the stepsize for
			double numScaleHeights=(pos).z/dustDens->scaleHeight((pos).r()*sin((pos).theta()));
			while (remainingDist >= stepsize/2) {
				double numScaleHeightsNext=(pos+dir*stepsize).z/dustDens->scaleHeight((pos+dir*stepsize).r()*sin((pos+dir*stepsize).theta()));
				if(!isShrunk && ((numScaleHeightsNext > -scaleHeightLimit && numScaleHeights < -scaleHeightLimit) || (numScaleHeightsNext < scaleHeightLimit && numScaleHeights > scaleHeightLimit))){
					//std::cout << "shrinking step size" << std::endl;
					isShrunk=true;
					stepsize /= contractionFactor;
					numScaleHeightsNext=(pos+dir*stepsize).z/dustDens->scaleHeight((pos+dir*stepsize).r()*sin((pos+dir*stepsize).theta()));
				}
				if(isShrunk && ((numScaleHeightsNext > scaleHeightLimit && numScaleHeights < scaleHeightLimit) || (numScaleHeightsNext < -scaleHeightLimit && numScaleHeights > -scaleHeightLimit))){
					//std::cout << "expanding step size" << std::endl;
					isShrunk=false;
					stepsize *= contractionFactor;
					numScaleHeightsNext=(pos+dir*stepsize).z/dustDens->scaleHeight((pos+dir*stepsize).r()*sin((pos+dir*stepsize).theta()));
				}
				numScaleHeights=numScaleHeightsNext;
				//std::cout << stepcount << "\t" << pos.z/AU << "\t" << numScaleHeights << "\t" << numScaleHeightsNext << std::endl;
				
				double d=(*dens)(pos.r(),pos.theta(),pos.phi());	//gas density
				double ddust=(*dustDens)(pos.r(),pos.theta(),pos.phi());	//seperate dust density
				//Previously, I just turned off the CO if the column (number) density of H2 
				//was less than 10^21. This is from C. Qi et al 2011, and Rosenfeld 2013
				//Now, it gets smoothly decreased over an order of magnitude around 10^21
				if(type!=continuum){
					double colNumDens=dens->colDensAbove(pos.r(),pos.theta(),pos.phi())/amu/2.0;
					//I'm further complicating this by smoothly phasing out the CO
					double factor=lineOpac->disassociationFactor(colNumDens);
					d*=factor;
				}
				
				double temp=(*diskTemp)(pos.r(),pos.theta(),pos.phi()); //starting temp
				temperatures.push_back(temp);
		
				double r_cyl=pos.r()*sin(pos.theta());	//cylindrical radius
				double azdif=(cameraPosition.phi()-pos.phi());//difference between the camera orientation and the point in question
				
				double velocity=sqrt(gravConst*starMass/r_cyl/r_cyl/r_cyl)*r_cyl*sin(cameraPosition.theta())*cos(-pi/2+azdif);//velocity with respect to LOS
				double soundspeed=sqrt(temp*kboltzmann/(3.819239518e-24));
			
				for(int i=0;i<freqsAVX.size();i++){
					Vec4d newfreq=dopplerAVX(freqsAVX[i],velocity);
					Vec4d dustOp=(*dustOpac)(newfreq);
					Vec4d BB=blackBodyAVX(temp,newfreq);
					Vec4d opacity;
					if(type==continuum){
						opacity=ddust*dustOp;	
						valNew[i]=(stepsize*opacity*BB+value[i]*(1-stepsize*opacity/2))/(1+stepsize*opacity/2);
						//std::cout << "dens: " << ddust << "\tTemp: " << temp << "\tOpac: " << dustOp <<	std::endl;	
					}else if(type==normal){
						opacity=d*(*lineOpac)(temp,newfreq,turbulence*soundspeed,freezeout)+ddust*dustOp;
						valNew[i]=(stepsize*opacity*BB+value[i]*(1-stepsize*opacity/2))/(1+stepsize*opacity/2);
						//std::cout << "dens: " << ddust << "\tTemp: " << temp << "\tOpac: " << opacity <<	std::endl;
					//std::cout << value[i] << "\t" << valueSlow[i] << "\t" << opticalDepth[i] << "\t" << stepsize*opacity*BB << "\t" << exp(-opticalDepth[i]) << std::endl;
					//valueSlow[i]+=stepsize*opacity*BB * (exp(-opticalDepth[i]));
					opticalDepth[i]+=(stepsize*opacity);
					}
					//std::cout << pos/AU << "\t" << opacity << "\t" << opticalDepth[i] << std::endl;
				}
				value=valNew;
				pos+=dir*stepsize;
				remainingDist-=stepsize;
				stepcount++;
			}
		}
		
		//optional conversion to brightness temperature
		/*assert(value.size()==frequencies.size());
		for(int f=0;f<frequencies.size();f++){
			value[f]=tempFromSurfaceBrightness(value[f], frequencies[f]);
		}*/
		
		//as a last step, we need to unpack the Vec4ds back into just a vector of doubles.
		std::vector<double> finalValues;
		for(int i=0;i<freqsAVX.size()-1;i++){
			finalValues.push_back(value[i][0]);
			finalValues.push_back(value[i][1]);
			finalValues.push_back(value[i][2]);
			finalValues.push_back(value[i][3]);
		}
		//and don't forget the last couple
		for(int i=0;i<4;i++){
			if(freqsAVX.back()[i]!=-1.0)//don't add the padding ones
				finalValues.push_back(value.back()[i]);
		}
		
		return(finalValues);
	}
	
	/*interpSurface makeVelInterpSurf(int nr, int ntheta){
		if(ntheta%2!=0)
			ntheta++;
		double thetastride=(pi/2-tmin)/ntheta;
		double rstride=(rmax-rmin)/nr;
		std::vector<interpCell> data;
		for(int i=0; i<ntheta; i++){
			double tCenter=tmin+thetastride*(i+0.5);
			for(int j=0; j<nr; j++){
				double rCenter=rmin+rstride*(j+0.5);
				double r_cyl=rCenter*sin(tCenter);
				double z=rCenter*cos(tCenter);
				integrand intgrnd{r_cyl,z, dens};
				double vel=velocity([&](double rp){return intgrnd.evaluate(rp);},r_cyl,z);
				data.push_back(interpCell(rCenter,tCenter,0,vel));
				//std::cout << r_cyl/1.496e13 << "\t" << z/1.496e13 << "\t" << sigma0FromMass() << "\t" << sqrt(vel) << std::endl;
			}
		}
		interpSurface surf(data,rmin,rmax,tmin,tmax,pmin,pmax,2, nr);
		return surf;
	}*/
	
	/*struct integrand{
		double r, z;
		densityPert& dens;
		
		double evaluate(double rprime){
			gsl_mode_t mode = GSL_PREC_DOUBLE;
			
			double sigma=dens.surfaceMassDensity(rprime);
			double k= sqrt(4*r*rprime / ( (r+rprime)*(r+rprime) + z*z));
			double KK, EE;
			if(k>0.7){
				KK = gsl_sf_ellint_Kcomp (k, mode);
				EE = gsl_sf_ellint_Ecomp (k, mode);
			}else {
				double point= k*k;
				KK=1;
				EE=1;
				double coef=0.25;
				for(int i=1;i<12;i++){
					KK+=(coef*pow(point,i));
					EE-=(coef*pow(point,i)/(2*i-1));
					coef *= (2*i+1)*(2*i+1);
					coef /= ((2*i+2)*(2*i+2));
				}
				KK *= pi/2; EE *= pi/2;
			}
			double part1=(1.0/4.0)*(k*k/(1.0-k*k))*(rprime/r - r/rprime + z*z/(r*rprime))*EE;
			//std::cout << k << "\t" << rprime << "\tSigma: " << sigma << "\t" << ((KK-part1)*sqrt(rprime/r)*k*sigma) << std::endl;
			return ((KK-part1)*sqrt(rprime/r)*k*sigma);
		}
	};*/
	
	//velocity for generalized function. This can be used with the integrand class above
	//to find the velocity with self gravity added. Note that in my testing I never found the
	//self gravity contribution to be significant or even anything close to what we could 
	//actually observe.
	/*template<typename FunctionType>
	double velocity(FunctionType f, double r, double z){
		double (*wrapper)(double,void*)=[](double rp, void* params){
			FunctionType& f=*static_cast<FunctionType*>(params);
			return(f(rp));
		};
		
		gsl_function F;
		F.function = wrapper;
		F.params   = &f;
		
		gsl_integration_workspace* w = gsl_integration_workspace_alloc (5000);
		double result, error;
		
		//double pts[] = {rmin, r, rmax};
		//std::cout << start << "\t" << r << "\t" << end << std::endl;
		
		//std::ofstream out("grav.txt");
		//double rprime=rmin;
		//double stepsize= (rmax-rmin)/5000;
		
		gsl_integration_qag(&F, rmin, rmax, 0, 1.e-2, 5000, 1, w, &result, &error);
		result*=(gravConst/r/r);
		return result;
	}*/
	
	//the pressure of gas in the disk alters the speed at which the gas rotates
	//The amount isn't always noticeable, so this is optional.
	/*double presVel(double r, double theta, double phi, double turbulenceParameter) const{
		double mmw=3.819239518e-24; //mean molecular weight: 2.3 amu (in grams)
		double density=dens(r,theta, phi);
		double temp=surf.interp2d(r, theta);
		double soundSpeed=sqrt(kboltzmann*temp/mmw);
		
		double v_turb = turbulenceParameter*soundSpeed;
		
		double p1_th = density*kboltzmann*temp/mmw;
		double p1_tu = density*v_turb*v_turb;
		double p1 = p1_th + p1_tu;
		
		double r2 = r*(1+1.e-7);
		double density2=dens(r2,theta, phi);
		double temp2=surf.interp2d(r2, theta);
		double p2_th = density2*kboltzmann*temp2/mmw;
		
		double p2_tu = density2*v_turb*v_turb;
		double p2 = p2_th + p2_tu;
		
		double dPsudr_t = (p2-p1)/(r2-r);
		
		dPsudr_t/=(r*density);
		return dPsudr_t;
	}*/
	
	//this is for the case where the disk is axisymmetric and imaged face on, so r is the only variable
	void radialImage(const std::string outFileName, const double centFreq, const double freqRange, const int freqBins, 
	unsigned int nsteps, prop_type type, double contractionFactor){
		std::vector<double> frequencies;
		double freq=centFreq;
		freq-=freqRange/2;
		for(int i=0;i<freqBins;i++){
			freq+=(freqRange/freqBins);
			frequencies.push_back(freq);
			//std::cout << freq << std::endl;
		}
		double radialDist=rmax-rmin;
		double radialStep=radialDist/nsteps;
		vect down(0,0,-1);
		vect pos(rmin+(radialStep/2),0,1e18); //start at rmin, and way above the disk
		std::vector<std::vector<double> > results;
		
		for(int i=0;i<nsteps;i++){
			line path(pos,down);
			results.push_back(propagateRay(path, frequencies, pos, type, nsteps, contractionFactor));
			pos.x+=radialStep;
		}
		//print to file if there is only one frequency.
		if(frequencies.size()==1){
			double rpos=rmin+(radialStep/2);
			std::ofstream outfile(outFileName.c_str());
			for(int i=0;i<nsteps;i++){
				outfile << rpos/AU << "\t" << results[i].front() << std::endl;
				rpos+=radialStep;
			}
			outfile.close();
		}
		
		//otherwise, calculate moment 0 and 8
		if(frequencies.size()>1){
			std::vector<double> moment0;
			std::vector<double> moment8;
			
			for(int i=0;i<nsteps;i++){
				double mom0=0;
				double mom8=0;
				for(int j=0;j<frequencies.size();j++){
					//std::cout << results[i][j] << std::endl;
					mom0+=results[i][j];
					if(results[i][j] > mom8)
						mom8=results[i][j];
				}
				//std::cout << "sum: " << mom0 << std::endl << std::endl;
				mom0/=frequencies.size();
				moment0.push_back(mom0);
				moment8.push_back(mom8);
			}
	
			std::ofstream outfile(outFileName.c_str());
			double rpos=rmin+(radialStep/2);
			for(int i=0;i<nsteps;i++){
				double mom0=moment0[i];
				double mom8=moment8[i];
				double tempMom0=tempFromSurfaceBrightness(mom0, centFreq);
				double tempMom8=tempFromSurfaceBrightness(mom8, centFreq);
				outfile << rpos/AU << "\t" << moment0[i] << "\t" << moment8[i] << "\t" << tempMom0 << "\t" << tempMom8 << std::endl;
				rpos+=radialStep;
			}
			outfile.close();	
		}
	}

	//special propagate for characterizing the emission region.
	//this works only for a single frequency, and only for "normal" propagation, since we have to do things the slow way
	std::vector<widthInfo> propagateRaySpecial(const line l, const double frequency, const vect cameraPosition) const{
		std::vector<std::pair<intersection,intersection> > diskPairs=entriesAndExits(l);

		vect dir=l.direction/-l.direction.mag(); //ensure unit vector
		//std::cout << "pos:" << l.position << "\tdir: " << l.direction << std::endl;
		//std::cout << "dir: " << dir << std::endl;
		double opticalDepth=0;
		std::vector<widthInfo> allThree2; //for each step, a pair of vect,taus for each freq
		std::vector<double> temperatures; //temps at each step through the disk
		double valueSlow=0;
		unsigned int nsteps=500;
	
		for(auto pair : diskPairs){
			int stepcount=0;
			vect pos=pair.first.location; //starting point
			double remainingDist=(pair.second.location-pair.first.location).mag(); //the total distance along which to integrate
			double stepsize = remainingDist/nsteps;	//currently nonadaptive stepsize
			pos+=dir*stepsize/2;
		
			while (remainingDist >= stepsize/2) {
				double d=(*dens)(pos.r(),pos.theta(),pos.phi());	//gas density
				double ddust=(*dustDens)(pos.r(),pos.theta(),pos.phi());	//seperate dust density
				double temp=(*diskTemp)(pos.r(),pos.theta(),pos.phi());	//starting temp
				temperatures.push_back(temp);
		
				double r_cyl=pos.r()*sin(pos.theta());	//cylindrical radius
				double azdif=(cameraPosition.phi()-pos.phi());//difference between the camera orientation and the point in question
				
				double velocity=sqrt(gravConst*starMass/r_cyl/r_cyl/r_cyl)*r_cyl*sin(cameraPosition.theta())*cos(-pi/2+azdif);//velocity with respect to LOS
				double soundspeed=sqrt(temp*kboltzmann/(3.819239518e-24));
				widthInfo posIntTau;//position,intensity, and optical depth at s
			
				double newfreq=doppler(frequency,velocity);
				double dustOp=(*dustOpac)(newfreq);
				double BB=blackBody(temp,newfreq);
				double opacity;
				opacity=d*(*lineOpac)(temp,newfreq,turbulence*soundspeed,freezeout)+ddust*dustOp;
				valueSlow+=stepsize*opacity*BB * (exp(-opticalDepth));
				opticalDepth+=(stepsize*opacity);
					
				posIntTau=widthInfo(pos,valueSlow,opticalDepth);
				allThree2.push_back(posIntTau);
				pos+=dir*stepsize;
				remainingDist-=stepsize;
				stepcount++;
			}
		}
				
		/*std::ofstream out("junk.txt");
		for(auto i : allThree2){
			double realTemp=temperature(i.pos.r(),i.pos.theta(),i.pos.phi());
			double bTemp=tempFromSurfaceBrightness(i.intensity, frequency);
			out << i.pos.z/AU << "\t" << i.intensity << "\t" << i.tau << "\t" << realTemp << "\t" << bTemp << std::endl;
		}
		out.close();*/
		
		//calculating the width
		//special case if the ray didn't intersect the disk. Just set the width to zero
		if(allThree2.size()==0){//empty
			//std::cout << "missed the disk." << std::endl;
			std::vector<widthInfo> empty;
			for(int i=0;i<5;i++)
				empty.push_back(widthInfo());
			return empty;
		}
		std::vector<double> widths;
		//because the optical depth is accrued as you progress through, it can only increase monotonically.
		//the range is therefore given by first and the last.
		double min=allThree2.front().intensity; double max=allThree2.back().intensity;
		//now, we define the width as the range where it changes rapidly, which I'm arbitrarily deciding 
		//is from 5% to 95%.
		double range=max-min;
		double lowerlimit = min + 0.05*range;
		double upperlimit = min + 0.95*range;
		double lowMid = min + 0.33*range;
		double highMid = min + 0.67*range;
		int lowerindex=0; int upperindex=0; int lowmidindex=0; int highmidindex=0;
		for(int j=0;j<allThree2.size();j++){
			if(allThree2[j].intensity > lowerlimit){
				break;
			}
			lowerindex++;
		}
		for(int j=0;j<allThree2.size();j++){
			upperindex++;
			if(allThree2[j].intensity > upperlimit){
				break;
			}
		}
		if(upperindex==allThree2.size())	//make sure not to go off the edge
			upperindex--;
			
		for(int j=0;j<allThree2.size();j++){
			if(allThree2[j].intensity > lowMid){
				break;
			}
			lowmidindex++;
		}
		for(int j=0;j<allThree2.size();j++){
			highmidindex++;
			if(allThree2[j].intensity > highMid){
				break;
			}
		}
		if(highmidindex==allThree2.size())	//make sure not to go off the edge
			highmidindex--;
			
		//the "centroid" of the emission is the step closest to the middle of the range
		//int centroid = (lowerindex+upperindex)/2; //index of the middle point
			
		//double bTemp=tempFromSurfaceBrightness(allThree2.back()[i].intensity,frequencies[i]);
		double bTemp=tempFromSurfaceBrightness(valueSlow,frequency);
		//find where the brightness temperature matches the physical temperature most closely
		int tempindex=0;
		double mindif=std::numeric_limits<double>::max(); //we want to minimize the difference between t_phys and b_temp
		//std::cout << "btemp: " << bTemp << std::endl;
		for(int j=0;j<allThree2.size();j++){
			double tPhysical=(*diskTemp)(allThree2[j].pos.r(),allThree2[j].pos.theta(),allThree2[j].pos.phi());
			double dif = std::abs(tPhysical-bTemp);
			//std::cout << tempindex << "\t" << tPhysical << "\t" << dif << "\t" << allThree2[j][i].pos.z/AU << std::endl;
			if(dif<mindif){
				mindif=dif;
				tempindex=j;
			}
		}
		if(tempindex==allThree2.size())	//make sure not to go off the edge
			tempindex--;
		//std::cout << "final position: " << allThree2[tempindex][i].pos.x/AU << "\t" << std::abs(allThree2[tempindex][i].pos.z/AU) << std::endl;
		//std::cout << lowerindex << "\t" << lowmidindex << "\t" << tempindex << "\t" << highmidindex << "\t" << upperindex << std::endl;
			
		//return a vector with three width infos in it.
		//The 5% point, the 95% point, and the match point
		std::vector<widthInfo> wvect;
		wvect.push_back(allThree2[lowerindex]);
		wvect.push_back(allThree2[lowmidindex]);
		wvect.push_back(allThree2[tempindex]);
		wvect.push_back(allThree2[highmidindex]);
		wvect.push_back(allThree2[upperindex]);
			
		//double width= (allThree2[upperindex][i].pos - allThree2[lowerindex][i].pos).mag();
		//double width= allThree2[upperindex][i].first.z - allThree2[lowerindex][i].first.z;
		//width/=AU;
		//std::cout << allThree2[upperindex][i].pos.z/AU << "\t" << allThree2[lowerindex][i].pos.z/AU << "\t" << width << std::endl;
		//widths.push_back(width);
		return(wvect);
	}
	
	std::vector<std::vector<double> > propagateRayDiagnostic(const line l, const double frequency, const vect cameraPosition) const{
		std::vector<std::pair<intersection,intersection> > diskPairs=entriesAndExits(l);

		vect dir=l.direction/-l.direction.mag(); //ensure unit vector
		//std::cout << "pos:" << l.position << "\tdir: " << l.direction << std::endl;
		//std::cout << "dir: " << dir << std::endl;
		double opticalDepth=0;
		std::vector<double> temperatures; //temps at each step through the disk
		double valueSlow=0;
		unsigned int nsteps=100;
		
		std::vector<std::vector<double> > results;
		
		for(auto pair : diskPairs){
			int stepcount=0;
			vect pos=pair.first.location; //starting point
			double remainingDist=(pair.second.location-pair.first.location).mag(); //the total distance along which to integrate
			double stepsize = remainingDist/nsteps;	//currently nonadaptive stepsize
			pos+=dir*stepsize/2;
			
			while (remainingDist >= stepsize/2) {
				std::vector<double> stats; //z position, temperature, density, intensity etc for this ray
				double d=(*dens)(pos.r(),pos.theta(),pos.phi());	//gas density
				double ddust=(*dustDens)(pos.r(),pos.theta(),pos.phi());	//seperate dust density
				double temp=(*diskTemp)(pos.r(),pos.theta(),pos.phi());	//starting temp
				temperatures.push_back(temp);
		
				double r_cyl=pos.r()*sin(pos.theta());	//cylindrical radius
				double azdif=(cameraPosition.phi()-pos.phi());//difference between the camera orientation and the point in question
				
				double velocity=sqrt(gravConst*starMass/r_cyl/r_cyl/r_cyl)*r_cyl*sin(cameraPosition.theta())*cos(-pi/2+azdif);//velocity with respect to LOS
				double soundspeed=sqrt(temp*kboltzmann/(3.819239518e-24));

				stats.push_back(pos.z);
				stats.push_back(temp);
				stats.push_back(ddust);
				stats.push_back(valueSlow);
			
				double newfreq=doppler(frequency,velocity);
				double dustOp=(*dustOpac)(newfreq);
				double BB=blackBody(temp,newfreq);
				double opacity;
				opacity=d*(*lineOpac)(temp,newfreq,turbulence*soundspeed,freezeout)+ddust*dustOp;
				valueSlow+=stepsize*opacity*BB * (exp(-opticalDepth));
				opticalDepth+=(stepsize*opacity);
				
				pos+=dir*stepsize;
				remainingDist-=stepsize;
				stepcount++;
				results.push_back(stats);
			}
		}
		
		return(results);
	}
	
};
#endif
