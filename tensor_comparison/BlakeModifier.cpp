#include "Modifier.h"
#include <iostream>
#include <math.h>

Modifier::Modifier(void)
{
}

Modifier::~Modifier(void)
{
}

PinBackbone::PinBackbone(vector<Polymer*> sys, double initpk, double initsk, double initstiffk, int numpin)
{
	this->pk = initpk;
	this->sk = initsk;
	this->stiffk = initstiffk;
        this->linksize = 1.0;
	this->numpin = numpin;
	for(unsigned int i=0; i<sys.size() ; i++){
		for(int j=0; j<numpin; j++){
			pins.push_back(sys[i]->Loc[j]);
		}
	}
}
void
PinBackbone::Act(vector<Polymer*> sys, Polymer* yolk)
{

	//do the chain forces for each polymer
	for(unsigned int i=0; i<sys.size() ; i++){

		//first the springforce backbone
		for(int j=1; j<sys[i]->Length-1; j++){
			f1 = sys[i]->Loc[j-1] - sys[i]->Loc[j];
			f2 = sys[i]->Loc[j+1] - sys[i]->Loc[j];
			mag1 = sqrt(f1*f1);
			mag2 = sqrt(f2*f2);
			sys[i]->Vel[j] += ((mag1-linksize)*sk/mag1)*f1 + ((mag2-linksize)*sk/mag2)*f2;
		}
		f1 = sys[i]->Loc[1] - sys[i]->Loc[0];
		f2 = sys[i]->Loc[sys[i]->Length-2] - sys[i]->Loc[sys[i]->Length-1];
		mag1 = sqrt(f1*f1);
		mag2 = sqrt(f2*f2);
		sys[i]->Vel[0] += ((mag1-linksize)*sk/mag1)*f1;
		sys[i]->Vel[sys[i]->Length-1] += ((mag2-linksize)*sk/mag2)*f2;

		//and now, the stiffness
		for(int j=2; j<sys[i]->Length-2; j++){
			sys[i]->Vel[j] += (sys[i]->Loc[j]*2.0 - sys[i]->Loc[j+2] - sys[i]->Loc[j-2])*stiffk;
		}
		sys[i]->Vel[1] += (sys[i]->Loc[1] - sys[i]->Loc[3])*stiffk;
		sys[i]->Vel[0] += (sys[i]->Loc[0] - sys[i]->Loc[2])*stiffk;
		sys[i]->Vel[sys[i]->Length-2] += (sys[i]->Loc[sys[i]->Length-2] - sys[i]->Loc[sys[i]->Length-4])*stiffk;
		sys[i]->Vel[sys[i]->Length-1] += (sys[i]->Loc[sys[i]->Length-1] - sys[i]->Loc[sys[i]->Length-3])*stiffk;

		//finally, pin the first point to its initial location
		for(int j=0; j<numpin; j++){
			sys[i]->Vel[j] += (pins[i*numpin + j] - sys[i]->Loc[j])*pk;
		}
		//cout << sys[i]->Loc[0]<<endl;

		//STEVE: add in force that orients first link so that it is normal to wall
		//sys[i]->Vel[1] += (vect_d(0.0,0.0,1.0) - f1/mag1)*3.0;//sloppily at first
		//NEVER MIND. Let basek = 100.0, this will simply need to be hard coded in for now.

		//cout<<basek;

		sys[i]->Vel[1] += -100.0*(f1 - vect_d(0.0,1.0)); 

		//wall repulsion
		for(int j=1; j<sys[i]->Length; j++){
			if(sys[i]->Loc[j].ycomp() < 0.5){
				sys[i]->Vel[j] += (-100.0 + 100.0/(sys[i]->Loc[j].ycomp()*sys[i]->Loc[j].ycomp()*sys[i]->Loc[j].ycomp()*sys[i]->Loc[j].ycomp()*16.0))*y_hat;
				if(sys[i]->Loc[j].ycomp() < 0){
					cout<<"wall broken!! " << "index:" << i <<"-"<< j << " "<< sys[i]->Loc[j].ycomp() <<endl;
					throw exception();
				}
			}
		}
	}
}

KinesinShift::KinesinShift(double ck){
	this->k = ck;
}

void
KinesinShift::Act(vector<Polymer*> sys, Polymer* yolk){

	for(unsigned int i=0; i<sys.size() ; i++){
		for(int j=1; j<sys[i]->Length-1; j++){
			sys[i]->Vel[j] += (sys[i]->Loc[j-1] - sys[i]->Loc[j+1])*k;
		}
		sys[i]->Vel[0] += (sys[i]->Loc[0] - sys[i]->Loc[1])*k;
		sys[i]->Vel[sys[i]->Length-1] += (sys[i]->Loc[sys[i]->Length-2] - sys[i]->Loc[sys[i]->Length-1])*k;
	}
}

OseenTensor::OseenTensor(vector<Polymer*> csys, double ck){
	this->k = ck;
	for(unsigned int i=0; i< csys.size(); i++){
		tempsys.push_back(new Polymer(csys[i]->Length));
	}
}
void
OseenTensor::Act(vector<Polymer*> sys, Polymer* yolk){
	if(1){


	for(unsigned int n=0; n < sys.size(); n++){
		for(int i=0; i<sys[n]->Length; i++){
			tempsys[n]->Vel[i].set(0,0);
			//double loop all the way...
			for(unsigned int m=0; m < sys.size(); m++){
				for(int j=0; j<sys[m]->Length; j++){
					//if looking at yourself

					/*
					if( i==j && m==n ){
						tempsys[n]->Vel[i] += 0.3*sys[m]->Vel[j];
						continue;
					}
					*/
					//calculate some values that will be needed in the oseen tensor, and cache them temporarily

					dif = sys[m]->Loc[j] - sys[n]->Loc[i];
					normsq = dif*dif;
					norm = sqrt(normsq);
					//wall calculation caching
					mirdif = dif;
					mirdif.setY(sys[m]->Loc[j].ycomp() + sys[n]->Loc[i].ycomp());
					mirnormsq = mirdif*mirdif;
					//mirnorm4 = mirnormsq*mirnormsq;
					mirnorm = sqrt(mirnormsq);
					h = sys[n]->Loc[i].ycomp();

					//hminusz = h - mirdif.ycomp();
					vect_d hardcore;


					/*
					if( norm == 0.0 ){
						continue;
					}
					*/
					if( i!=j || m!=n){
					//hardcore repulsion
					if( norm < .5 ){
						hardcore = (1.0 - 1.0/(normsq*normsq*16.0))*dif;//tempsys[n]->Vel[i] += (1.0 - 1.0/(normsq*normsq*16.0))*dif;
					}
					else{
						hardcore = dif*0.0;
					}
					if( norm < .05 ){
						cout << "close " << norm<< endl;
					}
					//calculate the oseen tensor

						tempsys[n]->Vel[i] += (sys[m]->Vel[j] + hardcore + (dif * (dif*(sys[m]->Vel[j]+hardcore)))*(1.0/normsq) )*(k/norm);
					}
					else{
						//hardcore = dif*0.0;
						tempsys[n]->Vel[i] += sys[m]->Vel[j];//+hardcore;
					}
					//now to include a wall
					
					if( i!=j || m!=n ){
						mirVel = sys[m]->Vel[j]+hardcore;
						mirVel.setY(-mirVel.ycomp());
						mVeldotmDif = mirVel*mirdif;
						tempsys[n]->Vel[i] += ( -1.0*(sys[m]->Vel[j]+hardcore +(mirdif * (mirdif*(sys[m]->Vel[j]+hardcore)))*(1.0/mirnormsq) ) //reflected oseen
												+ (2*h/mirnormsq)*(
													- (sys[n]->Loc[i].ycomp())*mirVel //mirrored force direction
													+ (3*mVeldotmDif*sys[n]->Loc[i].ycomp()/mirnormsq - mirVel.ycomp())*mirdif //mirrored position direction
													+ (mVeldotmDif)*y_hat //y_hat component correction
												)
											) * (k/mirnorm);
					}
					
											
				}
			} //end inner loop
			for(int j=0;j<yolk->Length;j++){
					dif = sys[n]->Loc[i] - yolk->Loc[j];
					normsq = dif*dif;
					norm = sqrt(normsq);
					//wall calculation caching
					mirdif = dif;
					mirdif.setY(sys[n]->Loc[i].ycomp() + yolk->Loc[j].ycomp());
					mirnormsq = mirdif*mirdif;
					//mirnorm4 = mirnormsq*mirnormsq;
					mirnorm = sqrt(mirnormsq);
					h = yolk->Loc[j].ycomp();

					vect_d hardcore;

					//hardcore repulsion STEVE:
					if( norm < .5 ){
						hardcore = (1.0 - 1.0/(normsq*normsq*16.0))*dif;//tempsys[n]->Vel[i] += (1.0 - 1.0/(normsq*normsq*16.0))*dif;
					}
					else{
						hardcore = dif*0.0;
					}
					//calculate the oseen tensor

					yolk->Vel[j] += (sys[n]->Vel[i] + hardcore + (dif * (dif*(sys[n]->Vel[i]+hardcore)))*(1.0/normsq) )*(k/norm);

					//now to include a wall
					
					mirVel = sys[n]->Vel[i]+hardcore;
					mirVel.setY(-mirVel.ycomp());
					mVeldotmDif = mirVel*mirdif;
					yolk->Vel[j] += ( -1.0*(sys[n]->Vel[i]+hardcore +(mirdif * (mirdif*(sys[n]->Vel[i]+hardcore)))*(1.0/mirnormsq) ) //reflected oseen
											+ (2*h/mirnormsq)*(
												- (yolk->Loc[j].ycomp())*mirVel //mirrored force direction
												+ (3*mVeldotmDif*yolk->Loc[j].ycomp()/mirnormsq - mirVel.ycomp())*mirdif //mirrored position direction
												+ (mVeldotmDif)*y_hat //y_hat component correction
											)
										) * (k/mirnorm);
					
			}
		}
	} //end outer loop
	for(unsigned int n=0; n < sys.size(); n++){
		for(int i=0; i<sys[n]->Length; i++){
			sys[n]->Vel[i] = tempsys[n]->Vel[i];
		}
	}
	} //end toggle switch

}
