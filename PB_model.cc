#include <iostream>
#include <fstream>
#include <math.h>       /* sqrt */
#include <cmath>        // std::abs
#include <stdlib.h>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMath.h"
#include "TH3D.h"
#include "TF1.h"
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "TObjString.h"
#include "TBrowser.h"
#include "TRandom.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>

using namespace std;

// A fuction for modellign the signal noise and resolution of a proton beam in a monolithic scintillator
//
// TODO:
//      - remove the object fuctnion and the variables for deffining it
//      - remove the p_biining histogram
//      - Change the file naming so that it no longer contains incorrect parameters
//      - Remove references to variable for lateral deviations in energy and scattering properties
//      - Add argument based options rather than having to manually change parameters as needed
//
//
//
//
//
//
//


float Epsilon = pow(13.6,2);
float m_rest, p0v0, Einit;
//set beam sigma parameters for initial conditions
float sigma_t_y_0, sigma_theta_y_0, sigma_t_theta_y_0;
float sigma_t_z_0, sigma_theta_z_0, sigma_t_theta_z_0;
//contain beam sigma parameters at each step (at x=0 sigma_i = simga_0)
float sigma_t_y_i, sigma_theta_y_i, sigma_t_theta_y_i;
float sigma_t_z_i, sigma_theta_z_i, sigma_t_theta_z_i;
//initial sigma parameters at the point for begining the PSF
float sigma_t_y_i_point, sigma_theta_y_i_point, sigma_t_theta_y_i_point;
float sigma_t_z_i_point, sigma_theta_z_i_point, sigma_t_theta_z_i_point;
// variance, probaility density for beam, probability density for PSF, initial probability density at start of PSF, variance squared
float var, p, p_point,p_point_0 , Var_sq;

//Number of particles just adds a scaling factor to the front
int N_particles;
//Number of pencil beams
int NPBY, NPBZ;
//Pencil beam spacing (calculated automatically from number and field size)
float PB_spacing_y, PB_spacing_z;
float N_point, N_pix, N_tot;
//position of the point to calcualte the PSF from
float point_x, point_y, point_z;
// toggle variables
int calc_point, PSF_fill;

//RSP value of the detector
float RSP_det;
//recording parameters
int NbinsX, NbinsY, NbinsZ;
float Xmax, Ymax, Zmax;
float Xmin, Ymin, Zmin;
// spotsizes taken from a clinical IBA system from Luis Perles
float spotsize[] = {6.8, 6.6,6.4,6.2,5.9,5.8,5.6,5.3,5.1,4.9,4.7,4.7,4.5,4.4,4.3,4.1,4,3.8,3.7,3.6,3.6,3.5,3.5,3.4,3.4,3.3,3.3,3.2,3.2,3,3};

// parameters from a defunct object funtion that allowed placing a object in front of the beam
int NbinsX_Obj, NbinsY_Obj, NbinsZ_Obj;
float Xmin_Obj, Xmax_Obj;
float Ymin_Obj, Ymax_Obj;
float Zmin_Obj, Zmax_Obj;
//Old parameters from an attempt to vary energy in the lateral axis to try and determine edges
float E_dep_vector[31];
float Einit_vector[31];
float Var_vector[31];
float Exp_E[31];
float Exp_E_sq[31];

// vectors for determining energy loss from G4 data
vector<double> Energy;
vector<double> dEdXBins;

//function to calculate beam sigmaat a depth in the scintillator
pair<float,pair<float,float>> sigma_cal(float x0, float x1,float p0v0, float sigma_t_i, float sigma_theta_i, float sigma_t_theta_i, float RSP);

//Gotschalk highland approcimation of fermi eyeges
TMatrixD Gott_High(float x0,float x1,float p0v0, float RSP);
// funtion to calculate the probability density at a depth in the scintillator
float calculate_contribution(float x0,float y0,float x1,float y1,float sigma_y,float sigma_z);
//defunct object function to be removed
TH3F* Generate_object(string object_type,int NbinsX_Obj,float Xmin_Obj,float Xmax_Obj,int NbinsY_Obj,float Ymin_Obj,float Ymax_Obj,int NbinsZ_Obj,float Zmin_Obj,float Zmax_Obj);
//storage of sigma parameters as calculated using the funtion above
pair<float,pair<float,float>>sigmas_z;
pair<float,pair<float,float>>sigmas_y;
pair<float,pair<float,float>>sigmas_z_point;
pair<float,pair<float,float>>sigmas_y_point;

int main(){
  string object_type = "WaterCylinder";
  //detector properties
  RSP_det = 1.0;
  NbinsX = 600;
  NbinsY = 200;
  NbinsZ = 200;

  Xmax = 300;
  Ymax = 50;
  Zmax = 50;

  Xmin = 0;
  Ymin = -50;
  Zmin = -50;

  //Object properties from depricated functionallity
  NbinsX_Obj = 300;
  NbinsY_Obj = 300;
  NbinsZ_Obj = 300;

  Xmin_Obj=0;
  Xmax_Obj=300;
  Ymin_Obj=-150;
  Ymax_Obj=150;
  Zmin_Obj=-150;
  Zmax_Obj=150;

  //positon variables (initial and current)
  float x0, y0, z0, x1, y1, z1;

  //set beam parameters - note these are reset later in the funtion so this could be removed if new file naming funtion was used
  N_particles = 1;
  Einit = 200;
  //initial y direction sigma
  sigma_t_y_0 = pow(3.2,2);
  sigma_theta_y_0 = 0;
  sigma_t_theta_y_0 = 0;
  //initial z direction sigmas
  sigma_t_z_0 = pow(3.2,2);
  sigma_theta_z_0 = 0;
  sigma_t_theta_z_0 =0;
  //Number of pencil beams in each direction
  NPBY = 1;
  NPBZ = 11;
  //Makes even symmetrical spacing
  PB_spacing_y = (Ymax-Ymin)/NPBY;
  PB_spacing_z = (Zmax-Zmin)/NPBZ;
  // depth of point for PSF determination
  point_x = 0;
  //toggle to store projections as individual beams or combined together
  int individual = 1;
  // proton rest mass
  m_rest= 938.27;
  //load G$ data for energy loss calculation
  std::string line;
  std::ifstream SPWater("dEdX/Water_Geant4.dat");
  double data[3];
  while(getline(SPWater, line)) {
    stringstream ss(line);
    for(int i=0;i<3;i++) ss >> data[i];
      Energy.push_back(data[0]);
      dEdXBins.push_back(data[1]);
  }

  // creates new file to save results in - note inaccurate file names from the previous defined parameters changing later in the function
  TFile *f1 = new TFile(Form("PB_model_%.0f_%.1f_%.1f_%i_%i_distal_noise_energy.root",Einit,sigma_t_y_0,sigma_t_z_0,NPBY,NPBZ),"recreate");

  //for quenching calculation
  float A    = 1.0;
  float kB   = 1.59E-2; // cm/MeV

  // 3D plots of the positiona, energy deposit and light yield within the scintillator
  TProfile3D *PB_model_scatter = new TProfile3D("PB_model_scatter","PB_model_scatter",NbinsX,Xmin,Xmax,NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
  TProfile3D *PB_model_Energy = new TProfile3D("PB_model_Energy","PB_model_Energy",NbinsX,Xmin,Xmax,NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
  TProfile3D *PB_model_Light = new TProfile3D("PB_model_Light","PB_model_Light",NbinsX,Xmin,Xmax,NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);


  // for saving the lateral projection
  TH2F *lateral_projection = new TH2F("lateral_projection", "lateral_projection",NbinsX,Xmin,Xmax,NbinsZ,Zmin,Zmax);

  //lateral variance projection is false and needs to be removed, lateral variance is determined by range straggling
  TH2F *lateral_variance_projection = new TH2F("lateral_variance_projection", "lateral_variance_projection",NbinsX,Xmin,Xmax,NbinsZ,Zmin,Zmax);
  // hist for storing the PSF, parameters defined later in the function
  TH2F *PSF_2D;

  TH2F *variance_projection_sq;
  TH2F *variance_projection;
  TH2F *distal_projection;
  // If indivifual is 0 then then only a single histogram of each of these projections is recorded
  if (individual==0){
    distal_projection = new TH2F("Distal_projection", "Distal_projection",NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
    variance_projection_sq = new TH2F("variance_projection_squared", "variance_projection_squared",NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
    variance_projection = new TH2F("variance_projection", "variance_projection",NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
  }
  // Deprecated histogram for removal
  TH2F *p_binning = new TH2F("p_bining", "p_binning",NbinsY,-7.,7.,NbinsZ,-7.,7.);

  // create the object from the unused object function
  TH3F *object = Generate_object(object_type,NbinsX_Obj,Xmin_Obj,Xmax_Obj,NbinsY_Obj,Ymin_Obj,Ymax_Obj,NbinsZ_Obj,Zmin_Obj,Zmax_Obj);
  //iterate of the number of defined pencil beams in the y direction
  for(int PBY = 0; PBY < NPBY; PBY++){
    // Position of current pencil beam in Y - set to 0 for central pencil beams
    y0 = 0;// Ymin + ((Ymax-Ymin)/(2*NPBY)) + (PBY*PB_spacing_y);

    // previous vectors from an attempt to vary properties laterally depending on the object chosen
    E_dep_vector[PBY]=0;
    Var_vector[PBY]=0;
    Exp_E[PBY]=0;
    Exp_E_sq[PBY]=0;

    int toggle = 0;
    // set thickness of material to pass through before the detector [mm of water]
    float WET = 50;
    //iterate of the number of pencil beams in the Z direction
    for(int PBZ = 0; PBZ < NPBZ; PBZ++){
      if(individual == 1){
        distal_projection = new TH2F(Form("Distal_projection_%i_%i",PBY,PBZ), "Distal_projection",NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
        variance_projection_sq = new TH2F(Form("variance_projection_squared_%i_%i",PBY,PBZ), "variance_projection_squared",NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);
        variance_projection = new TH2F(Form("variance_projection_%i_%i",PBY,PBZ), "variance_projection",NbinsY,Ymin,Ymax,NbinsZ,Zmin,Zmax);

      }
      //set initial energy (can be a function of pencil beam number in order to investigate this parameter)
      Einit = 100 + (10*PBZ);
      // set depth of point for determination of PSF
      point_x = 25;
      //set inital value of probability density at the point
      p_point_0 = 0;
      //create histogram for PSF on the first beam
      if(PBY == 0 &&PBZ ==0){
         PSF_2D = new TH2F(Form("PSF_%i_%i",PBZ,PBY), Form("PSF_%.1f",point_x),NbinsX,Xmin,Xmax,NbinsZ,Zmin,Zmax);
      }
      // Position of current pencil beam in Z - set to 0 for central pencil beams
      z0 = 0;//(Zmin + ((Zmax-Zmin)/(2*NPBZ))) + (PBZ*PB_spacing_z);

      // set lateral position of the point
      point_y = 0;
      point_z = 0.0;
      //determines whether the point ahs been reached and should be calcualted
      calc_point = 0;
      //determines whether the point has been calculated and the histogram needs to be filled
      PSF_fill =0;

      //initial y direction sigma
      sigma_t_y_i = pow(3,2);//sigma_t_y_0;
      sigma_theta_y_i =(sigma_theta_y_0);
      sigma_t_theta_y_i = (sigma_t_theta_y_0);

      x0 =0;
      //initial z direction sigmas
      sigma_t_z_i = pow(3,2);//sigma_t_z_0+PBY;
      sigma_theta_z_i = sigma_theta_z_0;
      sigma_t_theta_z_i = sigma_t_theta_z_0;
      float Energy_depo = 0;
      x0=0;
      //output beam parameters
      cout<<"Beam: "<<PBZ + (PBY*NPBZ)<<" Out of: "<<NPBY*NPBZ<<endl;
      cout<<"Spotsize: "<<TMath::Sqrt(sigma_t_y_i)<<" Energy: "<<Einit<<endl;
      //iterate of the X direction for the depth in teh scintillator
      for(float i=1; i<=Xmax+WET; i+=0.5){
          x1 = i;
          // fix for if the energy has become negative in the last iteration
          if(Einit < 0){
            if (toggle == 0){
            //Var_vector[PBY]=sigma_t_y_i;
            toggle = 1;
            }
            PSF_fill =1;
            Einit = 0;
         }
          //calculate momentum and velocity
          p0v0 = ((Einit + (2*m_rest))*Einit)/(Einit+m_rest);
          // get energy depositied in this step
          int idE    = lower_bound(Energy.begin(), Energy.end(), Einit)-Energy.begin();
          //calcualte the energy deposited
          Energy_depo = N_particles*(x1-x0)*RSP_det*dEdXBins[idE];

          //calcualte LET and quenching
          float LET = Energy_depo/(x1-x0);
          float L    = (x1-x0)*A*LET/(1+ kB*LET);
          // if the depth is greater than the depth of the point its PSF begins calculation
          if(x1>=point_x){
            if(calc_point == 0){
              //initial size of the point is very small and intiial divergence is taken from the main beam at this point only on first time it reaches this if
              sigma_t_y_i_point =0.001;
              sigma_t_z_i_point =0.001;
              sigma_theta_y_i_point =(sigma_theta_y_i);
              sigma_t_theta_y_i_point = (sigma_t_theta_y_i);
              sigma_theta_z_i_point =(sigma_theta_z_i);
              sigma_t_theta_z_i_point = (sigma_t_theta_z_i);
            }
            //calcualte sigma for points at this depth
            sigmas_y_point = sigma_cal(x0,x1,p0v0,sigma_t_y_i_point,sigma_theta_y_i_point,sigma_t_theta_y_i_point, RSP_det);
            sigmas_z_point = sigma_cal(x0,x1,p0v0,sigma_t_theta_z_i_point,sigma_theta_z_i_point,sigma_t_theta_z_i_point,RSP_det);

            //cout<<p<<endl;
            sigma_t_y_i_point = sigmas_y_point.first;
            sigma_theta_y_i_point = sigmas_y_point.second.first;
            sigma_t_theta_y_i_point = sigmas_y_point.second.second;
            //replace initial z direction sigmas
            sigma_t_z_i_point = sigmas_z_point.first;
            sigma_theta_z_i_point = sigmas_z_point.second.first;
            sigma_t_theta_z_i_point =sigmas_z_point.second.second;
            calc_point =1;
          }
          //calculate sigmas for main beam at this depth
          sigmas_y = sigma_cal(x0,x1,p0v0,sigma_t_y_i,sigma_theta_y_i,sigma_t_theta_y_i, RSP_det);
          sigmas_z = sigma_cal(x0,x1,p0v0,sigma_t_z_i,sigma_theta_z_i,sigma_t_theta_z_i,RSP_det);
          //cout<<x1<<endl;
          //iterate through lateral positions
          for(float j=Ymin; j<=Ymax; j+=Ymax/NbinsY){
              y1=j;
              for (float k=Zmin;k<=Zmax; k+=Zmax/NbinsZ){
                z1=k;
                N_point=0;
                N_tot=0;
                // calculate probability density at this depth for both the PSF and main beam
                p_point = calculate_contribution(point_y,point_z,y1,z1,sigmas_y_point.first,sigmas_z_point.first);
                p = calculate_contribution(y0,z0,y1,z1,sigmas_y.first,sigmas_z.first);
                //calculate the number of particles this corresponds to - for the point it is scaled by the initial number at this point
                N_point=p_point*N_particles*p_point_0;
                N_tot = p*N_particles;

                // sets initial probability at the start of the point depth
                if(x1==point_x && y1==point_y &&z1 ==point_z){
                  p_point_0 = p;
                  cout<<"p: "<<p<<endl;
                }
                //calcualte the variance and variance squared
                var = (p*pow(Energy_depo,2))-pow((p*Energy_depo),2);
                Var_sq = pow(var,2);

                //fills histograms if teh energy is realistic and the depth is greater than the thickness of material before the detector
                if (Einit > 0 && x1 >WET && Energy_depo>0){
                  PB_model_scatter->Fill(x1-WET,y1,z1,N_particles*p);
                  PB_model_Energy->Fill(x1-WET,y1,z1,p*Energy_depo);
                  PB_model_Light->Fill(x1-WET,y1,z1,p*L);
                  // old functionallity that needs to be removed
                  float p_y =(y1-y0)/sigma_t_y_i;
                  float p_z = (z1-z0)/sigma_t_z_i;
                  //if(p_y == 0) cout<<p<<endl;
                  p_binning->Fill(p_y, p_z, p*Energy_depo);
                  //fill PSF
                  if(p>0 && N_point>0 && Energy_depo>0){
                    if(PSF_fill == 0){
                      //cout<<N_point<<" "<<p<<" "<<(N_point/p)<<endl;
                      PSF_2D->Fill(x1-WET,z1,(N_point/N_tot));
                    }
                  }
                  distal_projection->Fill(y1,z1,p*Energy_depo);

                  lateral_projection->Fill(x1-WET,z1,p*Energy_depo);
                  if(var>0){
                    variance_projection_sq->Fill(y1,z1,Var_sq);
                    lateral_variance_projection->Fill(x1-WET,z1,Var_sq);
                  }
                  // from old attempt to add lateral deviation to energy and sigma - to be removed
                  if( y1== 0 && z1 == 0){
                    //cout<<Energy_depo<<endl;
                    E_dep_vector[PBY]+=Energy_depo;
                    Exp_E_sq[PBY] +=(p*pow(Energy_depo,2));
                    Exp_E[PBY] += (p*Energy_depo);
                    Var_vector[PBY]=sigmas_y.first;

                  }
                }
              }


          }
        // decrease energy by deposited amount
        Einit     -= Energy_depo/N_particles;
        //replace initial y direction sigma
        sigma_t_y_i = sigmas_y.first;
        sigma_theta_y_i = sigmas_y.second.first;
        sigma_t_theta_y_i = sigmas_y.second.second;
        //replace initial z direction sigmas
        x0 = x1;
        sigma_t_z_i = sigmas_z.first;
        sigma_theta_z_i = sigmas_z.second.first;
        sigma_t_theta_z_i =sigmas_z.second.second;
      }
      //old functionality attempt that needs to be removed
      Exp_E[PBY] += p*Energy_depo;
      //cout<<Exp_E_sq[PBY]<<" "<<Exp_E[PBY]<<" "<<TMath::Sqrt(Var_vector[PBY])<<" "<<Exp_E_sq[PBY]-pow(Exp_E[PBY],2)<<endl;
      //Write PSF hist after each beam ....
      PSF_2D->Write("",TObject::kOverwrite);
      // write projections if being recorded individually
      if (individual == 1){
        distal_projection->Write("",TObject::kOverwrite);
        // sums the variance in quadrature at each depth before writing to file
        for(int i=0;i<=NbinsY;i++){
          for(int j=0;j<=NbinsZ;j++){
              int var_sq_bin = variance_projection->GetBin(i,j);
            float var_sq_temp = variance_projection_sq->GetBinContent(var_sq_bin);
            variance_projection->SetBinContent(i,j,TMath::Sqrt(var_sq_temp));
          }
        }
        variance_projection->Write("",TObject::kOverwrite);
      }
    }

  }
  PB_model_scatter->Write("",TObject::kOverwrite);
  PB_model_Energy->Write("",TObject::kOverwrite);
  PB_model_Light->Write("",TObject::kOverwrite);
  p_binning->Write("",TObject::kOverwrite);
  distal_projection->Write("",TObject::kOverwrite);
  //sums variance at each depth in quadrature before writing
  for(int i=0;i<=NbinsY;i++){
    for(int j=0;j<=NbinsZ;j++){
      int var_sq_bin = variance_projection->GetBin(i,j);
      float var_sq_temp = variance_projection_sq->GetBinContent(var_sq_bin);
      variance_projection->SetBinContent(i,j,TMath::Sqrt(var_sq_temp));
    }
  }
  // if not writing individual projections then the total is written here
  if(individual == 0){
  variance_projection->Write("",TObject::kOverwrite);
  lateral_projection->Write("",TObject::kOverwrite);
  lateral_variance_projection->Write("",TObject::kOverwrite);
  }
  //object->Write("",TObject::kOverwrite);
  f1->Close();
}

TMatrixD Gott_High(float x0,float x1,float p0v0, float RSP){
//calculate the scattering power using Gottschalk's adaption of Highland's formula
  float step = x1-x0;

  float rad_length;
  if (RSP > 0.01) rad_length = 361;
  else rad_length = 3.5E4;

  float E0s0      = Epsilon*pow(1+0.038* TMath::Log(step/rad_length),2);
  float GH_t = E0s0*(pow(step,3)/(pow(p0v0,2)*rad_length));
  float GH_theta = E0s0*(step/(pow(p0v0,2)*rad_length));
  float GH_t_theta = E0s0*(pow(step,2)/(pow(p0v0,2)*rad_length));

  TMatrixD GH_matrix(3,1);
  GH_matrix(0,0) = GH_t;
  GH_matrix(1,0) = GH_theta;
  GH_matrix(2,0) = GH_t_theta;
  return GH_matrix;
}

pair<float,pair<float,float>>  sigma_cal(float x0, float x1,float p0v0, float sigma_t_i, float sigma_theta_i, float sigma_t_theta_i,float RSP){
//caluclate sigma at a depth x1 (sigma_f) given the sigma at a point x0
  TMatrixD GH_matrix(3,1);
  GH_matrix = Gott_High(x0,x1,p0v0,RSP);
  float step=x1-x0;
  float sigma_t_f = sigma_t_i + 2*sigma_t_theta_i*step + (sigma_theta_i*pow(step,2)) + GH_matrix(0,0);
  float sigma_theta_f = sigma_theta_i + GH_matrix(1,0);
  float sigma_t_theta_f = sigma_t_theta_i +(sigma_theta_i*step) + GH_matrix(2,0);
  //cout<<sigma_t_f<<endl;
  pair<float,float> first_pair = pair<float,float>(sigma_theta_f,sigma_t_theta_f);
  pair<float,pair<float,float>> sigmas = pair<float,pair<float,float>>(sigma_t_f, first_pair);

  return sigmas;

}

float calculate_contribution(float y0,float z0,float y1,float z1,float sigma_y,float sigma_z){
//calculate the lateral probability distribution for each x and y position
  //cout<<"other y1: "<<y1<<endl;
  float p = (1./(TMath::Pi()*(sigma_y + sigma_z)))*TMath::Exp(-1*((pow(y1-y0,2)+pow(z1-z0,2))/(sigma_y+sigma_z)));
  return p;
}

TH3F* Generate_object(string object_type, int NbinsX_Obj,float Xmin_Obj,float Xmax_Obj,int NbinsY_Obj,float Ymin_Obj,float Ymax_Obj,int NbinsZ_Obj,float Zmin_Obj,float Zmax_Obj){

  TH3F *object = new TH3F("object","object",NbinsX_Obj,Xmin_Obj,Xmax_Obj,NbinsY_Obj,Ymin_Obj,Ymax_Obj,NbinsZ_Obj,Zmin_Obj,Zmax_Obj);

  for(float x=Xmin_Obj; x<Xmax_Obj; x++){
    for(float y=Ymin_Obj; y<Ymax_Obj; y++){
      for(float z = Zmin_Obj; z<Zmax_Obj;z++){
        if(object_type == "WaterCylinder"){
          if( pow((x-150),2) + pow(y,2) < pow(100,2) && z < 150 && z > -150){
            object->Fill(x,y,z,1);
          }
          else{
            object->Fill(x,y,z,0.001);
          }
        }
      }
    }
  }
  return object;
}
