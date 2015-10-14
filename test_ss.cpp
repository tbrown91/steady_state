#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
using namespace std;

void steady_state(const double alpha,const double beta,const double gamma,const double delta,const double epsilon,const double phi,double* nuc_ss,double* cyt_ss){
  double lambda_1, lambda_2, lambda_3;
  //Calculate rate paramters for the Poisson distrubtions
  lambda_1 = (alpha*gamma)/(delta*(alpha+beta));
  lambda_2 = (alpha*gamma)/(epsilon*(alpha+beta));
  lambda_3 = (alpha*gamma)/(phi*(alpha+beta));
  //Calculate the nuclear distribution
  for (int i=0;i<50;++i){
    nuc_ss[i] = 0;
    for (int j=0;j<=i;++j){
      nuc_ss[i] += (pow(lambda_1,i-j)*pow(lambda_2,j)*exp(-(lambda_1+lambda_2)))/(tgamma(i+1-j)*tgamma(j+1));
    }
    nuc_ss[i] *= 1000;
  }
  //Calculate the cytoplasmic distribution
  for (int i=0;i<50;++i){
    cyt_ss[i] = ((pow(lambda_3,i)*exp(-lambda_3))/tgamma(i+1))*1000;
  }
}

void stochastic_simulation(const double alpha,const double beta,const double gamma,const double delta,const double epsilon,const double phi,double* nuc_sim,double* cyt_sim){
  for (int i=0;i<50;++i){
    nuc_sim[i] = 0;
    cyt_sim[i] = 0;
  }
  for (int it=0;it<1000;++it){
    double current_time = 0.0;
    double tau;
    double propensities[6] = {0};
    double tot_propensity;
    double concentrations[4] = {0};
    double time_rand, reac_rand;
    while (current_time < 200){
      //Calculate propensities
      //Promoter turns on
      propensities[0] = alpha * (1-concentrations[0]);
      //Promoter turns off
      propensities[1] = beta * concentrations[0];
      //New transcript initiated
      propensities[2] = gamma * concentrations[0];
      //Transcript is complete mRNA
      propensities[3] = delta * concentrations[1];
      //mRNA exported to the cytoplasm
      propensities[4] = epsilon * concentrations[2];
      //mRNA degraded in the cytoplasm
      propensities[5] = phi * concentrations[3];
      tot_propensity = 0.0;
      for (int i=0;i<6;++i){tot_propensity+=propensities[i];}
      //Calculate time to next reaction
      time_rand = (double)rand()/(RAND_MAX);
      tau = (-1/tot_propensity)*log(1-time_rand);
      current_time += tau;
      //Calculate which reaction occurs next
      reac_rand = ((double)rand()/(RAND_MAX))*tot_propensity;
      int reac_index=0;
      while (reac_rand > propensities[reac_index]){
        reac_rand -= propensities[reac_index];
        ++reac_index;
      }
      //Perform the relevant reaction
      switch (reac_index){
        case 0:
          concentrations[0] = 1;
          break;
        case 1:
          concentrations[0] = 0;
          break;
        case 2:
          concentrations[1] += 1;
          break;
        case 3:
          concentrations[1] -= 1;
          concentrations[2] += 1;
          break;
        case 4:
          concentrations[2] -= 1;
          concentrations[3] += 1;
          break;
        case 5:
          concentrations[3] -= 1;
          break;
        default:
          cout << "Reaction index error! " << endl;
          break;
      }
    }
    nuc_sim[(int)(concentrations[1]+concentrations[2])] += 1;
    cyt_sim[(int)concentrations[3]] += 1;
  }
}

void calculate_chiSq(const double* nuc_ss,const double* cyt_ss,const double* nuc_sim,const double* cyt_sim,double& nuc_chiSq,double& cyt_chiSq){
  nuc_chiSq = 0.0;
  cyt_chiSq = 0.0;
  for (int i=0;i<50;++i){
    if ((nuc_ss[i] >= 5) && (nuc_sim[i] >= 5)){
      nuc_chiSq += ((nuc_ss[i] - nuc_sim[i])*(nuc_ss[i] - nuc_sim[i]))/nuc_sim[i];
    }
    if ((cyt_ss[i] >= 5) && (cyt_sim[i] >= 5)){
      cyt_chiSq += ((cyt_ss[i] - cyt_sim[i])*(cyt_ss[i] - cyt_sim[i]))/cyt_sim[i];
    }
  }
}

int main(int argc, char* argv[]){

  double nuc_chiSq, cyt_chiSq;
  ofstream parameter_file;
  ofstream ss_nucFile, ss_cytFile;
  ofstream sim_nucFile, sim_cytFile;
  parameter_file.open("parameters.txt",ios_base::out);
  ss_nucFile.open("steady_statesNuc.txt",ios_base::out);
  sim_nucFile.open("simulations_nuc.txt",ios_base::out);
  ss_cytFile.open("steady_statesCyt.txt",ios_base::out);
  sim_cytFile.open("simulations_cyt.txt",ios_base::out);
  //Test range of paramter values
  int j = 1;
  for (double alpha=2;alpha<=2;alpha+=0.2){
    for (double beta=0.2;beta<=2;beta+=0.2){
      for (double gamma=0.2;gamma<=2;gamma+=0.2){
        for (double delta=0.2;delta<=2;delta+=0.2){
          for (double epsilon=0.2;epsilon<=2;epsilon+=0.2){
            for (double phi=0.2;phi<=2;phi+=0.2){
              double nuc_ss[50], cyt_ss[50];
              double nuc_sim[50], cyt_sim[50];
              steady_state(alpha,beta,gamma,delta,epsilon,phi,nuc_ss,cyt_ss);
              stochastic_simulation(alpha,beta,gamma,delta,epsilon,phi,nuc_sim,cyt_sim);
              //Write
              for (int i=0;i<50;++i){
                ss_nucFile << nuc_ss[i] << "\t";
                sim_nucFile << nuc_sim[i] << "\t";
              }
              ss_nucFile << endl;
              sim_nucFile << endl;
              for (int i=0;i<50;++i){
                ss_cytFile << cyt_ss[i] << "\t";
                sim_cytFile << cyt_sim[i] << "\t";
              }
              ss_cytFile << endl;
              sim_cytFile << endl;
              calculate_chiSq(nuc_ss,cyt_ss,nuc_sim,cyt_sim,nuc_chiSq,cyt_chiSq);
              parameter_file << nuc_chiSq << "\t" << cyt_chiSq << "\t" << alpha << "\t" << beta << "\t" << gamma << "\t" << delta << "\t" << epsilon << "\t" << phi << endl;
            }
          }
        }
      }
    }
  }
  parameter_file.close();
  ss_nucFile.close();
  ss_cytFile.close();
  sim_nucFile.close();
  sim_cytFile.close();
}
