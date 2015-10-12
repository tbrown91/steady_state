#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;

void calculate_cytoProb(const double lambda_1,const double lambda_2,const double lambda_3,const double mod_lambda,double* cyto_prob){
  //Calculate the probability of having x mRNA molecules in the cytoplasm for the given paramter values for x in [0,49]
  double decay_rate;
  decay_rate = exp(-mod_lambda + lambda_1 + lambda_2);
  for (int i=0;i<50;++i){
    cyto_prob[i] = (pow(lambda_3,i)*decay_rate)/tgamma(i+1);
  }
}

void calculate_nucProb(const double lambda_1,const double lambda_2,const double lambda_3,const double mod_lambda,double* nuc_prob){
  //Calculate the probability of having x mRNA molecules in the nucleus for the given paramter values for x in [0,49]
  double decay_rate;
  decay_rate = exp(-mod_lambda + lambda_3);
  for (int i=0;i<50;++i){
    nuc_prob[i]=0;
    for (int j=0;j<=i;++j){
      nuc_prob[i] += (pow(lambda_1,i-j)*pow(lambda_2,j)*decay_rate)/(tgamma(i+1-j)*tgamma(j+1));
    }
  }
}

double calculate_chiSq(double* nuc_prob,double* cyt_prob,int* nuc_data,int* cyt_data){
  double chi_sq;
  for (int i=0;i<50;++i){
    if ((nuc_prob[i] >= 5) && (nuc_data[5] > 0)){chi_sq += ((nuc_prob[i]-nuc_data[i])*(nuc_prob[i]-nuc_data[i]))/nuc_prob[i];}
    if ((cyt_prob[i] >= 5) && (cyt_prob[i] > 0)){chi_sq += ((cyt_prob[i]-cyt_data[i])*(cyt_prob[i]-cyt_data[i]))/cyt_prob[i];}
  }
  return chi_sq;
}

void test_parameters(const double alpha,const double beta,const double gamma,const double delta,int* nuc_data,int* cyt_data,ofstream& outFile,const int tot_data){
  //For a range of paramter values, calcuate the probability distribution of nucleear and cytoplasmic RNA and calculate
  //the metric determining how well the parameters fit the data

  double nuc_prob[50], cyt_prob[50];
  double chi_sq;
  double lambda_1, lambda_2,lambda_3, mod_lambda;
  lambda_1 = alpha/beta;
  lambda_2 = alpha/gamma;
  lambda_3 = alpha/delta;
  mod_lambda = sqrt(lambda_1*lambda_1 + lambda_2*lambda_2 + lambda_3*lambda_3);
  //Calculate the probability distribution for the nuclear spots
  calculate_nucProb(lambda_1,lambda_2,lambda_3,mod_lambda,nuc_prob);

  //Calculate the probability distribution for the cytoplasmic spots
  calculate_cytoProb(lambda_1,lambda_2,lambda_3,mod_lambda,nuc_prob);

  //Normalise probabilities to the raw data
  for (int j=0;j<50;++j){
    nuc_prob[j] *= tot_data;
    cyt_prob[j] *= tot_data;
  }

  //Calculate the Chi-Square statistic for the given set of paramters
  chi_sq = calculate_chiSq(nuc_prob,cyt_prob,nuc_data,cyt_data);

  //Write results to the output file
  outFile << chi_sq << "\t" << alpha << "\t" << beta << "\t" << gamma << "\t" << delta << endl;
}
/*
void read_data(string nuc_fileName, string cyt_fileName, ifstream nucFile, ifstream cytFile, int *nuc_data, int *cyt_data){
  //Read the data from file for nuclear and cytoplasmic dots
  nucFile.open(nuc_fileName,ios_base::in);
  cytFile.open(cyt_fileName,ios_base::in);
  size_t sz;
  if (nucFile.is_open()){
    string line_contents;
    int i = 0;
    while (nucFile.good()){
      getline(nucFile, line_contents);
      nuc_data[i] = stoi(line_contents,sz);
      ++i;
    }
    nucFile.close();
  }else{cout << "Nuclear data file given is incorrect" << endl;return}
  if (cytFile.is_open()){
    string line_contents;
    int i = 0;
    while (cytFile.good()){
      getline(cytFile, line_contents);
      cyt_data[i] = stoi(line_contents,sz);
      ++i;
    }
    cytFile.close();
  }else{cout << "Cytoplasmic data file given is incorrect" << endl;}
}
*/
int main(int argc, char* argv[]){

  cout << "Test statement" << endl;
  cout << "First argument: " << argv[1] << endl;
  cout << "Second argument: " << argv[2] << endl;
  cout << "Third argument: " << argv[3] << endl;

  ifstream nucFile, cytFile;
  int nuc_data[50], cyt_data[50];
  //read_data(argv[1],argv[2],nucFile,cytFile);

  //File to write results to
  ofstream outFile;
  outFile.open(argv[3]);
  //Write the file header
  outFile << "Chi_sq\talpha\tbeta\tgamma\tdelta" << endl;

  int tot_data = 0;
  for (int j=0;j<50;++j){
    tot_data += nuc_data[j];
  }

  //Introduce parameters, these correspond to:
  //alpha = init x on / (on + off), the rate of transcription initiation
  //beta = elong, the rate of elongation
  //gamme = exp, the rate of export from the nucleus to the cytoplasm
  //delta = deg, the rate of degradation in the cytoplasm
  for (double alpha=0.1;alpha<1;alpha+=0.1){
    for (double beta=0.1;beta<1;beta+=0.1){
      for (double gamma=0.1;gamma<1;gamma+=0.1){
        for (double delta=0.1;delta<1;delta+=0.1){
          test_parameters(alpha,beta,gamma,delta,nuc_data,cyt_data,outFile,tot_data);
        }
      }
    }
    cout << alpha << endl;
  }

  outFile.close();

  return 0;
}
