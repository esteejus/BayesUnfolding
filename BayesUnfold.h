#ifndef BAYESUNFOLD_H
#define BAYESUNFOLD_H

#include <iostream>
#include <vector>
#include <map>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TMatrixD.h"
#include "TF3.h"
#include "TMatrixDSparse.h"

using namespace std;

struct element {
  int    idx;    //index position either i or j depending on use
  double v; //value
};

struct neighbor {
  vector<int> n_bins; //neighbor bins in global index
  vector<double> weight; // weights of those global bins
};

class BayesUnfold
{
private:

  //  TMatrixD *response_m = nullptr;
  int dim_e = 0;//dimension of effect array
  int dim_c = 0;//dimension of cause array
  map< pair<int,int>, double> response_m;//response matrix

  map< pair<int,int>, double> dn_dn; //for error matrix of n(Ci) from n(E) contribution
  vector< map< pair<int,int>, double> > *dn_dP = nullptr; //for error matrix of n(Ci) from P(E|C) contribution

  map< pair<int,int>, double> dn_dn_prev;    //prev error matrix
  vector< map< pair<int,int>, double> > *dn_dP_prev = nullptr; //prev error matrix
  map< pair<int,int>, double> n_var;    //n_ci variance
  TMatrixDSparse var;

  vector<double> *e_vec = nullptr; //observed events
  vector<double> *f_vec = nullptr; //observed events
  vector<double> *c_vec = nullptr; //cause events 
  vector<double> *eff_vec = nullptr; //efficiency 
  vector<double> *sum_vec = nullptr; //sum of columns  

  //contains non zero elemnts of P(E|C) for each row i;  which has dimension of Effect
  vector< vector<element> > *P_EC_i = nullptr;
  //contains non zero elemnts of P(E|C) for each column j; which has dimension of Cause
  vector< vector<element> > *P_EC_j = nullptr; 

  //contains non zero elemnts of P(C|E) for each column which has dimension of Effect dimension
  // also known as M_ij which has dimension n_c x n_e 
  vector< vector<element> > *M_j = nullptr;
  vector< vector<element> > *dn_dn_prev_i = nullptr; //by cause dimension C
  vector< vector<element> > *dn_dn_prev_j = nullptr; //by cause dimension E
  
  vector<double> *p_vec        = nullptr;  //prior 
  vector<double> *p_vec_prev   = nullptr;  //previous itteration prior
  vector<double> *no_vec       = nullptr; //N_prime previous
  vector<double> *n_prime      = nullptr; //N_prime 

  vector<vector<int>> *e_gbin_map  = nullptr; //global bin map for observed events 
  vector<vector<int>> *c_gbin_map  = nullptr; //global bin map for cause events 
  vector<vector<int>> *c_gbin_next = nullptr; //global bin neighbor map

  vector<neighbor> *n_vect = nullptr; //global bin neighbor map

  TH3D *kernel = nullptr;

public:
  BayesUnfold()
  {
  }

  BayesUnfold(TH1 *obs_h, TH1 *true_h)
  {
    dim_e = obs_h -> GetNbinsX();
    dim_c = true_h -> GetNbinsX();

    P_EC_i = new vector< vector<element> >(dim_c,vector<element>(0));
    P_EC_j = new vector< vector<element> >(dim_e,vector<element>(0));  
    M_j = new vector< vector<element> >(dim_e,vector<element>(0));
    dn_dn_prev_i  = new vector< vector<element> >(dim_c,vector<element>(0));//by cause dimension C
    dn_dn_prev_j  = new vector< vector<element> >(dim_e,vector<element>(0));//by cause dimension C
    dn_dP_prev = new vector< map<pair<int,int>, double> >(dim_c,map<pair<int,int>,double>()); //prev error matrix
    dn_dP = new vector< map<pair<int,int>, double>>(dim_c,map<pair<int,int>,double>()); //prev error matrix
	  
    e_vec = new vector<double>(dim_e,0.);
    f_vec = new vector<double>(dim_e,0.);
    e_gbin_map = new vector<vector<int>>(dim_e,vector<int>(1,-1));
    
    c_vec = new vector<double>(dim_c,0.);
    n_vect = new vector<neighbor>(dim_c);
    c_gbin_map = new vector<vector<int>>(dim_c,vector<int>(1,-1));
    c_gbin_next = new vector<vector<int>>(dim_c,vector<int>(0));
    eff_vec = new vector<double>(dim_c,0.);
    sum_vec = new vector<double>(dim_c,0.);
  }

  BayesUnfold(TH2 *obs_h, TH2 *true_h)
  {
      dim_e = obs_h -> GetNbinsX() * obs_h -> GetNbinsY();
      dim_c = true_h -> GetNbinsX() * true_h -> GetNbinsY();

      P_EC_i = new vector< vector<element> >(dim_c,vector<element>(0));
      P_EC_j = new vector< vector<element> >(dim_e,vector<element>(0));  
      M_j = new vector< vector<element> >(dim_e,vector<element>(0));
      dn_dn_prev_i  = new vector< vector<element> >(dim_c,vector<element>(0));//by cause dimension C
      dn_dn_prev_j  = new vector< vector<element> >(dim_e,vector<element>(0));//by cause dimension C
      dn_dP_prev = new vector< map<pair<int,int>, double> >(dim_c,map<pair<int,int>,double>()); //prev error matrix
      dn_dP = new vector< map<pair<int,int>, double>>(dim_c,map<pair<int,int>,double>()); //prev error matrix

	
      e_vec = new vector<double>(dim_e,0.);
      f_vec = new vector<double>(dim_e,0.);
      e_gbin_map = new vector<vector<int>>(dim_e,vector<int>(2,-1));
      
      c_vec = new vector<double>(dim_c,0.);
      n_vect = new vector<neighbor>(dim_c);
      c_gbin_map = new vector<vector<int>>(dim_c,vector<int>(2,-1));
      c_gbin_next = new vector<vector<int>>(dim_c,vector<int>(0));
      eff_vec = new vector<double>(dim_c,0.);
      sum_vec = new vector<double>(dim_c,0.);
  }

  BayesUnfold(TH3 *obs_h, TH3 *true_h)
  {
    dim_e = obs_h -> GetNbinsX() * obs_h -> GetNbinsY() * obs_h -> GetNbinsZ();
    dim_c = true_h -> GetNbinsX() * true_h -> GetNbinsY() * true_h -> GetNbinsZ();

    P_EC_i = new vector< vector<element> >(dim_c,vector<element>(0));
    P_EC_j = new vector< vector<element> >(dim_e,vector<element>(0));  
    M_j = new vector< vector<element> >(dim_e,vector<element>(0));
    dn_dn_prev_i  = new vector< vector<element> >(dim_c,vector<element>(0));//by cause dimension C
    dn_dn_prev_j  = new vector< vector<element> >(dim_e,vector<element>(0));//by cause dimension E
    dn_dP_prev = new vector< map<pair<int,int>, double> >(dim_c,map<pair<int,int>,double>()); //prev error matrix
    dn_dP = new vector< map<pair<int,int>, double> >(dim_c,map<pair<int,int>,double>()); //prev error matr

    e_vec = new vector<double>(dim_e,0.);
    f_vec = new vector<double>(dim_e,0.);
    e_gbin_map = new vector<vector<int>>(dim_e,vector<int>(3,-1));
    
    c_vec = new vector<double>(dim_c,0.);
    n_vect = new vector<neighbor>(dim_c);
    c_gbin_map = new vector<vector<int>>(dim_c,vector<int>(3,-1));
    c_gbin_next = new vector<vector<int>>(dim_c,vector<int>(0));
    eff_vec = new vector<double>(dim_c,0.);
    sum_vec = new vector<double>(dim_c,0.);
  }

  void SetObservedData(TH1 *obs_h, TH1 *true_h);
  void SetObservedData(TH2 *obs_h, TH2 *true_h);
  void SetObservedData(TH3 *obs_h, TH3 *true_h);

  void SetNeighbors(TH1 *obs_h, TH1 *true_h);
  void SetNeighbors(TH2 *obs_h, TH2 *true_h);
  void SetNeighbors(TH3 *obs_h, TH3 *true_h);

  int GetGbin(int, int, int, TH3 *);
  int GetGbin(int, int, TH2 *);
    
  void FillMatrix(int, int, bool);
  //  void SetResponse(TMatrixD *mat);
  void FinalizeMatrix();
  void Finalize_P_EC();
  void SetResponse(map< pair<int,int>, double>);
    
  void InitPrior();
  void Init_No();
  void Get_Nprime();
  void UpdateNo();
  void UpdatePrior();
  void SmoothPrior();
  void SmoothPriorByNeighbor();
  void GetMeasured(TH2 *, TH2 *);
  void PlotReco(TH1 *);
  void PlotReco(TH2 *);
  void PlotReco(TH3 *);
  void GetFVec();
  void SmoothByKernel();


  void Get_dn_dn();
  void Get_dn_dP();
  void GetCovariance();
  double Get_dn_dn_prev(int,int);
  
  TGraph PlotNprime();
  //  double GetEntry(int row, int col){return response_m[row][col];};

  //  TMatrixD*  GetTMatrixD(){return response_m;};
  vector<double> * GetEffVec(){return eff_vec;};
  
};


#endif
