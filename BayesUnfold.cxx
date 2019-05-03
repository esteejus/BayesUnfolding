#include "BayesUnfold.h"





int BayesUnfold::GetGbin(int i , int j, int k, TH3 *hist)
{
  return hist -> GetNbinsZ() * hist -> GetNbinsY() * i + k * hist -> GetNbinsY() + j;
}

int BayesUnfold::GetGbin(int i , int j, TH2 *hist)
{
  return hist -> GetNbinsY() * i + j;
}

void BayesUnfold::SetObservedData(TH1 *obs_h,TH1 *true_h)
{
  
  for(int i = 0; i < obs_h -> GetNbinsX(); i++)
    {
      e_vec->at(i) = obs_h -> GetBinContent(i+1);
      e_gbin_map->at(i).at(0)  = i;//Observed index x
    }

  SetNeighbors(obs_h,true_h);
  
  return;
}

void BayesUnfold::SetObservedData(TH2 *obs_h, TH2 *true_h)
{
  
  for(int i = 0; i < obs_h -> GetNbinsX(); i++)
    {
      for(int j = 0; j < obs_h -> GetNbinsY(); j++)
	{
	  int gbin = GetGbin(i,j,obs_h);
	  e_vec->at(gbin) = obs_h -> GetBinContent(i+1,j+1);
	  cout<<"E_vec "<<gbin<<" "<<e_vec->at(gbin)<<endl;
	  e_gbin_map->at(gbin).at(0) = i;//Observed index x
	  e_gbin_map->at(gbin).at(1) = j;//Observed index y
	}
    }

  for(int i = 0; i < true_h -> GetNbinsX(); i++)
    {
      for(int j = 0; j < true_h -> GetNbinsY(); j++)
	{
	  int gbin = GetGbin(i,j,true_h);
	  c_gbin_map->at(gbin).at(0) = i;//Observed index x
	  c_gbin_map->at(gbin).at(1) = j;//Observed index y
	}
    }

  SetNeighbors(obs_h,true_h);
  
  return;
}

void BayesUnfold::SetObservedData(TH3 *obs_h, TH3 *true_h)
{

  for(int i = 0; i < obs_h -> GetNbinsX(); i++)
    {
      for(int j = 0; j < obs_h -> GetNbinsY(); j++)
	{
	  for(int k = 0; k < obs_h -> GetNbinsZ(); k++)
	    {
	      int gbin = GetGbin(i,j,k,obs_h);
	      e_vec->at(gbin) = obs_h -> GetBinContent(i+1,j+1,k+1);
	      e_gbin_map->at(gbin).at(0) = i;//Observed index x
	      e_gbin_map->at(gbin).at(1) = j;//Observed index y
	      e_gbin_map->at(gbin).at(2) = k;//Observed index z
	    }
	}
    }

  for(int i = 0; i < true_h -> GetNbinsX(); i++)
    {
      for(int j = 0; j < true_h -> GetNbinsY(); j++)
	{
	  for(int k = 0; k < true_h -> GetNbinsZ(); k++)
	    {
	      int gbin = GetGbin(i,j,k,true_h);
	      c_gbin_map->at(gbin).at(0) = i;//Observed index x
	      c_gbin_map->at(gbin).at(1) = j;//Observed index y
	      c_gbin_map->at(gbin).at(2) = k;//Observed index z
	    }
	}
    }

  SetNeighbors(obs_h,true_h);

  return;
}

void BayesUnfold::SetNeighbors(TH1 *obs_h, TH1 *true_h)
{
  double sigma = .5; //pixels
  int dim = 2*ceil(sigma) + 1;

  TH1D *kernel = new TH1D("kernel","kernel",dim,0,dim);
  TF1 *gaus = new TF1("gaus","[2]*TMath::Gaus(x,[0],[1],true)",0,dim);
  gaus->SetParameters(dim/2.,sigma,1.);
  double norm_v = 1./gaus->Integral(0,dim);
  gaus->SetParameters(dim/2.,sigma,norm_v);

  for(int i = 1; i <= dim; i++)
    {
      double int_v = gaus->Integral(i-1,i);
      kernel -> SetBinContent(i,int_v);
    }
  
  for(int i = 0; i < true_h -> GetNbinsX(); i++)
    {
      int gbin = i;
      neighbor gbin_neighbor;

      for(int l = 0; l < dim; l++)
	{
	  int xbin = i + l - ceil(sigma);
	  if(xbin < true_h->GetNbinsX() && xbin >= 0)
	    {
	      int gbin1 = xbin;
	      gbin_neighbor.n_bins.push_back(gbin1);
	      gbin_neighbor.weight.push_back(kernel->GetBinContent(l+1));
	    }
	      
	}
      n_vect->at(gbin) = gbin_neighbor;
    }

  /*
  for(int i = 0; i < true_h -> GetNbinsX(); i++)
    {
      c_gbin_map->at(i).at(0)  = i;//cause index x
      if( i-1 > 0  &&  i+1 < true_h->GetNbinsX() )
	{
	  c_gbin_next->at(i).push_back(i - 1);
	  c_gbin_next->at(i).push_back(i + 1);
	}
    }
  */
  return;
}

void BayesUnfold::SetNeighbors(TH2 *obs_h, TH2 *true_h)
{
  double sigma = .15; //pixels
  int dim = 2*ceil(sigma) + 1;

  TH2D *kernel = new TH2D("kernel","kernel",dim,0,dim,dim,0,dim);
  TF2 *gaus = new TF2("gaus","[4]*TMath::Gaus(x,[0],[2],true) * TMath::Gaus(y,[1],[3],true)",0,dim,0,dim);
  gaus->SetParameters(dim/2.,dim/2.,sigma,sigma,1.);
  double norm_v = 1./gaus->Integral(0,dim,0,dim);
  gaus->SetParameters(dim/2.,dim/2.,sigma,sigma,norm_v);

  for(int i = 1; i <= dim; i++)
    {
      for(int j = 1; j <= dim; j++)
	{
	  double int_v = gaus->Integral(i-1,i,j-1,j);
	  kernel -> SetBinContent(i,j,int_v);
	}
    }
  
  for(int i = 0; i < true_h -> GetNbinsX(); i++)
    {
      for(int j = 0; j < true_h -> GetNbinsY(); j++)
	{
	      int gbin = GetGbin(i,j,true_h);
	      neighbor gbin_neighbor;

	      for(int l = 0; l < dim; l++)
		{
		  for(int m = 0; m < dim ; m++)
		    {
			  int xbin = i + l - ceil(sigma);
			  int ybin = j + m - ceil(sigma);
							   
			  if(xbin < true_h->GetNbinsX() && ybin < true_h->GetNbinsY() && xbin >= 0 && ybin >= 0)
			    {
			      int gbin1 = GetGbin(xbin,ybin,true_h);
			      gbin_neighbor.n_bins.push_back(gbin1);
			      gbin_neighbor.weight.push_back(kernel->GetBinContent(l+1, m+1));
			    }

		    }
		}
	      n_vect->at(gbin) = gbin_neighbor;
	}
    }

  return;

}

void BayesUnfold::SetNeighbors(TH3 *obs_h, TH3 *true_h)
{
  double sigma = 2; //pixels
  int dim = 2*ceil(sigma) + 1;
  cout<<"Dimension of Kernel is "<<dim<<endl;
  TH3D *kernel = new TH3D("kernel","kernel",dim,0,dim,dim,0,dim,dim,0,dim);
  TF3 *gaus = new TF3("gaus","[6]*TMath::Gaus(x,[0],[3],true) * TMath::Gaus(y,[1],[4],true) * TMath::Gaus(z,[2],[5],true)",0,dim,0,dim,0,dim);
  gaus->SetParameters(dim/2.,dim/2.,dim/2.,sigma,sigma,sigma,1.);
  double norm_v = 1./gaus->Integral(0,dim,0,dim,0,dim);
  gaus->SetParameters(dim/2.,dim/2.,dim/2.,sigma,sigma,sigma,norm_v);

  for(int i = 1; i <= dim; i++)
    {
      for(int j = 1; j <= dim; j++)
	{
	  for(int k = 1; k <= dim; k++)
	    {
	      double int_v = gaus->Integral(i-1,i,j-1,j,k-1,k);
	      kernel -> SetBinContent(i,j,k,int_v);
	    }
	}
    }
  
  for(int i = 0; i < true_h -> GetNbinsX(); i++)
    {
      for(int j = 0; j < true_h -> GetNbinsY(); j++)
	{
	  for(int k = 0; k < true_h -> GetNbinsZ(); k++)
	    {
	      int gbin = GetGbin(i,j,k,true_h);
	      neighbor gbin_neighbor;
	      //	      cout<<"Bin "<<i<<" "<<j<<" "<<k<<endl;
	      if(gbin == 100)
		cout<<"i j k "<<i<<" "<<j<<" "<<k<<endl;
	      for(int l = 0; l < dim; l++)
		{
		  for(int m = 0; m < dim ; m++)
		    {
		      for(int n = 0; n < dim ; n++)
			{
			  int xbin = i + l - ceil(sigma);
			  int ybin = j + m - ceil(sigma);
			  int zbin = k + n - ceil(sigma);
							   
			  if(xbin < true_h->GetNbinsX() && ybin < true_h->GetNbinsY() && zbin < true_h->GetNbinsZ() && xbin >= 0 && ybin >= 0 && zbin >= 0)
			    {
			      int gbin1 = GetGbin(xbin,ybin,zbin,true_h);
			      gbin_neighbor.n_bins.push_back(gbin1);
			      gbin_neighbor.weight.push_back(kernel->GetBinContent(l+1, m+1, n+1));
			      if(gbin == 100)
				{

				  //  cout<<"Neighbor is "<<gbin1<<" "<<xbin<<" "<<ybin<<" "<<zbin<<" weight "<<kernel->GetBinContent(l+1,m+1,n+1)<<endl;
				}
			      //			      cout<<"l m n is "<<l<<" "<<m<<" "<<n<<endl;
			      //			      cout<<"Neighbor is "<<xbin<<" "<<ybin<<" "<<zbin<<" weight "<<kernel->GetBinContent(l,m,n)<<endl;
			    }

			}
		    }
		}//End of kernel loop
	      //	      if(gbin == 100)
	      //		cout<<"Gbin size "<<gbin_neighbor.n_bins.size()<<endl;
	      
	      n_vect->at(gbin) = gbin_neighbor;
	    }
	}
    }


  return;
}

/*
void BayesUnfold::SetResponse(TMatrixD *mat)
{
  response_m = mat;
  
  for(int j = 0; j < response_m -> GetNcols(); j++)
    {
      for(int i = 0; i < response_m -> GetNrows(); i++)
	eff_vec->at(j) += (*response_m)(i,j);
    }

  return;
}
*/

void BayesUnfold::InitPrior()
{
  if(c_vec == nullptr)
    {
      cout<<"Cause vector not properly set"<<endl;
      return;
    }
  else if(e_vec==nullptr)
    {
      cout<<"Effect vector not properly set"<<endl;
      return;
    }

  p_vec = new vector<double>(c_vec->size(),0.);
  for(int i = 0; i < p_vec->size(); i++)
    p_vec->at(i) = 1./p_vec->size();

  return;
}

void BayesUnfold::Init_No()
{
  if(c_vec == nullptr)
    {
      cout<<"Cause vector not properly set"<<endl;
      return;
    }
  else if(e_vec==nullptr)
    {
      cout<<"Effect vector not properly set"<<endl;
      return;
    }
  else if(p_vec==nullptr)
    {
      cout<<"Prior not initialized"<<endl;
      return;
    }

  no_vec = new vector<double>(c_vec->size(),0.);
  n_prime = new vector<double>(c_vec->size(),0.);
  
  double nobs = 0; //observed value
  for(int i = 0; i < e_vec->size(); i++)
    nobs += e_vec->at(i);

  for(int i = 0; i < p_vec->size(); i++)
    {
      no_vec->at(i) = nobs * p_vec->at(i);
    }
  return;
}



void BayesUnfold::Finalize_P_EC()
{
  double thresh = 1e-6; //threshold for zero element

  //Notation is P(E_j|C_i)
  //j rows and i columns
  //Get the non zero matrix elements
  //P_EC_i is the non zero elements of a given column i
  //P_EC_j is the non zero elements of a given row j


  //find position of matrix element
  for(auto const& x : response_m)
    {
      int j = x.first.first;  //row index j
      int i = x.first.second; //column index i
      double matrix_el = x.second; //matrix value

      if(matrix_el > thresh)
	{
	  element idx_j, idx_i;
	  idx_i.idx = i;
	  idx_i.v = matrix_el;
	  idx_j.idx = j;
	  idx_j.v = matrix_el;

	  P_EC_i ->at(i).push_back(idx_j);
	  P_EC_j ->at(j).push_back(idx_i);
	}
      
    }

  return;
}

void BayesUnfold::GetFVec()
{
  for(int i = 0; i < f_vec->size(); i++)//row index
    {
      double norm = 0;
      for(int l = 0; l < p_vec->size(); l++)//C_l
	{

	  map< pair<int,int>,double>::iterator it = response_m.find(make_pair(i,l));
	  if( it != response_m.end())
	    norm += it->second * p_vec->at(l);
	  //	  norm += response_m[i][l] * p_vec->at(l);
	}

      f_vec->at(i) = norm;
    }

  return;
}

void BayesUnfold::Get_dn_dP()
{



  return;
}

void BayesUnfold::Get_dn_dn()
{
  dn_dn.erase(dn_dn.begin(),dn_dn.end());
  //P_EC_i is the non zero elements of a given column i
  //P_EC_j is the non zero elements of a given row j

      
  for(int j = 0; j < dim_e; j++)  //M_ij->at(i) non zero j
    {
      for(auto i_el : M_j->at(j))
	{
	  int i = i_el.idx;
	  dn_dn[make_pair(i,j)] = i_el.v; //M_ij
	}
    }

  for(int i = 0; i < dim_c; i++)  //M_ij->at(i) non zero j
    {
      for(auto j_el : dn_dn_prev_i->at(i))
	{
	  int j = j_el.idx;
	  dn_dn[make_pair(i,j)] += j_el.v;
	}
    }

  for(int k = 0; k < dim_e; k++)  //M_ij->at(i) non zero j
    {
      for(auto i_el : M_j->at(k))
	{
	  int i = i_el.idx;
	  double M_ik = i_el.v;

	  for(auto l_el : M_j->at(k))
	    {
	      int l = l_el.idx;
	      double M_lk = l_el.v;

	      for(auto j_el : dn_dn_prev_i->at(l))
		{
		  int j = j_el.idx;
		  double dno_dn = j_el.v;
		  if(no_vec->at(l_el.idx) > 1e-12)
		    dn_dn[make_pair(i,j)] -= (e_vec->at(k) * eff_vec->at(l_el.idx)/no_vec->at(l_el.idx))*M_lk * M_ik * dno_dn;
		}
	    }

	}
    }

  for(int i = 0; i < dn_dn_prev_i->size(); i++)
      dn_dn_prev_i->at(i).clear(); 

  for(int i = 0; i < dn_dn_prev_j->size(); i++)
      dn_dn_prev_j->at(i).clear(); 


  for(auto const& x : dn_dn)
    {
      int i = x.first.first; 
      int j = x.first.second;

      element idx_i;
      idx_i.idx = j;
      idx_i.v = x.second; //matrix value
      dn_dn_prev_i ->at(i).push_back(idx_i);

      element idx_j;
      idx_j.idx = i;
      idx_j.v = x.second; //matrix value
      dn_dn_prev_j ->at(j).push_back(idx_j);
  }
  cout<<"Dn dn size "<<dn_dn.size()<<endl;

  return;
}

double BayesUnfold::Get_dn_dn_prev(int i, int j)
{
  map< pair<int,int>,double>::iterator it = dn_dn_prev.find(make_pair(i,j));
  if( it != dn_dn_prev.end())
    return it->second;
  else
    return 0.;
  }

void BayesUnfold::GetCovariance()
{ 

  double n_true = 0;
  for(int i = 0; i < dim_c; i++)
   n_true += n_prime->at(i);


  for(int i = 0; i < dim_c; i++)
    {
      for(auto k_el : dn_dn_prev_j->at(i))
	{
	  int k = k_el.idx;

	  for(auto l_el : dn_dn_prev_j->at(i))
	    {
	      int l = l_el.idx;
	      double matrix_el = k_el.v * l_el.v * e_vec->at(i);
	      n_var[make_pair(k,l)] += matrix_el;
	      //	      if(k==l)
	      //		cout<<"matrix el  l is "<<l<<" "<<matrix_el<<endl;
	    }
	}
    }

  //  for(int i = 0; i < dim_c; i++)
  //    cout<<"covar "<<n_var[make_pair(i,i)]<<endl;

  /*
  TMatrixDSparse nn_var(dim_e,dim_e);
  TMatrixDSparse dn_dn_m(dim_c,dim_c);
 
  int row[dn_dn.size()];
  int col[dn_dn.size()];
  double value[dn_dn.size()];
    
  int count = 0;
  for(auto const&x : dn_dn)
    {
      row[count] = x.first.first; //C 
      col[count] = x.first.second;//E
      value[count] = x.second;
      count++;
    }
  dn_dn_m.SetMatrixArray(dn_dn.size(),row,col,value);

  int row2[dim_e];
  int col2[dim_e];
  double value2[dim_e];
    
  for(int i = 0; i < dim_e; i++)
    {
      row2[i] = i;
      col2[i] = i;
      value2[i] = e_vec->at(i);
    }
  nn_var.SetMatrixArray(dim_e,row2,col2,value2);

  var = dn_dn_m * nn_var ;
  */
  return;
}


void BayesUnfold::Get_Nprime()
{
  std::fill(f_vec->begin(),f_vec->end(),0.);
  std::fill(n_prime->begin(),n_prime->end(),0.);

  for(int i = 0; i < M_j->size(); i++)
    M_j->at(i).clear();

  for(int j = 0; j < e_vec->size(); j++)
    {
      for(auto el : P_EC_j->at(j))
	{
	  f_vec->at(j) += el.v * p_vec->at(el.idx);
	}
    }

  for(int i = 0; i < c_vec->size(); i++)
    {
      //looping over P_EC by col over the Effect dimension
      for(auto el : P_EC_i->at(i))
	{
	  //el.idx  is the Effect index
	  if(eff_vec->at(i)* f_vec -> at(el.idx) > 1e-14)
	    {
	      n_prime->at(i) += e_vec->at(el.idx) * el.v * p_vec->at(i) / (eff_vec->at(i)*f_vec->at(el.idx));
	      element idx_i;
	      idx_i.idx = i;
	      idx_i.v = el.v * p_vec->at(i) / (eff_vec->at(i)*f_vec->at(el.idx));
	      M_j ->at(el.idx).push_back(idx_i);
	    }
	  else
	    n_prime->at(i) = 0.;
	}
    }

  // Get_dn_dn();
  UpdatePrior();
  UpdateNo();
  
  return;
}

TGraph BayesUnfold::PlotNprime()
{
  TGraph g(n_prime->size());

  for(int i = 0; i < n_prime->size(); i++)
    g.SetPoint(i,i,n_prime->at(i));

  return g;
}


void BayesUnfold::SmoothByKernel()
{
  for(int gbin = 0; gbin < p_vec->size(); gbin++)
    {
      //      neighbor gbin_n = n_vect->at(gbin);
      double p_value = 0;
      if(gbin == 100)
	{
	  //	  cout<<"Prior before "<<p_vec->at(gbin)<<" "<<gbin<<endl;
	  //	  cout<<"Size of neighbors "<<n_vect->at(gbin).weight.size()<<endl;
	}
      for(int idx = 0; idx < n_vect->at(gbin).weight.size(); idx++)
	{
	  double p_neigh_v = p_vec->at(n_vect->at(gbin).n_bins.at(idx));
	  double weight_v = n_vect->at(gbin).weight.at(idx);
	  p_value += p_neigh_v * weight_v;
	  //	  if(gbin == 100)
	  //	    cout<<"Inside loop  gbin "<<n_vect->at(gbin).n_bins.at(idx)<<" value "<<p_neigh_v<<" weight "<<weight_v<<endl;
	}
      //      if(gbin == 100)
      //	cout<<"Prior after "<<p_value<<" "<<gbin<<endl;

      p_vec->at(gbin) = p_value;
    }

  return;
}

void BayesUnfold::SmoothPriorByNeighbor()
{
  vector<double> new_prior(p_vec->size(),0);

  for(int i = 0; i < p_vec->size(); i++)
    {
      double sum = 0;

      if(c_gbin_next->at(i).size() == 0)
	cout<<"Smoothing algorithm will not smooth neighbor vector size is 0!!!"<<endl;

      for(int j = 0; j < c_gbin_next->at(i).size(); j++)
	{
	  int gbin = c_gbin_next->at(i).at(j);
	  sum += .3*p_vec->at(gbin);
	}
      sum += p_vec->at(i);
      new_prior.at(i) = sum/(c_gbin_next->at(i).size()+1);
    }
  
  for(int i = 0; i < p_vec->size(); i++)
    p_vec->at(i) = new_prior.at(i);
  
  return;
}


void BayesUnfold::SmoothPrior()
{
  if(p_vec->size() < 3)
    {
      cout<<"Vector too small to smooth" <<endl;
      return;
    }
  
  double sum = 0;
  for(int i = 0; i < p_vec->size(); i++)
    {
      if(i == 0)
	p_vec->at(0) = (p_vec->at(0) +  p_vec->at(1) + p_vec->at(2))/3.;
      else if( i == p_vec->size() - 1)
	p_vec->at(i) = (p_vec->at(i-2) + p_vec->at(i-1) + p_vec->at(i))/3.;
      else
	p_vec->at(i) = (p_vec->at(i-1) + p_vec->at(i) + p_vec->at(i+1))/3.;

      sum += p_vec->at(i);
    }

  for(int i = 0; i < p_vec->size(); i++)
    p_vec->at(i) /= sum;

  return;
}

void BayesUnfold::UpdateNo()
{
  for(int i = 0; i < no_vec->size(); i++)
    no_vec->at(i) = n_prime->at(i);

  return;
}

void BayesUnfold::UpdatePrior()
{

  double n_true = 0;
  for(int i = 0; i < n_prime->size(); i++)
    n_true += n_prime->at(i);

  cout<<"N_true "<<n_true<<endl;
  for(int i = 0; i < n_prime->size(); i++)
    p_vec->at(i) = n_prime->at(i)/n_true;


  SmoothByKernel();
  SmoothByKernel();
  SmoothByKernel();
  //   SmoothPriorByNeighbor();
 
  return;
}


void BayesUnfold::GetMeasured(TH2 *true_h, TH2 *measured)
{
  if(c_gbin_map == nullptr || e_gbin_map==nullptr)
    cout<<"Set up true and measured histograms first"<<endl;

  for(int i = 0; i < e_gbin_map->size(); i++)
    {

      for(int j = 0; j < c_gbin_map->size(); j++)
	{

	  int cbin_i = c_gbin_map->at(j).at(0) + 1;
	  int cbin_j = c_gbin_map->at(j).at(1) + 1;

	  int ebin_i = e_gbin_map->at(i).at(0) + 1;
	  int ebin_j = e_gbin_map->at(i).at(1) + 1;

	  double cause_v = true_h -> GetBinContent(cbin_i,cbin_j);
	  double measured_v =  0;//cause_v * response_m[i][j];
	  measured->SetBinContent(ebin_i,ebin_j,measured_v + measured->GetBinContent(ebin_i,ebin_j));
	}
    }
      
      
  return;
}

void BayesUnfold::PlotReco(TH1 *reco)
{
  map< pair<int,int>,double>::iterator it;
  double sum = 0;
  for(int j = 0; j < n_prime->size(); j++)
    {
      sum += n_prime->at(j);
      reco->SetBinContent(j+1,n_prime->at(j));
      it = n_var.find(make_pair(j,j));
      if( it != n_var.end())
      	{
      	  reco->SetBinError(j,sqrt(it->second));
	  //	  cout<<"error is "<<sqrt(it->second)<<endl;
	}

    }

  reco->SetEntries(sum);
  
return;
}

void BayesUnfold::PlotReco(TH2 *reco)
{
  map< pair<int,int>,double>::iterator it;
  double sum = 0;
  for(int j = 0; j < c_gbin_map->size(); j++)
    {

      int cbin_i = c_gbin_map->at(j).at(0) + 1;
      int cbin_j = c_gbin_map->at(j).at(1) + 1;
      sum += n_prime->at(j);
      cout<<"Reco "<<cbin_i<<" "<<cbin_j<<" "<<n_prime->at(j)<<endl;
      it = n_var.find(make_pair(j,j));
      if( it != n_var.end())
      	{
      	  reco->SetBinError(cbin_i,cbin_j,sqrt(it->second));
	  //	  cout<<"error is "<<sqrt(it->second)<<endl;
	}

      reco->SetBinContent(cbin_i,cbin_j,n_prime->at(j));
    }

  reco->SetEntries(sum);
  
return;
}


void BayesUnfold::PlotReco(TH3 *reco)
{
  map< pair<int,int>,double>::iterator it;

  double sum = 0;
  for(int j = 0; j < c_gbin_map->size(); j++)
    {

      int cbin_i = c_gbin_map->at(j).at(0) + 1;
      int cbin_j = c_gbin_map->at(j).at(1) + 1;
      int cbin_k = c_gbin_map->at(j).at(2) + 1;
      sum += n_prime->at(j);
      reco->SetBinContent(cbin_i,cbin_j,cbin_k,n_prime->at(j));
      //      cout<<j<<endl;
      //      cout<<"Error "<<n_var[make_pair(j,j)]<<endl;
      it = n_var.find(make_pair(j,j));
      if( it != n_var.end())
      	{
      	  reco->SetBinError(cbin_i,cbin_j,cbin_k,sqrt(it->second));
	  //	  cout<<"error is "<<sqrt(it->second)<<endl;
	}

    }
      

  reco->SetEntries(sum);
  
return;
}


void BayesUnfold::FillMatrix(int e_i, int c_j, bool measured)
{
  if(!(c_j < dim_c && e_i < dim_e && e_i >=0 && c_j >=0))
    {
      //      cout<<"Fill Element ("<<e_i<<","<<c_j<<") out of range. Size of matrix ";
      //      cout<< dim_e <<" X "<< dim_c <<endl;
      return;
    }

  if(measured)
    {
      response_m[make_pair(e_i,c_j)]++;
      //      response_m[e_i][c_j]++;
      sum_vec->at(c_j)++;
    }
  else
    sum_vec->at(c_j)++;

  return; 
}

void BayesUnfold::SetResponse(map< pair<int,int>, double> mat)
{
  response_m = mat;
  for(auto &x : response_m)
    {
      int i = x.first.first;  //row index E dim 
      int j = x.first.second; //column index C dim
      eff_vec->at(j) += x.second;
    }
  
  return;
}

void BayesUnfold::FinalizeMatrix()
{
  if(sum_vec->size() != dim_c)
    cout<<"Sum vector not the same dim as the response matrix "<<endl;

  for(auto &x : response_m)
    {

      int j = x.first.first;  //row index j
      int i = x.first.second; //column index i
      //      double matrix_el = x.second; //matrix value
      x.second /= sum_vec->at(i); 
      eff_vec->at(i) += x.second;
  }

  //  map< pair<int,int>,double>::iterator it;
  /*
  for(int j = 0; j < dim_c; j++)
    {
      for(int i = 0; i < dim_e; i++)
	{
	  if(sum_vec->at(j) > 0)
	    {
	      //	      cout<<i<<" "<<j<<endl;
	      //	      response_m[i][j] /= sum_vec->at(j);
	      //find position of matrix element
	      it = response_m.find(make_pair(i,j));
	      if( it != response_m.end())
		it -> second /= sum_vec->at(j); 
	    }
	  
	  it = response_m.find(make_pair(i,j));
	  if( it != response_m.end())
	    eff_vec->at(j) += it -> second;
	  //		eff_vec->at(j) += response_m[i][j];
	}
    }
  */
  
  return;
}

