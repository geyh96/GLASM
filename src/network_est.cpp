

#include "math.h"
#include <vector>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::plugins(cpp11)]]

double thresholdl1(double x, double thr);

//' calculate the conditional correlation among the X's
//' @param X The data matrix for the X
//' @param betaMat the beta matrix from lasso approach and its variantes. the betaMat can be estimated from the function SPMBgraphsqrtKeh: the sqrt lasso mathod is deployed.
//' @param n the sample nize. 
//' @param d the dimension. 
//' @export
//[[Rcpp::export]]
Eigen::ArrayXXd calculate_corrhat(Eigen::Map<Eigen::ArrayXXd> X, Eigen::Map<Eigen::ArrayXXd> betaMat, int n, int d)
{
    // ArrayXXd betaMat = ArrayXXd::Random(d,d);
    // ArrayXXd X = ArrayXXd::Random(n,d);
    Eigen::ArrayXXd ErrorMat;
    ErrorMat.resize(n,d);
    for(int m=0;m<d;m++){
        ErrorMat.col(m)= X.col(m).array() - (X.matrix()*betaMat.matrix().col(m)).array();
    }

    Eigen::ArrayXd diagCovhat;
    diagCovhat.resize(d);
    
    for(int k=0;k<d;k++){
            diagCovhat[k]=(ErrorMat.col(k).array() * ErrorMat.col(k).array()).mean(); 
    }

    Eigen::ArrayXXd Corrhat;
    Corrhat.resize(d,d);
    
    for(int k=0;k<d;k++){
        for(int l=0;l<d;l++){
            Corrhat(k,l)=-1*(                                                              \
                    ErrorMat.col(k).array() * ErrorMat.col(l).array() +                   \
                    ErrorMat.col(k).array()*ErrorMat.col(k).array()*betaMat(k,l) +          \
                    ErrorMat.col(l).array()*ErrorMat.col(l).array()*betaMat(l,k) ).mean()/ (sqrt(diagCovhat[k]*diagCovhat[l])); 
                    //for diagonal
            if(k==l){
                Corrhat(k,k)=(ErrorMat.col(k).array() * ErrorMat.col(l).array()).mean()/ (sqrt(diagCovhat[k]*diagCovhat[l]));
            }
        }
    }

    return Corrhat;
}


//' calculate the Test statistics for the independence of X.
//' @param X The data matrix for the X
//' @param betaMat the beta matrix from lasso approach and its variantes. the betaMat can be estimated from the function SPMBgraphsqrtKeh: the sqrt lasso mathod is deployed.
//' @param n the sample nize. dim(X)
//' @param d the dimension. dim(X)
//' @export
//[[Rcpp::export]]
Eigen::ArrayXd calculate_Testhat(Eigen::Map<Eigen::ArrayXXd> X, Eigen::Map<Eigen::ArrayXXd> betaMat, int n, int d)
{
    // ArrayXXd betaMat = ArrayXXd::Random(d,d);
    // ArrayXXd X = ArrayXXd::Random(n,d);
    Eigen::ArrayXXd ErrorMat;
    ErrorMat.resize(n,d);
    for(int m=0;m<d;m++){
        ErrorMat.col(m)= X.col(m).array() - (X.matrix()*betaMat.matrix().col(m)).array();
    }

    Eigen::ArrayXd diagCovhat;
    diagCovhat.resize(d);
    
    for(int k=0;k<d;k++){
            diagCovhat[k]=(ErrorMat.col(k).array() * ErrorMat.col(k).array()).mean(); 
    }

    // Eigen::ArrayXXd Corrhat;
    // Corrhat.resize(d,d);
    
    // for(int k=0;k<d;k++){
    //     for(int l=0;l<d;l++){
    //         Corrhat(k,l)=-1*(                                                              
    //                 ErrorMat.col(k).array() * ErrorMat.col(l).array() +                   
    //                 ErrorMat.col(k).array()*ErrorMat.col(k).array()*betaMat(k,l) +          
    //                 ErrorMat.col(l).array()*ErrorMat.col(l).array()*betaMat(l,k) ).mean()/ (sqrt(diagCovhat[k]*diagCovhat[l])); 
    //     }
    // }
    
    Eigen::ArrayXd TestAll;
    int numTest=d*(d-1)/2;
    TestAll.resize(numTest);
    int m=0;
    for(int k=1;k<d;k++){
        for(int l=0;l<k;l++){
            TestAll[m]=-1*(                                                              \
                    ErrorMat.col(k).array() * ErrorMat.col(l).array() +                   \
                    ErrorMat.col(k).array()*ErrorMat.col(k).array()*betaMat(k,l) +          \
                    ErrorMat.col(l).array()*ErrorMat.col(l).array()*betaMat(l,k)).mean()/   \
                    (sqrt(diagCovhat[k]*diagCovhat[l])); 
        m=m+1;
        }
    }

    // Eigen::ArrayXXd Corrhat;
    // Corrhat.resize(d,d);
    // for(int k=0;k<d;k++){
    //     for(int l=0;l<d;l++){
    //         Corrhat(k,l)=Covhat(k,l)/sqrt(Covhat(k,k)*Covhat(l,l));
    //     }
    // }

    return TestAll;
}





//' sqrt lasso to calculate the betaMat for the X. actually part of the TIGER method
//' @param data The data matrix for the X
//' @param lambda sequence of tuning parameter.
//' @param nlambda the number of lambda sequence.
//' @param d the dimension of X. 
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd SPMBgraphsqrtKeh(Eigen::Map<Eigen::MatrixXd> data, NumericVector lambda, int nlambda, int d)
{
  
  // Eigen::ArrayXd Xb, r, grad, w1, Y, XX, gr;
  Eigen::ArrayXd  XX;
  Eigen::ArrayXXd X;
  Eigen::MatrixXd tmp_icov,beta_matrix;
  tmp_icov.resize(d, d);
  tmp_icov.setZero();
  beta_matrix.resize(d,d);
  beta_matrix.setZero();
  std::vector<Eigen::MatrixXd > tmp_icov_p;
  tmp_icov_p.clear();
  for(int i = 0; i < nlambda; i++)
    tmp_icov_p.push_back(tmp_icov);
  int n = data.rows();

  int maxdf = 0;
  maxdf = (n < d ? n : d);
  NumericVector x(d*maxdf*nlambda);
  IntegerVector col_cnz(d+1);
  IntegerVector row_idx(d*maxdf*nlambda);
  X = data;
  XX.resize(d); 
  for (int j = 0; j < d; j++)
    XX[j] = (X.col(j)*X.col(j)).sum()/n;//XX record the variance of d covariates.
  
  double prec = 1e-5;
  int max_iter = 1000;
  int num_relaxation_round = 3;
 
  
  #pragma omp parallel for
  for(int m=0;m<d;m++)
  {
    //cout<<m<<endl;
     int cnz = 0;
    Eigen::ArrayXd Xb, r, w1, Y, gr;
     Xb.resize(n);
    Xb.setZero();

    gr.resize(d);
    gr.setZero();
    w1.resize(d);
    w1.setZero();
    r.resize(n);
    r.setZero();
    Y = X.col(m);
    
    Eigen::ArrayXd Xb_master(n);
    Eigen::ArrayXd w1_master(n);
    std::vector<int> actset_indcat(d, 0);
    std::vector<int> actset_indcat_master(d, 0);
    std::vector<int> actset_idx;
    std::vector<double> old_coef(d, 0);
    std::vector<double> grad(d, 0);
    std::vector<double> grad_master(d, 0);
    
    double a = 0, g = 0, L = 0, sum_r2 = 0;
    double tmp_change = 0, local_change = 0;
    
    r = Y - Xb;
    sum_r2 = r.matrix().dot(r.matrix());
    L = sqrt(sum_r2 / n); // sqrt part of loss function
    
    double dev_thr = fabs(L) * prec; 
    
    //cout<<dev_thr<<endl;
    
    
    for(int i = 0; i < d; i++)   //calculate the gradient of residuals and each covariate. 
    {
      grad[i] = (r * X.col(i)).sum() / (n*L);
    }
    for(int i = 0; i < d; i++) gr[i] = abs(grad[i]); // abs of the gradient for each covariate.
    w1_master = w1; //at the begining, the w1 is zeros vector
    Xb_master = Xb;
    for (int i = 0; i < d; i++) grad_master[i] = gr[i];
    
    std::vector<double> stage_lambdas(d, 0);
    
    for(int i=0;i<nlambda;i++)
    {
      w1 = w1_master; //so for each lambda .all start from the _master value of gradient
      Xb = Xb_master;
      
      for (int j = 0; j < d; j++)
      {
        gr[j] = grad_master[j];
        actset_indcat[j] = actset_indcat_master[j];
      }
      
      // init the active set
      //choose the covariate with gradient > 2*lambda
      double threshold;
      if (i > 0)
        threshold = 2 * lambda[i] - lambda[i - 1];
      else
        threshold = 2 * lambda[i];
      
      for (int j = 0; j < m; ++j) // 在Y那一列之前的类，干了这么个事情。主要是为了避开Y这个第m列
      {
        stage_lambdas[j] = lambda[i];  //what? 把参数付给一个向量，这个是干什么
        
        if (gr[j] > threshold) actset_indcat[j] = 1; //梯度大于threshold 的 就是active set。标成1
      }
      for (int j = m+1; j < d; ++j) //在Y那一列后面的列，就赶了个这么个事情。主要是为了避开Y这个第m列
      {
        stage_lambdas[j] = lambda[i];
        
        if (gr[j] > threshold) actset_indcat[j] = 1;//梯度大于threshold 的 就是active set。标成1
      }
      stage_lambdas[m] = lambda[i];  //Y这一列也初始化一下
      r = Y - Xb;  
      sum_r2 = r.matrix().dot(r.matrix());
      L = sqrt(sum_r2 / n);  //继续算个loss出来
      // loop level 0: multistage convex relaxation
      int loopcnt_level_0 = 0;
      int idx;
      double old_w1=0, updated_coord=0;
      while (loopcnt_level_0 < num_relaxation_round) //int num_relaxation_round = 3;
        {
        loopcnt_level_0++;
        
        // loop level 1: active set update
        int loopcnt_level_1 = 0;
        bool terminate_loop_level_1 = true;
        
        while (loopcnt_level_1 < max_iter)
        {
          loopcnt_level_1++;
          terminate_loop_level_1 = true;
          
          for (int j = 0; j < d; j++) old_coef[j] = w1[j];
          
          // initialize actset_idx
          actset_idx.clear();
          for (int j = 0; j < m; j++)
            if (actset_indcat[j]) //如果是空模型里面选出来的变量的话
            {
              g = 0.0;
              a = 0.0;
              
              double tmp;
              
              sum_r2 = r.matrix().dot(r.matrix()); //Eigen::ArrayXd 
              L = sqrt(sum_r2 / n);
              
              Eigen::ArrayXd wXX  = (1 - r*r/sum_r2) * X.col(j) * X.col(j);
              g = (wXX * w1[j] + r * X.col(j)).sum()/(n*L);
              a = wXX.sum()/(n*L);
              
              tmp = w1[j];
              w1[j] = thresholdl1(g, stage_lambdas[j]) / a;
              
              tmp = w1[j] - tmp;
              // Xb += delta*X[idx*n]
              Xb = Xb + tmp * X.col(j);
              
              sum_r2 = 0.0;
              // r -= delta*X
              r = r - tmp * X.col(j);  
              
              sum_r2 = r.matrix().dot(r.matrix());
              L = sqrt(sum_r2 / n);
              
              updated_coord = w1[j];
              
              if (fabs(updated_coord) > 0) actset_idx.push_back(j);
            }
            
            for (int j = m+1; j < d; j++)
              if (actset_indcat[j])
              {
                g = 0.0;
                a = 0.0;
                
                double tmp;
                
                sum_r2 = r.matrix().dot(r.matrix());
                L = sqrt(sum_r2 / n);
                //这就是在更新了，但是为啥这么更新我不明白，要看算法
                Eigen::ArrayXd wXX  = (1 - r*r/sum_r2) * X.col(j) * X.col(j);
                g = (wXX * w1[j] + r * X.col(j)).sum()/(n*L);
                a = wXX.sum()/(n*L);
                
                tmp = w1[j];
                w1[j] = thresholdl1(g, stage_lambdas[j]) / a;
                
                tmp = w1[j] - tmp;
                // Xb += delta*X[idx*n]
                Xb = Xb + tmp * X.col(j);
                
                sum_r2 = 0.0;
                // r -= delta*X
                r = r - tmp * X.col(j);
                
                sum_r2 = r.matrix().dot(r.matrix());
                L = sqrt(sum_r2 / n);
                
                updated_coord = w1[j];
                
                if (fabs(updated_coord) > 0) actset_idx.push_back(j);
              }
              
              // loop level 2: proximal newton on active set
              int loopcnt_level_2 = 0;
              bool terminate_loop_level_2 = true;
              while (loopcnt_level_2 < max_iter)
              {
                loopcnt_level_2++;
                terminate_loop_level_2 = true;
                
                for (unsigned int k = 0; k < actset_idx.size(); k++)
                {
                  idx = actset_idx[k];
                  
                  old_w1 = w1[idx];
                  g = 0.0;
                  a = 0.0;
                  
                  double tmp;
                  
                  sum_r2 = r.matrix().dot(r.matrix());
                  L = sqrt(sum_r2 / n);
                  
                  Eigen::ArrayXd wXX  = (1 - r*r/sum_r2) * X.col(idx) * X.col(idx);
                  g = (wXX * w1[idx] + r * X.col(idx)).sum()/(n*L);
                  a = wXX.sum()/(n*L);
                  
                  tmp = w1[idx];
                  w1[idx] = thresholdl1(g, stage_lambdas[idx]) / a;
                  
                  tmp = w1[idx] - tmp;
                  // Xb += delta*X[idx*n]
                  Xb = Xb + tmp * X.col(idx);
                  
                  sum_r2 = 0.0;
                  // r -= delta*X
                  r = r - tmp * X.col(idx);
                  
                  sum_r2 = r.matrix().dot(r.matrix());
                  L = sqrt(sum_r2 / n);
                  
                  updated_coord = w1[idx];
                  tmp_change = old_w1 - w1[idx];
                  double a =  (X.col(idx) * X.col(idx) * (1 - r * r/(L*L*n))).sum()/(n*L);
                  local_change = a * tmp_change * tmp_change / (2 * L * n);
                  if (local_change > dev_thr)
                    terminate_loop_level_2 = false;
                }
                if (terminate_loop_level_2)
                  break;
              }
              
              terminate_loop_level_1 = true;
              // check stopping criterion 1: fvalue change
              for (unsigned int k = 0; k < actset_idx.size(); ++k)
              {
                idx = actset_idx[k];
                tmp_change = old_w1 - w1[idx];
                double a =  (X.col(idx) * X.col(idx) * (1 - r * r/(L*L*n))).sum()/(n*L);
                local_change = a * tmp_change * tmp_change / (2 * L * n);
                if (local_change > dev_thr)
                  terminate_loop_level_1 = false;
              }
              
              r = Y - Xb;
              sum_r2 = r.matrix().dot(r.matrix());
              L = sqrt(sum_r2 / n);
              
              if (terminate_loop_level_1)
                break;
              
              
              // check stopping criterion 2: active set change
              bool new_active_idx = false;
              for (int k = 0; k < m; k++)
                if (actset_indcat[k] == 0)
                {
                  grad[idx] = (r * X.col(idx)).sum() / (n*L);
                  //cout<<grad[idx];
                  gr[k] = fabs(grad[k]);
                  if (gr[k] > stage_lambdas[k])
                  {
                    actset_indcat[k] = 1;
                    new_active_idx = true;
                  }
                }
                for (int k = m+1; k < d; k++)
                  if (actset_indcat[k] == 0)
                  {
                    grad[idx] = (r * X.col(idx)).sum() / (n*L);
                    //cout<<grad[idx]
                    gr[k] = fabs(grad[k]);
                    if (gr[k] > stage_lambdas[k])
                    {
                      actset_indcat[k] = 1;
                      new_active_idx = true;
                    }
                  }
                  // 20210102
                  if(!new_active_idx)
                    break;
        }
        if (loopcnt_level_0 == 1)
        {
          for (int j = 0; j < d; j++)
          {
            w1_master[j] = w1[j];
            
            grad_master[j] = gr[j];
            actset_indcat_master[j] = actset_indcat[j];
          }
          
          for (int j = 0; j < n; j++) Xb_master[j] = Xb[j];
        }
        }
      for(unsigned int j=0;j<actset_idx.size();j++)
      {
        int w_idx = actset_idx[j];
        x[cnz] = w1[w_idx];
        row_idx[cnz] = i*d+w_idx;
        cnz++;
        //cout<<cnz<<"    ";
      }
      double tal = 0;
      Eigen::MatrixXd temp;
      temp.resize(n, 1);
      for(int j = 0; j < n; j++)
      {
        temp(j, 0) = 0;
        for(int k = 0; k < d; k++)
          temp(j, 0) += X.matrix()(j, k)*w1[k];
        temp(j, 0) = Y[j] - temp(j, 0);
      }
      //temp = Y.matrix() - X.matrix().transpose()*w1.matrix();
      for(int j = 0; j < n; j++)
        tal += temp(j, 0)*temp(j, 0);
      tal = sqrt(tal)/sqrt(n);
//#pragma omp critical
//{   
      //tmp_icov = tmp_icov_p[i];
      // if(0){
      // tmp_icov(m, m) = pow(tal, -2);
      // for(int j = 0; j < m; j++)
      //   tmp_icov(j, m) = -tmp_icov(m, m)*w1[j];
      // for(int j = m+1; j < d; j++)
      //   tmp_icov(j, m) = -tmp_icov(m, m)*w1[j];
      //   }
      //tmp_icov_p[i] = tmp_icov;
//  }
    for(int j=0;j<d;j++){
      beta_matrix(j,m)=w1[j];
    }
    col_cnz[m+1]=cnz;
}
  }
  
  //   tmp_icov = (tmp_icov.transpose()+tmp_icov)/2;
  // return tmp_icov;
  return beta_matrix;
}



double thresholdl1(double x, double thr) {
  if (x > thr)
    return x - thr;
  else if (x < -thr)
    return x + thr;
  else
    return 0;
}
