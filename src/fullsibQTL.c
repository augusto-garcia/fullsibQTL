/**********************************************************************
 * 
 * fullsibQTL.c
 *
 * copyright (c) 2011, Rodrigo Gazaffi (fullsibQTL)
 *
 * C functions for the fullsibQTL package
 *
 * Licensed under the GNU General Public License version 2 (June, 1991)
 *
 * This files has a series of functions written in C to perform 
 * QTL mapping in outcross context. The considered approach is 
 * Interval Mapping considering mixture models, under EM algorithm 
 * presented in Gazaffi (2009, PHD thesis) and Gazaffi et al. (in prep)
 *
 * This files has some functions that was first present in R/qtl package
 * written by Karl Broman and Hao Wu:
 * reorg_genoprob, allocate_double, allocate_dmatrix, reorg_errlod

 * Function written by Rodrigo Gazaffi:
 * R_scan_qtl, scan_qtl, H0test, em, constrat used, restr_effects, 
 * R_char_qtl, char_qtl, qtl_3effects, qtl_2effects, qtl_1effects 
 *
 * 'Luiz de Queiroz' College of Agriculture 
 * Departament of Genetics - Sao Paulo - Brazil
 * Contact: rgazaffi@gmail.com
 * First version: Oct-06-2011
 * Last update: Oct-06-2011
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include "fullsibQTL.h"

#define colin0 0
#define colin1 1
#define colin2 2
#define colin3 3
#define colin4 4
#define colin5 5 
#define colin6 6 
#define colin7 7 
#define colin8 8
#define colin9 9
#define colin10 10 
#define colin11 11
#define colin12 12
#define restr0 0
#define restr1 1
#define restr2 2
#define restr3 3
#define restr123 123
#define restr12 12
#define means4 4

#define effects3 3
#define effects2 2
#define effects1 1

#define verb0 0

/**********************************************************************
 * 
 * R_scan_qtl
 * 
 * Wrapper function for call from R;
 * it calls scan_qtl
 *
 * function to scan the genome for QTL mapping. IM approach is used
 * under maximum likelihood approach - using EM algorithm - here is
 * possible to whether use or not covariables to control some fixed
 * effects. This function is also used in CIM approach, controlling
 * the markers used as cofactor.
 *
 **********************************************************************/

void R_scan_qtl(int *n_ind, int *n_pos, int *n_gen, int *colin_tp, 
		double *genoprob, double *addcov, int *nc_addcov,
		double *gamma, double *addcov_ls, double *pheno, 
		double *result, int *maxit, double *tol, int *verbose)
{
  double ***Genoprob, **Addcov, **Addcov_ls;

  /* reorganize vector for matrices */
  reorg_genoprob(*n_ind, *n_pos, means4, genoprob, &Genoprob); 
  reorg_errlod(*n_ind, *nc_addcov, addcov, &Addcov);
  reorg_errlod(*n_ind, *nc_addcov, addcov_ls, &Addcov_ls);

  /* do the QTL mapping scan */
  scan_qtl(*n_ind, *n_pos, n_gen, Genoprob, pheno, Addcov, *nc_addcov, 
	   gamma, Addcov_ls, result, *maxit, *tol,  colin_tp, *verbose);
}

/**********************************************************************
 * scan_qtl
 *
 * function to scan the genome using IM approach with or without 
 * considering covariables. The covariables are fixed effects in the 
 * linear model
 * 
 * n_ind      Number of individuals
 *
 * n_pos      Number of marker positions - position on the "genome walk"
 *
 * n_gen      Number of different genotypes (can be 2, 3 or 4)
 *            Here it is a vector containing the information along L.G.
 *            in other functions can be just a number instead of vector
 *
 * Genoprob   Array of conditional genotype probs [n_gen][n_pos][n_ind]
 *
 * pheno      Phenotype data, as a vector
 *
 * Addcov     Additive covariates - can be marker used as cofactor or 
 *            any other source of variation present in phenotype
 *
 * nc_Addcov  Number of columns in Addcov
 *
 * gamma      Coefficients of the multiple linear regression parameters
 *            used as starting values in EM - i do in R and import for C
 *            (for now just thinking in practice not in speed yet!)
 *
 * Addcov_ls  (Addcov'Addcov)^{-1}Addcov': least square matrix for Addcov
 *
 * result     Upon exit, the log_likelihood
 * 
 * maxit      Maximum number of iterations in the EM algorithm
 *
 * tol        Tolerance for determining convergence in EM
 * 
 * colin_tp   Colinearity Type for a a single position:
 *            Here it is a vector containing the information along L.G.
 *
 * verbose    if 1, it prints the number of interations,
 *            the genetic effects estimatives and sqrt(variance) 
 *            after EM has converged.
 * 
 **********************************************************************/

void scan_qtl(int n_ind, int n_pos, int *n_gen, double ***Genoprob,
	      double *pheno,  double **Addcov, int nc_Addcov, 
	      double *gamma, double **Addcov_ls, double *result, 
	      int maxit, double tol, int *colin_tp, int verbose)
{
  int i,j,k, hypot=0, restric;
  double **Geno1pos, *effects, H0, *gamma2work;
  int p;

  allocate_dmatrix(means4, n_ind, &Geno1pos);
  allocate_double(nc_Addcov, &gamma2work);  
  allocate_double(effects3 + 1, &effects);  

  /* effects can have until 3 effects, i allocate the maximum and 
     just control the used position with for()     */

  
  for(i=0; i<n_pos; i++) { /* loop over marker positions */
    
    /*setting conditional probs for i-th position*/
    for (j=0; j<n_ind; j++){
      for (k = 0; k < means4; k++){
	Geno1pos[k][j] = Genoprob[k][i][j];
      }
    }

    /*Attention: gamma and come from R but i need to update gamma in EM 
      I pass the values of gamma, so i keep R information safe. */
    for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
    
    
    if(verbose == 1) Rprintf("\nHa: details for %d\n", i+1);

    result[i] = em(n_ind, n_gen[i], Geno1pos, pheno, Addcov, nc_Addcov,
		   gamma2work, Addcov_ls, maxit, tol, effects, hypot, 
		   colin_tp[i], verbose);
    
    /*again updating for H0*/
    for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];

    
    /* model under H0 */
    restric = H0test(colin_tp[i]);
    H0 = em(n_ind, n_gen[i], Geno1pos, pheno, Addcov, nc_Addcov, 
	    gamma2work, Addcov_ls, maxit, tol, effects, restric, 
	    colin_tp[i], verb0);


    if(verbose == 1) {
      Rprintf("\nln-likel for Ha\t%4.8f", result[i]);
      Rprintf("\nln-likel for H0\t%4.8f\n", H0);
    }
    
    /* test */
    result[i] -= H0;
  } /* end for cM */
}

/**********************************************************************
 * H0test
 *
 * function to determine the restriction model will be considered in EM
 * if 3 effects, effects[0] = effects[1] = effects[2] = 0
 * if 2 effects, effects[0] = effects[1] = 0 
 * if 1 effect,  effects[0] = 0
 **********************************************************************/

int H0test (int colin_tp){
  int restric;
  switch(colin_tp){
  case 0: restric = restr123;break;
  case 1: restric = restr1;  break;
  case 2: restric = restr1;  break;
  case 3: restric = restr12; break;
  case 4: restric = restr12; break;
  case 5: restric = restr12; break;
  case 6: restric = restr12; break;
  case 7: restric = restr12; break;
  case 8: restric = restr12; break;
  case 9: restric = restr1;  break;
  case 10: restric= restr1;  break;
  case 11: restric= restr1;  break;
  case 12: restric= restr1;  break;
  }
  return restric;
}

/**********************************************************************
 * em
 *
 * function to perform EM algorithm
 **********************************************************************/

double em(int n_ind, int n_gen, double **Geno1pos, double *pheno, 
	  double **Addcov, int nc_Addcov, double *gamma, 
	  double **Addcov_ls, int maxit, double tol, double *effects,
	  int hypot, int colin_tp, int verbose) 
{
  int j, k, k1,k2, s, p, w, verb;
  double s1,s2, s3, old_lik, new_lik, conv=1000, sigma=0.0;
  double  *cont_ef, *r, *E_hat;
  double *temp, *temp2;
  double **contr, *den_post, *PIcontrE, *y_Ac_gam, **bigPI, **M, **V;
  double *aux1, *aux2, *aux3, *aux4, *aux5;
  
  /*alocating variables */
  allocate_double(means4, &cont_ef); /* receives result: Contr %*% Eff */
  allocate_double(n_ind, &aux1); 
  
  allocate_double(means4, &aux2); 
  allocate_double(means4, &aux3); 
  allocate_double(n_ind, &aux4); 
  allocate_double(means4, &aux5); 


 
  /* aux vectors with n x 1 dim*/
  allocate_double(n_ind, &den_post); /* For bigP11 + bigP12, bigP21, bigP22  */
  allocate_double(n_ind, &PIcontrE); /* For PI %*% contr %*% E  */
  allocate_double(n_ind, &y_Ac_gam); /* For y - Addcov %*% gamma*/

 
  allocate_dmatrix(means4, n_ind, &bigPI); /* where i place E pwjfwj / sum_w pwjfwj */

  allocate_double(means4, &temp); 
  allocate_double(means4, &temp2); 
  allocate_double(n_gen - 1, &r);
  allocate_double(n_gen - 1, &E_hat);

  /*matrix M and V are defined in my estimators, i know it is not smart using matrix in C, 
    but i prefer to do this, for didatic reason */
  allocate_dmatrix(3, 3, &M);
  allocate_dmatrix(3, 3, &V);  

  /* Defining the contrasts to be used in that position*/
  allocate_dmatrix(means4, n_gen - 1, &contr);
  contrast_used (colin_tp, contr);


  /*
  for(k = 0; k < (n_gen -1); k++){
    for(w = 0; w < means4; w++){
      Rprintf("%1.0f\t",contr[w][k]);
    }
    Rprintf("\n");
  }
  */


  /*starting values:*/
  /*variance*/
  s1 = s2 = 0.0;
  
  for(j=0; j<n_ind; j++) {
    s1 += pheno[j];
    s2 += pheno[j]*pheno[j];
  }
  s1 = (s1*s1) / ((double) n_ind);
  sigma = sqrt((s2 - s1)/ ((double) n_ind));

  /*y - Zgama*/
  for(j=0; j<n_ind; j++) {
    s1 = 0.0;
    for(p = 0; p < nc_Addcov; p++) s1 += Addcov[p][j] * gamma[p];
    y_Ac_gam[j] = pheno[j] -  s1;
  }

  /*effects*/
  for (k = 0; k < (n_gen - 1); k++) effects[k]=0.0;
  
  /* E step */

  for(w = 0; w < means4; w++){
    cont_ef[w] = 0.0;
    for (k = 0; k < (n_gen -1); k++) cont_ef[w] += contr[w][k] * effects[k]; 
    for (j = 0; j < n_ind; j++){
      aux1[j] = y_Ac_gam[j] - cont_ef[w];
      bigPI[w][j] = Geno1pos[w][j] * dnorm(aux1[j], (double) 0.0, sigma, 0);
    }
  }
  
  /*bigPI contains the posterior probability*/
  for(j = 0; j < n_ind; j++){
    den_post[j] = 0.0;
    for (w = 0; w < means4; w++) den_post[j] += bigPI[w][j]; 
    for (w = 0; w < means4; w++) bigPI[w][j] =  bigPI[w][j] / den_post[j]; 
  }
  
  old_lik = 0.0;
  for(j = 0; j < n_ind; j++){
    old_lik += log(den_post[j]);
  }
  
  /* verbose to print name of results that will be printed*/
  if(verbose == 1){
    Rprintf("\niter \t old_likeli \t new_likeli \t difference \t convergence \t mean\t");
    for(verb = 0; verb < (n_gen - 1); verb++) 
      Rprintf("\t Effect %d", verb+1);
    Rprintf("\t sqrt(sigma2)\n", sigma);
  }

  /*to count the number of iteration*/
  s = 0;
  
  while (conv > tol) {
    R_CheckUserInterrupt(); /* check for ^C */
    
    for(k = 0; k < (n_gen -1); k++){

      /* Starting M step */
      /*calculating r[k]*/
      s1 = 0.0;
      for (j = 0; j < n_ind; j++){
	aux1[j] = 0.0;
	for(w = 0; w < means4; w++) aux1[j] += (bigPI[w][j]*contr[w][k]);
	s1 += (y_Ac_gam[j] * aux1[j]);
      }
      
      s2 = 0.0;
      for(w = 0; w < means4; w++) aux2[w] = contr[w][k] * contr[w][k]; /*D_k # D_k */
      for (j = 0; j < n_ind; j++){
	aux1[j] = 0.0;
	for(w = 0; w < means4; w++) aux1[j] +=bigPI[w][j]*aux2[w]; /*\sum =PI * (D_k # D_k) */
	s2 +=aux1[j];
      }
      r[k] = s1 / s2;
    }
    
    /* calculating M */
    for(k1 = 0; k1 < (n_gen - 1); k1++) {
      for(k2 = 0; k2 < (n_gen - 1); k2++) {
	M[k1][k2] = 0.0;
	s1 = s2 = 0.0;
	if(k1 != k2){
	  for(w = 0; w < means4; w++){
	    temp[w] = contr[w][k1] * contr[w][k2];
	    temp2[w] = contr[w][k1] * contr[w][k1];
	  }
	  for (j = 0; j < n_ind; j++){
	    for(w = 0; w < means4; w++){
	      s1 +=  (bigPI[w][j] * temp[w]);
	      s2 +=  (bigPI[w][j] * temp2[w]);
	    }
	  }
	  M[k1][k2] = (s1 / s2);
	}
      }
    }
    
    /* new estimatives E(t+1) = r[k] - M%*%E[t] */
    for(k1 = 0; k1 < (n_gen - 1); k1++) {
      s1 = 0.0;
      for(k2 = 0; k2 < (n_gen - 1); k2++) {
	s1 += (M[k1][k2] * effects[k2]);
      }
      E_hat[k1] = r[k1] - s1;
    }
    
    for(k = 0; k < (n_gen -1); k++){
      effects[k] = E_hat[k];
    }
    
    /* restrictions for marginal tests*/
    if(hypot != 0)
      restr_effects(n_gen, hypot, effects);
   

    /* preliminary step to update gamma*/
    for (w = 0; w < means4; w++){
      cont_ef[w] = 0.0;
      for (k = 0; k < (n_gen -1); k++) cont_ef[w] += (contr[w][k] * effects[k]); 
    }

    for(j=0; j<n_ind; j++) {
      PIcontrE[j] = 0.0;
      for (w = 0; w < means4; w++){
	PIcontrE[j] += bigPI[w][j] * cont_ef[w];
      }
    }
  
    
    for(j=0; j<n_ind; j++) aux4[j] = pheno[j] - PIcontrE[j];    
    /* gamma updated*/
    for(p = 0; p < nc_Addcov; p++){    
      gamma[p] = 0.0;
      for(j=0; j<n_ind; j++) gamma[p] += Addcov_ls[p][j] * aux4[j];
    }

    
    /* sigma2 */
    /* y - Zgamma*/
    for(j=0; j<n_ind; j++) {
      s1 = 0.0;
      for(p = 0; p < nc_Addcov; p++) s1 += Addcov[p][j] * gamma[p]; 
      y_Ac_gam[j] = pheno[j] -  s1;
    }

    /* PI*D*theta */
    for(j=0; j<n_ind; j++) {
      PIcontrE[j] = 0.0;
      for (w = 0; w < means4; w++){
	PIcontrE[j] += bigPI[w][j] * cont_ef[w];
      }
    }
  
    /* 1'PI */
    for(w=0; w <means4; w++){
      aux3[w] = 0.0;
      for (j=0; j<n_ind; j++){
	aux3[w] += bigPI[w][j];
      }
    }

    s1 = s2 = 0.0;
    for(j = 0; j < n_ind; j++){
      s1 += y_Ac_gam[j]*y_Ac_gam[j]; /*(y - Z gamma)'(y - Z gamma)*/
      s2 += y_Ac_gam[j]*PIcontrE[j]; /*(y - Z gamma)'(PI * D * theta)*/
      y_Ac_gam[j] = pheno[j];
      for(p = 0; p < nc_Addcov; p++) y_Ac_gam[j] -= Addcov[p][j]*gamma[p];
    }
    
    /* obtaining V matrix and s3 (theta' * V * theta) */
    s3 = 0.0;
    for(k1 = 0; k1 < (n_gen - 1); k1++) {
      for(k2 = 0; k2 < (n_gen - 1); k2++) {
	V[k1][k2] = 0.0;
	for(w = 0; w < means4; w++) aux5[w] = (contr[w][k1]*contr[w][k2]);
	for(w = 0; w < means4; w++) V[k1][k2] += aux3[w] * aux5[w];
	s3 += effects[k1]*V[k1][k2]*effects[k2];
      }
    }

    /* sigma2 updated - M step finished*/
    sigma = sqrt( (s1 - (2*s2) + s3) / ((double) n_ind));
  
    /* E estep again */
    for(w = 0; w < means4; w++){
      cont_ef[w] = 0.0;
      for (k = 0; k < (n_gen -1); k++) cont_ef[w] += contr[w][k] * effects[k];
      for (j = 0; j < n_ind; j++){
	bigPI[w][j] = exp(log(Geno1pos[w][j]) + log(dnorm(y_Ac_gam[j] - cont_ef[w], 0, sigma, 0)));
      }
    }
    
    for(j = 0; j < n_ind; j++){
      den_post[j] = 0.0;
      for (w = 0; w < means4; w++) den_post[j] += bigPI[w][j]; 
      for (w = 0; w < means4; w++) bigPI[w][j] = exp(log(bigPI[w][j]) - log(den_post[j]));
    }

    new_lik = 0.0;
    for(j = 0; j < n_ind; j++){
      new_lik += log(den_post[j]);
    }

    conv = fabs((new_lik - old_lik) / old_lik);
    //conv = (new_lik - old_lik) / fabs(old_lik);
    // conv = (new_lik - old_lik);

    if(verbose == 1){
      Rprintf("%d\t%2.8f\t%2.8f\t%2.8f\t%2.8f", s+1, old_lik, new_lik, (new_lik - old_lik), conv);
      Rprintf("\t%f", gamma[0]); 
      for(verb = 0; verb < (n_gen - 1); verb++) 
	Rprintf("\t%f", effects[verb]);
      Rprintf("\t%f\n", sigma);
    }

   
    if(s > maxit){
      warning("Didn't converge!");
      conv = 0.0; /* instead of break... */
    }
    else
      s++;
    /* if likelihood ratio between old and new is lower than a threshold or
       iteration limit is achieved stops the algorithm */

    old_lik = new_lik;
  } /*end-while loop to obtain MLE*/


  effects[n_gen] = sigma; /* make sure that sigma is here */ 

  return new_lik;

}


//problemas nessa funcao, acho que o jeito de preencher que deve ser melhorado...
void contrast_used (int colin_tp, double **contr){

  switch (colin_tp){
  case 0:
    //alpha_p
    contr[0][0] = contr[1][0] = 1.0;
    contr[2][0] = contr[3][0] = -1.0; 
    //alpha_q
    contr[0][1] = contr[2][1] = 1.0;
    contr[1][1] = contr[3][1] = -1.0; 
    //delta_pq
    contr[0][2] = contr[3][2] = 1.0;
    contr[1][2] = contr[2][2] = -1.0;
    break;
  case 1:
    //alpha_p
    contr[0][0] = contr[1][0] = 1.0;
    contr[2][0] = contr[3][0] = -1.0; 
    break;

  case 2:
    //alpha_q
    contr[0][0] = contr[2][0] = 1.0;
    contr[1][0] = contr[3][0] = -1.0; 
    break;

  case 3:
    //alpha_p
    contr[0][0] = contr[1][0] = 1.0;
    contr[2][0] = contr[3][0] = -1.0; 
    //alpha_q
    contr[0][1] = contr[1][1] = 0.0;
    contr[2][1] = 1.0;
    contr[3][1] = -1.0; 
    
    break;
  case 4:
    //alpha_p
    contr[0][0] = contr[1][0] = 1.0;
    contr[2][0] = contr[3][0] = -1.0; 
    //alpha_q
    contr[0][1] = 1.0;
    contr[1][1] = -1.0;
    contr[2][1] = contr[3][1] = 0.0;
    
    break;
  case 5:
    //alpha_p
    contr[0][0] = contr[2][0] = 0.0;
    contr[1][0] = 1.0;
    contr[3][0] = -1.0;
    //alpha_q
    contr[0][1] = contr[2][1] = 1.0;
    contr[1][1] = contr[3][1] = -1.0; 
    
    break;
  case 6:
    //alpha_p
    contr[0][0] = 1.0;
    contr[1][0] = contr[3][0] = 0.0;
    contr[2][0] = -1.0;
    //alpha_q
    contr[0][1] = contr[2][1] = 1.0;
    contr[1][1] = contr[3][1] = -1.0; 
    
    break;
  case 7:
    //alpha_p
    contr[0][0] = contr[3][0] = 0.0;
    contr[1][0] =  1.0;
    contr[2][0] = -1.0;
    //delta_pq
    contr[0][1] = contr[3][1] = 1.0;
    contr[1][1] = contr[2][1] = -1.0; 
    
    break;
  case 8:
    //alpha_p
    contr[0][0] =  1.0;
    contr[1][0] = contr[2][0] = 0.0;
    contr[3][0] = -1.0;
    //delta_pq
    contr[0][1] = contr[3][1] = 1.0;
    contr[1][1] = contr[2][1] = -1.0; 
    break;
  case 9:
    //alpha_p
    contr[0][0] = contr[1][0] = contr[2][0] = 1/3;
    contr[3][0] = -1;
      break;
  case 10:
    //alpha_p
    contr[0][0] = contr[1][0] = contr[3][0] = 1/3;
    contr[2][0] = -1;
      break;
  case 11:
    //alpha_p
    contr[0][0] = contr[2][0] = contr[3][0] = 1/3;
    contr[1][0] = -1;
      break;
  case 12:
    //alpha_p
    contr[1][0] = contr[2][0] = contr[3][0] = 1/3;
    contr[0][0] = -1;
      break;
  }
}


void restr_effects(int n_gen, int hypot, double *effects)
{
  switch(n_gen){
  case 4: /* 4 differents genotypes */
    switch(hypot){
    case 1: 
      effects[0] = 0.0;
      break;
    case 2:
      effects[1] = 0.0;
      break;
    case 3:
      effects[2] = 0.0;
      break;
    case 123:
      effects[0] = effects[1] = effects[2] = 0.0;
      break;
    }
    break;
    
  case 3: /* 3 differents genotypes */
    switch(hypot){
    case 1: 
      effects[0] = 0.0;
      break;
    case 2:
      effects[1] = 0.0;
      break;
    case 12:
      effects[0] = effects[1] = 0.0;
      break;
    }
    break;
 
  case 2: /* 2 differents genotypes */
    switch(hypot){
    case 1: 
      effects[0] = 0.0;
      break;
    }
    break;
  }
}





/**********************************************************************
 * 
 * reorg_genoprob
 *
 * Reorganize the genotype probability data so that it is a triply 
 * indexed array rather than a single long vector
 *
 * Afterwards, genoprob indexed like Genoprob[gen][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_genoprob(int n_ind, int n_pos, int n_gen, 
		    double *genoprob, double ****Genoprob)
{
  int i, j;
  double **a;

  *Genoprob = (double ***)R_alloc(n_gen, sizeof(double **));

  a = (double **)R_alloc(n_pos*n_gen, sizeof(double *));

  (*Genoprob)[0] = a;
  for(i=1; i< n_gen; i++) 
    (*Genoprob)[i] = (*Genoprob)[i-1]+n_pos;
  
  for(i=0; i<n_gen; i++) 
    for(j=0; j<n_pos; j++) 
      (*Genoprob)[i][j] = genoprob + i*n_ind*n_pos + j*n_ind;
}

/**********************************************************************
 * 
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_double(int n, double **vector)
{
  *vector = (double *)R_alloc(n, sizeof(double));
}

/**********************************************************************
 * 
 * allocate_dmatrix
 *
 * Allocate space for a matrix of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_dmatrix(int n_row, int n_col, double ***matrix)
{
  int i;

  *matrix = (double **)R_alloc(n_row, sizeof(double *));
  
  (*matrix)[0] = (double *)R_alloc(n_col*n_row, sizeof(double));

  for(i=1; i<n_row; i++) 
    (*matrix)[i] = (*matrix)[i-1]+n_col;
}

/**********************************************************************
 * 
 * reorg_errlod
 *
 * Just like reorg_geno(), only for a matrix of doubles.
 *
 * Afterwards, errlod indexed like Errlod[mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void reorg_errlod(int n_ind, int n_mar, double *errlod, double ***Errlod)
{
  int i;

  *Errlod = (double **)R_alloc(n_mar, sizeof(double *));

  (*Errlod)[0] = errlod;
  for(i=1; i< n_mar; i++) 
    (*Errlod)[i] = (*Errlod)[i-1] + n_ind;
}



/**********************************************************************
 * 
 * R_char_qtl
 * 
 * Wrapper function for call from R;
 * it calls char_qtl to characterize the QTL...
 **********************************************************************/

void R_char_qtl(int *n_ind, int *n_gen, int *colin_tp, 
		double *genoprob, double *addcov, int *nc_Addcov, 
		double *gamma, double *addcov_ls, double *pheno,
		double *result, int *maxit, double *tol, int *verbose)
{

  double **Genoprob, **Addcov, **Addcov_ls;

  reorg_errlod(*n_ind, means4, genoprob, &Genoprob);
  reorg_errlod(*n_ind, *nc_Addcov, addcov, &Addcov);
  reorg_errlod(*n_ind, *nc_Addcov, addcov_ls, &Addcov_ls);

  char_qtl(*n_ind, *n_gen, Genoprob, pheno, Addcov, *nc_Addcov,
	   gamma, Addcov_ls, result, *maxit, *tol, *colin_tp, *verbose);

}


/**********************************************************************
 * char_qtl
 *
 * function to characterize the mapped QTL, e.g., for a given position 
 * genetic effects are obtained and tested, giving enough information 
 * infer QTL segregation and its linkage phase with flanking markers.
 * 
 * n_ind      Number of individuals
 *
 * n_gen      Number of different genotypes (can be 2, 3 or 4)
 *
 * colin_tp   Colinearity Type for a a single position:
 *
 * Genoprob   Array of conditional genotype probs [n_gen][n_pos][n_ind]
 *
 * pheno      Phenotype data, as a vector
 *
 * weights    Vector of positive weights, of length n_ind
 *
 * result     Upon exit, contains all QTL characterization:
 *            Global LOD, mu, alpha_p, LOD, alpha_q, LOD,
 *            delta_pq, LOD, H4_LOD, H5_LOD, H6_LOD, used model
 * 
 * maxit     Maximum number of iterations in the EM algorithm
 * 
 * tol       Tolerance for determining convergence in EM
 * 
 * verbose   if 1, it prints the number of interations,
 *           the genetic effects estimatives and sqrt(variance) 
 *           after EM has converged.
 * 
 **********************************************************************/

void char_qtl(int n_ind, int n_gen, double **Genoprob,
	      double *pheno,  double **Addcov, int nc_Addcov, 
	      double *gamma, double **Addcov_ls, double *result, 
	      int maxit, double tol, int colin_tp, int verbose) 
{
  switch(n_gen){
  case 4:
    qtl_3effects(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
		 gamma, Addcov_ls, result, maxit, tol, colin_tp, verbose);
    break;
  case 3:
    qtl_2effects(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
		 gamma, Addcov_ls, result, maxit, tol, colin_tp, verbose);
    break;   
  case 2:
    qtl_1effects(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
		 gamma, Addcov_ls, result, maxit, tol, colin_tp, verbose);
    break;
  }

}


/**********************************************************************
 * qtl_3effects
 *
 * function to characterize the QTL when 3 contrasts are estimables
 * 
 * n_ind      Number of individuals
 *
 * n_gen      Number of different genotypes (here is 4)
 *
 * colin_tp   Colinearity Type for a a single position:
 *
 * Genoprob   Array of conditional genotype probs [n_gen][n_pos][n_ind]
 *
 * pheno      Phenotype data, as a vector
 *
 * weights    Vector of positive weights, of length n_ind
 *
 * result     Upon exit, contains all QTL characterization:
 *            Global LOD, mu, alpha_p, LOD, alpha_q, LOD,
 *            delta_pq, LOD, H4_LOD, H5_LOD, H6_LOD, used model
 * 
 * maxit     Maximum number of iterations in the EM algorithm
 *
 * tol       Tolerance for determining convergence in EM
 * 
 * verbose   if 1, it prints the number of interations,
 *           the genetic effects estimatives and sqrt(variance) 
 *           after EM has converged.
 * 
 **********************************************************************/
 
void qtl_3effects(int n_ind, int n_gen, double **Genoprob,
		  double *pheno,  double **Addcov, int nc_Addcov, 
		  double *gamma, double **Addcov_ls, double *result, 
		  int maxit, double tol, int colin_tp, int verbose)
{
  double Ha=0.0, H0=0.0, H1=0.0, H2=0.0, H3=0.0, H4=0.0, H5=0.0, H6=0.0;
  double mu=0.0, alpha_p=0.0, alpha_q=0.0, delta_pq=0.0;
  double s1=0.0, s2=0.0;
  double *effects;
  double *gamma2work;
  int i, p, j;

  allocate_double(effects3 + 1, &effects); /* + 1 to storage sigma */
  allocate_double(nc_Addcov, &gamma2work);  

  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];


  /* model considering QTL without any restrictions */
  if(verbose == 1) Rprintf("\nHere, 8 Hypothesis are considered for Ha and H0; H1, H2 and H3; H4, H5 and H6\n");

  if(verbose == 1) Rprintf("\nDetails for Ha\n");
  Ha = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr0, colin_tp, verbose);
  

  /* getting the estimatives */
  mu = gamma[0];
  alpha_p = effects[0];
  alpha_q = effects[1];
  delta_pq = effects[2];

  if(verbose == 1) Rprintf("\nDetails for H0:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H0 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr123, colin_tp, verbose);
 
  /* testing models */
  if(verbose == 1) Rprintf("\nDetails for H1:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H1 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr1, colin_tp, verbose);

  if(verbose == 1) Rprintf("\nDetails for H2:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H2 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr2, colin_tp, verbose);

  if(verbose == 1) Rprintf("\nDetails for H3:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H3 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr3, colin_tp, verbose);
 
   /* testing if effects are similar or not */
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  if(verbose == 1) Rprintf("\nDetails for H4:\n");
  if(fabs(alpha_p)/alpha_p == fabs(alpha_q)/alpha_q){
    H4 = em(n_ind, n_gen - 1, Genoprob, pheno, Addcov, nc_Addcov,
	    gamma2work, Addcov_ls, maxit, tol, effects,
	    restr0, colin8, verbose);
  }
  else{
    H4 = em(n_ind, n_gen - 1, Genoprob, pheno, Addcov, nc_Addcov,
	    gamma2work, Addcov_ls, maxit, tol, effects,
	    restr0, colin7, verbose);
  }


  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  if(verbose == 1) Rprintf("\nDetails for H5:\n");
  if(fabs(alpha_p)/alpha_p == fabs(delta_pq)/delta_pq){
    H5 = em(n_ind, n_gen - 1, Genoprob, pheno, Addcov, nc_Addcov,
	    gamma2work, Addcov_ls, maxit, tol, effects,
	    restr0, colin6, verbose);
  }
  else{
    H5 = em(n_ind, n_gen - 1, Genoprob, pheno, Addcov, nc_Addcov,
	    gamma2work, Addcov_ls, maxit, tol, effects,
	    restr0, colin5, verbose);
  }

  if(verbose == 1) Rprintf("\nDetails for H6:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  if(fabs(alpha_q)/alpha_q == fabs(delta_pq)/delta_pq){
    H6 = em(n_ind, n_gen - 1, Genoprob, pheno, Addcov, nc_Addcov,
	    gamma2work, Addcov_ls, maxit, tol, effects,
	    restr0, colin4, verbose);
  }
  else{
    H6 = em(n_ind, n_gen - 1, Genoprob, pheno, Addcov, nc_Addcov,
	    gamma2work, Addcov_ls, maxit, tol, effects,
	    restr0, colin3, verbose);
  }
  
  /* Preparing the vector that contains all the results */
  result[ 0] = (Ha - H0)/log(10); /* Global LOD */
    
  result[ 1] = mu;
  result[ 2] = alpha_p;
  result[ 3] = (Ha - H1)/log(10); /* test for alpha_p == 0 */
  result[ 4] = alpha_q;
  result[ 5] = (Ha - H2)/log(10); /* test for alpha_q == 0 */
  result[ 6] = delta_pq;
  result[ 7] = (Ha - H3)/log(10); /* test for delta_pq == 0 */

  result[ 8] = (Ha - H4)/log(10);
  result[ 9] = (Ha - H5)/log(10);
  result[10] = (Ha - H6)/log(10);

  result[11] = colin_tp; /* show the used model */
}

/**********************************************************************
 * qtl_2effects
 *
 * function to characterize the QTL in cases with 3 genotypic classes,
 * e.g., in situations that is possible to estimate 2 contrasts.
 * 
 * n_ind      Number of individuals
 *
 * n_gen      Number of different genotypes (here is 3)
 *
 * colin_tp   Colinearity Type for a a single position:
 *
 * Genoprob   Array of conditional genotype probs [n_gen][n_pos][n_ind]
 *
 * pheno      Phenotype data, as a vector
 *
 * weights    Vector of positive weights, of length n_ind
 *
 * result     Upon exit, contains all QTL characterization:
 *            Global LOD,  mu, alpha_p, LOD, alpha_q, LOD, delta_pq, 
 *            LOD, H4_LOD (or H5_LOD, if markers are B3) and used model
 * 
 * maxit     Maximum number of iterations in the EM algorithm
 *
 * tol       Tolerance for determining convergence in EM
 * 
 * verbose   if 1, it prints the number of interations,
 *           the genetic effects estimatives and sqrt(variance) 
 *           after EM has converged.
 * 
 **********************************************************************/

void qtl_2effects(int n_ind, int n_gen, double **Genoprob,
		  double *pheno,  double **Addcov, int nc_Addcov, 
		  double *gamma, double **Addcov_ls, double *result, 
		  int maxit, double tol, int colin_tp, int verbose)
{
  int colin, p,j;
  double Ha, H0, H1, H2, H4, mu, alpha_p, alpha_q;
  double s1=0.0, s2=0.0;
  double *effects;
  double *gamma2work;

  allocate_double(effects2 + 1, &effects); /* + 1 to storage sigma */
  allocate_double(nc_Addcov, &gamma2work);  

  if(verbose == 1) Rprintf("\nHere, 5 Hypothesis are considered for Ha and H0; H1, H2 and H4:\n");

  if(verbose == 1) Rprintf("\nDetails for Ha:\n");

  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  Ha = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr0, colin_tp, verbose);
  
  /* getting the estimatives */
  mu = gamma[0];
  alpha_p = effects[0];
  alpha_q = effects[1];

  if(verbose == 1) Rprintf("\nDetails for H0:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H0 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr12, colin_tp, verbose);

  /* testing models */
  if(verbose == 1) Rprintf("\nDetails for H1:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H1 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr1, colin_tp, verbose);

  if(verbose == 1) Rprintf("\nDetails for H2:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H2 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr2, colin_tp, verbose);

  /* check colin type and storage the correct colin for H4 */
  switch(colin_tp){
  case colin3:
    if(fabs(alpha_p)/alpha_p == fabs(alpha_q)/alpha_q)
      colin = colin9;
    else 
      colin = colin10;
    break;
  case colin4:
    if(fabs(alpha_p)/alpha_p == fabs(alpha_q)/alpha_q)
      colin = colin12;
    else 
      colin = colin11;
    break;

  case colin5:
    if(fabs(alpha_p)/alpha_p == fabs(alpha_q)/alpha_q)
      colin = colin9;
    else 
      colin = colin11;
    break;
  case colin6:
    if(fabs(alpha_p)/alpha_p == fabs(alpha_q)/alpha_q)
      colin = colin12;
    else 
      colin = colin10;
    break;

  case colin7:
    if(fabs(alpha_p)/alpha_p == fabs(alpha_q)/alpha_q)
      colin = colin10;
    else 
      colin = colin11;
    break;
  case colin8:
    if(fabs(alpha_p)/alpha_p == fabs(alpha_q)/alpha_q)
      colin = colin12;
    else 
      colin = colin9;
    break;
  }

  /* testing if effects are similar or not */
  if(verbose == 1) Rprintf("\nDetails for H4:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H4 = em(n_ind, n_gen - 1, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr0, colin, verbose);

  /* Preparing the vector that contains all the results */
  if ((colin_tp == 7) || (colin_tp == 8)){ /* it is like a F2 */
    result[ 0] = (Ha - H0)/log(10);  /* Global LOD */

    result[ 1] = mu;
    result[ 2] = alpha_p;
    result[ 3] = (Ha - H1)/log(10); /* test for alpha_p == 0 */
    result[ 4] = NA_REAL;
    result[ 5] = NA_REAL; 
    result[ 6] = alpha_q;
    result[ 7] = (Ha - H2)/log(10); /* test for delta_pq == 0 */
    result[ 8] = NA_REAL;
    result[ 9] = (Ha - H4)/log(10); /* should be H5, but in just that case */
    result[10] = NA_REAL;
    result[11] = colin_tp; /* show the used model */  
  }
  else{
    result[ 0] = (Ha - H0)/log(10);  /* Global LOD */

    result[ 1] = mu;
    result[ 2] = alpha_p;
    result[ 3] = (Ha - H1)/log(10); /* test for alpha_p == 0 */
    result[ 4] = alpha_q;
    result[ 5] = (Ha - H2)/log(10); /* test for alpha_q == 0 */
    result[ 6] = NA_REAL;
    result[ 7] = NA_REAL; 
    result[ 8] = (Ha - H4)/log(10);
    result[ 9] = NA_REAL;
    result[10] = NA_REAL;
    result[11] = colin_tp; /* show the used model */
  }

}

/**********************************************************************
 * qtl_2means
 *
 * function to characterize the QTL in cases with 2 genotypic classes,
 * e.g., in situations that is possible to estimate 1 contrasts.
 * 
 * n_ind      Number of individuals
 *
 * n_gen      Number of different genotypes (here is 2)
 *
 * colin_tp   Colinearity Type for a a single position:
 *
 * Genoprob   Array of conditional genotype probs [n_gen][n_pos][n_ind]
 *
 * pheno      Phenotype data, as a vector
 *
 * Addcov   Additive covariates
 *
 * n_addcov Number of columns in Addcov
  * weights    Vector of positive weights, of length n_ind
 *
 * result     Upon exit, contains all QTL characterization:
 *            Global LOD, mu, alpha_p, LOD, alpha_q, LOD,
 *            delta_pq, LOD, H4_LOD, H5_LOD, H6_LOD, used model
 * 
 * maxit     Maximum number of iterations in the EM algorithm
 *
 * tol       Tolerance for determining convergence in EM
 * 
 * verbose   if 1, it prints the number of interations,
 *           the genetic effects estimatives and sqrt(variance) 
 *           after EM has converged.
 * 
 **********************************************************************/

void qtl_1effects(int n_ind, int n_gen, double **Genoprob,
		  double *pheno,  double **Addcov, int nc_Addcov, 
		  double *gamma, double **Addcov_ls, double *result, 
		  int maxit, double tol, int colin_tp, int verbose)
{
  int p, j;
  double Ha, H0, mu, alpha;
  double s1=0.0, s2=0.0;
  double *effects;
  double *gamma2work;
  
  allocate_double(effects2 + 1, &effects); /* + 1 to storage sigma */
  allocate_double(nc_Addcov, &gamma2work);  

  if(verbose == 1) Rprintf("\nHere, 2 Hypothesis are considered for Ha and H0:\n");

  if(verbose == 1) Rprintf("\nDetails for Ha:\n");

  /* model considering QTL without any restrictions */
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  Ha = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr0, colin_tp, verbose);
  
  /* getting the estimatives */
  mu = gamma[0];
  alpha = effects[0];

  if(verbose == 1) Rprintf("\nDetails for H0:\n");
  for (p = 0; p < nc_Addcov; p++) gamma2work[p] = gamma[p];
  H0 = em(n_ind, n_gen, Genoprob, pheno, Addcov, nc_Addcov,
	  gamma2work, Addcov_ls, maxit, tol, effects,
	  restr1, colin_tp, verbose);

  /* Preparing the vector that contains all the results */ 
  if (colin_tp == colin2){
    result[ 0] = (Ha - H0)/log(10); /* Global LOD */

    result[ 1] = mu;
    result[ 2] = NA_REAL;
    result[ 3] = NA_REAL;
    result[ 4] = alpha;
    result[ 5] = (Ha - H0)/log(10); /* test for alpha_q == 0 */
    result[ 6] = NA_REAL;
    result[ 7] = NA_REAL; 
    result[ 8] = NA_REAL;
    result[ 9] = NA_REAL;
    result[10] = NA_REAL;
    result[11] = colin_tp; /* show the used model */
  }
  else{
    result[ 0] = (Ha - H0)/log(10); /* Global LOD */

    result[ 1] = mu;
    result[ 2] = alpha;
    result[ 3] = (Ha - H0)/log(10); /* test for alpha_p == 0 */
    result[ 4] = NA_REAL;
    result[ 5] = NA_REAL;
    result[ 6] = NA_REAL;
    result[ 7] = NA_REAL; 
    result[ 8] = NA_REAL;
    result[ 9] = NA_REAL;
    result[10] = NA_REAL;
    result[11] = colin_tp; /* show the used model */
  }
}
