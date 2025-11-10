#ifndef FL_HELPERS_H
#define FL_HELPERS_H

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;
#define tolr 1e-5
#define crossprod(x) symmatu((x).t() * (x))
#define tcrossprod(x) symmatu((x) * (x).t())

/*

==================================================================================

This file contains all the general-purpose helper functions extracted from

your original 'DL_linear_split_merge_package.cpp' file.

These functions (math, random number generators, etc.) can be re-used by

both your original MCMC code and the new federated EM code.

==================================================================================
*/

// [[Rcpp::export]]
double c_lmvgamma (double x, int p) {
int i;
double ans = 0;
if (p < 1)
Rcpp::stop("p must be greater than or equal to 1.");
if (x <= 0)
Rcpp::stop("x must be greater than 0.");
ans =(p * (p - 1)/4.0) * M_LNPI;
for (i = 0; i < p; ++i){
ans +=  (lgamma(x  - (i/2.0) ));
}
return ans;
}

double log_marg_dens (const mat &chol_marginal_delta,int n,double t1ppp
,double  niw_nu,double  niw_kap,double  diag_psi_iw   ){
int k=chol_marginal_delta.n_cols;
double dens;
dens=-(k/2)(n M_LNPI+log1p(n/niw_kap)-niw_nu*log(diag_psi_iw)- (niw_nu+n)*log(t1ppp) )
+c_lmvgamma ( (niw_nu+n)/2 ,  k)- c_lmvgamma ( niw_nu/2 ,  k)- (niw_nu+n)*sum(log(chol_marginal_delta.diag()) );
return dens;
}

// [[Rcpp::export]]
double rig(double mu){
double y = randn<double>();
y = y;
double mu2 = gsl_pow_2(mu);
double quad = 4 * mu * y + mu2 * gsl_pow_2(y);
// double x = mu + y * mu2 / 2 - mu / 2  * sqrt(quad);
double  x_div_mu=(2+muy - sqrt(quad) )/2;
double  x = (mu* x_div_mu);
if(x<=0){
// cout<<"mu= "<<mu<<" y= "<<y<<endl;
// Rcpp::stop("Error at rig!!");
// x=0;
// cout<<"manual return at rig!!"<<endl;
return(1e-5);
}

double u = log (randu<double>());
// if(u <= (mu / (x + mu))) return x;
if(u <=  -log1p(x_div_mu) ) return x;
else return mu / x_div_mu;
}

double    _unur_bessel_k_nuasympt (double x, double nu, int islog, int expon_scaled)    {
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */

double z;                   /* rescaled argument for K_nu() /
double sz, t, t2, eta;      / auxiliary variables /
double d, u1t,u2t,u3t,u4t;  / (auxiliary) results for Debye polynomials /
double res;                 / value of log(K_nu(x)) [= result] */

/* rescale: we comute K_nu(z * nu) */
z = x / nu;

/* auxiliary variables /
sz = hypot(1,z);   / = sqrt(1+z^2) /
t = 1. / sz;
t2 = tt;

eta = (expon_scaled) ? (1./(z + sz)) : sz;
eta += log(z) - log1p(sz);                  /* = log(z/(1+sz)) */

/* evaluate Debye polynomials u_j(t) /
u1t = (t * (3. - 5.t2))/24.;
u2t = t2 * (81. + t2(-462. + t2 * 385.))/1152.;
u3t = tt2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
u4t = t2*t2 * (4465125.
+ t2 * (-94121676.
+ t2 * (349922430.
+ t2 * (-446185740.
+ t2 * 185910725.)))) / 39813120.;
d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;

               /* log(K_nu(x)) */
               res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);

               return (islog ? res : exp(res));


}

double _gig_mode(double lambda, double omega)
/---------------------------------------------------------------------------/
/* Compute mode of GIG distribution.                                         /
{
if (lambda >= 1.)
/ mode of fgig(x) /
return (sqrt((lambda-1.)(lambda-1.) + omegaomega)+(lambda-1.))/omega;
else
/ 0 <= lambda < 1: use mode of f(1/x) /
return omega / (sqrt((1.-lambda)(1.-lambda) + omegaomega)+(1.-lambda));
} / end of _gig_mode() */

void


_rgig_ROU_shift_alt (double res, int n, double lambda, double lambda_old, double omega, double alpha)
/---------------------------------------------------------------------------/
/ Type 8:                                                                   /
/ Ratio-of-uniforms with shift by 'mode', alternative implementation.       /
{
double xm, nc;     / location of mode; c=log(f(xm)) normalization constant /
double s, t;       / auxiliary variables /
double U, V, X;    / random variables */

int i;             /* loop variable (number of generated random variables) /
int count = 0;     / counter for total number of iterations */

double a, b, c;    /* coefficent of cubic /
double p, q;       / coefficents of depressed cubic /
double fi, fak;    / auxiliary results for Cardano's rule */

double y1, y2;     /* roots of (1/x)*sqrt(f((1/x)+m)) */

double uplus, uminus;  /* maximum and minimum of x*sqrt(f(x+m)) */

/* -- Setup -------------------------------------------------------------- */

/* shortcuts */
t = 0.5 * (lambda-1.);
s = 0.25 * omega;

/* mode = location of maximum of sqrt(f(x)) */
xm = _gig_mode(lambda, omega);

/* normalization constant: c = log(sqrt(f(xm))) /
nc = tlog(xm) - s*(xm + 1./xm);

/* location of minimum and maximum of (1/x)*sqrt(f(1/x+m)):  */

/* compute coeffients of cubic equation y^3+ay^2+by+c=0 /
a = -(2.(lambda+1.)/omega + xm);       /* < 0 /
b = (2.(lambda-1.)*xm/omega - 1.);
c = xm;

/* we need the roots in (0,xm) and (xm,inf) */

/* substitute y=z-a/3 for depressed cubic equation z^3+pz+q=0 /
p = b - aa/3.;
q = (2.aaa)/27. - (a*b)/3. + c;

/* use Cardano's rule /
fi = acos(-q/(2.sqrt(-(ppp)/27.)));
fak = 2.*sqrt(-p/3.);
y1 = fak * cos(fi/3.) - a/3.;
y2 = fak * cos(fi/3. + 4./3.*M_PI) - a/3.;

/* boundaries of minmal bounding rectangle:                  /
/ we us the "normalized" density f(x) / f(xm). hence        /
/ upper boundary: vmax = 1.                                 /
/ left hand boundary: uminus = (y2-xm) * sqrt(f(y2)) / sqrt(f(xm)) /
/ right hand boundary: uplus = (y1-xm) * sqrt(f(y1)) / sqrt(f(xm)) /
uplus  = (y1-xm) * exp(tlog(y1) - s*(y1 + 1./y1) - nc);
uminus = (y2-xm) * exp(tlog(y2) - s(y2 + 1./y2) - nc);

/* -- Generate sample ---------------------------------------------------- */

for (i=0; i<n; i++) {
do {
++count;
U = uminus + unif_rand() * (uplus - uminus);    /* U(u-,u+)  /
V = unif_rand();                                / U(0,vmax) /
X = U/V + xm;
}                                         / Acceptance/Rejection /
while ((X <= 0.) || ((log(V)) > (tlog(X) - s*(X + 1./X) - nc)));

/* store random point */
res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);


}

/* -- End ---------------------------------------------------------------- */

return;
} /* end of _rgig_ROU_shift_alt() */

void _rgig_newapproach1 (double res, int n, double lambda, double lambda_old, double omega, double alpha)
/---------------------------------------------------------------------------/
/ Type 4:                                                                   /
/ New approach, constant hat in log-concave part.                           /
{
/ parameters for hat function /
double A[3], Atot;  / area below hat /
double k0;          / maximum of PDF /
double k1, k2;      / multiplicative constant */

double xm;          /* location of mode */
double x0;          /* splitting point T-concave / T-convex */
double a;           /* auxiliary variable */

double U, V, X;     /* random numbers */
double hx;          /* hat at X */

int i;              /* loop variable (number of generated random variables) */
int count = 0;      /* counter for total number of iterations */

/* -- Check arguments ---------------------------------------------------- */

if (lambda >= 1. || omega >1.)
  Rcpp::stop ("invalid parameters");

/* -- Setup -------------------------------------------------------------- */

/* mode = location of maximum of sqrt(f(x)) */
xm = _gig_mode(lambda, omega);

/* splitting point */
x0 = omega/(1.-lambda);

/* domain [0, x_0] */
k0 = exp((lambda-1.)*log(xm) - 0.5*omega*(xm + 1./xm));     /* = f(xm) */
A[0] = k0 * x0;

/* domain [x_0, Infinity] */
if (x0 >= 2./omega) {
  k1 = 0.;
  A[1] = 0.;
  k2 = pow(x0, lambda-1.);
  A[2] = k2 * 2. * exp(-omega*x0/2.)/omega;
}

else {
  /* domain [x_0, 2/omega] */
  k1 = exp(-omega);
  A[1] = (lambda == 0.)
    ? k1 * log(2./(omega*omega))
      : k1 / lambda * ( pow(2./omega, lambda) - pow(x0, lambda) );

  /* domain [2/omega, Infinity] */
  k2 = pow(2/omega, lambda-1.);
  A[2] = k2 * 2 * exp(-1.)/omega;
}

/* total area */
Atot = A[0] + A[1] + A[2];

/* -- Generate sample ---------------------------------------------------- */

for (i=0; i<n; i++) {
  do {
    ++count;

    /* get uniform random number */
    V = Atot * unif_rand();

    do {

      /* domain [0, x_0] */
      if (V <= A[0]) {
        X = x0 * V / A[0];
        hx = k0;
        break;
      }

      /* domain [x_0, 2/omega] */
      V -= A[0];
      if (V <= A[1]) {
        if (lambda == 0.) {
          X = omega * exp(exp(omega)*V);
          hx = k1 / X;
        }
        else {
          X = pow(pow(x0, lambda) + (lambda / k1 * V), 1./lambda);
          hx = k1 * pow(X, lambda-1.);
        }
        break;
      }

      /* domain [max(x0,2/omega), Infinity] */
      V -= A[1];
      a = (x0 > 2./omega) ? x0 : 2./omega;
      X = -2./omega * log(exp(-omega/2. * a) - omega/(2.*k2) * V);
      hx = k2 * exp(-omega/2. * X);
      break;

    } while(0);

    /* accept or reject */
    U = unif_rand() * hx;

    if (log(U) <= (lambda-1.) * log(X) - omega/2. * (X+1./X)) {
      /* store random point */
      res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);
      break;
    }
  } while(1);

}

/* -- End ---------------------------------------------------------------- */

return;


} /* end of _rgig_newapproach1() */

void _rgig_ROU_noshift (double *res, int n, double lambda, double lambda_old, double omega, double alpha)


/---------------------------------------------------------------------------/
/* Tpye 1:                                                                   /
/ Ratio-of-uniforms without shift.                                          /
{
double xm, nc;     / location of mode; c=log(f(xm)) normalization constant /
double ym, um;     / location of maximum of x*sqrt(f(x)); umax of MBR /
double s, t;       / auxiliary variables /
double U, V, X;    / random variables */

int i;             /* loop variable (number of generated random variables) /
int count = 0;     / counter for total number of iterations */

/* -- Setup -------------------------------------------------------------- */

/* shortcuts */
t = 0.5 * (lambda-1.);
s = 0.25 * omega;

/* mode = location of maximum of sqrt(f(x)) */
xm = _gig_mode(lambda, omega);

/* normalization constant: c = log(sqrt(f(xm))) /
nc = tlog(xm) - s*(xm + 1./xm);

/* location of maximum of xsqrt(f(x)):           /
/ we need the positive root of                   /
/ omega/2y^2 - (lambda+1)y - omega/2 = 0    /
ym = ((lambda+1.) + sqrt((lambda+1.)(lambda+1.) + omegaomega))/omega;

/* boundaries of minmal bounding rectangle:                   /
/ we us the "normalized" density f(x) / f(xm). hence         /
/ upper boundary: vmax = 1.                                  /
/ left hand boundary: umin = 0.                              /
/ right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(xm)) /
um = exp(0.5(lambda+1.)log(ym) - s(ym + 1./ym) - nc);

/* -- Generate sample ---------------------------------------------------- */

for (i=0; i<n; i++) {
do {
++count;
U = um * unif_rand();        /* U(0,umax) /
V = unif_rand();             / U(0,vmax) /
X = U/V;
}                              / Acceptance/Rejection /
while (((log(V)) > (tlog(X) - s*(X + 1./X) - nc)));

/* store random point */
res[i] = (lambda_old < 0.) ? (alpha / X) : (alpha * X);


}

/* -- End ---------------------------------------------------------------- */

return;
} /* end of _rgig_ROU_noshift() */

#define ZTOL (DOUBLE_EPS10.0)
// [[Rcpp::export]]
/---------------------------------------------------------------------------/
double rgig( double lambda, double chi, double psi)
{
double omega, alpha;     / parameters of standard distribution /
// SEXP sexp_res;           / results */
double *res;
int i;
int n=1;

/* check sample size */
if (n<=0) {
Rcpp::stop("sample size 'n' must be positive integer.");
}

/* check GIG parameters: */
if ( !(R_FINITE(lambda) && R_FINITE(chi) && R_FINITE(psi)) ||
(chi <  0. || psi < 0)      ||
(chi == 0. && lambda <= 0.) ||
(psi == 0. && lambda >= 0.) ) {
// cout<<"lambda="<<lambda<<", chi="<<chi<<", psi"=psi<<endl;
printf("lambda= %lf, chi=%lf, psi= %lf\n",lambda, chi,psi);
Rcpp::stop("invalid parameters for GIG distribution!!");
}

/* allocate array for random sample */
// PROTECT(sexp_res = NEW_NUMERIC(n));
res = (double )malloc(nsizeof(double)) ;

if (chi < ZTOL) {
/* special cases which are basically Gamma and Inverse Gamma distribution */
if (lambda > 0.0) {
for (i=0; i<n; i++) res[i] = R::rgamma(lambda, 2.0/psi);
}
else {
for (i=0; i<n; i++) res[i] = 1.0/R::rgamma(-lambda, 2.0/psi);
}
}

else if (psi < ZTOL) {
/* special cases which are basically Gamma and Inverse Gamma distribution */
if (lambda > 0.0) {
for (i=0; i<n; i++) res[i] = 1.0/R::rgamma(lambda, 2.0/chi);
}
else {
for (i=0; i<n; i++) res[i] = R::rgamma(-lambda, 2.0/chi);
}

}

else {
double lambda_old = lambda;
if (lambda < 0.) lambda = -lambda;
alpha = sqrt(chi/psi);
omega = sqrt(psi*chi);

/* run generator */
do {
  if (lambda > 2. || omega > 3.) {
    /* Ratio-of-uniforms with shift by 'mode', alternative implementation */
    _rgig_ROU_shift_alt(res, n, lambda, lambda_old, omega, alpha);
    break;
  }

  if (lambda >= 1.-2.25*omega*omega || omega > 0.2) {
    /* Ratio-of-uniforms without shift */
    _rgig_ROU_noshift(res, n, lambda, lambda_old, omega, alpha);
    break;
  }

  if (lambda >= 0. && omega > 0.) {
    /* New approach, constant hat in log-concave part. */
    _rgig_newapproach1(res, n, lambda, lambda_old, omega, alpha);
    break;
  }

  /* else */
  Rcpp::stop("parameters must satisfy lambda>=0 and omega>0.");

} while (0);


}

/* return result /
// UNPROTECT(1);
double ret=res[0];
free(res);
return ret;
} / end of do_rgig() */

double log_sum_exp(const arma::vec &x) {
double maxVal= x.max();

double sum_exp=sum(exp(x-maxVal));
return log(sum_exp)+maxVal ;


}

inline double multiplier_fn(double niw_nu, double niw_kap, int k){
return (niw_kap +1) / ( niw_kap* (niw_nu-k+1 )  );
}

// [[Rcpp::export]]
inline double log_t_density(const double nu,const  vec &x,const vec &mu,const mat &lower_chol){
int k=x.n_elem;
double det_sig_half= sum(log(lower_chol.diag() ) );
vec resid=solve(trimatl ( lower_chol ) , x - mu );

double quad_form=dot(resid,resid)/nu;
// gsl_sf_lnpoch(a,b)=log (gamma(a+b)/gamma(a)  )
double density =gsl_sf_lnpoch(nu/2 , ((double)k)/2)-// lgamma((nu+k)/2 ) -lgamma(nu/2) -
  (k* log( (datum::pi)*nu) + (nu+k)*log1p(quad_form) )/2 -det_sig_half;

return density;


}

// [[Rcpp::export]]
inline double log_t_density_empty(const double nu,const double kap,const double diag_psi_iw, const vec &x){
int k=x.n_elem;
double t1ppp=(kap+1)/(kapnu);
/ adjustment_factor.diag()+= (diag_psi_iw* t1ppp);
mat lower_chol=chol(adjustment_factor,"lower");*/

double adjustment_factor = (diag_psi_iw* t1ppp);

double det_sig_half= k* log(adjustment_factor)/2;  //sum(log(lower_chol.diag() ) );
// vec resid= solve(trimatl ( lower_chol ) , x  );

double quad_form= dot(x,x)/  (adjustment_factor*nu);
// gsl_sf_lnpoch(a,b)=log (gamma(a+b)/gamma(a)  )
double density =gsl_sf_lnpoch(nu/2 , ((double)k)/2)-  //lgamma((nu+k)/2 ) -lgamma(nu/2) -
  (k* log((datum::pi)*nu) + (nu+k)*log1p(quad_form ) )/2 -det_sig_half;
// cout<<"t1ppp="<<t1ppp<<" det_sig_half"<<det_sig_half<<" quad_form="<<quad_form<<" density="<<density;
return density;


}

uvec std_setdiff(arma::uvec a, arma::uvec b) {

// std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
// std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
a=sort(a);b=sort(b);

std::vector<unsigned> out;

std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                    std::inserter(out, out.end()));

return arma::conv_to<arma::uvec>::from(out);


}

// Helper for log-likelihood of multivariate normal
inline double log_norm_pdf(const rowvec& x, const vec& mean, const mat& cov) {
int k = x.n_elem;
double log_det_cov = log(det(cov));
mat chol_cov = chol(cov, "lower");
vec resid = solve(trimatl(chol_cov), (x.t() - mean));
double quad_form = dot(resid, resid);

return -0.5 * (k * log(2 * datum::pi) + log_det_cov + quad_form);


}

#endif // FL_HELPERS_H