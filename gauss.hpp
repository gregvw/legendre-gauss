#include<vector>
#include<cmath>
#include<iostream>


#define pi 3.1415926535897932384 
 
typedef std::vector<double> dvec;

void printvector(const dvec& v) {
    int n = v.size();
    for(int i=0;i<n;++i) {
        std::cout << v[i] << std::endl;
    }  
} 



// Compute maximum pointwise difference
double dist_max(const dvec& a, const dvec& b){
    int n = a.size();
    double dist;
    double max_dist;

    for(int i=0; i<n; ++i) {
        dist = std::abs(a[i]-b[i]);
        if(dist>max_dist) {
            max_dist = dist; 
        }  
    }
    return max_dist;
}



// Compute Chebyshev Gauss Nodes
double chebnodes(dvec& x) {
    int n = x.size();

    for(int i=0;i<n;++i) {
        x[i] = cos(pi*double(2*i+1)/double(2*n));    
    } 
} 


double initguess(dvec& x) {
    int n = x.size();

    for(int i=0;i<n;++i) {
        x[i] = -cos(pi*double(i+.75)/double(n+0.5));    
    } 

    
}



// Compute Legendre polynomial recursion coefficients
void rec_legendre(dvec& a, dvec& b){
    int n = a.size();
    for(int i=0;i<n;++i) {
        a[i] = (2.0*i+1.0)/(i+1.0);
        b[i] = double(i)/(i+1.0);
    }
}

// Evaluate the nth order Legendre polynomial and its derivative
void legendre(const dvec& a, const dvec& b, const dvec& x, dvec& L, dvec& Lp) {
    int n = x.size();
    dvec L0(n,1.0);
    dvec L1(n,0.0);
    
    // Iterate over grid points
    for(int j=0;j<n;++j) {
        L1[j] = x[j];
        // Iterate over polynomials  
        for(int k=1;k<n;++k) {

            L[j] = a[k]*x[j]*L1[j]-b[k]*L0[j];
            L0[j] = L1[j];
            L1[j] = L[j];
        }
        Lp[j] = n*(L0[j]-x[j]*L[j])/(1.0-x[j]*x[j]); 
    } 
} 

// Update grid points as well as the nth Legendre polynomial and its derivative  
void newton_step(const dvec& a, const dvec& b, dvec& x, dvec& L, dvec& Lp) {
    int n = x.size();
    legendre(a,b,x,L,Lp);
    for(int i=0;i<n;++i) {
        x[i] -= L[i]/Lp[i];
    } 
}

int gauss(dvec& x, dvec& w) {

    int n = x.size();
    double dist = 1;
    double tol = 1e-15;
    int iter = 0;

    // Use Chebyshev-Gauss nodes and initial guess
    initguess(x);
//    chebnodes(x);

    dvec x0(x);
    dvec L(n,0.0);
    dvec Lp(n,0.0);
    dvec a(n,0.0);
    dvec b(n,0.0);

    rec_legendre(a,b);

    // Iteratively correct nodes with Newton's method
    while(dist>tol) {     
        newton_step(a,b,x,L,Lp);
        dist = dist_max(x,x0); 
        ++iter;
        x0 = x;
    } 

    // Compute Weights
    for(int i=0;i<n;++i){
        w[i] = 2.0/((1-x[i]*x[i])*(Lp[i]*Lp[i]));
    }

    return iter;
}
