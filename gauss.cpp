#include<iostream>
#include"gauss.hpp"




int main(int argc, char *argv[]) {

    // Number of quadrature points 
    int n = 3;
    int k = 1;
 
    if (argc > 1) {
        n = atoi( argv[1] );
    }
    if (argc > 2) {
        k = atoi( argv[2] );
    }


    // Quadrature points (allocate vector)
    dvec x(n,0.0);

    // Quadrature weights (allocate vector)
    dvec w(n,0.0);

    gauss(x,w); 

    double sum = 0;
    for(int i=0;i<n;++i) {
        sum += w[i]*std::cos(k*pi*x[i]/2.0);
    } 

    double exact = 4.0*std::sin(k*pi/2.0)/(pi*k); 

    std::cout << "Computed integral of cos(" << k << "*pi*x/2) on [-1,1]" << std::endl; 
    std::cout << "Q[f] = " << sum << std::endl;
    std::cout << "Exact value = 4*sin(k*pi/2)/(k*pi)" << std::endl;
    std::cout << "Error = " << std::abs(exact-sum) << std::endl;
  
    return 0;
}
