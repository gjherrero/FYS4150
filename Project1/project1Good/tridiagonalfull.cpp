#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
// use namespace for output and input
using namespace std;

// object for output files
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Begin main program
int main(){
    ofstream ofile;
    int exponent;
    string filename;
    int maxExponent;
    cin>>maxExponent;
    exponent = maxExponent;
    filename = to_string(maxExponent);
    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i-
      fileout.append(argument);
      fileout.append("full.txt");
      double h = 1.0/(n+1);
      double hh = h*h;
      // Set up arrays for the full case (question b)
      double *a = new double [n+1];double *b = new double [n+1]; double *c = new double [n+1];double *d = new double [n+1]; double *solution = new double [n+1];
      double *x = new double[n+1];
      // Quick setup of updated diagonal elements and value of
      b[0] = 0; b[1] = b[n] = 2; solution[0] = 0.0; // d es la diagonal y solo 1 elemento
      a[0] = a[1] = 0; a[n] = -1;
      c[0] = c[n] = 0; c[1] = -1;
      for (int i = 2; i < n; i++){
          a[i] = -1;
          b[i] = 2;
          c[i] = -1;
        }
      x[0]=d[0]=0;
      for (int i = 1; i <= (n); i++){
        x[i]= i*h;
        d[i] = hh*f(i*h);
      }

      // Start timing
      clock_t start, finish;
      start = clock();
      double btemp = b[1];
      double *ctemp = new double[n+1];
      solution[1] = d[1]/btemp;
      // Forward substitution
      for (int i = 2; i <= n; i++){
          ctemp[i] = c[i-1]/btemp;
          btemp = b[i]-a[i]*ctemp[i];
          solution[i] = (d[i] - a[i]*solution[i-1])/btemp;
      }
      // Backward substitution
      for (int i = n-1; i >= 1; i--){
          solution[i] -= ctemp[i+1]*solution[i+1];
      }
      //Stop clock and read
      finish = clock();
      double timeused = (double) (finish - start)/((double) CLOCKS_PER_SEC )*1000;
      cout << setiosflags(ios::showpoint | ios::uppercase);
      cout << setprecision(10) << setw(20) << "Time used  for  computation=" << timeused<<"ms"<<endl;

      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i <= n;i++) {
         double xval = x[i];
         double RelativeError = fabs((exact(xval)-solution[i])/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << solution[i];
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] a; delete [] b; delete [] c; delete [] d; delete [] solution; delete [] x; delete [] ctemp;
    }
    return 0;
}


