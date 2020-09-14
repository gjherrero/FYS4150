//#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <cmath>
//#include <string>
//// use namespace for output and input
//using namespace std;

//// object for output files
//// Functions used
//inline double f(double x){return 100.0*exp(-10.0*x);
//}
//inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

//// Begin main program
//int main(){
//    ofstream ofile;
//    int exponent;
//    string filename;
//    int maxExponent;
//    cin>>maxExponent;
//    exponent = maxExponent;
//    filename = to_string(maxExponent);
//    // We read also the basic name for the output file and the highest power of 10^n we want
//    // Loop over powers of 10
//    for (int i = 1; i <= exponent; i++){
//      int  n = (int) pow(10.0,i);
//      // Declare new file name
//      string fileout = filename;
//      // Convert the power 10^i to a string
//      string argument = to_string(i);
//      // Final filename as filename-i-
//      fileout.append(argument);
//      fileout.append("simp.txt");
//      double h = 1.0/(n+1);
//      double hh = h*h;
//      // Set up arrays for the simple case
//      double *d = new double [n+2]; double *b = new double [n+2]; double *solution = new double [n+2];
//      double *x = new double[n+2];
//      // Quick setup of updated diagonal elements and value of
//      d[0] = d[n+1] = 2; solution[0] = solution[n+1] = 0.0;
//      for (int i = 1; i <= n; i++) d[i] = (i+1.0)/( (double) i);
//      for (int i = 0; i <= (n+1); i++){
//        x[i]= i*h;
//        b[i] = hh*f(i*h);
//      }
//      // Start timing
//      clock_t start, finish;
//      start = clock();
//      // Forward substitution
//      for (int i = 2; i <= n; i++) b[i] = b[i] + b[i-1]/d[i-1];
//      // Backward substitution
//      solution[n] = b[n]/d[n];
//      for (int i = n-1; i > 0; i--) solution[i] = (b[i]+solution[i+1])/d[i];
//      //Stop clock and read
//      finish = clock();
//      double timeused = (double) (finish - start)/((double) CLOCKS_PER_SEC );
//      cout << setiosflags(ios::showpoint | ios::uppercase);
//      cout << setprecision(10) << setw(20) << "Time used  for  computation=" << timeused;

//              ofile.open(fileout);
//      ofile << setiosflags(ios::showpoint | ios::uppercase);
//      //      ofile << "       x:             approx:          exact:       relative error" << endl;
//      for (int i = 1; i <= n;i++) {
//    double xval = x[i];
//     double RelativeError = fabs((exact(xval)-solution[i])/exact(xval));
//         ofile << setw(15) << setprecision(8) << xval;
//         ofile << setw(15) << setprecision(8) << solution[i];
//         ofile << setw(15) << setprecision(8) << exact(xval);
//         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
//      }
//      ofile.close();
//      delete [] x; delete [] d; delete [] b; delete [] solution;
//    }
//    return 0;
//}


