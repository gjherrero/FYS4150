////  Solves linear equations for simple tridiagonal matrix using LU decomposition
////  This is armadillo version that calls the function solve.
//#include <iostream>
//#include <fstream>
//#include <iomanip>
//#include <cmath>
//#include <string>
//// use namespace for output and input
//using namespace std;
//#define   ZERO       1.0E-15

///* function declarations */
//double ** AllocateMatrix(int, int);
//void DeallocateMatrix(double **, int, int);
//void MatrixInverse(double **, int);
//void WriteMatrix(double **, int);
//void MatrixMultiplication(double **, double **, int);
//void LUDecomposition(double **, int, int *);
//void LUBackwardSubstitution(double **, int, int *, double *);

//// Begin main function, reads from terminal mode the dimension
//// Functions used
//double f(double x){return 100.0*exp(-10.0*x);
//}
//double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

//// Begin main program
//int main(int argc, char *argv[]){
//    ofstream ofile;
//    int exponent;
//    string filename;
//    int maxExponent;
//    cin>>maxExponent;
//    exponent = maxExponent;
//    filename = to_string(maxExponent);

//    // Loop over powers of 10
//    for (int i = 1; i <= exponent; i++){
//      int  n = (int) pow(10.0,i);
//      // Declare new file name
//      string fileout = filename;
//      // Convert the power 10^i to a string
//      string argument = to_string(i);
//      // Final filename as filename-i-
//      fileout.append(argument);
//      fileout.append("LU.txt");
//      double h = 1.0/(n);
//      double hh = h*h;
//      n = n-1;  //  shift so that only points between endpoints are studied

//      double **A = AllocateMatrix(n, n);
//      double *b = new double [n+1]; double *x = new double[n+1];
//      A[0][0] = 2.0;  A[0][1] = -1;
//      x[0] = h; b[0] =  hh*f(x[0]);
//      x[n-1] = x[0]+(n-1)*h; b[n-1] = hh*f(x[n-1]);

//      for (int i = 1; i < n-1; i++){
//        x[i] = x[i-1]+h;
//        b[i] = hh*f(x[i]);
//        A[i][i-1]  = -1.0;
//        A[i][i]    = 2.0;
//        A[i][i+1]  = -1.0;
//      }
//      A[n-1][n-1] = 2.0; A[n-2][n-1] = -1.0; A[n-1][n-2] = -1.0;
//  // solve Ax = b
//      // Start timing
//      clock_t start, finish;
//      start = clock();
//      MatrixInverse(A,n);
//      double *solution = new double[n+1];
//      for(int i=0;i < n;i++){
//           double sum = 0.0;
//           for (int k = 0; k < n; k++) sum += A[i][k]*b[k];
//           solution[i] = sum;
//      }
//      //Stop clock and read
//      finish = clock();
//      double timeused = (double) (finish - start)/((double) CLOCKS_PER_SEC )*1000;
//      cout << setiosflags(ios::showpoint | ios::uppercase);
//      cout << setprecision(10) << setw(20) << "Time used  for  computation=" << timeused<<"ms"<<endl;

//      ofile.open(fileout);
//      ofile << setiosflags(ios::showpoint | ios::uppercase);
//      for (int i = 0; i < n ;i++) {
//        double xval = x[i];
//        double RelativeError = fabs((exact(xval)-solution[i])/exact(xval));
//        ofile << setw(15) << setprecision(8) << xval;
//        ofile << setw(15) << setprecision(8) << solution[i];
//        ofile << setw(15) << setprecision(8) << exact(xval);
//        ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
//      }
//      ofile.close();
//      delete [] b; delete [] solution; delete [] x; delete [] A;
//    }
//    return 0;
//}


///*
//   The function MatrixInverse() performs a matrix inversion
//   of a square matrix a[][] of  dimension n.
//*/

//void MatrixInverse(double **A, int n)
//{
//  // allocate space in memory
//  int *indx;
//  double *column;
//  indx = new int[n];
//  column  = new double[n];
//  double **Y    = AllocateMatrix(n,n);
//  // Perform the LU decomposition
//  LUDecomposition(A, n, indx);   // LU decompose  a[][]
//  //cout << "LU decomposed matrix  A:" << endl;
//  //WriteMatrix(A,n);
//  // find inverse of a[][] by columns
//  for(int j = 0; j < n; j++) {
//    // initialize right-side of linear equations
//    for(int i = 0; i < n; i++) column[i] = 0.0;
//    column[j] = 1.0;
//    LUBackwardSubstitution(A, n, indx, column);
//    // save result in y[][]
//    for(int i = 0; i < n; i++) Y[i][j] = column[i];
//  }
//  // return the inverse matrix in A[][]
//  for(int i = 0; i < n; i++) {
//    for(int j = 0; j < n; j++) A[i][j] = Y[i][j];
//  }
//  DeallocateMatrix(Y, n, n);     // release local memory
//  delete [] column;
//  delete []indx;

//}  // End: function MatrixInverse()


//// Allocate memory for a matrix and initialize the elements to zero

//double ** AllocateMatrix(int m, int n){
//  double ** Matrix;
//  Matrix = new double*[m];
//  for(int i=0;i<m;i++){
//    Matrix[i] = new double[n];
//    for(int j=0;j<m;j++)
//      Matrix[i][j] = 0.0;
//  }
//  return Matrix;
//}

//// Free memory

//void DeallocateMatrix(double ** Matrix, int m, int n){
//  for(int i=0;i<m;i++)
//    delete[] Matrix[i];
//  delete[] Matrix;
//}

//// Write out a given matrix
//void WriteMatrix(double ** Matrix, int n){
//  for(int i=0;i < n;i++){
//    cout << endl;
//     for (int j=0 ; j < n;j++){
//        printf("  A[%2d][%2d] = %12.4E",i, j, Matrix[i][j]);
//     }
//  }
//    cout << endl;
//}

//// Straightforward matrix-matrix multiplication

//void MatrixMultiplication(double ** a, double **b, int n){
//  double **c = AllocateMatrix(n, n);
//  for(int i=0;i < n;i++){
//     for (int j=0 ; j < n;j++){
//       double sum = 0.0;
//       for (int k = 0; k < n; k++) sum += a[i][k]*b[k][j];
//       c[i][j] = sum;
//     }
//  }
//  WriteMatrix(c,n);
//}

///*
//    The function
//    void LUDecomposition(double **a, int n, int *indx)
//    takes as input a two-dimensional matrix a[][] of dimension n and
//    replaces it by the LU decomposition of a rowwise permutation of
//    itself. The results is stored in a[][]
//    The vector
//    indx[] records the row permutation effected by the partial pivoting.
//*/

//void LUDecomposition(double **a, int n, int *indx)
//{
//   int      i, imax, j, k;
//   double   big, dum, sum, temp, *vv;

//  vv = new double [n];
//   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
//      big = ZERO;
//      for(j = 0; j < n; j++) {
//         if((temp = fabs(a[i][j])) > big) big = temp;
//      }
//      if(big == ZERO) {
//         printf("\n\nSingular matrix in routine ludcmp()\n");
//         exit(1);
//      }
//      vv[i] = 1.0/big;
//   }

//   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
//      for(i = 0; i< j; i++) {   // not i = j
//         sum = a[i][j];
//     for(k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
//     a[i][j] = sum;
//      }
//      big = ZERO;   // initialization for search for largest pivot element
//      for(i = j; i< n; i++) {
//         sum = a[i][j];
//     for(k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
//     a[i][j] = sum;
//     if((dum = vv[i]*fabs(sum)) >= big) {
//        big = dum;
//        imax = i;
//     }
//      } // end i-loop
//      if(j != imax) {    // do we need to interchange rows ?
//         for(k = 0;k< n; k++) {       // yes
//        dum        = a[imax][k];
//        a[imax][k] = a[j][k];
//        a[j][k]    = dum;
//     }
//     vv[imax] = vv[j];         // also interchange scaling factor
//      }
//      indx[j] = imax;
//      if(fabs(a[j][j]) < ZERO)  a[j][j] = ZERO;
//      if(j < (n - 1)) {                   // divide by pivot element
//         dum = 1.0/a[j][j];
//     for(i=j+1;i < n; i++) a[i][j] *= dum;
//      }
//   } // end j-loop over columns
//   delete [] vv;   // release local memory
//}


///*
//     The function
//       void LUBackwardSubstitution(double **a, int n, int *indx, double *b)
//     solves the set of linear equations A X = B of dimension n.
//     a[][] is input, not as the matrix A[][] but rather as
//     its LU decomposition, indx[] is input as the permutation vector returned by
//     ludcmp(). b[] is input as the right-hand side vector B,
//     The solution X is returned in B. The input data a[][],
//     n and indx[] are not modified. This routine takes into
//     account the possibility that b[] will begin with many
//     zero elements, so it is efficient for use in matrix
//     inversion.
//*/

//void LUBackwardSubstitution(double **a, int n, int *indx, double *b)
//{
//   int        i, ii = -1, ip, j;
//   double     sum;

//   for(i = 0; i< n; i++) {
//      ip    = indx[i];
//      sum   = b[ip];
//      b[ip] = b[i];
//      if(ii > -1)   for(j = ii; j < i; j++) sum -= a[i][j] * b[j];
//      else if(sum) ii = i;
//      b[i] = sum;
//   }
//   for(i = n - 1; i >= 0; i--) {
//      sum = b[i];
//      for(j = i+1; j < n; j++) sum -= a[i][j] * b[j];
//      b[i] = sum/a[i][i];
//   }
//}








