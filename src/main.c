#include  <stdio.h>
#include  <stdlib.h>
#include  <time.h>
#include  "sparseMatrix.h"

int main(int argc, char *argv[])
{

  //seeding the random number generator
  srand(time(0));

  int choice;
  printf("\n1. Test assignment parts a & b together.");
  printf("\n2. Test assignment part c.");
  printf("\n3. Test assignment part d.\n");
  scanf("%d", &choice);

  switch (choice) {
    case 1: {

      //defining variables for various input parameters
      int m;
      int n;
      double density_nz;
      double density_negative;

      //pointer to COO struct which will be used by the generateCooMat function to fill in the matrix
      struct CooMat *cooMat;

      //pointer to CSR and CSC struct which will be used to store the results of converted matrices from COO
      struct CsrMat *csrMat;
      struct CscMat *cscMat;

      //taking inputs for generating the matrix
      printf("\nEnter the row(m), col(n), density of non-zero elements and negative elements repectively: ");
      scanf("%d %d %lf %lf", &m, &n, &density_nz, &density_negative);

      //generating the matrix with the specified inputs
      printf("\nGenerating random matrix in COO format for the entered values of m: %d, n: %d, density_nz: %lf, density_negative: %lf", m, n, density_nz, density_negative);
      generateCooMat(m,n,density_nz,density_negative,&cooMat);
      printCooMatFull(&cooMat);

      //converting the matrix to CSR and CSC formats
      printf("\nConverting the generated matrix to CSR and CSC formats:");
      convert_coo_csr(&cooMat, &csrMat);
      convert_coo_csc(&cooMat, &cscMat);

      //printing the matrices in the new formats
      printCsrMatFull(&csrMat);
      printCscMatFull(&cscMat);

      //freeing the memory allocated for the different matrices
      free(cooMat);
      free(csrMat);
      free(cscMat);

      break;
    }

    case 2: {

      //pointer to sturctures for storing the read matrices
      struct CsrMat *csrMat;
      struct CscMat *cscMat;

      //pointer to COO struct for storing the result of multiplication
      struct CooMat *cooMat;
      struct CooMat *cooMat2;

      //reading the two matrices
      printf("\nFirst enter the matrix A (CSR) below:\n");
      readCsrMat(&csrMat);
      printf("\nNow enter the matrix B (CSC) below:\n");
      readCscMat(&cscMat);

      //multiplying the two matrices
      printf("\nMultiplying the input matrix A & B and printing the resultant matrix AB in COO format\n");
      multiply_csr_csc(&csrMat, &cscMat, &cooMat);
      printCooMatFull(&cooMat);

      //multiplying the two matrices in reverse orientation
      printf("\nMultiplying the input matrix B & A and printing the resultant matrix BA in COO format\n");
      multiply_csc_csr(&cscMat, &csrMat, &cooMat2);
      printCooMatFull(&cooMat2);

      //freeing the memory allocated for the different matrices
      free(cooMat);
      free(cooMat2);
      free(csrMat);
      free(cscMat);

      break;
    }
    case 3: {

      //variables to store input parameters
      int r1, r2, k1, k2;

      //pointer to sturctures for storing the matrices
      struct CscMat *cscMat;
      struct CscMat *resCscMat;

      //reading the CSC matrix
      printf("\nFirst enter the matrix in CSC format below:\n");
      readCscMat(&cscMat);

      //reading row transformation parameters
      printf("\nEnter the inputs (all are zero-indexed and integers only) for row transformation r1, k1, r2 and k2 respectively: ");
      scanf("%d %d %d %d", &r1, &k1, &r2, &k2);

      //performing the row rowTransformation
      rowTransformationCscMat(r1, k1, r2, k2, &cscMat, &resCscMat);
      printCscMatFull(&resCscMat);

      //freeing the memory for matrices
      free(cscMat);
      free(resCscMat);

      break;
    }
  }

  return 0;
}
