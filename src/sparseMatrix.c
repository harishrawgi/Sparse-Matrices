#include  <stdio.h>
#include  <stdlib.h>
#include  <time.h>
#include  "sparseMatrix.h"
#include  "utils.h"

//Functions to print the matrices of different formats in full format
void printCooMatFull(struct CooMat ** cooMat){

  int nz = (*cooMat)->nz;
  int m = (*cooMat)->n_row;
  int n = (*cooMat)->n_col;


  if(m>20 || n>20){
    printf("\nMatrix too large to print.\n");
    return;
  }

  int fullMatrix[m][n];

  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      fullMatrix[i][j]=0;
    }
  }

  for(int i=0; i<nz; i++){
    fullMatrix[(*cooMat)->row_ind[i]][(*cooMat)->col_ind[i]] = (*cooMat)->val[i];
  }

  printf("\nPrinting COO matrix in full.\n");
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      printf("%10d\t", fullMatrix[i][j]);
    }
    printf("\n");
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
}

void printCsrMatFull(struct CsrMat ** csrMat){

  int nz = (*csrMat)->nz;
  int m = (*csrMat)->n_row;

  int max=0;
  for(int i=0; i<nz; i++){
    if(((*csrMat)->col_ind[i]) > max){
      max = (*csrMat)->col_ind[i];
    }
  }
  int n = max+1;

  if(m>20 || n>20){
    printf("\nMatrix too large to print.\n");
    return;
  }

  int fullMatrix[m][n];

  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      fullMatrix[i][j]=0;
    }
  }

  int num_row[m];
  for(int i=1; i<m+1; i++){
    num_row[i-1] = (*csrMat)->row_ptr[i] - (*csrMat)->row_ptr[i-1];
  }

  int k = 0;
  int left = num_row[k];

  for(int i=0; i<nz; i++){
    while(!left){
      k++;
      left = num_row[k];
    }
    fullMatrix[k][(*csrMat)->col_ind[i]] = (*csrMat)->val[i];
    left--;
  }

  printf("\nPrinting CSR matrix in full.\n");
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      printf("%10d\t", fullMatrix[i][j]);
    }
    printf("\n");
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
}

void printCscMatFull(struct CscMat ** cscMat){

  int nz = (*cscMat)->nz;
  int n = (*cscMat)->n_col;

  int max=0;
  for(int i=0; i<nz; i++){
    if(((*cscMat)->row_ind[i]) > max){
      max = (*cscMat)->row_ind[i];
    }
  }
  int m = max+1;

  if(m>20 || n>20){
    printf("\nMatrix too large to print.\n");
    return;
  }

  int fullMatrix[m][n];

  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      fullMatrix[i][j]=0;
    }
  }

  int num_col[n];
  for(int i=1; i<n+1; i++){
    num_col[i-1] = (*cscMat)->col_ptr[i] - (*cscMat)->col_ptr[i-1];
  }

  int k = 0;
  int left = num_col[k];

  for(int i=0; i<nz; i++){
    while(!left){
      k++;
      left = num_col[k];
    }
    fullMatrix[(*cscMat)->row_ind[i]][k] = (*cscMat)->val[i];
    left--;
  }

  printf("\nPrinting CSC matrix in full.\n");
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  for(int i=0; i<m; i++){
    for(int j=0; j<n; j++){
      printf("%10d\t", fullMatrix[i][j]);
    }
    printf("\n");
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
}

//functions to print the matrices in specified sparse formats
void printCooMat(struct CooMat ** cooMat){

  printf("\nPrinting matrix in COO Format\n\n");
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Val:\t");
  for(int i=0; i<((*cooMat)->nz); i++){
    printf("%10d\t", ((*cooMat)->val)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Row:\t");
  for(int i=0; i<((*cooMat)->nz); i++){
    printf("%10d\t", ((*cooMat)->row_ind)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Col:\t");
  for(int i=0; i<((*cooMat)->nz); i++){
    printf("%10d\t", ((*cooMat)->col_ind)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
}

void printCsrMat(struct CsrMat ** csrMat){
  printf("\nPrinting matrix in CSR Format...\n\n");
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Val:\t");
  for(int i=0; i<((*csrMat)->nz); i++){
    printf("%10d\t", ((*csrMat)->val)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Row:\t");
  for(int i=0; i<((*csrMat)->n_row)+1; i++){
    printf("%10d\t", ((*csrMat)->row_ptr)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Col:\t");
  for(int i=0; i<((*csrMat)->nz); i++){
    printf("%10d\t", ((*csrMat)->col_ind)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
}

void printCscMat(struct CscMat ** cscMat){
  printf("\nPrinting matrix in CSC Format\n\n");
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Val:\t");
  for(int i=0; i<((*cscMat)->nz); i++){
    printf("%10d\t", ((*cscMat)->val)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Row:\t");
  for(int i=0; i<((*cscMat)->nz); i++){
    printf("%10d\t", ((*cscMat)->row_ind)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
  printf("Col:\t");
  for(int i=0; i<((*cscMat)->n_col)+1; i++){
    printf("%10d\t", ((*cscMat)->col_ptr)[i]);
  }
  printf("\n-------------------------------------------------------------------------------------------\n\n");
}


//Functions to read in the matrices in various formats
void readCooMat(struct CooMat ** cooMat){

  int n_row, n_col, nz;

  printf("\nEnter the number of rows in the COO matrix: ");
  scanf("%d", &n_row);
  if(n_row > MAX_ROW_SIZE || n_row < 1){
    printf("\nInvalid Input");
    exit(0);
  }
  printf("\nEnter the number of cols in the COO matrix: ");
  scanf("%d", &n_col);
  if(n_col > MAX_COL_SIZE || n_col < 1){
    printf("\nInvalid Input");
    exit(0);
  }
  printf("\nEnter the number of non-zero entries in the COO matrix: ");
  scanf("%d", &nz);
  if(nz > MAX_NZ || nz < 0){
    printf("\nInvalid Input");
    exit(0);
  }

  //allocating memory for the COO structure; will be free back in main()
  *cooMat = (struct CooMat *) malloc(sizeof(struct CooMat));
  if(NULL == cooMat){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }

  (*cooMat)->val = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->val){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->row_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->row_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->col_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->col_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->nz = nz;
  (*cooMat)->n_row = n_row;
  (*cooMat)->n_col = n_col;

  printf("\nEnter the matrix in COO format:\n");
  int x, read;
  for(int i=0; i<nz; i++){
    read = scanf("%d", &x);
    if(x==0 || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*cooMat)->val[i] = x;
  }
  for(int i=0; i<nz; i++){
    read = scanf("%d", &x);
    if(x<0 || x>n_row-1 || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*cooMat)->row_ind[i] = x;
  }
  for(int i=0; i<nz; i++){
    read = scanf("%d", &x);
    if(x<0 || x>n_col-1 || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*cooMat)->col_ind[i] = x;
  }
}

void readCsrMat(struct CsrMat ** csrMat){

  int n_row, nz;

  printf("\nEnter the number of rows in the CSR matrix: ");
  scanf("%d", &n_row);
  if(n_row > MAX_ROW_SIZE || n_row < 1){
    printf("\nInvalid Input");
    exit(0);
  }

  printf("\nEnter the number of non-zero entries in the CSR matrix: ");
  scanf("%d", &nz);
  if(nz > MAX_NZ || nz < 0){
    printf("\nInvalid Input");
    exit(0);
  }

  //allocating memory for the CSR matrix; will be free back in main()
  *csrMat = (struct CsrMat *) malloc(sizeof(struct CsrMat));
  if(NULL == csrMat){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }

  (*csrMat)->val = (int *) malloc(nz*sizeof(int));
  if(NULL == (*csrMat)->val){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }
  (*csrMat)->row_ptr = (int *) malloc((n_row+1)*sizeof(int));
  if(NULL == (*csrMat)->row_ptr){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }
  (*csrMat)->col_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*csrMat)->col_ind){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }
  (*csrMat)->nz = nz;
  (*csrMat)->n_row = n_row;

  printf("\nEnter the matrix in CSR format:\n");
  int x, read;
  for(int i=0; i<nz; i++){
    read = scanf("%d", &x);
    if(x==0 || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*csrMat)->val[i] = x;
  }
  for(int i=0; i<(n_row+1); i++){
    read = scanf("%d", &x);
    if(x<0 || x>nz || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*csrMat)->row_ptr[i] = x;
  }
  for(int i=0; i<nz; i++){
    read = scanf("%d", &x);
    if(x<0 || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*csrMat)->col_ind[i] = x;
  }
}

void readCscMat(struct CscMat ** cscMat){

  int  n_col, nz;

  printf("\nEnter the number of cols in the Csc matrix: ");
  scanf("%d", &n_col);
  if(n_col > MAX_COL_SIZE || n_col < 1){
    printf("\nInvalid Input");
    exit(0);
  }
  printf("\nEnter the number of non-zero entries in the Csc matrix: ");
  scanf("%d", &nz);
  if(nz > MAX_NZ || nz < 0){
    printf("\nInvalid Input");
    exit(0);
  }

  //allocating memory for the CSC matrix; will be free back in main()
  *cscMat = (struct CscMat *) malloc(sizeof(struct CscMat));
  if(NULL == cscMat){
    printf("\nUnable to allocate memory for CSC struct");
    exit(0);
  }

  (*cscMat)->val = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cscMat)->val){
    printf("\nUnable to allocate memory for CSC struct");
    exit(0);
  }
  (*cscMat)->row_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cscMat)->row_ind){
    printf("\nUnable to allocate memory for CSC struct");
    exit(0);
  }
  (*cscMat)->col_ptr = (int *) malloc((n_col+1)*sizeof(int));
  if(NULL == (*cscMat)->col_ptr){
    printf("\nUnable to allocate memory for CSC struct");
    exit(0);
  }
  (*cscMat)->nz = nz;
  (*cscMat)->n_col = n_col;

  printf("\nEnter the matrix in CSC format:\n");
  int x, read;
  for(int i=0; i<nz; i++){
    read = scanf("%d", &x);
    if(x==0 || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*cscMat)->val[i] = x;
  }
  for(int i=0; i<nz; i++){
    read = scanf("%d", &x);
    if(x<0 || read!=1){
      printf("\nInvalid input\n");
      exit(0);
    }
    (*cscMat)->row_ind[i] = x;
  }
  for(int i=0; i<(n_col+1); i++){
    read = scanf("%d", &x);
    if(x<0 || x>nz || read!=1){
      printf("\nInvalid input.\n");
      exit(0);
    }
    (*cscMat)->col_ptr[i] = x;
  }
}


/*
  Function: Generates a random matrix in COO format based on the input parameters.
            The generated matrix is stored in the last input parameter.
  Inputs:
    int m:  Number of rows of the matrix to be generated
    int n:  Number of cols of the matrix to be generated
    double density_nz:  Density of the non-zero elements wrt total number of elements
    double density_negative:  Density of the negative elements wrt total number of elements

  Outputs:
    struct CooMat ** cooMat:  Address of the pointer to the object of COO Matrix structure
                              The generated matrix will be stored using this.
*/
void generateCooMat(int m, int n, double density_nz, double density_negative, struct CooMat ** cooMat){


    //checking if all inputs are valid and within limits
    if(m > MAX_ROW_SIZE || n > MAX_COL_SIZE || density_nz > 1 || density_negative > 1){
        printf("\nInvalid input(s)");
        exit(0);
    }
    if(m < 1 || n < 1 || density_nz <= 0 || density_negative < 0){
        printf("\nInvalid input(s)");
        exit(0);
    }


    //calculating the number of non-zero and negative elements
    int nz = m*n*density_nz;
    int neg = m*n*density_negative;
    int pos = nz-neg;

    //allocating memory for the COO structure; will be free back in main()
    *cooMat = (struct CooMat *) malloc(sizeof(struct CooMat));
    if(NULL == cooMat){
      printf("\nUnable to allocate memory for COO struct");
      exit(0);
    }

    (*cooMat)->val = (int *) malloc(nz*sizeof(int));
    if(NULL == (*cooMat)->val){
      printf("\nUnable to allocate memory for COO struct");
      exit(0);
    }
    (*cooMat)->row_ind = (int *) malloc(nz*sizeof(int));
    if(NULL == (*cooMat)->row_ind){
      printf("\nUnable to allocate memory for COO struct");
      exit(0);
    }
    (*cooMat)->col_ind = (int *) malloc(nz*sizeof(int));
    if(NULL == (*cooMat)->col_ind){
      printf("\nUnable to allocate memory for COO struct");
      exit(0);
    }

    (*cooMat)->nz = nz;
    (*cooMat)->n_row = m;
    (*cooMat)->n_col = n;

    //generating +ve non-zero numbers
    for(int i=0; i<pos;i++){
      ((*cooMat)->val)[i] = 1 + rand() % RAND_MAX;
    }

    //generating remaining -ve non-zero numbers
    for(int i=pos; i<nz; i++){
      ((*cooMat)->val)[i] = (-1)*(1 + rand() % RAND_MAX);
    }

    //shuffling the val array to make sure the positions of +ve and -ve nos. are random
    shuffle(((*cooMat)->val), nz);

    //now we need to distribute these random numbers randomly between the row and col indices
    int upper_limit = m*n-1;
    int lower_limit = 0;
    int index_array[nz];

    //initializing the array to zero
    for(int i=0; i<nz; i++)
      index_array[i] = 0;

    //generating a random 1D index for the matrix which we will later convert to 2D index
    for(int i=0; i<nz; i++){
      int index = (rand() % (upper_limit - lower_limit + 1)) + lower_limit;

      //flag to check whether this random number has been drawn before
      int flag=0;
      //iterate over all previous 1D indices to make sure they are not repeated
      for(int k=0; k<nz; k++){
          if(index_array[k] == index){
            flag = 1;
            break;
          }
      }
      if(flag){
        //repeated random number, need to redraw another one
        i--;
      }
      else{
        //assigning the random number as one of the 1D index
        index_array[i] = index;
      }

    }

    //sorting the 1D array which will result in final sorting by rows and then cols
    bubbleSort(index_array, nz);

    //converting the 1D index to 2D index and thus generating the row_ind and col_ind arrays
    for(int i=0; i<nz; i++){
      int index = index_array[i];
      ((*cooMat)->row_ind)[i] = index/n;
      ((*cooMat)->col_ind)[i] = index%n;
    }

}


//functions to convert a COO matrix to col major and row major format
void convertCooToColumnMajor(struct CooMat **cooMat, int n){
  for (int i = 0; i < n-1; i++)
    for (int j = 0; j < n-i-1; j++)
        if ( (*cooMat)->col_ind[j] >  (*cooMat)->col_ind[j+1]){

          int t = (*cooMat)->col_ind[j];
          (*cooMat)->col_ind[j] = (*cooMat)->col_ind[j+1];
          (*cooMat)->col_ind[j+1] = t;

          t = (*cooMat)->val[j];
          (*cooMat)->val[j] = (*cooMat)->val[j+1];
          (*cooMat)->val[j+1] = t;

          t = (*cooMat)->row_ind[j];
          (*cooMat)->row_ind[j] = (*cooMat)->row_ind[j+1];
          (*cooMat)->row_ind[j+1] = t;
        }
}
void convertCooToRowMajor(struct CooMat **cooMat, int n){
  for (int i = 0; i < n-1; i++)
    for (int j = 0; j < n-i-1; j++)
        if ( (*cooMat)->row_ind[j] >  (*cooMat)->row_ind[j+1]){

          int t = (*cooMat)->col_ind[j];
          (*cooMat)->col_ind[j] = (*cooMat)->col_ind[j+1];
          (*cooMat)->col_ind[j+1] = t;

          t = (*cooMat)->val[j];
          (*cooMat)->val[j] = (*cooMat)->val[j+1];
          (*cooMat)->val[j+1] = t;

          t = (*cooMat)->row_ind[j];
          (*cooMat)->row_ind[j] = (*cooMat)->row_ind[j+1];
          (*cooMat)->row_ind[j+1] = t;
        }
}

//function to convert a matrix from COO format to CSR format
void convert_coo_csr(struct CooMat **cooMat, struct CsrMat **csrMat){


    int nz = (*cooMat)->nz;
    int m = (*cooMat)->n_row;

    //allocating memory for the CSR matrix; will be free back in main()
    *csrMat = (struct CsrMat *) malloc(sizeof(struct CsrMat));
    if(NULL == csrMat){
      printf("\nUnable to allocate memory for CSR struct");
      exit(0);
    }

    (*csrMat)->nz = (*cooMat)->nz;
    (*csrMat)->n_row = (*cooMat)->n_row;
    //(*csrMat)->n_col = (*cooMat)->n_col;



    //sharing the same memory for val and col_ind
    (*csrMat)->val = (int *) malloc(nz*sizeof(int));
    if(NULL == (*csrMat)->val){
      printf("\nUnable to allocate memory for CSR struct");
      exit(0);
    }
    (*csrMat)->val = (*cooMat)->val;
    (*csrMat)->col_ind = (int *) malloc(nz*sizeof(int));
    if(NULL == (*csrMat)->col_ind){
      printf("\nUnable to allocate memory for CSR struct");
      exit(0);
    }
    (*csrMat)->col_ind = (*cooMat)->col_ind;

    //size of row_ptr is fixed at number of rows + 1
    (*csrMat)->row_ptr = (int *) malloc((m+1)*sizeof(int));
    if(NULL == (*csrMat)->row_ptr){
      printf("\nUnable to allocate memory for CSR struct");
      exit(0);
    }

    //count_ar is used to count the number of non-zero entries in each row
    int count_ar[m];

    //initializing count_ar to zero
    for(int i=0; i<m; i++){
        count_ar[i] = 0;
    }

    //counting the number of non-zero entries in each row
    for(int i=0; i<nz; i++){
        count_ar[(*cooMat)->row_ind[i]]++;
    }

    //generating the row_ptr with the help of count array
    (*csrMat)->row_ptr[0] = 0;
    for(int i=1; i<m+1; i++){
        (*csrMat)->row_ptr[i] = (*csrMat)->row_ptr[i-1] + count_ar[i-1];
    }

}

//function to convert a matrix from COO format to CSC format
void convert_coo_csc(struct CooMat **cooMat, struct CscMat **cscMat){

    int nz = (*cooMat)->nz;
    int n = (*cooMat)->n_col;

    //converting COO matrix from row major to column major
    convertCooToColumnMajor(cooMat, nz);

    //allocating memory for the CSC matrix; will be free back in main()
    *cscMat = (struct CscMat *) malloc(sizeof(struct CscMat));
    if(NULL == cscMat){
      printf("\nUnable to allocate memory for CSC struct");
      exit(0);
    }

    (*cscMat)->nz = (*cooMat)->nz;
    //(*cscMat)->n_row = (*cooMat)->n_row;
    (*cscMat)->n_col = (*cooMat)->n_col;

    //can't share the same memory because CSC must be in Column Major format
    (*cscMat)->val = (int *) malloc(nz*sizeof(int));
    if(NULL == (*cscMat)->val){
      printf("\nUnable to allocate memory for CSC struct");
      exit(0);
    }
    (*cscMat)->row_ind = (int *) malloc(nz*sizeof(int));
    if(NULL == (*cscMat)->row_ind){
      printf("\nUnable to allocate memory for CSC struct");
      exit(0);
    }

    for(int i=0; i<nz; i++){
        (*cscMat)->val[i] = (*cooMat)->val[i];
        (*cscMat)->row_ind[i] = (*cooMat)->row_ind[i];
    }

    //size of col_ptr is fixed at number of cols + 1
    (*cscMat)->col_ptr = (int *) malloc((n+1)*sizeof(int));
    if(NULL == (*cscMat)->col_ptr){
      printf("\nUnable to allocate memory for CSC struct");
      exit(0);
    }

    //count_ar is used to count the number of non-zero entries in each col
    int count_ar[n];

    //initializing count_ar to zero
    for(int i=0; i<n; i++){
        count_ar[i] = 0;
    }

    //counting the number of non-zero entries in each col
    for(int i=0; i<nz; i++){
        count_ar[(*cooMat)->col_ind[i]]++;
    }

    //generating the col_ptr with the help of count array
    (*cscMat)->col_ptr[0] = 0;
    for(int i=1; i<n+1; i++){
        (*cscMat)->col_ptr[i] = (*cscMat)->col_ptr[i-1] + count_ar[i-1];
    }

    //converting COO matrix back to row major
    convertCooToRowMajor(cooMat, nz);


}

//function to convert a matrix from CSC format to COO format
void convert_csc_coo(struct CscMat **cscMat, struct CooMat **cooMat){

  int nz = (*cscMat)->nz;
  int n = (*cscMat)->n_col;

  //allocating memory for the COO structure; will be free back in main()
  *cooMat = (struct CooMat *) malloc(sizeof(struct CooMat));
  if(NULL == cooMat){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->val = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->val){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->row_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->row_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->col_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->col_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->nz = nz;
  (*cooMat)->n_col = n;

  int max=0;
  for(int i=0; i<nz; i++){
    if((*cscMat)->row_ind[i] > max){
      max = (*cscMat)->row_ind[i];
    }
  }
  (*cooMat)->n_row = max+1;

  int num_col[n];
  for(int i=1; i<n+1; i++){
    num_col[i-1] = (*cscMat)->col_ptr[i] - (*cscMat)->col_ptr[i-1];
  }

  int k = 0;
  int left = num_col[k];

  for(int i=0; i<nz; i++){
    while(!left){
      k++;
      left = num_col[k];
    }

    (*cooMat)->val[i] = (*cscMat)->val[i];
    (*cooMat)->row_ind[i] = (*cscMat)->row_ind[i];
    (*cooMat)->col_ind[i] = k;
    left--;
  }

}

//function to convert a matrix from CSR format to COO format
void convert_csr_coo(struct CsrMat **csrMat, struct CooMat **cooMat){

  int nz = (*csrMat)->nz;
  int m = (*csrMat)->n_row;

  //allocating memory for the COO structure; will be free back in main()
  *cooMat = (struct CooMat *) malloc(sizeof(struct CooMat));
  if(NULL == cooMat){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->val = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->val){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->row_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->row_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->col_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->col_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->nz = nz;
  (*cooMat)->n_row = m;

  int max=0;
  for(int i=0; i<nz; i++){
    if((*csrMat)->col_ind[i] > max){
      max = (*csrMat)->col_ind[i];
    }
  }
  (*cooMat)->n_col = max+1;

  int num_row[m];
  for(int i=1; i<m+1; i++){
    num_row[i-1] = (*csrMat)->row_ptr[i] - (*csrMat)->row_ptr[i-1];
  }

  int k = 0;
  int left = num_row[k];

  for(int i=0; i<nz; i++){
    while(!left){
      k++;
      left = num_row[k];
    }

    (*cooMat)->val[i] = (*csrMat)->val[i];
    (*cooMat)->col_ind[i] = (*csrMat)->col_ind[i];

    (*cooMat)->row_ind[i] = k;
    left--;
  }

}

//function to reallocate the size of COO mat to new size
void reallocCooMat(struct CooMat **cooMat, int size){



  int *ptr = realloc((*cooMat)->val, size * sizeof(int));
  if (ptr == NULL) // reallocated pointer ptr1
  {
      printf("\nInsufficient memory; Exiting!!");
      free(ptr);
      exit(0);
  }
  else
  {
      (*cooMat)->val = ptr;           // the reallocation succeeded, we can overwrite our original pointer now
  }

  ptr = realloc((*cooMat)->row_ind, size * sizeof(int));
  if (ptr == NULL) // reallocated pointer ptr1
  {
      printf("\nInsufficient memory; Exiting!!");
      free(ptr);
      exit(0);
  }
  else
  {
      (*cooMat)->row_ind = ptr;           // the reallocation succeeded, we can overwrite our original pointer now
  }

  ptr = realloc((*cooMat)->col_ind, size * sizeof(int));
  if (ptr == NULL) // reallocated pointer ptr1
  {
      printf("\nInsufficient memory; Exiting!!");
      free(ptr);
      exit(0);
  }
  else
  {
      (*cooMat)->col_ind = ptr;           // the reallocation succeeded, we can overwrite our original pointer now
  }

  (*cooMat)->nz = size;

}


/*
  Function: Multiplies two matrices A (in CSR format) and B (in CSC format)
            The resultant matrix is stored in the last input parameter in form of COO.
  Inputs:
    struct CsrMat **csrMat:  Input matrix A in CSR format
    struct CscMat **cscMat:  Input matrix B in CSC format

  Outputs:
    struct CooMat ** cooMat:  Address of the pointer to the object of COO Matrix structure
                              The resultant matrix will be stored using this.
*/
void multiply_csr_csc(struct CsrMat **csrMat, struct CscMat **cscMat, struct CooMat **cooMat){

  int m1 = (*csrMat)->n_row;
  int n2 = (*cscMat)->n_col;

  //initializing the no. of non-zero entries; will be adjusted later in the function
  int nz = ((*csrMat)->nz) + ((*cscMat)->nz);

  //the no. of rows of B (CSC) or the no. of columns of A(CSR) is unknown.
  //We can take it to be the maximum index of col in A (CSR)
  int m2 = 0;
  for(int i=0; i<((*csrMat)->nz); i++){
    if(((*csrMat)->col_ind[i]) > m2)
      m2 = ((*csrMat)->col_ind[i]);
  }
  m2++;

  //allocating memory for the COO structure; will be free back in main()
  *cooMat = (struct CooMat *) malloc(sizeof(struct CooMat));
  if(NULL == cooMat){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->val = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->val){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->row_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->row_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->col_ind = (int *) malloc(nz*sizeof(int));
  if(NULL == (*cooMat)->col_ind){
    printf("\nUnable to allocate memory for COO struct");
    exit(0);
  }
  (*cooMat)->nz = nz;
  (*cooMat)->n_row = m1;
  (*cooMat)->n_col = n2;

  //allocating memory for storing some intermediate values; will be made free at the end
  int *x = (int *)malloc(sizeof(int));
  int *y = (int *)malloc(sizeof(int));

  int coo_ind=0;

  //iterating each column of the CSC matrix
  for(int current_col=0; current_col<n2;current_col++){

    //analysing where the currect col lies in the CSC matrix format
    int cur_col_ptr = (*cscMat)->col_ptr[current_col];
    int size = (*cscMat)->col_ptr[current_col+1] - (*cscMat)->col_ptr[current_col];
    if(size != 0){ //this means there is atleast one non-zero entry in current col

      //cvec will hold the current col
      int *cvec = (int *)malloc(m2*sizeof(int));
      for(int z=0; z<m2;z++){
        cvec[z] = 0;
      }

      //filling the entries of the current col in cvec
      for(int k=cur_col_ptr; k<cur_col_ptr+size; k++){
        cvec[(*cscMat)->row_ind[k]] = (*cscMat)->val[k];
      }

      //now after we have the currect col in cvec, we will iterate over all rows of CSR and perform multiplication
      for(int current_row=0; current_row<m1; current_row++){

        //analysing where the current row lies in the CSR matrix
        int cur_row_ptr = (*csrMat)->row_ptr[current_row];
        int size2 = (*csrMat)->row_ptr[current_row+1] - (*csrMat)->row_ptr[current_row];
        if(size2 != 0){ //we have non-zero entries in this row

          //res will store the result of multiplication of current row with cvec
          //it will be one of the entries in the resultant matrix
          int res = 0;
          for(int k = cur_row_ptr; k<cur_row_ptr+size2; k++){
            *x = 0;
            *y = 0;
            if(multiplyOverflow(((*csrMat)->val[k]), (cvec[(*csrMat)->col_ind[k]]), x)){
              printf("\nInteger overflow due to multiplication");
              exit(0);
            }
            if(additionOverflow(*x, res, y)){
              printf("\nInteger overflow due to addition");
              exit(0);
            }

            res = *y;
          }

          //reallocate more memory if we need to expand the COO matrix
          if(coo_ind >= ((*cooMat)->nz)){
            printf("\nCurrent size of COO MAT is not enough, doubling current size from %d to %d", ((*cooMat)->nz), 2*((*cooMat)->nz));
            reallocCooMat(cooMat, ((*cooMat)->nz)*2);
          }

          //insert the res value only if it is non-zero
          if(res!=0){
            ((*cooMat)->val)[coo_ind] = res;
            ((*cooMat)->row_ind)[coo_ind] = current_row;
            ((*cooMat)->col_ind)[coo_ind] = current_col;
            coo_ind++;
          }

        }

      }

      printf("\nValue of coo_ind after iterating through one cvec: %d", coo_ind);

    }

  }

  //final reallocation to remove any additional memory which may have been allocated
  printf("\nChanging current size of Coo Mat from %d to %d", ((*cooMat)->nz), coo_ind);
  reallocCooMat(cooMat, coo_ind);

  //converting the resultant COO matrix to row major
  convertCooToRowMajor(cooMat, (*cooMat)->nz);


  free(x);
  free(y);

}


/*
  Function: Multiplies two matrices A (in CSC format) and B (in CSR format)
            The resultant matrix is stored in the last input parameter in form of COO.
  Inputs:
    struct CsrMat **csrMat:  Input matrix A in CSC format
    struct CscMat **cscMat:  Input matrix B in CSR format

  Outputs:
    struct CooMat ** cooMat:  Address of the pointer to the object of COO Matrix structure
                              The resultant matrix will be stored using this.
*/
void multiply_csc_csr(struct CscMat **cscMat, struct CsrMat **csrMat, struct CooMat **cooMat){

  //when multiplying CSC with CSR, we must check for matching sizes
  if((*cscMat)->n_col != (*csrMat)->n_row){
    printf("\n Cannot multiply matrices as the number of cols in first matrix don't match number of rows in second matrix.\n");
    exit(0);
  }

  //pointers to various structures to hold intermediated forms of matrices
  struct CooMat *intermediateCooMat1;
  struct CooMat *intermediateCooMat2;
  struct CsrMat *csrMatA;
  struct CscMat *cscMatB;

  //convert the given CSC matrix to CSR
  convert_csc_coo(cscMat, &intermediateCooMat1);
  convertCooToRowMajor(&intermediateCooMat1, intermediateCooMat1->nz);
  convert_coo_csr(&intermediateCooMat1, &csrMatA);

  //convert the given CSR matrix to CSC
  convert_csr_coo(csrMat, &intermediateCooMat2);
  convert_coo_csc(&intermediateCooMat2, &cscMatB);

  //now we can use the same function to multiply CSC with CSR
  multiply_csr_csc(&csrMatA, &cscMatB, cooMat);

  //free the local pointers
  free(intermediateCooMat1);
  free(intermediateCooMat2);
  free(csrMatA);
  free(cscMatB);

}


/*
  Function: Creates a elementary row operation matrix in CSR format
            The generated matrix is stored in the last input parameter.
  Inputs:
    int r1: Represents Row1 of the row transformation operation
    int k1: Represents the coefficient of Row1 in the row transformation operation
    int r2: Represents Row2 of the row transformation operation
    int k2: Represents the coefficient of Row2 in the row transformation operation

  Outputs:
    struct CsrMat ** csrMat:  Address of the pointer to the object of CSR Matrix structure
                              The created matrix will be stored using this.
*/
void createElementaryCsrMat(int r1, int k1, int r2, int k2, int n, struct CsrMat **csrMat){

  //allocating memory for the CSR matrix; will be free back in calling function
  *csrMat = (struct CsrMat *) malloc(sizeof(struct CsrMat));
  if(NULL == csrMat){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }

  (*csrMat)->val = (int *) malloc((n+1)*sizeof(int));
  if(NULL == (*csrMat)->val){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }
  (*csrMat)->row_ptr = (int *) malloc((n+1)*sizeof(int));
  if(NULL == (*csrMat)->row_ptr){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }
  (*csrMat)->col_ind = (int *) malloc((n+1)*sizeof(int));
  if(NULL == (*csrMat)->col_ind){
    printf("\nUnable to allocate memory for CSR struct");
    exit(0);
  }
  (*csrMat)->nz = n+1;
  (*csrMat)->n_row = n;
  (*csrMat)->n_col = n;


  //filling all values upto row r2; all such rows will have value 1 in diagonal entries
  int row_ptr_ct=0;
  int i=0;
  (*csrMat)->row_ptr[i] = 0;
  for(i; i<r2; i++){
    row_ptr_ct++;
    (*csrMat)->val[i] = 1;
    (*csrMat)->row_ptr[i+1] = row_ptr_ct;
    (*csrMat)->col_ind[i] = i;
  }

  //for row r2, we will have two values one in r2 and one in r1; both values will be k2 and k1 respectively
  if(i==r2){
    row_ptr_ct+=2;
    (*csrMat)->row_ptr[i+1] = row_ptr_ct;
    if(r1<r2){ //if r1<r2 meaninng k1 value will be first because we need to sort COO in row major order
      (*csrMat)->val[i] = k1;
      (*csrMat)->val[i+1] = k2;
      (*csrMat)->col_ind[i] = r1;
      (*csrMat)->col_ind[i+1] = r2;
    }
    else{
      (*csrMat)->val[i] = k2;
      (*csrMat)->val[i+1] = k1;
      (*csrMat)->col_ind[i] = r2;
      (*csrMat)->col_ind[i+1] = r1;
    }
  }

  //filling in remaing diagonal values as 1 in the matrix
  for(i=r2+2; i<n+1; i++){
    row_ptr_ct++;
    (*csrMat)->val[i] = 1;
    (*csrMat)->row_ptr[i] = row_ptr_ct;
    (*csrMat)->col_ind[i] = i-1;
  }

  //printf("\nCreated the following elementary CSR Matrix:\n");
  //printCsrMatFull(csrMat);
}


/*
  Function: Performs a row transformation on a matrix in CSC format
            The new matrix is stored in the last input parameter.
  Inputs:
    int r1: Represents Row1 of the row transformation operation
    int k1: Represents the coefficient of Row1 in the row transformation operation
    int r2: Represents Row2 of the row transformation operation
    int k2: Represents the coefficient of Row2 in the row transformation operation
    struct CscMat ** cscMat:  Address of the pointer to the object of CSC Matrix structure
                              This is the matrix on which the row operation will be performed

  Outputs:
    struct CscMat ** resCscMat: Address of the pointer to the object of CSC Matrix structure
                                The new matrix will be stored using this.
*/
void rowTransformationCscMat(int r1, int k1, int r2, int k2, struct CscMat **cscMat, struct CscMat **resCscMat){

  //checking if the row transformation parameters are valid
  if(r1<0 || r2<0 || r1>MAX_ROW_SIZE || r2>MAX_ROW_SIZE){
    printf("\nInvalid input for row transformation parameters.\n");
    exit(0);
  }

  //the no. of rows of CSC is unknown.
  //We can take it to be the maximum index of row in the matrix
  int m = 0;
  for(int i=0; i<((*cscMat)->nz); i++){
    if(((*cscMat)->row_ind[i]) > m)
      m = ((*cscMat)->row_ind[i]);
  }
  m++;

  if(m>MAX_ROW_SIZE){
    printf("\nNo. of rows in matrix exceed max row size allowed.\n");
    exit(0);
  }


  //now if r2 or r1 are greater than m, then we will assign m to be the max of them
  if((r1+1)>m){
    m = r1+1;
  }
  if((r2+1)>m){
    m = r2+1;
  }

  //creating a elementary matrix for the given row transformation in CSR format
  struct CsrMat *eleCsrMat;
  createElementaryCsrMat(r1, k1, r2, k2, m, &eleCsrMat);

  //multiplying the elementary matrix by the given CSC matrix to get the transformed matrix in COO format
  struct CooMat *cooMat;
  multiply_csr_csc(&eleCsrMat, cscMat, &cooMat);

  //converting the transformed matrix into CSC form
  convert_coo_csc(&cooMat, resCscMat);

  //freeing the memory for intermediate matrices
  free(eleCsrMat);
  free(cooMat);
}
