#define MAX_ROW_SIZE 150
#define MAX_COL_SIZE 100
#define MAX_NZ 15000

//structure defining the COO format of sparse matrix
struct CooMat{
  int * val;
  int * row_ind;
  int * col_ind;
  int nz;
  int n_row;
  int n_col;
};

//structure defining the CSR format of sparse matrix
struct CsrMat{
  int * val;
  int * row_ptr;
  int * col_ind;
  int nz;
  int n_row;
  int n_col;
};

//structure defining the CSC format of sparse matrix
struct CscMat{
  int * val;
  int * row_ind;
  int * col_ptr;
  int nz;
  int n_row;
  int n_col;
};

//Functions to print the matrices of different formats in full format
void printCooMatFull(struct CooMat ** cooMat);
void printCsrMatFull(struct CsrMat ** csrMat);
void printCscMatFull(struct CscMat ** cscMat);

//functions to print the matrices in specified sparse formats
void printCooMat(struct CooMat ** cooMat);
void printCsrMat(struct CsrMat ** csrMat);
void printCscMat(struct CscMat ** cscMat);

//functions to read the sparse matrices in different formats
void readCooMat(struct CooMat ** cooMat);
void readCsrMat(struct CsrMat ** csrMat);
void readCscMat(struct CscMat ** cscMat);


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
void generateCooMat(int m, int n, double density_nz, double density_negative, struct CooMat ** cooMat);


//functions to convert a COO matrix to col major and row major format
void convertCooToColumnMajor(struct CooMat **cooMat, int n);
void convertCooToRowMajor(struct CooMat **cooMat, int n);


//functions to convert a matrix from COO format to CSR and CSC format
void convert_coo_csr(struct CooMat **cooMat, struct CsrMat **csrMat);
void convert_coo_csc(struct CooMat **cooMat, struct CscMat **cscMat);


//functions to convert a matrix from CSC/CSR format to COO format
void convert_csc_coo(struct CscMat **cscMat, struct CooMat **cooMat);
void convert_csr_coo(struct CsrMat **csrMat, struct CooMat **cooMat);


//function to reallocate the size of COO mat to a new size
void reallocCooMat(struct CooMat **cooMat, int size);


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
void multiply_csr_csc(struct CsrMat **csrMat, struct CscMat **cscMat, struct CooMat **cooMat);


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
void multiply_csc_csr(struct CscMat **cscMat, struct CsrMat **csrMat, struct CooMat **cooMat);


/*
  Function: Creates a elementary row operation matrix in CSR format
            The generated matrix is stored in the last input parameter.
  Inputs:
    int r1: Represents Row1 of the row transformation operation
    double k1: Represents the coefficient of Row1 in the row transformation operation
    int r2: Represents Row2 of the row transformation operation
    double k2: Represents the coefficient of Row2 in the row transformation operation

  Outputs:
    struct CsrMat ** csrMat:  Address of the pointer to the object of CSR Matrix structure
                              The created matrix will be stored using this.
*/
void createElementaryCsrMat(int r1, int k1, int r2, int k2, int n, struct CsrMat **csrMat);


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
void rowTransformationCscMat(int r1, int k1, int r2, int k2, struct CscMat **cscMat, struct CscMat **resCscMat);
