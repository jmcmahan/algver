#include "linalgwrapper.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <cstring>



/////////////////////////////////////////////////// 
///////////////// LAW_Vec Members ///////////////// 

/* Allocate space for the vector */
inline void LAW_Vec::alloc_space()
{
    m_x = new double[m_N];
}


/* Free space for the vector on destruction */
inline void LAW_Vec::free_space()
{
    delete m_x;
}

double LAW_Vec::inprod(LAW_Vec *y)
{
    double r[1];
    gsl_vector_view xv = gsl_vector_view_array(m_x, m_N);
    gsl_vector_view yv = gsl_vector_view_array(y->m_x, y->m_N);
    gsl_blas_ddot(&xv.vector, &yv.vector, r);

    return r[0];
}

/* Compute z = x - y where y is this vector, x,z are the arguments */
void LAW_Vec::subtract(LAW_Vec *x, LAW_Vec *z)
{

    gsl_vector_view xv = gsl_vector_view_array(x->m_x, x->m_N);
    gsl_vector_view yv = gsl_vector_view_array(m_x, m_N);

    gsl_vector_sub(&xv.vector, &yv.vector);
    //gsl_blas_daxpy(-1.0, &xv.vector, &yv.vector);
    memcpy(z->m_x, gsl_vector_ptr(&(xv.vector), 0), sizeof(double)*(m_N));
}


/* Compute z = x + y where y is this vector, x,z are the arguments */
void LAW_Vec::add(LAW_Vec *x, LAW_Vec *z)
{

    gsl_vector_view xv = gsl_vector_view_array(x->m_x, x->m_N);
    gsl_vector_view yv = gsl_vector_view_array(m_x, m_N);

    gsl_vector_add(&xv.vector, &yv.vector);
    //gsl_blas_daxpy(-1.0, &xv.vector, &yv.vector);
    memcpy(z->m_x, gsl_vector_ptr(&(xv.vector), 0), sizeof(double)*(m_N));
}



/* Print MATLAB code containing the vector for verification */
void LAW_Vec::Print(const char *varname, std::ostream &os)
{
    int i, j;


    os << varname << "=[";
    //printf("%s = [", varname);
    for (i=0; i<m_N; i++) {
        //printf("%1.6e;", m_x[i]);
        os << m_x[i] << ";";
    }
    //printf("];\n");
    os << "];\n";
}

// Constructors / Destructors 

LAW_Vec::LAW_Vec(unsigned int N)
{
    m_N = N;
    alloc_space();
}


LAW_Vec::LAW_Vec(unsigned int N, double *x)
{
    m_N = N;
    alloc_space();
    memcpy(m_x, x, sizeof(double) * m_N);
}


LAW_Vec::~LAW_Vec()
{
    free_space();
}


/////////////////////////////////////////////////// 
///////////////// LAW_Mat Members ///////////////// 

 
/* Allocate space for the matrix */
inline void LAW_Mat::alloc_space()
{
    m_a = new double[m_M * m_N];
}


/* Free space for the matrix on destruction */
inline void LAW_Mat::free_space()
{
    delete m_a;
}

/* Return pointer to a particular index (starting from 0) */
double *LAW_Mat::ind(unsigned int i, unsigned int j)
{
    return &m_a[i*m_N + j];
}

/* Return pointer to a particular index (starting from 0) for the transpose
   of the matrix. */
double *LAW_Mat::trind(unsigned int i, unsigned int j)
{

    return &m_a[j*m_N + i];
}

/* Compute C := a*A*B + b*C  where A is the matrix instance  */
void LAW_Mat::lmult(LAW_Mat *B, LAW_Mat *C, double a, double b, bool transA, bool transB)
{
    
    CBLAS_TRANSPOSE_t At;
    CBLAS_TRANSPOSE_t Bt;

    At = (transA ? CblasTrans : CblasNoTrans);
    Bt = (transB ? CblasTrans : CblasNoTrans);

#if 1

    gsl_matrix_view Av = gsl_matrix_view_array(m_a, m_M, m_N);
    gsl_matrix_view Bv = gsl_matrix_view_array(B->m_a, B->m_M, B->m_N);
    gsl_matrix_view Cv = gsl_matrix_view_array(C->m_a, C->m_M, C->m_N);
    gsl_blas_dgemm(At, Bt, a, &Av.matrix, &Bv.matrix, b, &Cv.matrix);
    memcpy(C->m_a, gsl_matrix_ptr(&(Cv.matrix), 0, 0), sizeof(double)*(C->m_M)*(C->m_N));
#else
    gsl_matrix *Am;
    gsl_matrix *Bm;
    gsl_matrix *Cm;

    Am = gsl_matrix_alloc(m_M, m_N);
    Bm = gsl_matrix_alloc(B->m_M, B->m_N);
    Cm = gsl_matrix_alloc(C->m_M, C->m_N);

    memcpy(gsl_matrix_ptr(Am, 0, 0), m_a, sizeof(double)*m_M*m_N);
    memcpy(gsl_matrix_ptr(Bm, 0, 0), B->m_a, sizeof(double)*(B->m_M)*(B->m_N));

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
                    Am, Bm, 0.0, Cm);

    memcpy(C->m_a, gsl_matrix_ptr(Cm, 0, 0), sizeof(double)*(C->m_M)*(C->m_N));
    gsl_matrix_free(Am);
    gsl_matrix_free(Bm);
    gsl_matrix_free(Cm);
#endif
}

/* Compute C := A*B where A is the matrix instance*/
void LAW_Mat::lmult(LAW_Mat *B, LAW_Mat *C)
{
    lmult(B, C, 1.0, 0.0, false, false);
}

/* Compute y = A*x where A is the matrix instance */
void LAW_Mat::matvec(LAW_Vec *x, LAW_Vec *y)
{
    gsl_matrix_view Av = gsl_matrix_view_array(m_a, m_M, m_N);
    gsl_vector_view xv = gsl_vector_view_array(x->m_x, x->m_N);
    gsl_vector_view yv = gsl_vector_view_array(y->m_x, y->m_N);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &Av.matrix, &xv.vector, 0.0, &yv.vector);
    memcpy(y->m_x, gsl_vector_ptr(&(yv.vector), 0), sizeof(double)*(y->m_N));
}

/* Compute y = A^T*x where A is the matrix instance */
void LAW_Mat::matvec(LAW_Vec *x, LAW_Vec *y, bool transA)
{
    CBLAS_TRANSPOSE_t At;
    At = (transA ? CblasTrans : CblasNoTrans);

    gsl_matrix_view Av = gsl_matrix_view_array(m_a, m_M, m_N);
    gsl_vector_view xv = gsl_vector_view_array(x->m_x, x->m_N);
    gsl_vector_view yv = gsl_vector_view_array(y->m_x, y->m_N);
    gsl_blas_dgemv(At, 1.0, &Av.matrix, &xv.vector, 0.0, &yv.vector);
    memcpy(y->m_x, gsl_vector_ptr(&(yv.vector), 0), sizeof(double)*(y->m_N));
}

/* Invert the matrix with the Cholesky decomposition. Need matrix to be spd */
void LAW_Mat::cholinv()
{
    gsl_matrix_view Av = gsl_matrix_view_array(m_a, m_M, m_N);
    // Cholesky decomposition is needed before inversion
    gsl_linalg_cholesky_decomp(&Av.matrix);
    gsl_linalg_cholesky_invert(&Av.matrix);
    memcpy(m_a, gsl_matrix_ptr(&(Av.matrix), 0, 0), sizeof(double)*(m_M)*(m_N));

}

/* Solve A*x = b, resulting in x = inv(A)*b where A is the matrix instance.
   To be clear, b should be the given value here and the result is returned
   in x.  */
void LAW_Mat::cholsolve(LAW_Vec *b, LAW_Vec *x)
{
    gsl_matrix_view Av = gsl_matrix_view_array(m_a, m_M, m_N);
    gsl_vector_view xv = gsl_vector_view_array(x->m_x, x->m_N);
    gsl_vector_view bv = gsl_vector_view_array(b->m_x, b->m_N);
    // Cholesky decomposition is needed before inversion
    gsl_linalg_cholesky_decomp(&Av.matrix);
    gsl_linalg_cholesky_solve(&Av.matrix, &bv.vector, &xv.vector);
    memcpy(x->m_x, gsl_vector_ptr(&(xv.vector), 0), sizeof(double)*(x->m_N));
}


/* Print MATLAB code containing the matrix for verification. varname is
   the name of the matlab variable.  */
void LAW_Mat::Print(const char *varname, std::ostream &os)
{
    int i, j;


    //printf("%s = [", varname);
    os << varname << "=[";
    for (i=0; i<m_M; i++) {
        for (j=0; j<m_N; j++) {
            //printf("%1.6e ", m_a[i*m_N + j]);
            os << m_a[i*m_N + j] << " ";
        }
        os << ";\n";
        //printf(";\n");
    }
    os << "];\n";
    //printf("];\n");
}
// Constructors / Destructors

LAW_Mat::LAW_Mat(unsigned int M, unsigned int N)
{
    m_M = M;
    m_N = N;
    alloc_space();
}


LAW_Mat::LAW_Mat(unsigned int M, unsigned int N, double *a)
{
    m_M = M;
    m_N = N;
    alloc_space();
    memcpy(m_a, a, sizeof(double) * m_M * m_N);
}


LAW_Mat::~LAW_Mat()
{
    free_space();
}

