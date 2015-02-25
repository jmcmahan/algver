#ifndef _LINALGWRAPPER_HPP_
#define _LINALGWRAPPER_HPP_

#include <iostream>
#include <fstream>

/* Minimal code to represent matrices, vectors, and their operations 
   generically */


class LAW_Vec {
    private:
        void alloc_space();
        void free_space();

    public:
        double *m_x;            // Array of vector components
        unsigned int m_N;       // Size of the vector

        double inprod(LAW_Vec*);
        void subtract(LAW_Vec*, LAW_Vec*);
        void add(LAW_Vec*, LAW_Vec*);

        void Print(const char *varname, std::ostream &os);

        // Constructors / Destructors
        LAW_Vec(unsigned int);
        LAW_Vec(unsigned int, double*);
        ~LAW_Vec();
};

class LAW_Mat {
    private:
        void alloc_space();
        void free_space();

    public:
        double *m_a;            // Array of matrix values (column index
                                // increments first).
        unsigned int m_M;       // Number of rows
        unsigned int m_N;       // Number of columns

        double *ind(unsigned int i, unsigned int j);
        double *trind(unsigned int i, unsigned int j);
        void lmult(LAW_Mat *B, LAW_Mat *C, double a, double b, bool transA, bool transB);
        void lmult(LAW_Mat *B, LAW_Mat *C);
        void matvec(LAW_Vec*, LAW_Vec*);
        void matvec(LAW_Vec*, LAW_Vec*, bool transA);
        void cholinv();
        void cholsolve(LAW_Vec *b, LAW_Vec *x);

        void Print(const char *varname, std::ostream &os);

        // Constructors / Destructors
        LAW_Mat(unsigned int, unsigned int);
        LAW_Mat(unsigned int, unsigned int, double*);
        ~LAW_Mat();
};

#endif /* _LINALGWRAPPER_HPP_ */
