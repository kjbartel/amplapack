#pragma once
#ifndef AMPLAPACK_TEST_H
#define AMPLAPACK_TEST_H

#include <limits>
#include <memory>

// test listing
void potrf_test();
void getrf_test();
void geqrf_test();

// matrix-matrix multiply is used in a number of reconstruction algorithms
template <typename value_type>
void gemm(char transa, char transb, int m, int n, int k, value_type alpha, const value_type* a, int lda, const value_type* b, int ldb, value_type beta, value_type* c, int ldc);

// row swap is used in the GETRF reconstruction
template <typename value_type>
void laswp(int n, value_type* a, int lda, int k1, int k2, int* ipiv, int incx);

// returns the one-norm of an m by n matrix
template <typename value_type>
value_type one_norm(int m, int n, value_type* a, int lda)
{
    value_type norm = value_type();

    for (int j = 0; j < n; j++)
    {
        value_type sum = 0;
        for (int i = 0; i < m; i++)
            sum += abs( a[j*lda+i] );

        norm = std::max(norm, sum);
    }

    return norm/n;
}

class high_resolution_timer
{
public:
    high_resolution_timer();
    ~high_resolution_timer();

    void restart();
    double elapsed();

private:
    struct impl;
    std::unique_ptr<impl> pimpl;
};

#endif // AMPLAPACK_TEST_H
