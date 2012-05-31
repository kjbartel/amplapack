#pragma once
#ifndef AMPLAPACK_TEST_H
#define AMPLAPACK_TEST_H

#include <limits>

// matrix-matrix multiply is used in a number of reconstruction algorithms
template <typename value_type>
void gemm(char transa, char transb, int m, int n, int k, value_type alpha, const value_type* a, int lda, const value_type* b, int ldb, value_type beta, value_type* c, int ldc);

// returns the one-norm of an m by n matrix
template <typename value_type>
value_type one_norm(int m, int n, value_type* a, int lda)
{
    value_type norm = 0;
    value_type epsilon = std::numeric_limits<value_type>::epsilon();

    for (int j = 0; j < n; j++)
    {
        value_type sum = 0;
        for (int i = 0; i < m; i++)
            sum += abs( a[i+j*lda] );

        norm = std::max(norm, sum);
    }

    return (norm*epsilon) / n;
}

#endif // AMPLAPACK_TEST_H
