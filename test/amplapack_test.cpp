#include "amplapack_test.h"
#include "lapack_host.h"

template <>
void gemm(char transa, char transb, int m, int n, int k, float alpha, const float* a, int lda, const float* b, int ldb, float beta, float* c, int ldc)
{
    LAPACK_SGEMM(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

template <>
void laswp(int n, float* a, int lda, int k1, int k2, int* ipiv, int incx)
{
    LAPACK_SLASWP(&n, a, &lda, &k1, &k2, ipiv, &incx);
}

int main()
{
    potrf_test();
    getrf_test();
    //geqrf_test();
}