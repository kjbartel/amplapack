#include <vector>
#include <algorithm>
#include <iostream>

#include "amplapack_test.h"
#include "ampxlapack.h"

// host GEMM used for reconstruction
#include "lapack_host.h"

template <typename value_type>
void potrf_test(char uplo, int n, int lda_offset = 0)
{
    // header
    std::cout << "Testing xPOTRF for UPLO=" << uplo << " N=" << n << " LDA=" << n+lda_offset << std::endl;

    // create data
    int lda = n + lda_offset;
    std::vector<value_type> a(lda*n, value_type(1));

    // scale diagonal
    for (int i = 0; i < (lda*n); i += (lda+1))
        a[i] = value_type(n);

    // obliterate unused half
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            if (i < j && uplo == 'L')
                a[j*lda+i] = value_type();
            else if (i > j && uplo == 'U')
                a[j*lda+i] = value_type();

    // mark data outside of the leading dimension for debugging purposes
    for (int j = 0; j < n; j++)
        for (int i = 0; i < lda; i++)
            if (i >= n)
                a[j*lda+i] = value_type(-1);

    // backup a for reconstruction purposes
    std::vector<value_type> a_in(a);
    
    int info;
    amplapack_status status = amplapack_potrf(uplo, n, &a.at(0), lda, &info);

    switch(status)
    {
    case amplapack_success:
        std::cout << "Success!";
        break;
    case amplapack_data_error:
        std::cout << "Date Error @ " << info << std::endl;
        break;
    case amplapack_argument_error:
        std::cout << "Argument Error @ " << -info << std::endl;
        break;
    case amplapack_runtime_error:
        std::cout << "Runtime Error" << std::endl;
        break;
    case amplapack_memory_error:
        std::cout << "Insuffecient Memory" << std::endl;
        break;
    default:
        std::cout << "Unexpected Error?" << std::endl;
        break;
    }

    if ( status == amplapack_success )
    {
        // reconstruction
        if (uplo == 'L' || uplo == 'l')
        {
            // reflect Ain
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++)
                    if (i < j)
                        a_in[j*lda+i] = a_in[i*lda+j];

            // a = a - l*l'
            gemm('n', 't', n, n, n, value_type(1), a.data(), lda, a.data(), lda, value_type(-1), a_in.data(), lda);
        }
        else
        {   
            // reflect Ain
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++)
                    if (i > j)
                        a_in[i*lda+j] = a_in[j*lda+i];

            // a = a - u'*u
            gemm('t', 'n', n, n, n, value_type(1), a.data(), lda, a.data(), lda, value_type(-1), a_in.data(), lda);
        }

        // norm
        std::cout << " Error = " << one_norm(n, n, a_in.data(), lda) << std::endl;
    }
}

int main()
{
    // quick tests
    potrf_test<float>('L', 256);
    potrf_test<float>('U', 256);
    potrf_test<float>('L', 512, 2);
    potrf_test<float>('U', 512, 2);

    // performance tests
    potrf_test<float>('L', 4096);
    
    // requires full doubles!
    // potrf_test<double>('U', 256);
}