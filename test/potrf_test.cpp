#include <vector>
#include <algorithm>
#include <iostream>

#include "amplapack_test.h"
#include "ampxlapack.h"

// host GEMM used for reconstruction
#include "lapack_host.h"

template <typename value_type>
double gflops(double sec, double n)
{
    double operation_count = double(1)/double(3)*n*n*n;

    return operation_count / (sec * double(1e9));
}

template <typename value_type>
void do_potrf_test(char uplo, int n, int lda_offset = 0)
{
    // header
    std::cout << "Testing xPOTRF for UPLO=" << uplo << " N=" << n << " LDA=" << n+lda_offset << "... ";

    // performance timer
    high_resolution_timer timer;

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

    timer.restart();
    amplapack_status status = amplapack_potrf(uplo, n, &a.at(0), lda, &info);
    double sec = timer.elapsed();

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

    if (status == amplapack_success)
    {
        // reconstruction
        if (uplo == 'L' || uplo == 'l')
        {
            // reflect Ain
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    if (j > i)
                    {
                        a_in[j*lda+i] = a_in[i*lda+j];
                    }
                }
            }

            // a = a - l*l'
            gemm('n', 't', n, n, n, value_type(1), a.data(), lda, a.data(), lda, value_type(-1), a_in.data(), lda);
        }
        else
        {   
            // reflect Ain
            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    if (i > j)
                    {
                        a_in[j*lda+i] = a_in[i*lda+j];
                    }
                }
            }

            // a = a - u'*u
            gemm('t', 'n', n, n, n, value_type(1), a.data(), lda, a.data(), lda, value_type(-1), a_in.data(), lda);
        }

        // norm
        std::cout << " Error = " << one_norm(n, n, a_in.data(), lda) << " GLFOPs = " << gflops<value_type>(sec, n) << std::endl;
    }
}

void potrf_test()
{
    // performance tests
    do_potrf_test<float>('L', 1024);
    do_potrf_test<float>('U', 1024);
}