#include <vector>
#include <algorithm>
#include <iostream>

#include "amplapack_test.h"
#include "ampxlapack.h"

// host GEMM used for reconstruction
#include "lapack_host.h"

template <typename value_type>
value_type random_value(value_type min, value_type max)
{
    value_type val = value_type(rand()) / value_type(RAND_MAX);
    val *= (max-min);
    val += min;
    return val;
}

template <typename value_type>
double gflops(double sec, double m, double n)
{
    const double k = std::min(m,n);
    const double l = std::max(m,n);

    double operation_count = double(2)*l*k*k - double(2)/double(3)*k*k*k;

    return operation_count / (sec * double(1e9));
}

template <typename value_type>
void do_geqrf_test(int m, int n, int lda_offset = 0)
{
    // header
    std::cout << "Testing xGEQRF for M=" << m << " N=" << n << " LDA=" << m+lda_offset << "... ";

    // performance timer
    high_resolution_timer timer;

    // create data
    int k = std::min(m,n);
    int lda = m + lda_offset;
    std::vector<value_type> a(lda*n, value_type(1));
    std::vector<int> ipiv(k);

    // fill with random values
    std::for_each(a.begin(), a.end(), [&](value_type& val) {
        val = random_value(value_type(-1), value_type(1));
    });

    // mark data outside of the leading dimension for debugging purposes
    for (int j = 0; j < n; j++)
        for (int i = 0; i < lda; i++)
            if (i >= n)
                a[j*lda+i] = value_type(-1);

    // backup a for reconstruction purposes
    std::vector<value_type> a_in(a);

    // tau input
    std::vector<value_type> tau(k); 

    // run routine
    int info;

    timer.restart();
    amplapack_status status = amplapack_geqrf(m, n, a.data(), lda, tau.data(), &info);
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
        // extract v/q and r 
        std::vector<value_type> r = a; 
        int ldr = lda;
        std::vector<value_type> q(m*m, value_type());
        int ldq = m;
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                if (i > j)
                {
                    q[ldq*j+i] = a[lda*j+i];
                    r[lda*j+i] = value_type();
                }
            }
        }

        // work query
        int info;
        int lwork = -1;
        value_type work_size;
        LAPACK_SORGQR(&m, &m, &k, q.data(), &ldq, tau.data(), &work_size, &lwork, &info);
        lwork = int(work_size);
        std::vector<value_type> work(lwork);

        // generate q
        LAPACK_SORGQR(&m, &m, &k, q.data(), &ldq, tau.data(), work.data(), &lwork, &info);

        // a = a - qr
        gemm('n', 'n', m, k, m, value_type(1), q.data(), ldq, r.data(), ldr, value_type(-1), a_in.data(), lda);

        // norm
        std::cout << " Error = " << one_norm(n, n, a_in.data(), lda);

        // gflops
        std::cout << " GFLOPs = " << gflops<value_type>(1.0, m,n) << std::endl; 
    }
}

void geqrf_test()
{
    // quick tests
    do_geqrf_test<float>(1024, 1024); 
}