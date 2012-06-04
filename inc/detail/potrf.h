/*----------------------------------------------------------------------------
* Copyright © Microsoft Corp.
*
* Licensed under the Apache License, Version 2.0 (the "License"); you may not 
* use this file except in compliance with the License.  You may obtain a copy 
* of the License at http://www.apache.org/licenses/LICENSE-2.0  
* 
* THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
* KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED 
* WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE, 
* MERCHANTABLITY OR NON-INFRINGEMENT. 
*
* See the Apache Version 2.0 License for specific language governing 
* permissions and limitations under the License.
*---------------------------------------------------------------------------
* 
* potrf.h
*
*---------------------------------------------------------------------------*/

#include <amp.h>

#include "ampblas.h"
#include "amplapack.h"

#include "amplapack_runtime.h"
#include "lapack_host.h"

// external lapack functions

namespace amplapack {
namespace _detail {

//
// External LAPACK Wrappers
// 

namespace lapack {

template <typename value_type>
void potrf(char uplo, int n, value_type* a, int lda, int& info);

template <>
void potrf(char uplo, int n, float* a, int lda, int& info)
{ 
    LAPACK_SPOTRF(&uplo, &n, a, &lda, &info); 
}

template <>
void potrf(char uplo, int n, double* a, int lda, int& info)
{ 
    LAPACK_DPOTRF(&uplo, &n, a, &lda, &info); 
}

template <>
void potrf(char uplo, int n, ampblas::complex<float>* a, int lda, int& info)
{ 
    LAPACK_CPOTRF(&uplo, &n, a, &lda, &info); 
}

template <>
void potrf(char uplo, int n, ampblas::complex<double>* a, int lda, int& info)
{ 
    LAPACK_ZPOTRF(&uplo, &n, a, &lda, &info); 
}

} // namespace lapack

//
// Host Wrapper
//

namespace host {

template <enum class option::ordering storage_type, typename value_type>
void potrf(const concurrency::accelerator_view& /*av*/, enum class option::uplo uplo, concurrency::array_view<value_type,2>& a)
{
    static_assert(storage_type == option::ordering::column_major, "hybrid functionality requires column major ordering");

    const int n = require_square(a);
   
    // TODO: take from a pool
    const int lda = n;
    std::vector<value_type> hostVector(lda*n);

    // copy from acclerator to host
    concurrency::copy(a, hostVector.begin());

    // run host function
    int info = 0;
    lapack::potrf(to_char(uplo), n, hostVector.data(), lda, info);

    // check for errors
    info_check(info);

    // copy from host to accelerator
    // requires -D_SCL_SECURE_NO_WARNINGS
    concurrency::copy(hostVector.begin(), hostVector.end(), a);
}

} // namespace host

//
// Unblocked Factorization
//

template <enum class option::ordering storage_type, typename value_type>
void potf2(const concurrency::accelerator_view& av, enum class option::uplo uplo, concurrency::array_view<value_type,2>& a)
{
    using concurrency::array_view;
    using concurrency::index;
    using concurrency::extent;

    typedef ampblas::real_type<value_type>::type real_type;

    const int n = require_square(a);

    if (uplo == option::uplo::upper)
    {
        for (int j = 0; j < n; j++)
        {
            // compute u(j,j)
            real_type a_jj = ampblas::real(a[concurrency::index<2>(j,j)]);
            value_type dot = value_type();

            if (j != 0)
            {
                array_view<const value_type,2> a_vector_2d = get_sub_matrix<storage_type>(a, index<2>(0,j), extent<2>(j,1));
                auto a_vector_1d = make_subvector_view(a_vector_2d);
                dot = ampblas::dot<value_type>(av, a_vector_1d, a_vector_1d);
            }

            a_jj -= ampblas::real(dot);

            // check for non-positive-definiteness
            // +1 for Fortran indexing
            if (a_jj <= real_type())
                data_error(j+1);

            a_jj = sqrt(a_jj);
            a[concurrency::index<2>(j,j)] = value_type(a_jj);

            // compute elements j+1:n of column j
            if (j < n)
            {
                int m_ = j;
                int n_ = n-j-1;

                if (n_ > 0)
                {
                    array_view<value_type,2> y_vector_2d = get_sub_matrix<storage_type>(a, index<2>(j,j+1), extent<2>(1,n_));
                    auto y_vector_1d = make_subvector_view(y_vector_2d);
                    
                    if (m_ > 0)
                    {
                        array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(0,j+1), extent<2>(m_,n_));
                        array_view<const value_type,2> x_vector_2d = get_sub_matrix<storage_type>(a, index<2>(0,j), extent<2>(m_,1));
                        auto x_vector_1d = make_subvector_view(x_vector_2d);
                        
                        ampblas::gemv(av, AmpblasTrans, value_type(-1), a_sub, x_vector_1d, value_type(1), y_vector_1d);
                    }

                    ampblas::scal(av, value_type(1)/a_jj, y_vector_1d);
                }
            }
        }
    }
    else if (uplo == option::uplo::lower)
    {
        for (int j = 0; j < n; j++)
        {
            // compute l(j,j)
            real_type a_jj = ampblas::real(a[concurrency::index<2>(j,j)]);
            value_type dot = value_type();

            if (j != 0)
            {
                array_view<const value_type,2> a_vector_2d = get_sub_matrix<storage_type>(a, index<2>(j,0), extent<2>(1,j));
                auto a_vector_1d = make_subvector_view(a_vector_2d);
                dot = ampblas::dot<value_type>(av, a_vector_1d, a_vector_1d);
            }

            a_jj -= ampblas::real(dot);

            // check for non-positive-definiteness
            // +1 for Fortran indexing
            if (a_jj <= real_type())
                data_error(j+1);

            a_jj = sqrt(a_jj);
            a[concurrency::index<2>(j,j)] = value_type(a_jj);

            // compute elements j+1:n of column j
            if (j < n)
            {
                int m_ = n-j-1;
                int n_ = j;

                if (m_ > 0)
                {
                    array_view<value_type,2> y_vector_2d = get_sub_matrix<storage_type>(a, index<2>(j+1,j), extent<2>(m_,1));
                    auto y_vector_1d = make_subvector_view(y_vector_2d);
                    
                    if (n_ > 0)
                    {
                        array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(j+1,0), extent<2>(m_,n_));
                        array_view<const value_type,2> x_vector_2d = get_sub_matrix<storage_type>(a, index<2>(j,0), extent<2>(1,n_));
                        auto x_vector_1d = make_subvector_view(x_vector_2d);
                        
                        ampblas::gemv(av, AmpblasNoTrans, value_type(-1), a_sub, x_vector_1d, value_type(1), y_vector_1d);
                    }

                    ampblas::scal(av, value_type(1)/a_jj, y_vector_1d);
                }
            }
        }
    }
}

//
// Blocked Factorization
//

template <int block_size, int look_ahead_depth, enum class option::ordering storage_type, enum class block_factor_location location, typename value_type>
void potrf(const concurrency::accelerator_view& av, enum class option::uplo uplo, const concurrency::array_view<value_type,2>& a)
{
    using concurrency::array_view;
    using concurrency::index;
    using concurrency::extent;

    // matrix size
    const int n = require_square(a);

    // stored int lower triangular matrix
    if (uplo == option::uplo::upper)
    {
        // block stepping
        for (int j = 0; j < n; j += block_size)
        {
            // current block size
            int jb = std::min(block_size, n-j);

            // update diagonal block
            {
                int n_ = jb;
                int k_ = j;

                if (k_ > 0 && n_ > 0)
                {
                    array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(0,j), extent<2>(k_,n_));
                    array_view<value_type,2> c_sub = get_sub_matrix<storage_type>(a, index<2>(j,j), extent<2>(n_,n_));
                    ampblas::syrk(av, AmpblasUpper, AmpblasTrans, value_type(-1), a_sub, value_type(1), c_sub);
                }
            }

            // factorize current block
            try
            {
                int n_ = jb;
                array_view<value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(j,j), extent<2>(n_,n_));

                if (location == block_factor_location::accelerator)
                    potf2<storage_type>(av, uplo, a_sub);
                else if (location == block_factor_location::host)
                    host::potrf<storage_type>(av, uplo, a_sub);
            }
            catch(const data_error_exception& e)
            {
                // offset local block error
                data_error(e.get() + j);
            }

            // this currently has no look ahead optimizations
            if (j+jb < n)
            {
                // compute the current block column
                {
                    int m_ = jb;
                    int n_ = n-j-jb;
                    int k_ = j;
                    if (m_ > 0 && n_ > 0 && k_ > 0)
                    {
                        array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(0,j), extent<2>(k_,m_));
                        array_view<const value_type,2> b_sub = get_sub_matrix<storage_type>(a, index<2>(0,j+jb), extent<2>(k_,n_));
                        array_view<value_type,2> c_sub = get_sub_matrix<storage_type>(a, index<2>(j,j+jb), extent<2>(m_,n_));
                        ampblas::gemm(av, AmpblasTrans, AmpblasNoTrans, value_type(-1), a_sub, b_sub, value_type(1), c_sub);
                    }
                }

                // solve for this row
                {
                    int m_ = jb;
                    int n_ = n-j-jb;
                    if (m_ > 0 && n_ > 0)
                    {
                        array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(j,j), extent<2>(m_,m_));
                        array_view<value_type,2> b_sub = get_sub_matrix<storage_type>(a, index<2>(j,j+jb), extent<2>(m_,n_));
                        ampblas::trsm(av, AmpblasLeft, AmpblasUpper, AmpblasTrans, AmpblasNonUnit, value_type(1), a_sub, b_sub);
                    }
                }
            }
        }
    }
    else if (uplo == option::uplo::lower)
    {
        // block stepping
        for (int j = 0; j < n; j += block_size)
        {
            // current block size
            int jb = std::min(block_size, n-j);

            // update diagonal block
            {
                int n_ = jb;
                int k_ = j;

                if (k_ > 0 && n_ > 0)
                {
                    array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(j,0), extent<2>(n_,k_));
                    array_view<value_type,2> c_sub = get_sub_matrix<storage_type>(a, index<2>(j,j), extent<2>(n_,n_));
                    ampblas::syrk(av, AmpblasLower, AmpblasNoTrans, value_type(-1), a_sub, value_type(1), c_sub);
                }
            }

            // factorize current block
            try
            {
                int n_ = jb;
                array_view<value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(j,j), extent<2>(n_,n_));

                if (location == block_factor_location::accelerator)
                    potf2<storage_type>(av, uplo, a_sub);
                else if (location == block_factor_location::host)
                    host::potrf<storage_type>(av, uplo, a_sub);
            }
            catch(const data_error_exception& e)
            {
                // offset local block error
                data_error(e.get() + j);
            }

            // this currently has no look ahead optimizations
            if (j+jb < n)
            {
                // compute the current block column
                {
                    int m_ = n-j-jb;
                    int n_ = jb;
                    int k_ = j;
                    if (m_ > 0 && n_ > 0 && k_ > 0)
                    {
                        array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(j+jb,0), extent<2>(m_,k_));
                        array_view<const value_type,2> b_sub = get_sub_matrix<storage_type>(a, index<2>(j,0), extent<2>(n_,k_));
                        array_view<value_type,2> c_sub = get_sub_matrix<storage_type>(a, index<2>(j+jb,j), extent<2>(m_,n_));
                        ampblas::gemm(av, AmpblasNoTrans, AmpblasTrans, value_type(-1), a_sub, b_sub, value_type(1), c_sub);
                    }
                }

                // solve for this row
                {
                    int m_ = n-j-jb;
                    int n_ = jb;
                    if (m_ > 0 && n_ > 0)
                    {
                        array_view<const value_type,2> a_sub = get_sub_matrix<storage_type>(a, index<2>(j,j), extent<2>(n_,n_));
                        array_view<value_type,2> b_sub = get_sub_matrix<storage_type>(a, index<2>(j+jb,j), extent<2>(m_,n_));
                        ampblas::trsm(av, AmpblasRight, AmpblasLower, AmpblasTrans, AmpblasNonUnit, value_type(1), a_sub, b_sub);
                    }
                }
            }
        }
    }
}

//
// Forwarding Function
//

template <enum class option::ordering storage_type, typename value_type>
void potrf(const concurrency::accelerator_view& av, enum class option::uplo uplo, const concurrency::array_view<value_type,2>& a)
{
    // TODO: a tuning framework
    const int block_size = 4;
    const int look_ahead_depth = 1;

    _detail::potrf<block_size, look_ahead_depth, storage_type, block_factor_location::host>(av, uplo, a);
}

} // namespace _detail

//
// Host Interface Function
//

template <typename value_type>
void potrf(concurrency::accelerator_view& av, char uplo, int n, value_type* a, int lda)
{
    // quick return
    if (n == 0)
        return;
    
    // error checking
    uplo = static_cast<char>(toupper(uplo));

    if (uplo != 'L' && uplo != 'U')
        argument_error(2);
    if (n < 0)
        argument_error(3);
    if (a == nullptr)
        argument_error(4);
    if (lda < n)
        argument_error(5);

    // host views
    concurrency::array_view<value_type,2> host_view_a(n, lda, a);
    concurrency::array_view<value_type,2> host_view_a_sub = host_view_a.section(concurrency::index<2>(0,0), concurrency::extent<2>(n,n));

    // accelerator array (allocation and copy)
    concurrency::array<value_type,2> accl_a(host_view_a_sub);

    // accelerator view
    concurrency::array_view<value_type,2> accl_view_a(accl_a);

    // forwarding function
    _detail::potrf<option::ordering::column_major>(av, to_option(uplo), accl_view_a);

    // copy back to host
    concurrency::copy(accl_view_a, host_view_a_sub);
}

} // namespace amplapack
