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

namespace amplapack {
namespace _detail {

template <typename value_type>
void require_square(const concurrency::array_view<value_type,2>& a)
{
    if (a.extent[0] != a.extent[1])
        throw runtime_error();
}

//
// Unblocked Cholesky Factorization
//

template <typename value_type>
void potf2(const concurrency::accelerator_view& av, enum class Uplo uplo, const concurrency::array_view<value_type,2>& a)
{
    require_square(a);
    const int n = a.extent[0];

    if (uplo == Uplo::upper)
    {
        for (int j = 0; j < n; j++)
        {
        }
    }
    else if (uplo == Uplo::lower)
    {
        for (int j = 0; j < n; j++)
        {
            // compute l(j,j)
            value_type a_jj = a[concurrency::index<2>(j,j)];
            value_type dot = value_type();

            if (j != 0)
            {
                concurrency::array_view<const value_type,2> a_vector_2d = a.section(concurrency::index<2>(0,j), concurrency::extent<2>(j,1));
                auto a_vector_1d = make_subvector_view(a_vector_2d);
                dot = ampblas::dot<value_type>(av, a_vector_1d, a_vector_1d);
            }

            a_jj -= dot;

            // check for non-positive-definiteness
            if (a_jj <= value_type())
                throw data_error(j);

            a_jj = sqrt(a_jj);
            a[concurrency::index<2>(j,j)] = a_jj;

            // compute elements j+1:n of column j
            if (j < n)
            {
                int m_ = n-j-1;
                int n_ = j;

                if (m_ > 0)
                {
                    concurrency::array_view<value_type,2> y_vector_2d = a.section(concurrency::index<2>(j,j+1), concurrency::extent<2>(1,m_));
                    auto y_vector_1d = make_subvector_view(y_vector_2d);
                    
                    if (n_ > 0)
                    {
                        concurrency::array_view<const value_type,2> a_sub = a.section(concurrency::index<2>(0,j+1), concurrency::extent<2>(n_,m_));
                        concurrency::array_view<const value_type,2> x_vector_2d = a.section(concurrency::index<2>(0,j), concurrency::extent<2>(n_,1));
                        auto x_vector_1d = make_subvector_view(x_vector_2d);
                        
                        ampblas::gemv(av, AmpblasNoTrans, value_type(-1), a_sub, x_vector_1d, value_type(1), y_vector_1d);
                    }

                    ampblas::scal(av, value_type(1)/a_jj, y_vector_1d);
                }
            }
        }
    }
}

template <int block_size, int look_ahead_depth, typename value_type>
void potrf(const concurrency::accelerator_view& av, enum class Uplo uplo, const concurrency::array_view<value_type,2>& a)
{
    // matrix size
    const int n = a.extent[0];

    // left looking cholesky factorization
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
                concurrency::array_view<const value_type,2> a_sub = a.section(concurrency::index<2>(0,j), concurrency::extent<2>(k_,n_));
                concurrency::array_view<value_type,2> c_sub = a.section(concurrency::index<2>(j,j), concurrency::extent<2>(n_,n_));
                ampblas::syrk(av, AmpblasLower, AmpblasNoTrans, value_type(-1), a_sub, value_type(1), c_sub);
            }
        }

        // factorize current block
        try
        {
            int n_ = jb;
            concurrency::array_view<value_type,2> a_sub = a.section(concurrency::index<2>(j,j), concurrency::extent<2>(n_,n_));
            _detail::potf2(av, uplo, a_sub);
        }
        catch(const data_error& e)
        {
            // offset block error
            throw data_error(e.get() + jb);
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
                    concurrency::array_view<const value_type,2> a_sub = a.section(concurrency::index<2>(0,j+jb), concurrency::extent<2>(k_,m_));
                    concurrency::array_view<const value_type,2> b_sub = a.section(concurrency::index<2>(0,j), concurrency::extent<2>(k_,n_));
                    concurrency::array_view<value_type,2> c_sub = a.section(concurrency::index<2>(j,j+jb), concurrency::extent<2>(n_,m_));
                    ampblas::gemm(av, AmpblasNoTrans, AmpblasTrans, value_type(-1), a_sub, b_sub, value_type(1), c_sub);
                }
            }

            // 
            {
                int m_ = n-j-jb;
                int n_ = jb;
                if (m_ > 0 && n_ > 0)
                {
                    concurrency::array_view<const value_type,2> a_sub = a.section(concurrency::index<2>(j,j), concurrency::extent<2>(n_,n_));
                    concurrency::array_view<value_type,2> c_sub = a.section(concurrency::index<2>(j,j+jb), concurrency::extent<2>(n_,m_));
                    ampblas::trsm(av, AmpblasRight, AmpblasLower, AmpblasTrans, AmpblasNonUnit, value_type(1), a_sub, c_sub);
                }
            }
        }
    }
}

template <typename value_type>
void potrf(const concurrency::accelerator_view& av, enum class Uplo uplo, const concurrency::array_view<value_type,2>& a)
{
    // TODO: tuning framework
    const int block_size = 64;
    const int look_ahead_depth = 1;

    _detail::potrf<block_size, look_ahead_depth>(av, uplo, a);
}

} // namespace _detail

template <typename value_type>
void potrf(const concurrency::accelerator_view& av, char uplo, int n, value_type* a, int lda)
{
    // TODO: wrap into class
    concurrency::array_view<value_type,2> array_view_a(n, lda, a);
    concurrency::array_view<value_type,2> array_view_a_sub = array_view_a.section(concurrency::extent<2>(n, n));

    // forwarding function
    _detail::potrf(av, Uplo::lower, array_view_a_sub);

    // copy back to host
    array_view_a_sub.synchronize();
}

} // namespace amplapack