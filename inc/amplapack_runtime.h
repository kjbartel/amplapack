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
 * amplapack_config.h
 *
 *---------------------------------------------------------------------------*/

#pragma once
#ifndef AMPLAPACK_RUNTIME_H
#define AMPLAPACK_RUNTIME_H

#include <functional>
#include <amp.h>

#include "ampblas_runtime.h"

#include "amplapack.h"
#include "ampclapack.h"
#include "ampblas_complex.h"

namespace amplapack {

// casts to internal complex types
ampblas::complex<float>* amplapack_cast(amplapack_fcomplex* ptr) 
{ 
    return reinterpret_cast<ampblas::complex<float>*>(ptr); 
}

ampblas::complex<double>* amplapack_cast(amplapack_dcomplex* ptr)
{ 
    return reinterpret_cast<ampblas::complex<double>*>(ptr); 
}

// exception safe execution wrapper
amplapack_status safe_call_interface(std::function<void(concurrency::accelerator_view& av)>& functor, int& info)
{
    try
    {
        // for now, simply use the default view
        concurrency::accelerator_view av(concurrency::accelerator().default_view);

        // safely call the amplack routine with the specific view
        functor(av);
    }
    catch(const data_error_exception& e)
    {
        // in traditional LAPACK interface a positive info represents a data error
        // this data error differs from routine to routine, but is typically used to indicate a specific point of failure
        info = e.get();
        return amplapack_data_error;
    }
    catch(const argument_error_exception& e)
    {
        // in traditional LAPACK interfaces a negative info represents a parameter error
        info = -static_cast<int>(e.get());
        return amplapack_argument_error;
    }
    catch(const runtime_error_exception&)
    {
        // this typically indicates an algorithmic error
        return amplapack_internal_error;
    }
    catch (const std::bad_alloc&)
    {
        return amplapack_memory_error;
    }
    catch(const concurrency::runtime_exception& e)
    {
        // 
        info = e.get_error_code();
        return amplapack_runtime_error;
    }
    catch(...)
    {
        // this should not be encountered under normal operation
        return amplapack_unknown_error;
    }

    // no error!
    info = 0;
    return amplapack_success;
}

// creates a row or column vector from a 2d array with either the 1st or 2nd dimension being 1 
template <typename value_type>
class subvector_view
{
public:

    // match some of the interface of array_view
    typedef typename value_type value_type;
    concurrency::extent<1> extent;

    subvector_view(const concurrency::array_view<value_type,2>& base_view) restrict(cpu,amp)
    : base_view(base_view)
	{
        //
        if (base_view.extent[1] == 1)
        {
            direction = 0;
        }
        else if (base_view.extent[0] == 1)
        {
            direction = 1;            
        }
        else
        {
            // must be a logical 2d vector 
            runtime_error();
        }

        extent = concurrency::extent<1>(base_view.extent[direction]);
	}

	~subvector_view() restrict(cpu,amp) {}

	value_type& operator[] (const concurrency::index<1>& idx) const restrict(cpu,amp)
	{
        concurrency::index<2> idx_2d = (direction == 0 ? concurrency::index<2>(idx[0],0) : concurrency::index<2>(0,idx[0]));
		return base_view[idx_2d];
	}

private:

    // the dimension of the data:
    //   0 is a row vector
    //   1 is a column vector
    unsigned int direction;
    concurrency::array_view<value_type,2> base_view;
};

template <typename value_type>
subvector_view<value_type> make_subvector_view(const concurrency::array_view<value_type,2>& base_view) restrict(cpu, amp)
{
    return subvector_view<value_type>(base_view);
}

// option used to specify where factorization takes place
enum class block_factor_location { host, accelerator };

// LAPACK character option casting
inline char to_char(enum class option::uplo uplo)
{
    switch (uplo)
    {
    case option::uplo::upper:
        return 'u';
    case option::uplo::lower:
    default:
        return 'l';
    }
}

inline enum class option::uplo to_option(const char& uplo)
{
    switch (uplo)
    {
    case 'L':
    case 'l':
        return option::uplo::lower;
    case 'U':
    case 'u':
    default:
        return option::uplo::upper;
    }
}

// runtime checks
template <typename value_type>
int require_square(const concurrency::array_view<value_type,2>& a)
{
    if (a.extent[0] != a.extent[1])
        runtime_error();

    return a.extent[0];
}

// checks returns from a host LAPACK call
inline void info_check(int info)
{
    if (info > 0)
        data_error(info);
    else if (info < 0)
        argument_error(-info);
}

// data layout aware sub-matrix extraction method
template <enum class option::ordering storage_type, typename value_type>
concurrency::array_view<value_type,2> get_sub_matrix(const concurrency::array_view<value_type,2>& matrix, const concurrency::index<2>& location, const concurrency::extent<2>& size)
{
    if (storage_type == option::ordering::row_major)
    {
        return matrix.section(location, size);
    }
    else
    {
        return matrix.section(concurrency::index<2>(location[1], location[0]), concurrency::extent<2>(size[1], size[0]));
    }
}

} // namesapce amplapack

#endif AMPLAPACK_RUNTIME_H


