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

#include "amplapack.h"

namespace amplapack {

int safe_call_interface(std::function<void(const concurrency::accelerator_view& av)>& functor)
{
    // TODO: fully flesh this out
    try
    {
        // TODO: thread safety
        concurrency::accelerator_view av(concurrency::accelerator().default_view);
        functor(av);
    }
    catch(const data_error& e)
    {
        // positive info is a data error code
        return e.get();
    }
    catch(const paramater_error& e)
    {
        // negative info is a parameter error code
        return -static_cast<int>(e.get());
    }
    catch(const runtime_error& e)
    {
        // TODO: what to do with these...
    }
    catch(...)
    {
        printf("uncaught exception!\n");
    }

    // no error
    return 0;
}

template <typename value_type>
class subvector_view
{
public:

    // match some of the interface of array_view
    typedef typename value_type value_type;
    concurrency::extent<1> extent;

	// The stride provided may be negative, in which case elements are retrieved in reverse order.
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
            throw runtime_error();
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

} // namesapce amplapack

#endif AMPLAPACK_RUNTIME_H


