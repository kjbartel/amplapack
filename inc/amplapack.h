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
 * amplapack.h
 *
 *---------------------------------------------------------------------------*/

#pragma once
#ifndef AMPLAPACK_H
#define AMPLAPACK_H

namespace amplapack {

// internal options
namespace option {
enum class ordering {row_major, column_major};
enum class transpose {no_trans, trans, conj_trans};
enum class uplo {upper, lower};
enum class diag {non_unit, unit};
enum class side {left, right};
} // namespace option

// exeptions
class argument_error_exception
{
public:
    argument_error_exception(unsigned int index)
        : index(index) 
    {}

    // return index of invalid parameter
    unsigned int get() const
    {
        return index;
    }

private:
    unsigned int index;
};

void argument_error(unsigned int index)
{
    throw argument_error_exception(index);
}

class data_error_exception
{
public:
    data_error_exception(unsigned int index)
        : index(index) 
    {}

    // return location of data error
    unsigned int get() const
    {
        return index;
    }

private:
    unsigned int index;
};

void data_error(unsigned int index)
{
    throw data_error_exception(index);
}

class runtime_error_exception
{
    // TODO: more information?
};

void runtime_error()
{
    throw runtime_error_exception();
}

} // namespace amplapack

#endif // AMPLAPACK_H
