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
enum class Order {row_major, col_major};
enum class Transpose {no_trans, trans, conj_trans};
enum class Uplo {upper, lower};
enum class Diag {non_unit, unit};
enum class Side {left, right};

// exeptions
class paramater_error
{
public:
    paramater_error(unsigned int index)
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

class data_error
{
public:
    data_error(unsigned int index)
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

class runtime_error {};

} // namespace amplapack

#endif // AMPLAPACK_H
