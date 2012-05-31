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
 * Overloaded C++ headers
 *
 *---------------------------------------------------------------------------*/

#pragma once
#ifndef AMPXLAPACK_H
#define AMPXLAPACK_H

#include "ampclapack.h"

template <typename value_type>
inline amplapack_status amplapack_potrf(int m, int n, value_type* a, int lda, int* info);

template <>
inline amplapack_status amplapack_potrf(int m, int n, float* a, int lda, int* info) 
{
    return amplapack_spotrf(m, n, a, lda, info);
}

template <>
inline amplapack_status amplapack_potrf(int m, int n, double* a, int lda, int* info) 
{
    return amplapack_dpotrf(m, n, a, lda, info);
}

template <>
inline amplapack_status amplapack_potrf(int m, int n, amplapack_fcomplex* a, int lda, int* info) 
{
    return amplapack_cpotrf(m, n, a, lda, info);
}

template <>
inline amplapack_status amplapack_potrf(int m, int n, amplapack_dcomplex* a, int lda, int* info) 
{
    return amplapack_zpotrf(m, n, a, lda, info);
}

#endif // AMPXLAPACK_H
