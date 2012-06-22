/* *
 * The C Protein Folding Library.
 * Copyright (C) 2009 Andres Colubri.
 * Contact: andres.colubri 'AT' gmail.com
 *
 * This library was written at the Institute of Biophysical Dynamics at the University of Chicago.
 * Gordon Center for Integrated Science, W101. 929 East 57th Street, Chicago, IL 60637, USA.
 * Homepage: http://ibd.uchicago.edu/
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

/**
 * This file defines some basic vector operations (addition, subtraction, dot product, etc).
 *
 */

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <math.h>
#include "errorcodes.h"

// 180 / PI
#define RADIANS_TO_DEGREES 57.295779513082320876798154814105

// PI / 180
#define DEGREES_TO_RADIANS 0.017453292519943295769236907684886

/**
 * Adds vectors (vx, vy, vz) and (wx, wy, wz) into (rx, ry, rz).
 * @param rx FloatValue*
 * @param ry FloatValue*
 * @param rz FloatValue*
 * @param vx FloatValue
 * @param vy FloatValue
 * @param vz FloatValue
 * @param wx FloatValue
 * @param wy FloatValue
 * @param wz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode vector_add(FloatValue *rx, FloatValue *ry, FloatValue *rz, FloatValue vx, FloatValue vy, FloatValue vz, FloatValue wx, FloatValue wy, FloatValue wz)
{
    *rx = vx + wx;
    *ry = vy + wy;
    *rz = vz + wz;
    return NO_ERROR; 
}

/**
 * Subtract vector (wx, wy, wz)  from (vx, vy, vz) and puts the result into (rx, ry, rz).
 * @param rx FloatValue*
 * @param ry FloatValue*
 * @param rz FloatValue*
 * @param vx FloatValue
 * @param vy FloatValue
 * @param vz FloatValue
 * @param wx FloatValue
 * @param wy FloatValue
 * @param wz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode vector_subtract(FloatValue *rx, FloatValue *ry, FloatValue *rz, FloatValue vx, FloatValue vy, FloatValue vz, FloatValue wx, FloatValue wy, FloatValue wz)
{
    *rx = vx - wx;
    *ry = vy - wy;
    *rz = vz - wz;
    return NO_ERROR; 
}

/**
 * Calculate the square of the euclidian distance from vector (wx, wy, wz) to (vx, vy, vz) and puts the result into dist.
 * @param sqrdist FloatValue*
 * @param vx FloatValue
 * @param vy FloatValue
 * @param vz FloatValue
 * @param wx FloatValue
 * @param wy FloatValue
 * @param wz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode sqr_vector_distance(FloatValue *sqrdist, FloatValue vx, FloatValue vy, FloatValue vz, FloatValue wx, FloatValue wy, FloatValue wz)
{
    FloatValue dx,dy,dz;
    dx=vx-wx;
    dy=vy-wy;
    dz=vz-wz;
    *sqrdist=dx*dx + dy*dy + dz*dz;
    return NO_ERROR; 
}

/**
 * Stores in r the dot product between vectors (vx, vy, vz) and (wx, wy, wz).
 * @param r FloatValue*
 * @param vx FloatValue
 * @param vy FloatValue
 * @param vz FloatValue
 * @param wx FloatValue
 * @param wy FloatValue
 * @param wz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode vector_dotprod(FloatValue *r, FloatValue vx, FloatValue vy, FloatValue vz, FloatValue wx, FloatValue wy, FloatValue wz)
{
    *r = vx * wx + vy * wy + vz * wz;
    return NO_ERROR; 
}

/**
 * Stores in r the angle (in degrees) between vectors (vx, vy, vz) and (wx, wy, wz).
 * @param r FloatValue*
 * @param vx FloatValue
 * @param vy FloatValue
 * @param vz FloatValue
 * @param wx FloatValue
 * @param wy FloatValue
 * @param wz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode vector_angle(FloatValue *r, FloatValue vx, FloatValue vy, FloatValue vz, FloatValue wx, FloatValue wy, FloatValue wz)
{
    
    FloatValue dot,norm,normsq;

    //normalizing vector (vx,vy,vz)
    normsq=vx*vx+vy*vy+vz*vz;
    norm = sqrt(normsq);
    vx/=norm;
    vy/=norm;
    vz/=norm;

    //normalizing vector (wx,wy,wz)
    normsq=wx*wx+wy*wy+wz*wz;
    norm = sqrt(normsq);
    wx/=norm;
    wy/=norm;
    wz/=norm;
   
    dot = vx * wx + vy * wy + vz * wz; // this is the dot product of the normalized vectors

    *r = RADIANS_TO_DEGREES * acos ( dot );

   return NO_ERROR; 
}

/**
 * Stores in (rx, ry, rz) the cross product between (vx, vy, vz) and (wx, wy, wz).
 * @param rx FloatValue*
 * @param ry FloatValue*
 * @param rz FloatValue*
 * @param vx FloatValue
 * @param vy FloatValue
 * @param vz FloatValue
 * @param wx FloatValue
 * @param wy FloatValue
 * @param wz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode vector_crossprod(FloatValue *rx, FloatValue *ry, FloatValue *rz, FloatValue vx, FloatValue vy, FloatValue vz, FloatValue wx, FloatValue wy, FloatValue wz)
{
    *rx = vy * wz - vz * wy;
    *ry = vz * wx - vx * wz;
    *rz = vx * wy - vy * wx;
    return NO_ERROR;
}

/**
 * Calculates the length of the vector (rx, ry, rz) and stores it in length.
 * @param length FloatValue*
 * @param rx FloatValue
 * @param ry FloatValue
 * @param rz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode vector_length(FloatValue *length, FloatValue rx, FloatValue ry, FloatValue rz)
{
    *length = sqrt(rx * rx + ry * ry + rz * rz);
    return NO_ERROR;
}


/**
 * Normalizes the vector (rx, ry, rz).
 * @param rx FloatValue*
 * @param ry FloatValue*
 * @param rz FloatValue*
 * @return ErrorCode
 */
static inline ErrorCode vector_normalize(FloatValue *rx, FloatValue *ry, FloatValue *rz)
{
    FloatValue normsq = *rx * *rx + *ry * *ry + *rz * *rz;
    if (fabs(normsq) < MIN_FLOAT_DIFF)
    {
        return FLOAT_PRECISION_ERROR;
    }
    FloatValue norm = sqrt(normsq);
    *rx /= norm;
    *ry /= norm;
    *rz /= norm;
    return NO_ERROR;
}

/**
 * Generates the components of the matrix M =((matxx, matxy, matxz), (matyx, matyy, matyz), (matzx, matzy, matzz)) representing a rotation angle around 
 * an axis defined by the vector (axisx, axisy, axisz). Angles must be given in degrees, and the rotation axis must be normalized.
 * @param matxx FloatValue*
 * @param matxy FloatValue*
 * @param matxz FloatValue*
 * @param matyx FloatValue*
 * @param matyy FloatValue*
 * @param matyz FloatValue*
 * @param matzx FloatValue*
 * @param matzy FloatValue*
 * @param matzz FloatValue*
 * @param angle FloatValue
 * @param axisx FloatValue
 * @param axisy FloatValue
 * @param axisz FloatValue
 * @return ErrorCode
 */
static inline ErrorCode create_rotation_matrix(FloatValue *matxx, FloatValue *matxy, FloatValue *matxz, 
                                               FloatValue *matyx, FloatValue *matyy, FloatValue *matyz, 
                                               FloatValue *matzx, FloatValue *matzy, FloatValue *matzz, 
                                               FloatValue angle, FloatValue axisx, FloatValue axisy, FloatValue axisz)
{
    FloatValue cosine = cos(DEGREES_TO_RADIANS* angle);
    FloatValue sine = sin(DEGREES_TO_RADIANS * angle);
    FloatValue oneminuscosine = 1.0 - cosine;

    *matxx = axisx * axisx + (1 - axisx * axisx) * cosine;
    *matxy = axisx * axisy * oneminuscosine + axisz * sine;
    *matxz = axisx * axisz * oneminuscosine - axisy * sine;

    *matyx = axisx * axisy * oneminuscosine - axisz * sine;
    *matyy = axisy * axisy + (1 - axisy * axisy) * cosine;
//    *matyz = axisx * axisz * oneminuscosine + axisx * sine; // very likely error in previous version
    *matyz = axisy * axisz * oneminuscosine + axisx * sine; // very likely error in previous version

    *matzx = axisx * axisz * oneminuscosine + axisy * sine;
    *matzy = axisy * axisz * oneminuscosine - axisx * sine;
    *matzz = axisz * axisz + (1 - axisz * axisz) * cosine;

    return NO_ERROR;
}

#endif
