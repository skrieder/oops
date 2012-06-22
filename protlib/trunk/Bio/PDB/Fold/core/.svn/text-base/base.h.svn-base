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
 * This file contains basic preprocessor defines used to set the numerical precision and range,
 *.as well as char constants used by Python to identify C types.
 *
 */

#ifndef __BASE_H__
#define __BASE_H__

#include <stddef.h>
#include <stdlib.h>
#include <float.h>

// Comment out the following line to use single precision.
#define __DOUBLE_PRECISION__

// Floating-point type.
#ifdef __DOUBLE_PRECISION__
    #define FloatValue double

    #ifndef DBL_EPSILON
        #define DBL_EPSILON 2.2204460492503131e-16
    #endif

    // Minimum difference between two doubles to be considered different numbers.
    static const double MIN_FLOAT_DIFF = DBL_EPSILON;

    // This char represents a double argument in the CPython functions.
    static const char FloatFormatUnit = 'd';
#else
    #define FloatValue float

    #ifndef FLT_EPSILON
        #define FLT_EPSILON 1.19209290e-7
    #endif

    // Minimum difference between two floats to be considered different numbers.
    static const float MIN_FLOAT_DIFF = FLT_EPSILON;

    // This char represents a float argument in the CPython functions.
    static const char FloatFormatUnit = 'f';
#endif

// Comment out the following line to use int as integer values.
//#define __LONG_INT_VALUES__

// Integer type.
#ifdef __LONG_INT_VALUES__
    #define IntValue long

    // This char represents a long argument in the CPython functions.
    static const char IntFormatUnit = 'l';
#else
    #define IntValue int

    // This char represents an int argument in the CPython functions.
    static const char IntFormatUnit = 'i';
#endif

// Comment out the following line to use unsigned int as index type.
//#define __LONG_INDEXES__

// Index type.
#ifdef __LONG_INDEXES__
    #define IndexValue unsigned long

    // This char represents an unsigend long argument in the CPython functions.
    static const char IndexFormatUnit = 'L';
#else
    #define IndexValue unsigned int

    // This char represents an unsigend int argument in the CPython functions.
    static const char IndexFormatUnit = 'I';
#endif

// Uncomment next line to avoid calculating square root of distance
//#define __USE_SQUARED_DISTANCES__

// Uncomment next line to avoid calculating distances between atoms within the same residue.
//#define __MASK_INTER_RESIDUE_PAIRS__

// Type-code type: just tiny number used to indicate the type of a certain entity.
#define TypeCode char

// This char represents a char argument in the CPython functions.
static const char TypeFormatUnit = 'b';

// Error-code type: 
#define ErrorCode char

// This char represents a char argument in the CPython functions.
static const char ErrorFormatUnit = 'b';

// Boolean type: C99 has a bool type, but just in case we are using a non-C99 compliant compiler.
#define BoolValue char
#define FALSE 0
#define TRUE  1

// This char represents a char argument in the CPython functions.
static const char BoolFormatUnit = 'b';

static const char StringFormatUnit = 's';

#endif