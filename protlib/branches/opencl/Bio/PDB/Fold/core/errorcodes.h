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
 * This file defines the error codes used to identify problems while running the library.
 *
 */

#ifndef __ERRORCODES_H__
#define __ERRORCODES_H__

#include "base.h"

static const ErrorCode NO_ERROR                   = 0;
static const ErrorCode MEMORY_ALLOCATION_ERROR    = 1;
static const ErrorCode ARRAY_INCONSISTENCY_ERROR  = 2; 
static const ErrorCode INDEX_OUT_OF_BOUNDS_ERROR  = 3;
static const ErrorCode ITEM_NOT_FOUND_ERROR       = 4;
static const ErrorCode BUILDING_BOND_ERROR        = 5;
static const ErrorCode ACCESSING_UNSET_BOND_ERROR = 6;
static const ErrorCode FLOAT_PRECISION_ERROR      = 7;
static const ErrorCode DRIVER_NOT_SET_ERROR       = 8;
static const ErrorCode SETTING_SYSTEM_ERROR       = 9;
static const ErrorCode UNKNOWN_PARAMETER_ERROR    = 10;
static const ErrorCode CANNOT_OPEN_FILE_ERROR     = 11;
static const ErrorCode MISMATCH                   = 12;

#endif
