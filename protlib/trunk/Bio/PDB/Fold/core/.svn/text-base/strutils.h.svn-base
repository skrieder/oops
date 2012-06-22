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
 * This file defines some basic C-string operations.
 *
 */

#ifndef __STRUTILS_H__
#define __STRUTILS_H__

#include <stdio.h>
#include "errorcodes.h"

/**
 * Appends string src2 to src1 and puts the result into dest, first allocating the memory needed by dest in the case it doesn't have the correct size.
 * The char pointer *dest should be freed after its use in the invoking code to avoid memory leaks.
 * @param dest char **
 * @param src const char*
 * @return ErrorCode
 */
static inline ErrorCode append_string(char **dest, const char *src1, const char * src2)
{
//    char *ptr = *dest;
    char *ptr = NULL;
    char *tmpstring1 = NULL;
    char *tmpstring2 = NULL;
    int l = strlen(src1)+strlen(src2);
    if (ptr == NULL || strlen(ptr) != l) 
    {
        // Allocating the correct amount of memory space needed by the destination string.
        if (ptr != NULL) free(ptr);
        ptr = malloc((l + 1) * sizeof(char));
        if (ptr == NULL) 
        {
            return MEMORY_ALLOCATION_ERROR;
        }
    }

    tmpstring1=malloc( (strlen(src1)+1)*sizeof(char)); 
    strcpy(tmpstring1,src1);
    tmpstring2=malloc( (strlen(src2)+1)*sizeof(char)); 
    strcpy(tmpstring2,src2);
    sprintf(ptr,"%s%s",tmpstring1,tmpstring2);
    *dest = ptr;

    free(tmpstring1);
    free(tmpstring2);
    return NO_ERROR;
}

/**
 * Copies string src into dest, first allocating the memory needed by dest in the case it doesn't have the correct size.
 * The char pointer *dest should be freed after its use in the invoking code to avoid memory leaks.
 * @param dest char **
 * @param src const char*
 * @return ErrorCode
 */
static inline ErrorCode copy_string(char **dest, const char *src)
{
    char *ptr = *dest;
    int l = strlen(src);
    if (ptr == NULL || strlen(ptr) != l) 
    {
        // Allocating the correct amount of memory space needed by the destination string.
        if (ptr != NULL) free(ptr);
        ptr = malloc((l + 1) * sizeof(char));
        if (ptr == NULL) 
        {
            return MEMORY_ALLOCATION_ERROR;
        }
    }
    strcpy(ptr, src);
    *dest = ptr;
    return NO_ERROR;
}

/**
 * Trims the character c from the right of string str.
 * @param str char *
 * @param c char
 * @return ErrorCode
 */
static inline ErrorCode trim_right(char *str, char c)
{
    char *end;
    int len;

    len = strlen(str);
    while (*str && len)
    {
        end = str + len - 1;
        if (c == *end)
            *end = '\0';
        else
            break;
        len = strlen(str);
    }
    
    return NO_ERROR;
}

/**
 * Trims the character c from the left of string str.
 * @param str char **
 * @param c char
 * @return ErrorCode
 */
static inline ErrorCode trim_left(char *str, char c)
{
    char *end;
    int len,i;

    len = strlen(str);
    while (*str && len)
    {
        end = str + len - 1;
        if (c == *str){
            for(i=0;i<len-1;i++) str[i]=str[i+1];
            *end = '\0';
           }
        else
            break;
        len = strlen(str);
    }
    
    return NO_ERROR;
}

/**
 * Trims not only spaces, but also tabs and line feeds, from the left and right of string str.
 * @param str char **
 * @return ErrorCode
 */
static inline ErrorCode trim_all(char *str)
{
    //NOTE: trim_left must be passed by reference and trim_right must be passed by value
    trim_left(str, 0x09);   // horizontal tab      
    trim_left(str, 0x0b);   // vertical tab
    trim_left(str, 0x0c);   // form feed
    trim_left(str, 0x20);   // space

    trim_right(str, 0x09);   // horizontal tab
    trim_right(str, 0x0a);   // linefeed
    trim_right(str, 0x0b);   // vertical tab
    trim_right(str, 0x0c);   // form feed
    trim_right(str, 0x0d);   // carriage return
    trim_right(str, 0x20);   // space

    return NO_ERROR;
}

#endif
