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
 * This file implements error handling routines.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "errorhandler.h"

void print_error(ErrorCode error, const char *location, int ncustomstr, ...)
{
    char *customstr;

    if (0 < ncustomstr)
    {
        // Concatenating all the custom strings into customstr.
        int i, n, size;
        va_list args;
        char *str;
        char *tmp;
        char *str0 = ": ";
 
        // First, calculating the total memory needed for customstr.
        size = 0;  
        va_start(args, ncustomstr);
        for(i = 3; i < 3 + ncustomstr; i++) 
        {
           str = (char *)va_arg(args, char *);
           size += strlen(str);
        }
        va_end(args);

        // 2 for the initial ": ", size for the concatenation fo the custom strings, 1 for the final "." 
        // and 1 for final null character.
        customstr = malloc((2 + size + 1 + 1) * sizeof(char));

        strcpy(customstr, str0);
        customstr[2] = '\0';
        // Concatenating the custom strings into customstr.
        va_start(args, ncustomstr);
        n = 0;
        for(i = 3; i < 3 + ncustomstr; i++) 
        {
           str = (char *)va_arg(args, char *);

           // Using temporal variable to store the result of the partial concatenation for iteration i:
           tmp = malloc((strlen(customstr) + strlen(str) + 1) * sizeof(char));
           strcpy(tmp, customstr);
           strcat(tmp, str);
           strcpy(customstr, tmp);
           free(tmp);

           // Marking the current end of string at this iteration.
           n += strlen(str);
           customstr[2 + n] = '\0';
        }
        customstr[2 + size] = '.';
        customstr[2 + size + 1] = '\0';

        va_end(args);
    }
    else
    {
        customstr = malloc((1 + 1) * sizeof(char));
        strcpy(customstr, ".");
    }

    if (error == MEMORY_ALLOCATION_ERROR)    printf("Memory allocation error at %s%s\n", location, customstr);
    if (error == ARRAY_INCONSISTENCY_ERROR)  printf("Array inconsistency error at %s%s\n", location, customstr);
    if (error == INDEX_OUT_OF_BOUNDS_ERROR)  printf("Index out of bounds error at %s%s\n", location, customstr); 
    if (error == ITEM_NOT_FOUND_ERROR)       printf("Item not found error at %s%s\n", location, customstr); 
    if (error == BUILDING_BOND_ERROR)        printf("Building bond error at %s%s\n", location, customstr);
    if (error == ACCESSING_UNSET_BOND_ERROR) printf("Accessing unset bond error at %s%s\n", location, customstr);
    if (error == FLOAT_PRECISION_ERROR)      printf("Float precision error at %s%s\n", location, customstr);
    if (error == DRIVER_NOT_SET_ERROR)       printf("Driver not set error at %s%s\n", location, customstr);
    if (error == SETTING_SYSTEM_ERROR)       printf("Setting system error at %s%s\n", location, customstr);
    if (error == UNKNOWN_PARAMETER_ERROR)    printf("Unknown parameter error at %s%s\n", location, customstr);
    if (error == CANNOT_OPEN_FILE_ERROR)     printf("Cannot open file error at %s%s\n", location, customstr);
    if (error == MISMATCH)     printf("Mismatch at %s%s\n", location, customstr);


    free(customstr);
}
