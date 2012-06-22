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
 * This file defines various functions to read and write files. This mechanism allows for a single file to be
 * open at the time.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "strutils.h"
#include "fileutils.h"
#include "errorcodes.h"
#include "errorhandler.h"

IndexValue MAX_LINE_LENGTH = 1024;
BoolValue INPUT_FILE_OPEN = 0;
BoolValue INPUT_BUFFER_LOADED = 0;

#define CR 13            // Decimal code of Carriage Return char
#define LF 10            // Decimal code of Line Feed char
#define EOF_MARKER 26    // Decimal code of DOS end-of-file marker

IndexValue get_max_line_length(){return MAX_LINE_LENGTH;}

// Local variables used to read an input file.
FILE *infile = 0;

// Local variables used to read from memory buffer.
char *inbuffer = 0;
char *buffptr = 0;
IntValue lLineCount;             // Current line number 
IntValue lTotalChars;            // Total characters read 

// Local variables that hold current line.
char *instr = 0;
char *tokenstr = 0;
int tokenidx = -1;

ErrorCode update_input_token(const char *sep, IndexValue idx)
{
    if (tokenidx == -1 || idx != tokenidx + 1)
    {
        // Looks for the first token separated by sep in instr.
        tokenstr = strtok(instr, sep);
    }
    else 
    {
        // Moves to the next token in instr.
        tokenstr = strtok(NULL, sep);
    }
    tokenidx = idx;
    return NO_ERROR;
}

ErrorCode set_max_line_length(IndexValue maxlen)
{
    MAX_LINE_LENGTH = maxlen;
    return NO_ERROR;
}

ErrorCode open_input_text_file(const char* name)
{
    if (!INPUT_FILE_OPEN)
    {
        infile = fopen(name, "r");
        if (infile == NULL)
        {
            print_error(CANNOT_OPEN_FILE_ERROR, "open_input_text_file", 4, "opening ", name, " with ", strerror(errno));
            return CANNOT_OPEN_FILE_ERROR;
        }

        if (instr != NULL)
        {
            free(instr);
        }
        instr = malloc(MAX_LINE_LENGTH * sizeof(char));

        INPUT_FILE_OPEN = TRUE;
    }
    return NO_ERROR;
}

ErrorCode close_input_file()
{
    if (INPUT_FILE_OPEN)
    {
        fclose(infile);
        infile = NULL;
        free(instr);
        instr = NULL;
        INPUT_FILE_OPEN = FALSE;
    }
    return NO_ERROR;
}

BoolValue is_input_file_open()
{
    return INPUT_FILE_OPEN;
}

BoolValue read_input_line_from_file()
{
    if (!INPUT_FILE_OPEN) return FALSE;

    BoolValue res = fgets(instr, MAX_LINE_LENGTH, infile) != NULL;
    tokenidx = -1;
    return res;
}

ErrorCode get_input_line(char **str)
{
    return copy_string(str, instr);
}

ErrorCode get_input_strval(char **val, const char *sep, IndexValue idx)
{
    update_input_token(sep, idx);
    copy_string(val, tokenstr);
    trim_all(*val);
    return NO_ERROR;
}

ErrorCode get_input_ival(IntValue *val, const char *sep, IndexValue idx)
{
    char *str = NULL;
    update_input_token(sep, idx);
    if (tokenstr!=NULL){
        copy_string(&str, tokenstr);
        trim_all(str);
        *val = (IntValue) strtol(str, (char **)NULL, 10);
        free(str);
        return NO_ERROR;
    }
    else{
        return ITEM_NOT_FOUND_ERROR;
    }
}

ErrorCode get_input_fval(FloatValue *val, const char *sep, IndexValue idx)
{
    char *str = NULL;
    update_input_token(sep, idx);
    if (tokenstr!=NULL){
        copy_string(&str, tokenstr);
        trim_all(str);
        *val = (FloatValue) strtod(str, (char **)NULL);
        free(str);
        return NO_ERROR;
    }
    else{
        return ITEM_NOT_FOUND_ERROR;
    }
 
}

ErrorCode load_text_file_into_buffer(const char* name)
{
    if (!INPUT_BUFFER_LOADED)
    {
        // Based on code from http://www.mrx.net/c/source.html
        FILE *fhandle = fopen(name,"r");

        if (fhandle == NULL)
        {
            print_error(CANNOT_OPEN_FILE_ERROR, "load_text_file_into_buffer", 4, "opening ", name, " with ", strerror(errno));
            return CANNOT_OPEN_FILE_ERROR;
        }

        IntValue  flen;                // Length of file
        fseek(fhandle, 0L, SEEK_END);  // Position to end of file 
        flen = ftell(fhandle);         // Get file length
        rewind(fhandle);               // Back to start of file

        if (inbuffer == NULL || strlen(inbuffer) != flen) 
        {
            // Allocating the correct amount of memory space needed by the destination string.
            if (inbuffer != NULL) free(inbuffer);
            inbuffer = malloc((flen + 1) * sizeof(char));
            if (inbuffer == NULL) 
            {
                return MEMORY_ALLOCATION_ERROR;
            }
        }

        fread(inbuffer, flen, 1, fhandle); // Read the entire file into inbuffer.
        inbuffer[flen] = '\0';
        fclose(fhandle);

        if (instr != NULL)
        {
            free(instr);
        }
        instr = malloc(MAX_LINE_LENGTH * sizeof(char));

        buffptr = inbuffer;
        lLineCount  = 0L;
        lTotalChars = 0L;
        INPUT_BUFFER_LOADED = TRUE;
    }
    return NO_ERROR;
}

ErrorCode delete_input_buffer()
{
    if (INPUT_BUFFER_LOADED)
    {
        free(inbuffer);
        inbuffer = NULL;
        free(instr);
        instr = NULL;
        INPUT_BUFFER_LOADED = FALSE;
    }
    return NO_ERROR;
}

BoolValue is_input_buffer_loaded()
{
    return INPUT_BUFFER_LOADED;
}

BoolValue read_input_line_from_buffer()
{
    if (!INPUT_BUFFER_LOADED) return FALSE;

    BoolValue morelines = (*buffptr != 0);
    tokenidx = -1;
    if (morelines)
    { 
        IntValue  isNewline;              // Boolean indicating we've read a CR or LF    
        IntValue  lIndex;                 // Index into instr array
        IntValue  lLineLen;               // Current line length 
        IntValue  lStartPos;              // Offset of start of current line 

        lIndex    = 0L;                 // Reset counters and flags
        isNewline = 0;
        lStartPos = lTotalChars;
        while (*buffptr)  // Read until reaching null char
        {
            if (!isNewline)               // Haven't read a CR or LF yet 
            {
                if (*buffptr == CR || *buffptr == LF) // This char IS a CR or LF
                    isNewline = 1;                        // Set flag
            }
            else if (*buffptr != CR && *buffptr != LF) // Already found CR or LF
                break;                                     // Done with line

            instr[lIndex++] = *buffptr++; // Add char to output and increment
            ++lTotalChars;

        } // end while (*buffptr) 

        instr[lIndex] = '\0';     // Terminate the string 
        ++lLineCount;                 // Increment the line counter 
        lLineLen = strlen(instr); // Get length of line
    }
    return morelines;
}
