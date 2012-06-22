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
#include "fileio.h"
#include "errorcodes.h"
#include "errorhandler.h"

IndexValue MAX_LINE_LENGHT = 1024;
BoolValue INPUT_FILE_OPEN = 0;

// Local variables used to read an input file.
FILE *infile = 0;
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
    MAX_LINE_LENGHT = maxlen;
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
        instr = malloc(MAX_LINE_LENGHT * sizeof(char));

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

BoolValue read_input_line()
{
    BoolValue res = fgets(instr, MAX_LINE_LENGHT, infile) != NULL;
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
    copy_string(&str, tokenstr);
    trim_all(str);
    *val = (IntValue) strtol(str, (char **)NULL, 10);
    free(str);
    return NO_ERROR;
}

ErrorCode get_input_fval(FloatValue *val, const char *sep, IndexValue idx)
{
    char *str = NULL;
    update_input_token(sep, idx);
    copy_string(&str, tokenstr);
    trim_all(str);
    *val = (FloatValue) strtod(str, (char **)NULL);
    free(str);
    return NO_ERROR;
}
