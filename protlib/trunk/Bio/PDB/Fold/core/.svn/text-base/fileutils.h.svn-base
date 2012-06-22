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

#ifndef __FILEUTILS_H__
#define __FILEUTILS_H__

#include "base.h"

/**
 * Sets the maximum length of the lines that are read in the input files.
 * @param maxlen IndexValue
 * @return ErrorCode
 */
ErrorCode set_max_line_length(IndexValue maxlen);

/**
 * Opens the text file name for reading.
 * @param name const char*
 * @return ErrorCode
 */
ErrorCode open_input_text_file(const char* name);

/**
 * Closes the currently open input file.
 * @return ErrorCode
 */
ErrorCode close_input_file();

/**
 * Returns TRUE if there is an input file open. FALSE otherwise.
 * @return BoolValue
 */
BoolValue is_input_file_open();

/**
 * Tries to read a new line from the currently open input file. Returns TRUE is succesful,
 * FALSE otherwise (which probably just means that the end of the file has been reached).
 * @return BoolValue
 */
BoolValue read_input_line_from_file();

/**
 * Returns in *str the contents of the current input line string. The pointer *str should be freed after its use in the invoking code.
 * to avoid memory leaks in the program.
 * @param str char**
 * @return ErrorCode
 */
ErrorCode get_input_line(char **str);

/**
 * Returns in *val the contents of the idx-th substring in the current input line that is delimited from the next substring by 
 * the string sep.
 * @param val char**
 * @param sep const char*
 * @param idx IndexValue
 * @return ErrorCode
 */
ErrorCode get_input_strval(char **val, const char *sep, IndexValue idx);

/**
 * Returns in *val the integer contents of the idx-th substring in the current input line that is delimited from the next substring by 
 * the string sep.
 * @param val IntValue*
 * @param sep const char*
 * @param idx IndexValue
 * @return ErrorCode
 */
ErrorCode get_input_ival(IntValue *val, const char *sep, IndexValue idx);

/**
 * Returns in *val the floating-point contents of the idx-th substring in the current input line that is delimited from the next substring by 
 * the string sep.
 * @param val FloatValue*
 * @param sep const char*
 * @param idx IndexValue
 * @return ErrorCode
 */
ErrorCode get_input_fval(FloatValue *val, const char *sep, IndexValue idx);

/**
 * Reads the entire file name into a memory buffer.
 * @param name const char*
 * @return ErrorCode
 */
ErrorCode load_text_file_into_buffer(const char* name);

/**
 * Deletes from memory the currently loaded memory buffer.
 * @return ErrorCode
 */
ErrorCode delete_input_buffer();

/**
 * Returns TRUE if there is a loaded memory buffer. FALSE otherwise.
 * @return BoolValue
 */
BoolValue is_input_buffer_loaded();

/**
 * Tries to read a new line from the currently loaded memory buffer. Returns TRUE is succesful,
 * FALSE otherwise (which probably just means that the end of the buffer has been reached).
 * @return BoolValue
 */
BoolValue read_input_line_from_buffer();

#endif