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
 * This file implements the funtions to select and run a simulation driver.
 *
 */

#include "errorhandler.h"
#include "simulation.h"

// Code use to indicate if a driver has been selected already.
TypeCode driverstatus = 0;   // 0 means NOT_SET

// Index of the selected driver module.
IndexValue seldriver;

ErrorCode select_driver_by_name(char *name)
{
    IndexValue mod;
    ErrorCode error = get_driver_module(&mod, name);
    if (error == NO_ERROR)
    {
        driverstatus = READY;
        seldriver = mod;
    }
    print_error(error, "select_driver_by_name", 0);
    return error;
}

ErrorCode select_driver_by_index(IndexValue driver)
{
    ErrorCode error = NO_ERROR;
    IndexValue num;
    get_num_driver_modules(&num);
    if (0 <= driver && driver < num)
    {
        driverstatus = READY;
        seldriver = driver;
    }
    else error = INDEX_OUT_OF_BOUNDS_ERROR;
    print_error(error, "select_driver_by_index", 0);
    return error;
}

ErrorCode run_simulation(const char* cfgname)
{
    ErrorCode error0, error1, error2;

    if (driverstatus == NOT_SET)
    {
        error0 = DRIVER_NOT_SET_ERROR;
        print_error(error0, "run_simulation", 0);
        return error0;
    }

    error0 = execute_init_func_driver_module(seldriver, cfgname);
    if (error0 != NO_ERROR)
    {
        print_error(error0, "execute_init_func_driver_module", 0);
        return error0;
    }

    error1 = execute_run_func_driver_module(seldriver);
    if (error1 != NO_ERROR)
    {
        print_error(error1, "execute_run_func_driver_module", 0);
        return error1;
    }

    error2 = execute_finish_func_driver_module(seldriver);
    if (error2 != NO_ERROR)
    {
        print_error(error2, "execute_finish_func_driver_module", 0);
        return error2;
    }

    return NO_ERROR;
}
