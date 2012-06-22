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
 * Default simulation driver.
 *
 */
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "../../core/modules.h"
#include "../../core/errorhandler.h"
#include "../../core/fileio.h"
#include "../../core/random.h"
#include "oops.h"

// Local copy of system variables.
struct AtomData *atoms;          // Array of atom data.
FloatValue *coords;              // Array of atom coordinates.
FloatValue *distances;           // Array of atom pair-wise distances.
FloatValue *invdistances;        // Array of atom pair-wise 1/distances.
BoolValue *atommask;             // Array of atom masking values.
BoolValue *distmask;             // Array of distance masking values.  
IndexValue natoms;               // Total number of atoms.

struct ResidueData *residues;    // Array of residue data.
IndexValue nres;                 // Total number of residues.
 
struct ChainData *chains;        // Array of chain data.
IndexValue nchains;              // Total number of chains.

IndexValue samplingmod;
IntValue nsteps;

#ifdef __USE_OPENCL_CODEPATH__
// Local copy of OpenCL variables.
cl_program *clprogram;
cl_kernel *cldistkernel;
cl_kernel *cldistsqkernel;
cl_kernel *clinvdistkernel;
cl_context *clcontext;
cl_command_queue *clqueue;
cl_mem *clcoords;
cl_mem *cldistmask;
cl_mem *cldistances;
#endif

ErrorCode set_system_for_oops(struct AtomData *_atoms, FloatValue *_coords, 
                                                       FloatValue *_distances, FloatValue *_invdistances, 
                                                       BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                              struct ResidueData *_residues, IndexValue _nres, 
                              struct ChainData *_chains, IndexValue _nchains)
{
    atoms = _atoms;
    coords = _coords;
    distances = _distances;
    invdistances = _invdistances;
    atommask = _atommask;
    distmask = _distmask;
    natoms = _natoms;
    residues = _residues;
    nres = _nres;
    chains = _chains;
    nchains = _nchains;

    return NO_ERROR;
}

#ifdef __USE_OPENCL_CODEPATH__
ErrorCode set_opencl_for_oops(cl_program *_clprogram, cl_kernel *_cldistkernel, cl_kernel *_cldistsqkernel, cl_kernel *_clinvdistkernel, cl_context *_clcontext, cl_command_queue *_clqueue, cl_mem *_clcoords, cl_mem *_cldistmask, cl_mem *_cldistances)
{
    clprogram = _clprogram;
    cldistkernel = _cldistkernel;
    cldistsqkernel = _cldistsqkernel;
    clinvdistkernel = _clinvdistkernel;
    clcontext = _clcontext;
    clqueue = _clqueue;

    clcoords = _clcoords;
    cldistmask = _cldistmask;
    cldistances = _cldistances;

    return NO_ERROR;
}
#endif

ErrorCode init_oops(const char* cfg)
{
    ErrorCode error = NO_ERROR;
    char *linehead = NULL;
    char *samplingname = NULL;
    char *samplingcfg = NULL;

    printf("Initializing OOPS driver...\n");

    if (open_input_text_file(cfg) == NO_ERROR)
    {
        while (read_input_line())
        {
            get_input_strval(&linehead, "=", 0);
            if (strcmp(linehead, "SAMPLING MODULE") == 0)
            {
                get_input_strval(&samplingname, "=", 1);
            }
            else if (strcmp(linehead, "SAMPLING CONFIG") == 0)
            {
                get_input_strval(&samplingcfg, "=", 1);
            }
            else if (strcmp(linehead, "SIMULATION STEPS") == 0)
            {
                get_input_ival(&nsteps, "=", 1);
            }
            else
            {
                print_error(UNKNOWN_PARAMETER_ERROR, "init_oops", 0);
            }
            free(linehead); 
            linehead = NULL; // Make sure of set this string pointer to NULL, beacause it will be used again when reading the next line. 
        }
        close_input_file();    
    }

    error = get_sampling_module(&samplingmod, samplingname);
    if (error != NO_ERROR)
    {
        print_error(error, "init_oops", 0);
    }
    free(samplingname);

    // Initializing sampling module.
    execute_init_func_sampling_module(samplingmod, samplingcfg);
    free(samplingcfg);

    // Initializing the Mersenne random number generator using a random seed.
    srand(time(NULL));
    int seed = rand();
    set_seed(seed);

    printf("Done.\n");

    return error;
}

ErrorCode finish_oops(void)
{
    execute_finish_func_sampling_module(samplingmod);

    return NO_ERROR;
}

ErrorCode run_oops(void)
{
    IndexValue step;
    for (step = 0; step < nsteps; step++)
    {
        // Running simulation step.
        printf("Running simulation step %d\n", step);
        execute_pre_func_sampling_module(samplingmod);
        execute_step_func_sampling_module(samplingmod);
        execute_post_func_sampling_module(samplingmod);
    }
    return NO_ERROR;
}
