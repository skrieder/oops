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
 * Driver for debugging purposes.
 * Author: Glen Hocky
 *
 */

#include <stdio.h>
#include <string.h>
#include "../../core/modules.h"
#include "../../core/rotation.h"
#include "../../core/distance.h"
#include "../../core/datautils.h"
#include "../../core/errorhandler.h"
#include "../../core/fileio.h"
#include "../../core/debugutils.h"
#include "debug.h"

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

IndexValue samplingmod;
IndexValue energymod;
IndexValue custommod;
IntValue nsteps=1;

BoolValue printdebuginfo=0;

ErrorCode set_system_for_debug(struct AtomData *_atoms, FloatValue *_coords, 
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
ErrorCode set_opencl_for_debug(cl_program *_clprogram, cl_kernel *_cldistkernel, cl_kernel *_cldistsqkernel, cl_kernel *_clinvdistkernel, cl_context *_clcontext, cl_command_queue *_clqueue, cl_mem *_clcoords, cl_mem *_cldistmask, cl_mem *_cldistances)
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

    printf("Got OpenCL from the core\n");

    return NO_ERROR;
}
#endif

ErrorCode init_debug(const char* cfg)
{
    ErrorCode error = NO_ERROR;
    char *linehead = NULL;
    char *energyname = NULL;
    char *energycfg = NULL;
    char *samplingname = NULL;
    char *samplingcfg = NULL;
    char *customname = NULL;
    char *customcfg = NULL;

    if (open_input_text_file(cfg) == NO_ERROR)
    {
        while (read_input_line())
        {
            get_input_strval(&linehead, "=", 0);

            char firstchar=linehead[0];
            if (firstchar=='#') {
		free(linehead);
	        linehead = NULL;
                continue;
	    }

            if (strcmp(linehead, "ENERGY MODULE") == 0)
            {
                get_input_strval(&energyname, "=,", 1);
                get_input_strval(&energycfg, "=,", 2);
                numemodules++;
            }
            else if (strcmp(linehead, "SAMPLING MODULE") == 0)
            {
                get_input_strval(&samplingname, "=,", 1);
                get_input_strval(&samplingcfg, "=,", 2);
                numsampmodules++;
            }
            else if (strcmp(linehead, "CUSTOM MODULE") == 0)
            {
                get_input_strval(&customname, "=,", 1);
                get_input_strval(&customcfg, "=,", 2);
                numcustmodules++;
            }
            else if (strcmp(linehead, "CALCULATION STEPS") == 0)
            {
                get_input_ival(&nsteps, "=", 1);
            }
            else if (strcmp(linehead, "MOLECULE INFO") == 0)
            {
                get_input_ival(&printdebuginfo, "=", 1);
            }
            else
            {
                print_error(UNKNOWN_PARAMETER_ERROR, "init_debug", 0);
                printf("...line: %s\n",linehead);
            }
            free(linehead); 
            linehead = NULL; // Make sure of set this string pointer to NULL, beacause it will be used again when reading the next line. 
        }
        close_input_file();    
    }

    // Calculating initial assignement of torsional angles for the whole system.
    update_all_torsional_angles(residues, chains, nchains, coords);
    set_distmask_intervals_product(1,distmask,atoms,natoms,0,natoms-1,0,natoms-1);

    if (printdebuginfo) print_molecule_data(atoms, natoms, residues, nres, chains, nchains);


    if (numsampmodules>0) {
        error = get_sampling_module(&samplingmod, samplingname);
        if (error != NO_ERROR)
        {
            print_error(error, "init_debug", 0);
        }
        free(samplingname);
        execute_init_func_sampling_module(samplingmod, samplingcfg);
    }
    free(samplingcfg);

    // Initializing energy module.
    if (numemodules>0) {
        error = get_energy_module(&energymod, energyname);
        if (error != NO_ERROR)
        {
            print_error(error, "init_debug", 0);
        }
        free(energyname);
        execute_init_func_energy_module(energymod, energycfg);
    }
    free(energycfg);

    // Initializing custom module.
    if (numcustmodules>0) {
        error = get_custom_module(&custommod, customname);
        if (error != NO_ERROR)
        {
            print_error(error, "init_debug", 0);
        }
        free(customname);
        execute_init_func_custom_module(custommod, customcfg);
    }
    free(customcfg);

    return error;
}

ErrorCode finish_debug(void)
{
    // Finishing sampling module.
    execute_finish_func_sampling_module(samplingmod);

    return NO_ERROR;
}

ErrorCode run_debug(void)
{
    IndexValue step;
    FloatValue energy;
    if (numsampmodules>0) {
        for (step = 0; step < nsteps; step++)
        {
        // Running simulation step.
        printf("Running simulation step %d\n", step); 
        execute_pre_func_sampling_module(samplingmod);
        execute_step_func_sampling_module(samplingmod);
        execute_post_func_sampling_module(samplingmod);
        }
    }
    if (numemodules>0) {
        // Running simulation step.
        printf("Running energy calucation\n"); 
        execute_pre_func_energy_module(energymod);
        for (step = 0; step < nsteps; step++){
            execute_calc_energy_func_energy_module(energymod,&energy);
            printf("        Energy: %f\n",energy);
            execute_post_func_energy_module(energymod);
            }
        }
    if (numcustmodules>0) {
        FloatValue result;
        for (step = 0; step < nsteps; step++)
        {


        #ifdef __USE_OPENCL_CODEPATH__
        printf("OpenCL data: %d %d %d %d %d %d\n", *clqueue, *cldistkernel, *clcoords, *cldistmask, *cldistances, natoms);
        //printf("OpenCL data: %d %d %d %d %d\n", *clqueue, *cldistkernel, *clcoords, *cldistmask, natoms);

        run_opencl_distkernel(*clqueue, *cldistkernel, *clcoords, *cldistmask, *cldistances, natoms, distances);

       int i;
       for (i = 0; i < 10; i++)
           printf("%i %f\n", i, distances[i]);
       printf("\n");

        #else
        calculate_all_distances(distances, coords, natoms);
        #endif

        // Running simulation step.
        printf("Running custom calculation step %d\n", step); 
        execute_pre_func_custom_module(custommod);
        execute_run_funcf_custom_module(custommod, &result, 1);
        //res = (FloatValue*)result;

        printf("        Result: %f\n",result);
        execute_post_func_custom_module(custommod);
        }
   }
   return NO_ERROR;
}
