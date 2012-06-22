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
 * Sample sampling module.
 *
 */

#include <stdio.h>
#include <string.h>
#include "../../core/memory.h"
#include "../../core/modules.h"
#include "../../core/fileio.h"
#include "../../core/errorhandler.h"
#include "../../core/rotation.h"
#include "../../core/random.h"
#include "dsa.h"

#include "../shared/dsadope.h"   // Include declaring variables and constants shared between DSA and DOPE modules.
extern BoolValue *flags;         // Flags array shared with the DOPE module (has to be declared extern to be shared).
extern FloatValue T;             // Temperature variable shared with the DOPE module (has to be declared extern to be shared).

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

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

IndexValue energymod;
IndexValue rbasinsmod;

ErrorCode set_system_for_dsa(struct AtomData *_atoms, FloatValue *_coords, 
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

ErrorCode init_dsa(const char* cfg)
{
    ErrorCode error = NO_ERROR;
    char *linehead = NULL;
    char *energyname = NULL;
    char *energycfg = NULL;
    char *rbasinsname = NULL;
    char *rbasinscfg = NULL;

    printf("    Initializing DSA sampler...\n");

    if (open_input_text_file(cfg) == NO_ERROR)
    {
        while (read_input_line())
        {
            get_input_strval(&linehead, "=", 0);
            if (strcmp(linehead, "ENERGY FUNCTION") == 0)
            {
                get_input_strval(&energyname, "=", 1);
            }
            else if (strcmp(linehead, "ENERGY CONFIG") == 0)
            {
                get_input_strval(&energycfg, "=", 1);
            }
            else if (strcmp(linehead, "RB CALCULATOR") == 0)
            {
                get_input_strval(&rbasinsname, "=", 1);
            }
            else if (strcmp(linehead, "RB CONFIG") == 0)
            {
                get_input_strval(&rbasinscfg, "=", 1);
            }
            else
            {
                print_error(UNKNOWN_PARAMETER_ERROR, "init_oops", 2, "header ", linehead);
            }
            free(linehead); 
            linehead = NULL; // Make sure of set this string pointer to NULL, beacause it will be used again when reading the next line. 
        }
        close_input_file();    
    }

    error = get_energy_module(&energymod, energyname);
    if (error != NO_ERROR)
    {
        print_error(error, "init_dsa", 1, " when getting energy module");
    }
    free(energyname);

    error = get_custom_module(&rbasinsmod, rbasinsname);
    if (error != NO_ERROR)
    {
        print_error(error, "init_dsa", 1, " when getting rbasins module");
    }
    free(rbasinsname);

    // Initializing energy module.
    execute_init_func_energy_module(energymod, energycfg);
    free(energycfg);

    // Initializing rbasins module.
    execute_init_func_custom_module(rbasinsmod, rbasinscfg);
    free(rbasinscfg);

    // Calculating initial assignement of torsional angles for the whole system.
    update_all_torsional_angles(residues, chains, nchains, coords);

    // Initializing global variables..
    T = 350.0;
    create_bool_array(&flags, nres);

    printf("    Done.\n");

    return error;
}

ErrorCode finish_dsa(void)
{
    delete_bool_array(&flags, nres);

    execute_finish_func_custom_module(rbasinsmod);
    execute_finish_func_energy_module(energymod);
    return NO_ERROR;
}

ErrorCode pre_dsa_step(void)
{
    execute_pre_func_energy_module(energymod);
    execute_pre_func_custom_module(rbasinsmod);
    return NO_ERROR;
}

ErrorCode post_dsa_step(void)
{
    execute_post_func_custom_module(rbasinsmod);
    execute_post_func_energy_module(energymod);
    return NO_ERROR;
}

ErrorCode run_dsa_step(void)
{
    FloatValue energy;
    execute_calc_energy_func_energy_module(energymod, &energy);
    printf("    DOPE energy: %f\n", energy);
    T -= 10.0;

    printf("    Flags array: ");
    IndexValue i;
    for (i = 0; i < min(62, nres); i++)
        printf("%d", flags[i]);
    printf("\n");

    printf("    Just a random number from the Mersenne RNG: %f\n", get_frandom());

    execute_run_funcf_custom_module(rbasinsmod, NULL, 0);
    return NO_ERROR;
}
