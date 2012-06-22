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
 * Contains the registry arrays for each module type and the implementation of the calling functions.
 *
 */

#include <string.h>
#include "strutils.h"
#include "modules.h"
#include "errorcodes.h"
#include "errorhandler.h"

/******* THE CODE AFTER THIS LINE CONTAINS THE MODULE REGISTRATION MECHANISM. MODIFY TO ADD YOUR MODULES *******/

// Add here the header files of your modules:
#include "../modules/sampling/dsa.h"
#include "../modules/energy/dopecb.h"
#include "../modules/custom/rbasins.h"
#include "../modules/drivers/oops.h"
#include "../modules/drivers/debug.h"
#include "../modules/custom/carg.h"

// Add to this array of SamplingModuleData structures the functions of your sampling module, plus the name and doc strings of the module.
static struct SamplingModuleData samplingmodules[] = 
{
    {"DSA", (SetSystemModuleFunction)    set_system_for_dsa,
            (InitModuleFunction)         init_dsa, 
            (FinishModuleFunction)       finish_dsa, 
            (PreModuleFunction)          pre_dsa_step,
            (PostModuleFunction)         post_dsa_step, 
            (StepSamplingModuleFunction) run_dsa_step, 
     "Discrete Simulated Annealing"},
    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} // Sentinel array element. Don't remove, it is needed to know where the array ends.
};

// Add to this array of EnergyModuleData structures the functions of your energy module, plus the name and doc strings of the module.
static struct EnergyModuleData energymodules[] = 
{
    {"DOPE-CB", (SetSystemModuleFunction)  set_system_for_dopecb,
                (InitModuleFunction)       init_dopecb, 
                (FinishModuleFunction)     finish_dopecb, 
                (PreModuleFunction)        pre_dopecb_calc,
                (PostModuleFunction)       post_dopecb_calc, 
                (CalcEnergyModuleFunction) calc_dopecb_energy, 
                (GradEnergyModuleFunction) calc_dopecb_gradient, 
     "Discrete Optimized Protein Energy function (main-chain + CB atoms)"},
    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} // Sentinel array element. Don't remove, it is needed to know where the array ends.
};

// Add to this array of CustomModuleData structures the functions of your custom module, plus the name and doc strings of the module.
static struct CustomModuleData custommodules[] = 
{
    {"RB", (SetSystemModuleFunction)  set_system_for_rbasins,
           (InitModuleFunction)       init_rbasins, 
           (FinishModuleFunction)     finish_rbasins, 
           (PreModuleFunction)        pre_rbasins,
           (PostModuleFunction)       post_rbasins,
           (RunCustomModuleFunctionb) calc_rbasinsb,
           (RunCustomModuleFunctioni) calc_rbasinsi,
           (RunCustomModuleFunctionf) calc_rbasinsf, 
     "Ramachandran Basin calculator"},
    {"CARG", (SetSystemModuleFunction)  set_system_for_carg,
             (InitModuleFunction)       init_carg, 
             (FinishModuleFunction)     finish_carg, 
             (PreModuleFunction)        pre_carg,
             (PostModuleFunction)       post_carg, 
             (RunCustomModuleFunctionb) calc_cargi,
             (RunCustomModuleFunctioni) calc_cargi,
             (RunCustomModuleFunctionf) calc_cargf,
     "C-alpha radius of gyration calculator"},
    {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} // Sentinel array element. Don't remove, it is needed to know where the array ends.
};

// Add to this array of DriverModuleData structures the functions of your driver module, plus the name and doc strings of the module.
static struct DriverModuleData drivermodules[] = 
{
    {"OOPS", (SetSystemModuleFunction) set_system_for_oops,
             (SetOpenCLModuleFunction) set_opencl_for_oops,
             (InitModuleFunction)      init_oops, 
             (FinishModuleFunction)    finish_oops, 
             (RunDriverModuleFunction) run_oops, 
     "Default simulation driver"},
    {"DEBUG", (SetSystemModuleFunction)  set_system_for_debug,
              (SetOpenCLModuleFunction)  set_opencl_for_debug,
              (InitModuleFunction)       init_debug,
              (FinishModuleFunction)     finish_debug,
              (RunDriverModuleFunction)  run_debug,
     "Debug driver"},
    {NULL, NULL, NULL, NULL, NULL, NULL, NULL} // Sentinel array element. Don't remove, it is needed to know where the array ends.
};

/********************* DON'T EDIT THE CODE BELOW THIS LINE (UNLESS YOU KNOW WHAT YOU ARE DOING) *********************/

ErrorCode send_system_to_modules(struct AtomData *atoms, FloatValue *coords, 
                                                         FloatValue *distances, FloatValue *invdistances, 
                                                         BoolValue *atommask, BoolValue *distmask, IndexValue natoms, 
                                 struct ResidueData *residues, IndexValue nres, 
                                 struct ChainData *chains, IndexValue nchains)
{
    IndexValue i;
    BoolValue sentinelfound;
    ErrorCode error;

    // Setting system for all sampling modules.
    i = 0;
    sentinelfound = FALSE;
    do
    {
        sentinelfound = samplingmodules[i].setsysfunc == NULL;
        if (!sentinelfound) 
        {
            error = (*(samplingmodules[i].setsysfunc))(atoms, coords, distances, invdistances, atommask, distmask, natoms, 
                                                       residues, nres, chains, nchains);
            if (error != NO_ERROR)
            {
                print_error(error, "send_system_to_modules", 0); 
                return error; 
            }
            i = i + 1;
        }
    } while (!sentinelfound);

    // Setting system for all energy modules.
    i = 0;
    sentinelfound = FALSE;
    do
    {
        sentinelfound = energymodules[i].setsysfunc == NULL;
        if (!sentinelfound) 
        {
            error = (*(energymodules[i].setsysfunc))(atoms, coords, distances, invdistances, atommask, distmask, natoms, 
                                                     residues, nres, chains, nchains);
            if (error != NO_ERROR)
            {
                print_error(error, "send_system_to_modules", 0); 
                return error; 
            }
            i = i + 1;
        }
    } while (!sentinelfound);

    // Setting system for all custom modules.
    i = 0;
    sentinelfound = FALSE;
    do
    {
        sentinelfound = custommodules[i].setsysfunc == NULL;
        if (!sentinelfound) 
        {
            error = (*(custommodules[i].setsysfunc))(atoms, coords, distances, invdistances, atommask, distmask, natoms, 
                                                     residues, nres, chains, nchains);
            if (error != NO_ERROR)
            {
                print_error(error, "send_system_to_modules", 0); 
                return error; 
            }
            i = i + 1;
        }
    } while (!sentinelfound);

    // Setting system for all driver modules.
    i = 0;
    sentinelfound = FALSE;
    do
    {
        sentinelfound = drivermodules[i].setsysfunc == NULL;
        if (!sentinelfound) 
        {
            error = (*(drivermodules[i].setsysfunc))(atoms, coords, distances, invdistances, atommask, distmask, natoms, 
                                                     residues, nres, chains, nchains);
            if (error != NO_ERROR)
            {
                print_error(error, "send_system_to_modules", 0); 
                return error; 
            }
            i = i + 1;
        }
    } while (!sentinelfound);

    return NO_ERROR;
}


ErrorCode send_opencl_to_modules(cl_program *clprogram, cl_kernel *cldistkernel, cl_kernel *cldistsqkernel, cl_kernel *clinvdistkernel, cl_context *clcontext, cl_command_queue *clqueue, cl_mem *clcoords, cl_mem *cldistmask, cl_mem *cldistances)
{
    IndexValue i;
    BoolValue sentinelfound;
    ErrorCode error;

    i = 0;
    sentinelfound = FALSE;
    do
    {
        sentinelfound = drivermodules[i].setoclfunc == NULL;
        if (!sentinelfound) 
        {
            error = (*(drivermodules[i].setoclfunc))(clprogram, cldistkernel, cldistsqkernel, clinvdistkernel, clcontext, clqueue, clcoords, cldistmask, cldistances);
            if (error != NO_ERROR)
            {
                print_error(error, "send_opencl_to_modules", 0); 
                return error; 
            }
            i = i + 1;
        }
    } while (!sentinelfound);
    return NO_ERROR;
}

ErrorCode get_num_sampling_modules(IndexValue *num)
{
    *num = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = samplingmodules[*num].name == NULL &&
                        samplingmodules[*num].setsysfunc == NULL &&  
                        samplingmodules[*num].initfunc == NULL &&
                        samplingmodules[*num].finishfunc == NULL &&
                        samplingmodules[*num].prefunc == NULL &&
                        samplingmodules[*num].postfunc == NULL && 
                        samplingmodules[*num].stepfunc == NULL && 
                        samplingmodules[*num].doc == NULL;
        if (!sentinelfound) *num = *num + 1;

    } while (!sentinelfound);
    return NO_ERROR;
}

ErrorCode get_sampling_module(IndexValue *mod, char *name)
{
    IndexValue i = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = samplingmodules[i].name == NULL &&
                        samplingmodules[i].setsysfunc == NULL && 
                        samplingmodules[i].initfunc == NULL &&
                        samplingmodules[i].finishfunc == NULL &&
                        samplingmodules[i].prefunc == NULL &&
                        samplingmodules[i].postfunc == NULL && 
                        samplingmodules[i].stepfunc == NULL && 
                        samplingmodules[i].doc == NULL;

        if (!sentinelfound) 
        {
            if (strcmp(samplingmodules[i].name, name) == 0)
            {
                *mod = i;
                return NO_ERROR;
            }
            i = i + 1;
        }

    } while (!sentinelfound);
    return ITEM_NOT_FOUND_ERROR;
}

ErrorCode get_sampling_module_name(char **name, IndexValue mod)
{
    return copy_string(name, samplingmodules[mod].name);
}

ErrorCode get_sampling_module_doc(char **doc, IndexValue mod)
{
    return copy_string(doc, samplingmodules[mod].doc);
}

ErrorCode execute_init_func_sampling_module(IndexValue mod, const char *cfgname)
{
    return (*(samplingmodules[mod].initfunc))(cfgname);
}

ErrorCode execute_finish_func_sampling_module(IndexValue mod)
{
    return (*(samplingmodules[mod].finishfunc))();
}

ErrorCode execute_pre_func_sampling_module(IndexValue mod)
{
    return (*(samplingmodules[mod].prefunc))();
}

ErrorCode execute_post_func_sampling_module(IndexValue mod)
{
    return (*(samplingmodules[mod].postfunc))();
}

ErrorCode execute_step_func_sampling_module(IndexValue mod)
{
    return (*(samplingmodules[mod].stepfunc))();
}

ErrorCode get_num_energy_modules(IndexValue *num)
{
    *num = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = energymodules[*num].name == NULL &&
                        energymodules[*num].setsysfunc == NULL && 
                        energymodules[*num].initfunc == NULL &&
                        energymodules[*num].finishfunc == NULL &&
                        energymodules[*num].prefunc == NULL &&
                        energymodules[*num].postfunc == NULL && 
                        energymodules[*num].energyfunc == NULL &&
                        energymodules[*num].gradfunc == NULL && 
                        energymodules[*num].doc == NULL;
        if (!sentinelfound) *num = *num + 1;

    } while (!sentinelfound);
    return NO_ERROR;
}

ErrorCode get_energy_module(IndexValue *mod, char *name)
{
    IndexValue i = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = energymodules[i].name == NULL &&
                        energymodules[i].setsysfunc == NULL && 
                        energymodules[i].initfunc == NULL &&
                        energymodules[i].finishfunc == NULL &&
                        energymodules[i].prefunc == NULL &&
                        energymodules[i].postfunc == NULL && 
                        energymodules[i].energyfunc == NULL &&
                        energymodules[i].gradfunc == NULL && 
                        energymodules[i].doc == NULL;

        if (!sentinelfound) 
        {
            if (strcmp(energymodules[i].name, name) == 0)
            {
                *mod = i;
                return NO_ERROR;
            }
            i = i + 1;
        }

    } while (!sentinelfound);
    return ITEM_NOT_FOUND_ERROR;
}

ErrorCode get_energy_module_name(char **name, IndexValue mod)
{
    return copy_string(name, energymodules[mod].name);
}

ErrorCode get_energy_module_doc(char **doc, IndexValue mod)
{
    return copy_string(doc, energymodules[mod].doc);
}

ErrorCode execute_init_func_energy_module(IndexValue mod, const char *cfgname)
{
    return (*(energymodules[mod].initfunc))(cfgname);
}

ErrorCode execute_finish_func_energy_module(IndexValue mod)
{
    return (*(energymodules[mod].finishfunc))();
}

ErrorCode execute_pre_func_energy_module(IndexValue mod)
{
    return (*(energymodules[mod].prefunc))();
}

ErrorCode execute_post_func_energy_module(IndexValue mod)
{
    return (*(energymodules[mod].postfunc))();
}

ErrorCode execute_calc_energy_func_energy_module(IndexValue mod, FloatValue *energy)
{
    return (*(energymodules[mod].energyfunc))(energy);
}

ErrorCode execute_calc_gradient_func_energy_module(IndexValue mod, FloatValue *grad)
{
    return (*(energymodules[mod].gradfunc))(grad);
}

ErrorCode get_num_custom_modules(IndexValue *num)
{
    *num = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = custommodules[*num].name == NULL &&
                        custommodules[*num].setsysfunc == NULL &&
                        custommodules[*num].initfunc == NULL &&
                        custommodules[*num].finishfunc == NULL &&
                        custommodules[*num].prefunc == NULL &&
                        custommodules[*num].postfunc == NULL && 
                        custommodules[*num].runfuncb == NULL && 
                        custommodules[*num].runfunci == NULL && 
                        custommodules[*num].runfuncf == NULL && 
                        custommodules[*num].doc == NULL;
        if (!sentinelfound) *num = *num + 1;

    } while (!sentinelfound);
    return NO_ERROR;
}

ErrorCode get_custom_module(IndexValue *mod, char *name)
{
    IndexValue i = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = custommodules[i].name == NULL &&
                        custommodules[i].setsysfunc == NULL &&
                        custommodules[i].initfunc == NULL &&
                        custommodules[i].finishfunc == NULL &&
                        custommodules[i].prefunc == NULL &&
                        custommodules[i].postfunc == NULL && 
                        custommodules[i].runfuncb == NULL && 
                        custommodules[i].runfunci == NULL && 
                        custommodules[i].runfuncf == NULL && 
                        custommodules[i].doc == NULL;

        if (!sentinelfound) 
        {
            if (strcmp(custommodules[i].name, name) == 0)
            {
                *mod = i;
                return NO_ERROR;
            }
            i = i + 1;
        }

    } while (!sentinelfound);
    return ITEM_NOT_FOUND_ERROR;
}

ErrorCode get_custom_module_name(char **name, IndexValue mod)
{
    return copy_string(name, custommodules[mod].name);
}

ErrorCode get_custom_module_doc(char **doc, IndexValue mod)
{
    return copy_string(doc, custommodules[mod].doc);
}

ErrorCode execute_init_func_custom_module(IndexValue mod, const char *cfgname)
{
    return (*(custommodules[mod].initfunc))(cfgname);
}

ErrorCode execute_finish_func_custom_module(IndexValue mod)
{
    return (*(custommodules[mod].finishfunc))();
}

ErrorCode execute_pre_func_custom_module(IndexValue mod)
{
    return (*(custommodules[mod].prefunc))();
}

ErrorCode execute_post_func_custom_module(IndexValue mod)
{
    return (*(custommodules[mod].postfunc))();
}

ErrorCode execute_run_funcb_custom_module(IndexValue mod, BoolValue *data, IndexValue len)
{
    return (*(custommodules[mod].runfuncb))(data, len);
}

ErrorCode execute_run_funci_custom_module(IndexValue mod, IntValue *data, IndexValue len)
{
    return (*(custommodules[mod].runfunci))(data, len);
}

ErrorCode execute_run_funcf_custom_module(IndexValue mod, FloatValue *data, IndexValue len)
{
    return (*(custommodules[mod].runfuncf))(data, len);
}

ErrorCode get_num_driver_modules(IndexValue *num)
{
    *num = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = drivermodules[*num].name == NULL &&
                        drivermodules[*num].setsysfunc == NULL &&
                        drivermodules[*num].setoclfunc == NULL &&
                        drivermodules[*num].initfunc == NULL &&
                        drivermodules[*num].finishfunc == NULL &&
                        drivermodules[*num].runfunc == NULL && 
                        drivermodules[*num].doc == NULL;
        if (!sentinelfound) *num = *num + 1;

    } while (!sentinelfound);
    return NO_ERROR;
}

ErrorCode get_driver_module(IndexValue *mod, char *name)
{
    IndexValue i = 0;
    BoolValue sentinelfound = FALSE;
    do
    {
        sentinelfound = drivermodules[i].name == NULL &&
                        drivermodules[i].setsysfunc == NULL &&
                        drivermodules[i].setoclfunc == NULL &&
                        drivermodules[i].initfunc == NULL &&
                        drivermodules[i].finishfunc == NULL &&
                        drivermodules[i].runfunc == NULL && 
                        drivermodules[i].doc == NULL;

        if (!sentinelfound) 
        {
            if (strcmp(drivermodules[i].name, name) == 0)
            {
                *mod = i;
                return NO_ERROR;
            }
            i = i + 1;
        }

    } while (!sentinelfound);
    return ITEM_NOT_FOUND_ERROR;
}

ErrorCode get_driver_module_name(char **name, IndexValue mod)
{
    return copy_string(name, drivermodules[mod].name);
}

ErrorCode get_driver_module_doc(char **doc, IndexValue mod)
{
    return copy_string(doc, drivermodules[mod].doc);
}

ErrorCode execute_init_func_driver_module(IndexValue mod, const char *cfgname)
{
    return (*(drivermodules[mod].initfunc))(cfgname);
}

ErrorCode execute_finish_func_driver_module(IndexValue mod)
{
    return (*(drivermodules[mod].finishfunc))();
}

ErrorCode execute_run_func_driver_module(IndexValue mod)
{
    return (*(drivermodules[mod].runfunc))();
}
