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
 * This file contains the definition of the funtion pointer types and the module data structures
 * for the all types of modules supported by the library: sampling, energy, custom and drivers.
 * It also defines the calling function that access the function pointers of each registered module.
 *
 */

#ifndef __MODULES_H__
#define __MODULES_H__

#include "base.h"
#include "datatypes.h"

// Function pointer type for the system setter function in each module (called once when the memory for the system is allocated).
typedef ErrorCode (*SetSystemModuleFunction)(struct AtomData *atoms, FloatValue *coords, 
                                                                     FloatValue *distances, FloatValue *invdistances, 
                                                                     BoolValue *atommask, BoolValue *distmask, IndexValue natom, 
                                             struct ResidueData *residues, IndexValue nres, 
                                             struct ChainData *chains, IndexValue nchains);

// Function pointer type for module initialization functions (called once at the beginning of the simulation). 
// The parameter is the file name of the configuration file for the module.
typedef ErrorCode (*InitModuleFunction)(const char* cfgname);
// Function pointer type for module finalization functions (called once at the end of the simulation).
typedef ErrorCode (*FinishModuleFunction)(void);

// Function pointer type for module pre-calculation functions (called before performing a calculation).
typedef ErrorCode (*PreModuleFunction)(void);
// Function pointer type for module post-calculation functions (called after performing a calculation).
typedef ErrorCode (*PostModuleFunction)(void);

// Function pointer type for sampling step calculation functions (called to execute a single simulation step).
typedef ErrorCode (*StepSamplingModuleFunction)(void);
// Function pointer type for energy calculation functions (called to calculate the current value of an energy function).
typedef ErrorCode (*CalcEnergyModuleFunction)(FloatValue *energy);
// Function pointer type for gradient calculation functions (called to calculate the current value of an energy gradient).
typedef ErrorCode (*GradEnergyModuleFunction)(FloatValue *gradient);
// Function pointer type for custom module calculation functions that operate on bool arrays(called to execute a custom calculation).
typedef ErrorCode (*RunCustomModuleFunctionb)(BoolValue *data, IndexValue len);
// Function pointer type for custom module calculation functions that operate on int arrays(called to execute a custom calculation).
typedef ErrorCode (*RunCustomModuleFunctioni)(IntValue *data, IndexValue len);
// Function pointer type for custom module calculation functions that operate on float arrays(called to execute a custom calculation).
typedef ErrorCode (*RunCustomModuleFunctionf)(FloatValue *data, IndexValue len);

// Function pointer type for run funtions of driver modules (called to execute an entire simulation).
typedef ErrorCode (*RunDriverModuleFunction)(void);

/**
 * Record of sampling module data.
 */
struct SamplingModuleData
{
    /**
          * Module name.
          */
    char* name;

    /**
          * Module system setter function,
          */
    SetSystemModuleFunction setsysfunc;

    /**
          * Module initialization function,
          */
    InitModuleFunction initfunc;
    /**
          * Module finalization function,
          */
    FinishModuleFunction finishfunc;

    /**
          * Module pre-calculation function,
          */
    PreModuleFunction prefunc;
    /**
          * Module post-calculation function,
          */
    PostModuleFunction postfunc;

    /**
          * Module sampling step function,
          */
    StepSamplingModuleFunction stepfunc;

    /**
          * Module documentation.
          */
    char* doc; 
};

/**
 * Record of energy module data.
 */
struct EnergyModuleData
{
    /**
          * Module name.
          */
    char* name;

    /**
          * Module system setter function,
          */
    SetSystemModuleFunction setsysfunc;

    /**
          * Module initialization function,
          */
    InitModuleFunction initfunc;
    /**
          * Module finalization function,
          */
    FinishModuleFunction finishfunc;

    /**
          * Module pre-calculation function,
          */
    PreModuleFunction prefunc;
    /**
          * Module post-calculation function,
          */
    PostModuleFunction postfunc;

    /**
          * Module energy calculation function,
          */
    CalcEnergyModuleFunction energyfunc;
    /**
          * Module gradient calculation function,
          */
    GradEnergyModuleFunction gradfunc;

    /**
          * Module documentation.
          */
    char* doc; 
};

/**
 * Record of custom module data.
 */
struct CustomModuleData
{
    /**
          * Module name.
          */
    char* name;

    /**
          * Module system setter function,
          */
    SetSystemModuleFunction setsysfunc;

    /**
          * Module initialization function,
          */
    InitModuleFunction initfunc;
    /**
          * Module finalization function,
          */
    FinishModuleFunction finishfunc;

    /**
          * Module pre-calculation function,
          */
    PreModuleFunction prefunc;
    /**
          * Module post-calculation function,
          */
    PostModuleFunction postfunc;

    /**
           * Module custom calculation function for int arrays,
           */
    RunCustomModuleFunctionb runfuncb;
    /**
          * Module custom calculation function for int arrays,
          */
    RunCustomModuleFunctioni runfunci;
    /**
          * Module custom calculation function for float arrays,
          */
    RunCustomModuleFunctionf runfuncf;

    /**
          * Module documentation.
          */
    char* doc; 
};

/**
 * Record of driver module data.
 */
struct DriverModuleData
{
    /**
            * Module name.
            */
    char* name;

    /**
            * Module system setter function,
            */
    SetSystemModuleFunction setsysfunc;

    /**
            * Module initialization function,
            */
    InitModuleFunction initfunc;
    /**
            * Module finalization function,
            */
    FinishModuleFunction finishfunc;

    /**
            * Module main run function,
            */
    RunDriverModuleFunction runfunc;

    /**
            * Module documentation.
            */
    char* doc; 
};

/**
 * Sends the system data arrays to all the registered modules.
 * @param atoms struct AtomData*
 * @param coords FloatValue*
 * @param distances FloatValue*
 * @param invdistances FloatValue*
 * @param atommask BoolValue*
 * @param distmask BoolValue*
 * @param natoms IndexValue
 * @param residues struct ResidueData*
 * @param nres IndexValue
 * @param chains struct ChainData*
 * @param nchains IndexValue
 * @return ErrorCode
 */
ErrorCode send_system_to_modules(struct AtomData *atoms, FloatValue *coords, 
                                                         FloatValue *distances, FloatValue *invdistances, 
                                                         BoolValue *atommask, BoolValue *distmask, IndexValue natoms, 
                                 struct ResidueData *residues, IndexValue nres, 
                                 struct ChainData *chains, IndexValue nchains);

/**
 * Returns in *num the number of registered sampling modules.
 * @param num IndexValue*
 * @return ErrorCode
 */
ErrorCode get_num_sampling_modules(IndexValue *num);
/**
 * Returns in *mod the index of the first sampling module in the registry array with the specified name.
 * @param IndexValue mod*
 * @param name char*
 * @return ErrorCode
 */
ErrorCode get_sampling_module(IndexValue *mod, char *name);
/**
 * Returns in *name the name of sampling module mod. The pointer *name should be freed after its use in the invoking code.
 * @param name char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_sampling_module_name(char **name, IndexValue mod);
/**
 * Returns in *doc the documentation string of sampling module mod.
 * @param doc char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_sampling_module_doc(char **doc, IndexValue mod);
/**
 * Executes the registered initialization custom for driver module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_init_func_sampling_module(IndexValue mod, const char *cfgname);
/**
 * Executes the registered finalization custom for driver module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_finish_func_sampling_module(IndexValue mod);
/**
 * Executes the registered pre-calculation function for sampling module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_pre_func_sampling_module(IndexValue mod);
/**
 * Executes the registered post-calculation function for sampling module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_post_func_sampling_module(IndexValue mod);
/**
 * Executes the registered step calculation function for sampling module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_step_func_sampling_module(IndexValue mod);

/**
 * Returns in *num the number of registered energy modules.
 * @param num IndexValue*
 * @return ErrorCode
 */
ErrorCode get_num_energy_modules(IndexValue *num);
/**
 * Returns in *mod the index of the first energy module in the registry array with the specified name.
 * @param IndexValue mod*
 * @param name char*
 * @return ErrorCode
 */
ErrorCode get_energy_module(IndexValue *mod, char *name);
/**
 * Returns in *name the name of energy module mod. The pointer *name should be freed after its use in the invoking code.
 * @param name char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_energy_module_name(char **name, IndexValue mod);
/**
 * Returns in *doc the documentation string of energy module mod.
 * @param doc char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_energy_module_doc(char **doc, IndexValue mod);
/**
 * Executes the registered initialization function for custom module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_init_func_energy_module(IndexValue mod, const char *cfgname);
/**
 * Executes the registered finalization custom for driver module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_finish_func_energy_module(IndexValue mod);
/**
 * Executes the registered pre-calculation function for energy module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_pre_func_energy_module(IndexValue mod);
/**
 * Executes the registered post-calculation function for energy module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_post_func_energy_module(IndexValue mod);
/**
 * Executes the registered energy calculation function for energy module mod.
 * @param mod IndexValue
 * @param energy FloatValue*
 * @return ErrorCode
 */
ErrorCode execute_calc_energy_func_energy_module(IndexValue mod, FloatValue *energy);
/**
 * Executes the registered gradient calculation function for energy module mod.
 * @param mod IndexValue
 * @param grad FloatValue*
 * @return ErrorCode
 */
ErrorCode execute_calc_gradient_func_energy_module(IndexValue mod, FloatValue *grad);

/**
 * Returns in *num the number of registered custom modules.
 * @param num IndexValue*
 * @return ErrorCode
 */
ErrorCode get_num_custom_modules(IndexValue *num);
/**
 * Returns in *mod the index of the first custom module in the registry array with the specified name.
 * @param IndexValue mod*
 * @param name char*
 * @return ErrorCode
 */
ErrorCode get_custom_module(IndexValue *mod, char *name);
/**
 * Returns in *name the name of custom module mod. The pointer *name should be freed after its use in the invoking code.
 * @param name char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_custom_module_name(char **name, IndexValue mod);
/**
 * Returns in *doc the documentation string of custom module mod.
 * @param doc char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_custom_module_doc(char **doc, IndexValue mod);
/**
 * Executes the registered initialization function for custom module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_init_func_custom_module(IndexValue mod, const char *cfgname);
/**
 * Executes the registered finalization custom for driver module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_finish_func_custom_module(IndexValue mod);
/**
 * Executes the registered pre-calculation function for custom module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_pre_func_custom_module(IndexValue mod);
/**
 * Executes the registered post-calculation function for custom module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_post_func_custom_module(IndexValue mod);
/**
 * Executes the registered custom calculation function for bool arrays for custom module mod.
 * @param mod IndexValue
 * @param data BoolValue*
 * @param IndexValue len
 * @return ErrorCode
 */
ErrorCode execute_run_funcb_custom_module(IndexValue mod, BoolValue *data, IndexValue len);
/**
 * Executes the registered custom calculation function for int arrays for custom module mod.
 * @param mod IndexValue
 * @param data IntValue*
 * @param IndexValue len
 * @return ErrorCode
 */
ErrorCode execute_run_funci_custom_module(IndexValue mod, IntValue *data, IndexValue len);
/**
 * Executes the registered custom calculation function for float arrays for custom module mod.
 * @param mod IndexValue
 * @param data FloatValue*
 * @param IndexValue len
 * @return ErrorCode
 */
ErrorCode execute_run_funcf_custom_module(IndexValue mod, FloatValue *data, IndexValue len);

/**
 * Returns in *num the number of registered driver modules.
 * @param num IndexValue*
 * @return ErrorCode
 */
ErrorCode get_num_driver_modules(IndexValue *num);

/**
 * Returns in *mod the index of the first driver module in the registry array with the specified name.
 * @param IndexValue mod*
 * @param name char*
 * @return ErrorCode
 */
ErrorCode get_driver_module(IndexValue *mod, char *name);
/**
 * Returns in *name the name of driver module mod. The pointer *name should be freed after its use in the invoking code.
 * @param name char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_driver_module_name(char **name, IndexValue mod);
/**
 * Returns in *doc the documentation string of driver module mod.
 * @param doc char**
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode get_driver_module_doc(char **doc, IndexValue mod);
/**
 * Executes the registered initialization function for driver module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_init_func_driver_module(IndexValue mod, const char *cfgname);
/**
 * Executes the registered finalization custom for driver module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_finish_func_driver_module(IndexValue mod);
/**
 * Executes the registered main run function for driver module mod.
 * @param mod IndexValue
 * @return ErrorCode
 */
ErrorCode execute_run_func_driver_module(IndexValue mod);

#endif