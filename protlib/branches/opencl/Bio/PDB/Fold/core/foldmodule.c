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
 * This file contains the Python interface functions .
 *
 */

#include <Python.h>
#include "foldmodule.h"
#include "datatypes.h"
#include "memory.h"
#include "datautils.h"
#include "protein.h"
#include "system.h"
#include "modules.h"
#include "simulation.h"

#ifdef __USE_OPENCL_CODEPATH__
#include "clutils.h"
#endif

static PyObject *
pl_allocate_memory(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    char format[4];
    format[0] = IndexFormatUnit;
    format[1] = IndexFormatUnit;
    format[2] = IndexFormatUnit;
    format[3] = '\0';

    if (PyArg_ParseTuple(args, format, &natoms, &nres, &nchains))
    {
        create_atoms_array(&atoms, natoms);
        create_coords_array(&coords, natoms);
        create_distances_array(&distances, natoms);
        create_distances_array(&invdistances, natoms);
        create_bool_array(&atommask, natoms);
        create_distmask_array(&distmask, natoms);
        set_atom_coords(atoms, coords, natoms);

        create_residues_array(&residues, nres);

        create_chains_array(&chains, nchains);

        send_system_to_modules(atoms, coords, distances, invdistances, atommask, distmask, natoms, residues, nres, chains, nchains);

        #ifdef __USE_OPENCL_CODEPATH__
        init_opencl(&clcontext, &cldevices, &clqueue); 
        set_opencl_memory(&clcoords, &cldistmask, &cldistances, clcontext, coords, distmask, distances, natoms);
        compile_opencl_kernels(&clprogram, &cldistkernel, &cldistsqkernel, &clinvdistkernel, clcontext);
        send_opencl_to_modules(&clprogram, &cldistkernel, &cldistsqkernel, &clinvdistkernel, &clcontext, &clqueue, &clcoords, &cldistmask, &cldistances);
        #endif

        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_set_atom(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    TypeCode moltype;
    IndexValue atom;
    FloatValue x;
    FloatValue y;
    FloatValue z;
    IndexValue res;
    IndexValue chain;
    char *name = NULL;

    char format[9];
    format[0] = TypeFormatUnit;
    format[1] = IndexFormatUnit;
    format[2] = FloatFormatUnit;
    format[3] = FloatFormatUnit;
    format[4] = FloatFormatUnit;
    format[5] = IndexFormatUnit;
    format[6] = IndexFormatUnit;
    format[7] = StringFormatUnit;
    format[8] = '\0';

    if (PyArg_ParseTuple(args, format, &moltype, &atom, &x, &y, &z, &res, &chain, &name))
    {
        if (moltype == PROTEIN)
        {
            set_protein_atom(atoms, atom, x, y, z, res, chain, name);
        }
        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_set_residue(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    TypeCode moltype;
    IndexValue res;
    IndexValue firstatom;
    IndexValue lastatom; 
    IndexValue chain; 
    char *name = NULL;

    char format[7];
    format[0] = TypeFormatUnit;
    format[1] = IndexFormatUnit;
    format[2] = IndexFormatUnit;
    format[3] = IndexFormatUnit;
    format[4] = IndexFormatUnit;
    format[5] = StringFormatUnit;
    format[6] = '\0';

    if (PyArg_ParseTuple(args, format, &moltype, &res, &firstatom, &lastatom, &chain, &name))
    {
        if (moltype == PROTEIN)
        {
            set_protein_residue(residues, res, firstatom, lastatom, chain, name);
        }
        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_set_chain(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    TypeCode moltype;
    IndexValue chain;
    IndexValue firstres;
    IndexValue lastres; 

    char format[5];
    format[0] = TypeFormatUnit;
    format[1] = IndexFormatUnit;
    format[2] = IndexFormatUnit;
    format[3] = IndexFormatUnit;
    format[4] = '\0';

    if (PyArg_ParseTuple(args, format, &moltype, &chain, &firstres, &lastres))
    {
        if (moltype == PROTEIN)
        {
            set_protein_chain(chains, chain, firstres, lastres);
        }
        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_build_connectivity(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    TypeCode moltype;
    BoolValue createsc;

    char format[3];
    format[0] = TypeFormatUnit;
    format[1] = BoolFormatUnit;
    format[2] = '\0';

    if (PyArg_ParseTuple(args, format, &moltype, &createsc))
    {
        if (moltype == PROTEIN)
        {
            create_protein_bonds(atoms, natoms, residues, nres, chains, nchains, createsc);
        }
        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_get_num_drivers(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    IndexValue numdrivers;
    get_num_driver_modules(&numdrivers);
    result = Py_BuildValue("i", numdrivers);

    return result;
}

static PyObject *
pl_get_driver_name(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    IndexValue idx;
    char *name = NULL;

    char format[2];
    format[0] = IndexFormatUnit;
    format[1] = '\0';

    if (PyArg_ParseTuple(args, format, &idx))
    {
        get_driver_module_name(&name, idx);
        result = Py_BuildValue("s", name);
        free(name);
    }

    return result;
}

static PyObject *
pl_get_driver_doc(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    IndexValue idx;
    char *name = NULL;

    char format[2];
    format[0] = IndexFormatUnit;
    format[1] = '\0';

    if (PyArg_ParseTuple(args, format, &idx))
    {
        get_driver_module_doc(&name, idx);
        result = Py_BuildValue("s", name);
    }

    return result;
}

static PyObject *
pl_set_driver_by_name(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    char *name = NULL;

    char format[2];
    format[0] = StringFormatUnit;
    format[1] = '\0';

    if (PyArg_ParseTuple(args, format, &name))
    {
        select_driver_by_name(name);
        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_set_driver_by_index(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    IndexValue idx;

    char format[2];
    format[0] = IndexFormatUnit;
    format[1] = '\0';

    if (PyArg_ParseTuple(args, format, &idx))
    {
        select_driver_by_index(idx);
        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_run_simulation(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    char *name = NULL;

    char format[2];
    format[0] = StringFormatUnit;
    format[1] = '\0';

    if (PyArg_ParseTuple(args, format, &name))
    {
        run_simulation(name);
        result = Py_BuildValue("i", 0);
    }

    return result;
}

static PyObject *
pl_get_atom_coords(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    IndexValue atom;

    char format[2];
    format[0] = IndexFormatUnit;
    format[1] = '\0';

    FloatValue x = 0.0;
    FloatValue y = 0.0;
    FloatValue z = 0.0;

    if (PyArg_ParseTuple(args, format, &atom))
    {
        result = PyList_New(0);
        if (0 <= atom && atom < natoms)
        {
            get_atom_coords(&x, &y, &z, atoms, atom);

            PyList_Append(result, PyFloat_FromDouble(x));
            PyList_Append(result, PyFloat_FromDouble(y));
            PyList_Append(result, PyFloat_FromDouble(z));
        }
    }

    return result;
}

static PyObject *
pl_free_memory(PyObject *self, PyObject *args)
{
    PyObject *result = NULL;

    #ifdef __USE_OPENCL_CODEPATH__
    finish_opencl(&cldevices, &clcoords, &cldistmask, &cldistances);
    #endif

    delete_chains_array(&chains, &nchains);

    delete_residues_array(&residues, &nres);

    delete_bool_array(&atommask, natoms);
    delete_distmask_array(&distmask, natoms);
    delete_distances_array(&invdistances, natoms);
    delete_distances_array(&distances, natoms);
    delete_coords_array(&coords, natoms);
    delete_atoms_array(&atoms, &natoms);

    result = Py_BuildValue("i", 0);

    return result;
}

static PyMethodDef methods[] = {
    {"pl_allocate_memory",     (PyCFunction) pl_allocate_memory,     METH_VARARGS, pl_allocate_memory__doc__},
    {"pl_set_atom",            (PyCFunction) pl_set_atom,            METH_VARARGS, pl_set_atom__doc__},
    {"pl_set_residue",         (PyCFunction) pl_set_residue,         METH_VARARGS, pl_set_residue__doc__},
    {"pl_set_chain",           (PyCFunction) pl_set_chain,           METH_VARARGS, pl_set_chain__doc__},
    {"pl_build_connectivity",  (PyCFunction) pl_build_connectivity,  METH_VARARGS, pl_build_connectivity__doc__},
    {"pl_get_num_drivers",     (PyCFunction) pl_get_num_drivers,     METH_VARARGS, pl_get_num_drivers__doc__},
    {"pl_get_driver_name",     (PyCFunction) pl_get_driver_name,     METH_VARARGS, pl_get_driver_name__doc__},
    {"pl_get_driver_doc",      (PyCFunction) pl_get_driver_doc,      METH_VARARGS, pl_get_driver_doc__doc__},
    {"pl_set_driver_by_name",  (PyCFunction) pl_set_driver_by_name,  METH_VARARGS, pl_set_driver_by_name__doc__},
    {"pl_set_driver_by_index", (PyCFunction) pl_set_driver_by_index, METH_VARARGS, pl_set_driver_by_index__doc__},
    {"pl_run_simulation",      (PyCFunction) pl_run_simulation,      METH_VARARGS, pl_run_simulation__doc__},
    {"pl_get_atom_coords",     (PyCFunction) pl_get_atom_coords,     METH_VARARGS, pl_get_atom_coords__doc__},
    {"pl_free_memory",         (PyCFunction) pl_free_memory,         METH_VARARGS, pl_free_memory__doc__},
    {NULL, NULL, 0, NULL}
};

void initfold(void)
{
  PyObject *m;

  m = Py_InitModule4("fold",
                     methods,
                     "C Protein Folding Library",
                     NULL,
                     PYTHON_API_VERSION);
  if (m==NULL) return;

  if (PyErr_Occurred()) Py_FatalError("can't initialize module fold");
}
