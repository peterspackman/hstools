#include <Python.h>
#ifndef __APPLE__
    #pragma GCC diagnostic ignored "-Wcpp"
#endif
#define PY_ARRAY_UNIQUE_SYMBOL cio_ARRAY_API
#include <numpy/arrayobject.h>
#ifndef __APPLE__
  #pragma GCC diagnostic pop
#endif
#define MAIN_FILE
#include "cio.h"
static char module_docstring[] =
    "This module provides an interface for reading the data values from a .cxs file";
static char readcxsfile_docstring[] =
    "read data values from a given .cxs file";
static char cross3D_docstring[] =
    "compute the cross product of 2 3D vectors";
static char normal3D_docstring[] =
    "compute the normal of a 3D vector";


static PyObject * cio_readcxsfile(PyObject *self, PyObject * args);
static PyObject * cio_cross3D(PyObject * self, PyObject * args);
static PyObject * cio_normal3D(PyObject * sefl, PyObject * args);
static PyObject * CXSError;

static PyMethodDef module_methods[] = {
    {"readcxsfile", cio_readcxsfile, METH_VARARGS, readcxsfile_docstring},
    {"cross3D", cio_cross3D, METH_VARARGS, cross3D_docstring},
    {"normal3D", cio_normal3D, METH_VARARGS, normal3D_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initcio(void)
{
    PyObject *m = Py_InitModule3("cio", module_methods, module_docstring);
    if (m == NULL)
        return;
    CXSError = PyErr_NewException("cio.error", NULL, NULL);
    Py_INCREF(CXSError);
    PyModule_AddObject(m, "error", CXSError);
    import_array();
}

static PyObject * cio_readcxsfile(PyObject *self, PyObject * args)
{
    char * fname;
    if (!PyArg_ParseTuple(args, "s", &fname))
        return NULL;

    PyObject * ret = readcxsfile(fname);
    if(ret == NULL) {
        PyErr_SetString(CXSError,
                        error_string);
        return NULL; //this returns None to Python
    }

    return ret;
}

static PyObject * cio_cross3D(PyObject * self, PyObject * args)
{
    double v1[3];
    double v2[3];
    double r[3];

    if(! PyArg_ParseTuple(args, "dddddd", &v1[0], &v1[1], &v1[2], &v2[0], &v2[1], &v2[2])) return NULL;

    cross3D(v1,v2,r);

    PyObject * ret = Py_BuildValue("ddd", r[0], r[1], r[2]);
    return ret;
}

static PyObject * cio_normal3D(PyObject * self, PyObject * args)
{
    double v[3];

    if(!PyArg_ParseTuple(args, "ddd", &v[0],&v[1],&v[2])) return NULL;
    double r = normal3D(v);

    PyObject * ret = Py_BuildValue("d", r);
    return ret;
}

