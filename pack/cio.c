#include <Python.h>
#include <numpy/arrayobject.h>
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

PyMODINIT_FUNC initcio(void) {
  PyObject *m = Py_InitModule3("cio", module_methods, module_docstring);
  if (m == NULL)
    return;
  CXSError = PyErr_NewException("cio.error", NULL, NULL);
  Py_INCREF(CXSError);
  PyModule_AddObject(m, "error", CXSError);
  import_array();
}

static PyObject * cio_readcxsfile(PyObject *self, PyObject * args) {
  char * fname;
  if (!PyArg_ParseTuple(args, "s", &fname))
    return NULL;
  
  CXS_DATA * cxs = readcxsfile(fname);
  if(cxs == NULL) {
    PyErr_SetString(CXSError,
        error_string);
    return NULL; //this returns None to Python
  }
  //UNPACK THE VALUES FROM CXS_DATA
  int nfaces = cxs->nfaces;
  int nvertices = cxs->nvertices;
  //this is ugly but string manipulation in python is MUCH easier
  PyObject * formula =  PyString_FromString(cxs->formula);
  //dimensions etc
  npy_intp didims[1] = {nvertices};
  npy_intp vdims[2] = {nvertices, 3};
  npy_intp idims[2] = {nfaces, 3};
  npy_intp exdims[1] = {nfaces};
  npy_intp stride[1] = {2*sizeof(char)};
  //Only way to create an array of c strings with length 2 like we have
  PyArray_Descr * desc = PyArray_DescrNewFromType(NPY_STRING);
  desc->elsize = 2;
  const int FLAGS = NPY_CARRAY | NPY_OWNDATA;
  //Construct our numpy arrays
  PyObject * divals = PyArray_SimpleNewFromData(1, didims, NPY_FLOAT, cxs->divals);
  PyObject * devals = PyArray_SimpleNewFromData(1, didims, NPY_FLOAT, cxs->devals);
  PyObject * vertices = PyArray_SimpleNewFromData(2, vdims, NPY_FLOAT, cxs->vertices);
  PyObject * indices = PyArray_SimpleNewFromData(2, idims, NPY_INT, cxs->indices);
  PyObject * internal = PyArray_NewFromDescr(&PyArray_Type, desc,
                        1, exdims, stride, cxs->internal, FLAGS, NULL);
  PyObject * external = PyArray_NewFromDescr(&PyArray_Type, desc,
                        1, exdims, stride, cxs->external, FLAGS, NULL);
  //UPDATE THE FLAGS ON SIMPLE ARRAYS
  PyArray_UpdateFlags((PyArrayObject * ) vertices, FLAGS);
  PyArray_UpdateFlags((PyArrayObject * ) indices, FLAGS);
  PyArray_UpdateFlags((PyArrayObject * ) divals, FLAGS);
  PyArray_UpdateFlags((PyArrayObject * ) devals, FLAGS);

  free(cxs);
  PyObject * x = Py_BuildValue("OOOOO",formula, vertices, indices, internal, external);
  PyObject * ret = Py_BuildValue("OOO", divals, devals, x);
  return ret;
}

static PyObject * cio_cross3D(PyObject * self, PyObject * args) {
  double v1[3];
  double v2[3];
  double r[3];

  if(! PyArg_ParseTuple(args, "dddddd", &v1[0], &v1[1], &v1[2], &v2[0], &v2[1], &v2[2])) return NULL;

  cross3D(v1,v2,r);

  PyObject * ret = Py_BuildValue("ddd", r[0], r[1], r[2]);
  return ret;
}

static PyObject * cio_normal3D(PyObject * self, PyObject * args) {
  double v[3];

  if(!PyArg_ParseTuple(args, "ddd", &v[0],&v[1],&v[2])) return NULL;
  double r = normal3D(v);

  PyObject * ret = Py_BuildValue("d", r);
  return ret;
}

