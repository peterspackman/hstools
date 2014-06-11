#include <Python.h>
#include <numpy/arrayobject.h>
#include "readcxsfile.h"
static char module_docstring[] =
  "This module provides an interface for reading the data values from a .cxs file";
static char readcxsfile_docstring[] =
  "read data values from a given .cxs file";


static PyObject * cio_readcxsfile(PyObject *self, PyObject * args);

static PyMethodDef module_methods[] = {
  {"readcxsfile", cio_readcxsfile, METH_VARARGS, readcxsfile_docstring},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initcio(void) {
  PyObject *m = Py_InitModule3("cio", module_methods, module_docstring);
  if (m == NULL)
    return;
  import_array();
}

static PyObject * cio_readcxsfile(PyObject *self, PyObject * args) {
  char * fname;
  if (!PyArg_ParseTuple(args, "s", &fname))
    return NULL;

  CXS_DATA * cxs = readcxsfile(fname);
  if(cxs == NULL) {
    PyErr_SetString(PyExc_RuntimeError,
        "Problem reading cxs file.");
    return NULL;
  }
  //UNPACK THE VALUES FROM CXS_DATA
  int nfaces = cxs->nfaces;
  int nvertices = cxs->nvertices;
  PyObject * formula =  PyString_FromString(cxs->formula);
  //dimensions etc
  npy_intp didims[1] = {nvertices};
  npy_intp vdims[2] = {nvertices, 3};
  npy_intp idims[2] = {nfaces, 3};
  npy_intp exdims[1] = {nfaces};
  npy_intp stride[1] = {2*sizeof(char)};
  PyArray_Descr * desc = PyArray_DescrNewFromType(NPY_STRING);
  desc->elsize = 2;
  //Construct our numpy arrays
  PyObject * divals = PyArray_SimpleNewFromData(1, didims, NPY_FLOAT, cxs->divals);
  PyObject * devals = PyArray_SimpleNewFromData(1, didims, NPY_FLOAT, cxs->devals);
  PyObject * vertices = PyArray_SimpleNewFromData(2, vdims, NPY_FLOAT, cxs->vertices);
  PyObject * indices = PyArray_SimpleNewFromData(2, idims, NPY_INT, cxs->indices);
  PyObject * internal = PyArray_NewFromDescr(&PyArray_Type, desc,
                        1, exdims, stride, cxs->internal,NPY_C_CONTIGUOUS, NULL);
  PyObject * external = PyArray_NewFromDescr(&PyArray_Type, desc,
                        1, exdims, stride, cxs->external,NPY_C_CONTIGUOUS, NULL);

  free(cxs);
  PyObject * x = Py_BuildValue("OOOOO",formula, vertices, indices, internal, external);
  PyObject * ret = Py_BuildValue("OOO", divals, devals, x);
  return ret;
}
