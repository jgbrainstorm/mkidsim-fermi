/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\
 * Copyright (C) 2012 by Lukas Gamper <gamperl@gmail.com>  *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define PY_ARRAY_UNIQUE_SYMBOL pycosmo_c_ARRAY_API

#include "k_binspec.h"

#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <numpy/arrayobject.h>

// args: numpy.array lambda, numpy.array spectrum, numpy.array newlambda, numpy.array newspectrum
PyObject * k_binspec_py(PyObject *self, PyObject *args) {

	PyObject *py_lambda, *py_spectrum, *py_newlambda, *py_newspectrum;
	npy_intp len, new_len;

    if (!PyArg_ParseTuple(args, "OOOO", &py_lambda, &py_spectrum, &py_newlambda, &py_newspectrum)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
	}

	if (
		   PyArray_DESCR(py_lambda)->type_num != PyArray_FLOAT
		|| PyArray_DESCR(py_spectrum)->type_num != PyArray_FLOAT
		|| PyArray_DESCR(py_newlambda)->type_num != PyArray_FLOAT
		|| PyArray_DESCR(py_newspectrum)->type_num != PyArray_FLOAT
	) {
        PyErr_SetString(PyExc_ValueError, "The arguments needs to be a numpy.array([..], dtype=float32)");
        return NULL;
	}

	if (
		   PyArray_NDIM(py_lambda) != 1
		|| PyArray_NDIM(py_spectrum) != 1
		|| PyArray_NDIM(py_newlambda) != 1
		|| PyArray_NDIM(py_newspectrum) != 1
	) {
        PyErr_SetString(PyExc_ValueError, "The arguments needs to be a one dimensional numpy.array");
        return NULL;
	}

	if (
		   PyArray_DIMS(py_lambda)[0] != PyArray_DIMS(py_spectrum)[0]
		|| PyArray_DIMS(py_newlambda)[0] != PyArray_DIMS(py_newspectrum)[0]
	) {
        PyErr_SetString(PyExc_ValueError, "(lambda and spectrum) and (newlambda and newspectrum) needs to have the same length");
        return NULL;
	}

	len = PyArray_DIMS(py_lambda)[0];
	new_len = PyArray_DIMS(py_lambda)[0];

	float * lambda = (float *)PyArray_DATA(py_lambda);
	float * spectrum = (float *)PyArray_DATA(py_spectrum);
	float * newlambda = (float *)PyArray_DATA(py_newlambda);
	float * newspectrum = (float *)PyArray_DATA(py_newspectrum);

	k_binspec(lambda, spectrum, newlambda, newspectrum, len, new_len);

	Py_INCREF(Py_None);
	return Py_None;
}

PyMethodDef moduleMethods[] = {
	  {"k_binspec", k_binspec_py, METH_VARARGS, "TBD: add documentation."}
	, {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initkcorrect(void) {
	import_array();
	(void)Py_InitModule("kcorrect", moduleMethods);
}