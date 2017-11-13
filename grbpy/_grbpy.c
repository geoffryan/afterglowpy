#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#include <numpy/arrayobject.h>
#include "offaxis_struct.h"

static char grbpy_docstring[] = 
    "This module calculates emission from a semi-analytic GRB afterglow model.";
static char fluxDensity_docstring[] = 
    "Calculate the flux density at several times and frequencies";

static PyObject *error_out(PyObject *m);
static PyObject *grbpy_fluxDensity(PyObject *self, PyObject *args);

struct module_state
{
    PyObject *error;
};
#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state *) PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyMethodDef grbpy_methods[] = {
    {"fluxDensity", grbpy_fluxDensity, METH_VARARGS, fluxDensity_docstring},
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL}};

#if PY_MAJOR_VERSION >= 3

static int grbpy_traverse(PyObject *m, visitproc visit, void *arg)
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int grbpy_clear(PyObject *m)
{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef grbpy = {
    PyModuleDef_HEAD_INIT,
    "_grbpy", /* Module Name */
    grbpy_docstring,
    sizeof(struct module_state),
    grbpy_methods,
    NULL,
    grbpy_traverse,
    grbpy_clear,
    NULL
};
#define INITERROR return NULL

PyMODINIT_FUNC PyInit__grbpy(void)
#else
#define INITERROR return

void init_grbpy(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&grbpy);
#else
    PyObject *module = Py_InitModule("_grbpy", grbpy_methods);
#endif
    if(module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);
    st->error = PyErr_NewException("_grbpy.Error", NULL, NULL);
    if(st->error == NULL)
    {
        Py_DECREF(module);
        INITERROR;
    }

    //Load numpy stuff!
    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

static PyObject *error_out(PyObject *m)
{
    struct module_state *st = GETSTATE(m);
    PyErr_SetString(st->error, "something bad happened");
    return NULL;
}

static PyObject *grbpy_fluxDensity(PyObject *self, PyObject *args)
{
    Py_RETURN_NONE;
}
