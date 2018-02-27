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
    PyObject *module = Py_InitModule3("_grbpy", grbpy_methods, grbpy_docstring);
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
    PyObject *t_obj = NULL;
    PyObject *nu_obj = NULL;

    int jet_type, spec_type;
    double theta_obs, E_iso_core, theta_h_core, theta_h_wing, n_0,
           p, epsilon_E, epsilon_B, ksi_N, d_L;

    //Parse Arguments
    if(!PyArg_ParseTuple(args, "OOiidddddddddd", &t_obj, &nu_obj, &jet_type,
                &spec_type, &theta_obs, &E_iso_core, &theta_h_core, &theta_h_wing,
                &n_0, &p, &epsilon_E, &epsilon_B, &ksi_N, &d_L))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    //Grab NUMPY arrays
    PyArrayObject *t_arr;
    PyArrayObject *nu_arr;

    t_arr = (PyArrayObject *) PyArray_FROM_OTF(t_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    nu_arr = (PyArrayObject *) PyArray_FROM_OTF(nu_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);

    if(t_arr == NULL || nu_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not read input arrays.");
        Py_XDECREF(t_arr);
        Py_XDECREF(nu_arr);
        return NULL;
    }

    int t_ndim = (int) PyArray_NDIM(t_arr);
    int nu_ndim = (int) PyArray_NDIM(nu_arr);

    if(t_ndim != 1 || nu_ndim != 1)
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be 1-D.");
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(t_arr, 0);
    int Nnu = (int)PyArray_DIM(nu_arr, 0);

    if(N != Nnu)
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be same size.");
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        return NULL;
    }

    double *t = (double *)PyArray_DATA(t_arr);
    double *nu = (double *)PyArray_DATA(nu_arr);

    //Allocate output array

    npy_intp dims[1] = {N};
    PyObject *Fnu_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    if(Fnu_obj == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not make flux array.");
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        return NULL;
    }
    double *Fnu = PyArray_DATA((PyArrayObject *) Fnu_obj);

    // Calculate the flux!
    calc_flux_density(jet_type, spec_type, t, nu, Fnu, N, theta_obs, 
                        E_iso_core, theta_h_core, theta_h_wing, n_0, p, 
                        epsilon_E, epsilon_B, ksi_N, d_L);

    // Clean up!
    Py_DECREF(t_arr);
    Py_DECREF(nu_arr);

    //Build output
    PyObject *ret = Py_BuildValue("N", Fnu_obj);
    
    return ret;
}
