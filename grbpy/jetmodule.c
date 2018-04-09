#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#include <numpy/arrayobject.h>
#include "offaxis_struct.h"

static char jet_docstring[] = 
    "This module calculates emission from a semi-analytic GRB afterglow model.";
static char fluxDensity_docstring[] = 
    "Calculate the flux density at several times and frequencies";
static char emissivity_docstring[] = 
    "Calculate the instantaneous emissivity of a sector of a blastwave.";

static PyObject *error_out(PyObject *m);
static PyObject *jet_fluxDensity(PyObject *self, PyObject *args, 
                                    PyObject *kwargs);
static PyObject *jet_emissivity(PyObject *self, PyObject *args);

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

static PyMethodDef jetMethods[] = {
    {"fluxDensity", (PyCFunction)jet_fluxDensity, METH_VARARGS|METH_KEYWORDS,
        fluxDensity_docstring},
    {"emissivity", jet_emissivity, METH_VARARGS, emissivity_docstring},
    {"error_out", (PyCFunction)error_out, METH_NOARGS, NULL},
    {NULL, NULL, 0, NULL}};

#if PY_MAJOR_VERSION >= 3

static int jet_traverse(PyObject *m, visitproc visit, void *arg)
{
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int jet_clear(PyObject *m)
{
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef jetModule = {
    PyModuleDef_HEAD_INIT,
    "jet", /* Module Name */
    jet_docstring,
    sizeof(struct module_state),
    jetMethods,
    NULL,
    jet_traverse,
    jet_clear,
    NULL
};
#define INITERROR return NULL

PyMODINIT_FUNC PyInit_jet(void)
#else
#define INITERROR return

void initjet(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&jetModule);
#else
    PyObject *module = Py_InitModule3("jet", jetMethods, jet_docstring);
#endif
    if(module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);
    st->error = PyErr_NewException("jet.Error", NULL, NULL);
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

static PyObject *jet_fluxDensity(PyObject *self, PyObject *args, 
                                    PyObject *kwargs)
{
    PyObject *t_obj = NULL;
    PyObject *nu_obj = NULL;
    PyObject *mask_obj = NULL;

    int jet_type, spec_type;
    double theta_obs, E_iso_core, theta_h_core, theta_h_wing, n_0,
           p, epsilon_E, epsilon_B, ksi_N, d_L;

    int latRes = 5;
    double rtol = 1.0e-4;
    int tRes = 1000;
    static char *kwlist[] = {"t", "nu", "jetType", "specType", "thetaObs", 
                                "E0", "thetaCore", "thetaWing", "n0", "p",
                                "epsilon_e", "epsilon_B", "ksiN", "dL", 
                                "tRes", "latRes", "rtol", "mask", NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "OOiidddddddddd|iidO", kwlist,
                &t_obj, &nu_obj, &jet_type, &spec_type, &theta_obs, &E_iso_core,
                &theta_h_core, &theta_h_wing, &n_0, &p, &epsilon_E, &epsilon_B, 
                &ksi_N, &d_L, &tRes, &latRes, &rtol, &mask_obj))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    //Grab NUMPY arrays
    PyArrayObject *t_arr;
    PyArrayObject *nu_arr;
    PyArrayObject *mask_arr = NULL;

    t_arr = (PyArrayObject *) PyArray_FROM_OTF(t_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    nu_arr = (PyArrayObject *) PyArray_FROM_OTF(nu_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    if(mask_obj != NULL)
        mask_arr = (PyArrayObject *) PyArray_FROM_OTF(mask_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);

    if(t_arr == NULL || nu_arr == NULL || (mask_obj != NULL
                                            && mask_arr == NULL))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not read input arrays.");
        Py_XDECREF(t_arr);
        Py_XDECREF(nu_arr);
        Py_XDECREF(mask_arr);
        return NULL;
    }

    int t_ndim = (int) PyArray_NDIM(t_arr);
    int nu_ndim = (int) PyArray_NDIM(nu_arr);
    int mask_ndim = 0;
    if(mask_obj != NULL)
        mask_ndim = (int) PyArray_NDIM(mask_arr);

    if(t_ndim != 1 || nu_ndim != 1 || (mask_obj != NULL && mask_ndim != 1))
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be 1-D.");
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        if(mask_obj != NULL)
            Py_DECREF(mask_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(t_arr, 0);
    int Nnu = (int)PyArray_DIM(nu_arr, 0);
    int Nmask = 0;
    if(mask_obj != NULL)
        Nmask = (int)PyArray_DIM(mask_arr, 0);

    if(N != Nnu)
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be same size.");
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        if(mask_obj != NULL)
            Py_DECREF(mask_arr);
        return NULL;
    }
    if(mask_obj != NULL && Nmask%9 != 0)
    {
        PyErr_SetString(PyExc_RuntimeError, 
                            "Mask length must be multiple of 9.");
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        Py_DECREF(mask_arr);
        return NULL;
    }

    double *t = (double *)PyArray_DATA(t_arr);
    double *nu = (double *)PyArray_DATA(nu_arr);
    double *mask = NULL;
    if(mask_obj != NULL)
        mask = (double *)PyArray_DATA(mask_arr);
    int masklen = Nmask/9;

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
                        epsilon_E, epsilon_B, ksi_N, d_L, tRes, latRes, rtol,
                        mask, masklen);

    // Clean up!
    Py_DECREF(t_arr);
    Py_DECREF(nu_arr);
    if(mask_obj != NULL)
        Py_DECREF(mask_arr);

    //Build output
    PyObject *ret = Py_BuildValue("N", Fnu_obj);
    
    return ret;
}

static PyObject *jet_emissivity(PyObject *self, PyObject *args)
{
    int spec_type = 0;
    double nu, R, sinTheta, mu, te, u, us, n0, p, epse, epsB, ksiN;


    //Parse Arguments
    if(!PyArg_ParseTuple(args, "dddddddddddd|i", &nu, &R, &sinTheta, &mu, &te,
                            &u, &us, &n0, &p, &epse, &epsB, &ksiN, &spec_type))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    // Calculate it!
    double em = emissivity(nu, R, sinTheta, mu, te, u, us, n0, p, epse, epsB, 
                            ksiN, spec_type);

    //Build output
    PyObject *ret = Py_BuildValue("d", em);
    
    return ret;
}
