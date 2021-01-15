#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_11_API_VERSION
#include <numpy/arrayobject.h>
#include <time.h>
#include "offaxis_struct.h"

#define PROFILE
#define PROFILE1
#define PROFILE2
#define PROFILEOUTA

static const int jet_type_default = -1;
static const int spec_type_default = 0;
static const double theta_obs_default = 0.0;
static const double E_iso_core_default = 1.0e53;
static const double theta_h_core_default = 0.1;
static const double theta_h_wing_default = 0.4;
static const double b_default = 0.0;
static const double L0_default = 0.0;
static const double q_default = 0.0;
static const double ts_default = 0.0; 
static const double n_0_default = 1.0;
static const double p_default = 2.2;
static const double epsilon_E_default = 0.1;
static const double epsilon_B_default = 0.01;
static const double ksi_N_default = 1.0; 
static const double d_L_default = 1.0e27;

static const int latRes_default = 5;
static const int tRes_default = 1000;
static const int spread_default = 7;
static const int counterjet_default = 0;
static const int gamma_type_default = GAMMA_INF;
static const double g0_default = -1.0;
static const double E_core_global_default = 0.0;
static const double theta_h_core_global_default = 0.0;

static const double rtol_struct_default = 1.0e-2;
static const double rtol_theta_default = 1.0e-2;
static const double rtol_phi_default = 1.0e-2;
static const int int_type_default = INT_CADRE;
static const int nmax_phi_default = 1000;
static const int nmax_theta_default = 1000;

static char jet_docstring[] = 
    "This module calculates emission from a semi-analytic GRB afterglow model.";
static char fluxDensity_docstring[] = 
    "Calculate the flux density at several times and frequencies";
static char emissivity_docstring[] = 
    "Calculate the instantaneous emissivity of a sector of a blastwave.";
static char intensity_docstring[] = 
    "Calculate the position dependent intensity of a blastwave.";
static char shockVals_docstring[] = 
    "Calculate the shock values of the blastwave.";
static char shock_docstring[] = 
    "Calculate the evolution of a tophat shock.";
static char shockObs_docstring[] = 
    "Calculate the evolution of a tophat shock with reference to observer time.";
static char find_jet_edge_docstring[] = 
    "Find jet edge at given observer time, phi, viewing angle.";

static PyObject *error_out(PyObject *m);
static PyObject *jet_fluxDensity(PyObject *self, PyObject *args, 
                                    PyObject *kwargs);
static PyObject *jet_emissivity(PyObject *self, PyObject *args);
static PyObject *jet_intensity(PyObject *self, PyObject *args, 
                                    PyObject *kwargs);
static PyObject *jet_shockVals(PyObject *self, PyObject *args, 
                                    PyObject *kwargs);
static PyObject *jet_shock(PyObject *self, PyObject *args);
static PyObject *jet_shockObs(PyObject *self, PyObject *args);
static PyObject *jet_find_jet_edge(PyObject *self, PyObject *args);

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
    {"intensity", (PyCFunction)jet_intensity, METH_VARARGS|METH_KEYWORDS,
        intensity_docstring},
    {"shockVals", (PyCFunction)jet_shockVals, METH_VARARGS|METH_KEYWORDS,
        shockVals_docstring},
    {"shock", jet_shock, METH_VARARGS, shock_docstring},
    {"shockObs", jet_shockObs, METH_VARARGS, shockObs_docstring},
    {"find_jet_edge", jet_find_jet_edge, METH_VARARGS, 
        find_jet_edge_docstring},
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

    PyModule_AddIntConstant(module, "Cone", _cone);
    PyModule_AddIntConstant(module, "TopHat", _tophat);
    PyModule_AddIntConstant(module, "Gaussian", _Gaussian);
    PyModule_AddIntConstant(module, "PowerLaw", _powerlaw);
    PyModule_AddIntConstant(module, "PowerLawCore", _powerlaw_core);
    PyModule_AddIntConstant(module, "GaussianCore", _Gaussian_core);
    PyModule_AddIntConstant(module, "Exponential", _exponential);
    PyModule_AddIntConstant(module, "TwoComponent", _twocomponent);
    PyModule_AddIntConstant(module, "Spherical", _spherical);
    PyModule_AddIntConstant(module, "TrapFixed", INT_TRAP_FIXED);
    PyModule_AddIntConstant(module, "TrapAdapt", INT_TRAP_ADAPT);
    PyModule_AddIntConstant(module, "SimpFixed", INT_SIMP_FIXED);
    PyModule_AddIntConstant(module, "SimpAdapt", INT_SIMP_ADAPT);
    PyModule_AddIntConstant(module, "RombAdapt", INT_ROMB_ADAPT);
    PyModule_AddIntConstant(module, "TrapNL", INT_TRAP_NL);
    PyModule_AddIntConstant(module, "Hybrid", INT_HYBRID);
    PyModule_AddIntConstant(module, "Cadre", INT_CADRE);
    PyModule_AddIntConstant(module, "GK49Adapt", INT_GK49_ADAPT);
    PyModule_AddIntConstant(module, "GK715Adapt", INT_GK715_ADAPT);
    PyModule_AddIntConstant(module, "GK1021Adapt", INT_GK1021_ADAPT);
    PyModule_AddIntConstant(module, "GammaInf", GAMMA_INF);
    PyModule_AddIntConstant(module, "GammaFlat", GAMMA_FLAT);
    PyModule_AddIntConstant(module, "GammaEvenMass", GAMMA_EVENMASS);
    PyModule_AddIntConstant(module, "GammaStruct", GAMMA_STRUCT);
    PyModule_AddIntConstant(module, "ICCooling", IC_COOLING_FLAG);
    PyModule_AddIntConstant(module, "EpsEBar", EPS_E_BAR_FLAG);
    PyModule_AddIntConstant(module, "SSASmooth", SSA_SMOOTH_FLAG);
    PyModule_AddIntConstant(module, "SSASharp", SSA_SHARP_FLAG);

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

#ifdef PROFILE
    clock_t profClock1A, profClock1B, profClock2A, profClock2B;
#endif

#ifdef PROFILE1
    //Profile 1
    profClock1A = clock();
#endif

    int jet_type = jet_type_default;
    int spec_type = spec_type_default;
    double theta_obs = theta_obs_default;
    double E_iso_core = E_iso_core_default;
    double theta_h_core = theta_h_core_default;
    double theta_h_wing = theta_h_wing_default;
    double b = b_default;
    double L0 = L0_default;
    double q = q_default;
    double ts = ts_default; 
    double n_0 = n_0_default;
    double p = p_default;
    double epsilon_E = epsilon_E_default;
    double epsilon_B = epsilon_B_default;
    double ksi_N = ksi_N_default; 
    double d_L = d_L_default;

    int latRes = latRes_default;
    int tRes = tRes_default;
    double g0 = g0_default;
    double E_core_global = E_core_global_default;
    double theta_h_core_global = theta_h_core_global_default;

    double rtol_struct = rtol_struct_default;
    double rtol_theta = rtol_theta_default;
    double rtol_phi = rtol_phi_default;
    int int_type = int_type_default;
    int nmax_phi = nmax_phi_default;
    int nmax_theta = nmax_theta_default;

    int spread = spread_default;
    int counterjet = counterjet_default;
    int gamma_type = gamma_type_default;

    static char *kwlist[] = {"t", "nu", "jetType", "specType",
                                "thetaObs", "E0", "thetaCore", "thetaWing",
                                    "b", "L0", "q", "ts", "n0", "p",
                                    "epsilon_e", "epsilon_B", "xi_N", "d_L",
                                    "g0",
                                "E0Global", "thetaCoreGlobal",
                                "tRes", "latRes", "intType", "rtolStruct",
                                    "rtolPhi", "rtolTheta", "NPhi", "NTheta",
                                "mask",
                                "spread", "counterjet", "gammaType",
                                NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs,
                "OO|ii""ddddddddddddddd""dd""iiidddii""O""iii",
                //"OO|ii ddddddddddddddd dd iiidddii O iii",
                kwlist,
                &t_obj, &nu_obj, &jet_type, &spec_type,
                &theta_obs, &E_iso_core, &theta_h_core, &theta_h_wing, &b, &L0,
                    &q, &ts, &n_0, &p, &epsilon_E, &epsilon_B, &ksi_N, &d_L,
                    &g0,
                &E_core_global, &theta_h_core_global,
                &tRes, &latRes, &int_type, &rtol_struct, &rtol_phi,
                    &rtol_theta, &nmax_phi, &nmax_theta,
                &mask_obj,
                &spread, &counterjet, &gamma_type))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    if(int_type < 0 || int_type >= INT_UNDEFINED)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "intType out of range, unknown integrator");
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

#ifdef PROFILE2
    //Profile 2
    profClock2A = clock();
#endif

    //Set Up The Parameters!
    
    double ta = t[0];
    double tb = t[0];
    int i;
    for(i=0; i<N; i++)
    {
        if(t[i] < ta)
            ta = t[i];
        else if(t[i] > tb)
            tb = t[i];
    }
    
    struct fluxParams fp;
    setup_fluxParams(&fp, d_L, theta_obs, E_iso_core, theta_h_core,
                        theta_h_wing, b,
                        L0, q, ts,
                        n_0, p, epsilon_E, epsilon_B, ksi_N, g0, 
                        E_core_global, theta_h_core_global, ta, tb,
                        tRes, latRes, int_type,
                        rtol_struct, rtol_phi, rtol_theta,
                        nmax_phi, nmax_theta,
                        spec_type, mask, masklen,
                        spread, counterjet, gamma_type);

    // Calculate the flux!
    calc_flux_density(jet_type, spec_type, t, nu, Fnu, N, &fp);
   
    if(fp.error)
    {
        PyErr_SetString(PyExc_RuntimeError, fp.error_msg);
        free_fluxParams(&fp);
        return NULL;
    }

    //Free the parameters!
    free_fluxParams(&fp);

#ifdef PROFILE2
    //Profile 2
    profClock2B = clock();
#endif


    // Clean up!
    Py_DECREF(t_arr);
    Py_DECREF(nu_arr);
    if(mask_obj != NULL)
        Py_DECREF(mask_arr);

    //Build output
    PyObject *ret = Py_BuildValue("N", Fnu_obj);
    
#ifdef PROFILE1
    //Profile 1 and output
    profClock1B = clock();
#endif

#ifdef PROFILEOUT
#ifdef PROFILE2
    printf("C Eval Inner: %lf s\n",
            ((double) profClock2B-profClock2A)/CLOCKS_PER_SEC);
#endif
#ifdef PROFILE1
    printf("C Eval Outer: %lf s\n",
            ((double) profClock1B-profClock1A)/CLOCKS_PER_SEC);
#endif
#endif
    
    return ret;
}

static PyObject *jet_emissivity(PyObject *self, PyObject *args)
{
    int spec_type = 0;
    double nu, R, mu, te, u, us, n0, p, epse, epsB, xi_N;


    //Parse Arguments
    if(!PyArg_ParseTuple(args, "ddddddddddd|i", &nu, &R, &mu, &te,
                            &u, &us, &n0, &p, &epse, &epsB, &xi_N, &spec_type))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    // Calculate it!
    double em = emissivity(nu, R, mu, te, u, us, n0, p, epse, epsB, 
                            xi_N, spec_type);

    //Build output
    PyObject *ret = Py_BuildValue("d", em);
    
    return ret;
}

static PyObject *jet_intensity(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *theta_obj = NULL;
    PyObject *phi_obj = NULL;
    PyObject *t_obj = NULL;
    PyObject *nu_obj = NULL;
    PyObject *mask_obj = NULL;

    int jet_type = jet_type_default;
    int spec_type = spec_type_default;
    double theta_obs = theta_obs_default;
    double E_iso_core = E_iso_core_default;
    double theta_h_core = theta_h_core_default;
    double theta_h_wing = theta_h_wing_default;
    double b = b_default;
    double L0 = L0_default;
    double q = q_default;
    double ts = ts_default; 
    double n_0 = n_0_default;
    double p = p_default;
    double epsilon_E = epsilon_E_default;
    double epsilon_B = epsilon_B_default;
    double ksi_N = ksi_N_default; 
    double d_L = d_L_default;

    int latRes = latRes_default;
    int tRes = tRes_default;
    double g0 = g0_default;
    double E_core_global = E_core_global_default;
    double theta_h_core_global = theta_h_core_global_default;

    double rtol_struct = rtol_struct_default;
    double rtol_theta = rtol_theta_default;
    double rtol_phi = rtol_phi_default;
    int int_type = int_type_default;
    int nmax_phi = nmax_phi_default;
    int nmax_theta = nmax_theta_default;

    int spread = spread_default;
    int counterjet = counterjet_default;
    int gamma_type = gamma_type_default;

    static char *kwlist[] = {"theta", "phi", "t", "nu", "jetType", "specType",
                                "thetaObs", "E0", "thetaCore", "thetaWing",
                                    "b", "L0", "q", "ts", "n0", "p",
                                    "epsilon_e", "epsilon_B", "xi_N", "d_L",
                                    "g0",
                                "E0Global", "thetaCoreGlobal",
                                "tRes", "latRes", "intType", "rtolStruct",
                                    "rtolPhi", "rtolTheta", "NPhi", "NTheta",
                                "mask",
                                "spread", "counterjet", "gammaType",
                                NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs,
                "OOOO|ii""ddddddddddddddd""dd""iiidddii""O""iii",
                kwlist,
                &theta_obj, &phi_obj, &t_obj, &nu_obj, &jet_type, &spec_type,
                &theta_obs, &E_iso_core, &theta_h_core, &theta_h_wing, &b, &L0,
                    &q, &ts, &n_0, &p, &epsilon_E, &epsilon_B, &ksi_N, &d_L,
                    &g0,
                &E_core_global, &theta_h_core_global,
                &tRes, &latRes, &int_type, &rtol_struct, &rtol_phi,
                    &rtol_theta, &nmax_phi, &nmax_theta,
                &mask_obj,
                &spread, &counterjet, &gamma_type))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    if(int_type < 0 || int_type >= INT_UNDEFINED)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "intType out of range, unknown integrator");
        return NULL;
    }

    //Grab NUMPY arrays
    PyArrayObject *theta_arr;
    PyArrayObject *phi_arr;
    PyArrayObject *t_arr;
    PyArrayObject *nu_arr;
    PyArrayObject *mask_arr = NULL;

    theta_arr = (PyArrayObject *) PyArray_FROM_OTF(theta_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    phi_arr = (PyArrayObject *) PyArray_FROM_OTF(phi_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    t_arr = (PyArrayObject *) PyArray_FROM_OTF(t_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    nu_arr = (PyArrayObject *) PyArray_FROM_OTF(nu_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    if(mask_obj != NULL)
        mask_arr = (PyArrayObject *) PyArray_FROM_OTF(mask_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);

    if(theta_arr == NULL || phi_arr == NULL || t_arr == NULL || nu_arr == NULL
            || (mask_obj != NULL && mask_arr == NULL))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not read input arrays.");
        Py_XDECREF(theta_arr);
        Py_XDECREF(phi_arr);
        Py_XDECREF(t_arr);
        Py_XDECREF(nu_arr);
        Py_XDECREF(mask_arr);
        return NULL;
    }

    int theta_ndim = (int) PyArray_NDIM(theta_arr);
    int phi_ndim = (int) PyArray_NDIM(phi_arr);
    int t_ndim = (int) PyArray_NDIM(t_arr);
    int nu_ndim = (int) PyArray_NDIM(nu_arr);
    int mask_ndim = 0;
    if(mask_obj != NULL)
        mask_ndim = (int) PyArray_NDIM(mask_arr);

    if(theta_ndim != 1 || phi_ndim != 1 || t_ndim != 1 || nu_ndim != 1
            || (mask_obj != NULL && mask_ndim != 1))
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be 1-D.");
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        if(mask_obj != NULL)
            Py_DECREF(mask_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(theta_arr, 0);
    int Np = (int)PyArray_DIM(phi_arr, 0);
    int Nt = (int)PyArray_DIM(t_arr, 0);
    int Nnu = (int)PyArray_DIM(nu_arr, 0);
    int Nmask = 0;
    if(mask_obj != NULL)
        Nmask = (int)PyArray_DIM(mask_arr, 0);

    if(N!=Np || N!=Nt || N!=Nnu || Np!=Nt || Np!=Nnu || Nt!=Nnu)
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be same size.");
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
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
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        Py_DECREF(mask_arr);
        return NULL;
    }

    double *theta = (double *)PyArray_DATA(theta_arr);
    double *phi = (double *)PyArray_DATA(phi_arr);
    double *t = (double *)PyArray_DATA(t_arr);
    double *nu = (double *)PyArray_DATA(nu_arr);
    double *mask = NULL;
    if(mask_obj != NULL)
        mask = (double *)PyArray_DATA(mask_arr);
    int masklen = Nmask/9;

    //Allocate output array

    npy_intp dims[1] = {N};
    PyObject *Inu_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    if(Inu_obj == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not make intensity array.");
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
        Py_DECREF(t_arr);
        Py_DECREF(nu_arr);
        return NULL;
    }
    double *Inu = PyArray_DATA((PyArrayObject *) Inu_obj);
    
    //Set Up The Parameters!
    
    double ta = t[0];
    double tb = t[0];
    int i;
    for(i=0; i<N; i++)
    {
        if(t[i] < ta)
            ta = t[i];
        else if(t[i] > tb)
            tb = t[i];
    }
    
    struct fluxParams fp;
    setup_fluxParams(&fp, d_L, theta_obs, E_iso_core, theta_h_core,
                        theta_h_wing, b,
                        L0, q, ts,
                        n_0, p, epsilon_E, epsilon_B, ksi_N, g0, 
                        E_core_global, theta_h_core_global, ta, tb,
                        tRes, latRes, int_type,
                        rtol_struct, rtol_phi, rtol_theta,
                        nmax_phi, nmax_theta,
                        spec_type, mask, masklen,
                        spread, counterjet, gamma_type);

    // Calculate the intensity!
    calc_intensity(jet_type, spec_type, theta, phi, t, nu, Inu, N, &fp);
   
    if(fp.error)
    {
        PyErr_SetString(PyExc_RuntimeError, fp.error_msg);
        free_fluxParams(&fp);
        return NULL;
    }
    
    // Free the parameters!
    free_fluxParams(&fp);

    // Clean up!
    Py_DECREF(theta_arr);
    Py_DECREF(phi_arr);
    Py_DECREF(t_arr);
    Py_DECREF(nu_arr);
    if(mask_obj != NULL)
        Py_DECREF(mask_arr);

    //Build output
    PyObject *ret = Py_BuildValue("N", Inu_obj);
    
    return ret;
}

static PyObject *jet_shockVals(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyObject *theta_obj = NULL;
    PyObject *phi_obj = NULL;
    PyObject *t_obj = NULL;
    PyObject *mask_obj = NULL;

    int jet_type = jet_type_default;
    double theta_obs = theta_obs_default;
    double E_iso_core = E_iso_core_default;
    double theta_h_core = theta_h_core_default;
    double theta_h_wing = theta_h_wing_default;
    double b = b_default;
    double L0 = L0_default;
    double q = q_default;
    double ts = ts_default; 
    double n_0 = n_0_default;
    double p = p_default;
    double epsilon_E = epsilon_E_default;
    double epsilon_B = epsilon_B_default;
    double ksi_N = ksi_N_default; 
    double d_L = d_L_default;

    int latRes = latRes_default;
    int tRes = tRes_default;
    int spread = spread_default;
    int counterjet = counterjet_default;
    int gamma_type = gamma_type_default;
    int spec_type = spec_type_default;
    double g0 = g0_default;
    double E_core_global = E_core_global_default;
    double theta_h_core_global = theta_h_core_global_default;

    double rtol_struct = rtol_struct_default;
    double rtol_theta = rtol_theta_default;
    double rtol_phi = rtol_phi_default;
    int int_type = int_type_default;
    int nmax_phi = nmax_phi_default;
    int nmax_theta = nmax_theta_default;

    static char *kwlist[] = {"theta", "phi", "t", "jetType", "specType",
                                "thetaObs", "E0", "thetaCore", "thetaWing",
                                    "b", "L0", "q", "ts", "n0", "p",
                                    "epsilon_e", "epsilon_B", "xi_N", "d_L",
                                    "g0",
                                "E0Global", "thetaCoreGlobal",
                                "tRes", "latRes", "intType", "rtolStruct",
                                    "rtolPhi", "rtolTheta", "NPhi", "NTheta",
                                "mask",
                                "spread", "counterjet", "gammaType",
                                NULL};

    //Parse Arguments
    if(!PyArg_ParseTupleAndKeywords(args, kwargs,
                "OOO|ii""ddddddddddddddd""dd""iiidddii""O""iii",
                kwlist,
                &theta_obj, &phi_obj, &t_obj, &jet_type, &spec_type,
                &theta_obs, &E_iso_core, &theta_h_core, &theta_h_wing, &b, &L0,
                    &q, &ts, &n_0, &p, &epsilon_E, &epsilon_B, &ksi_N, &d_L,
                    &g0,
                &E_core_global, &theta_h_core_global,
                &tRes, &latRes, &int_type, &rtol_struct, &rtol_phi,
                    &rtol_theta, &nmax_phi, &nmax_theta,
                &mask_obj,
                &spread, &counterjet, &gamma_type))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    if(int_type < 0 || int_type >= INT_UNDEFINED)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "intType out of range, unknown integrator");
        return NULL;
    }

    //Grab NUMPY arrays
    PyArrayObject *theta_arr;
    PyArrayObject *phi_arr;
    PyArrayObject *t_arr;
    PyArrayObject *mask_arr = NULL;

    theta_arr = (PyArrayObject *) PyArray_FROM_OTF(theta_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    phi_arr = (PyArrayObject *) PyArray_FROM_OTF(phi_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    t_arr = (PyArrayObject *) PyArray_FROM_OTF(t_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    if(mask_obj != NULL)
        mask_arr = (PyArrayObject *) PyArray_FROM_OTF(mask_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);

    if(theta_arr == NULL || phi_arr == NULL || t_arr == NULL
            || (mask_obj != NULL && mask_arr == NULL))
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not read input arrays.");
        Py_XDECREF(theta_arr);
        Py_XDECREF(phi_arr);
        Py_XDECREF(t_arr);
        Py_XDECREF(mask_arr);
        return NULL;
    }

    int theta_ndim = (int) PyArray_NDIM(theta_arr);
    int phi_ndim = (int) PyArray_NDIM(phi_arr);
    int t_ndim = (int) PyArray_NDIM(t_arr);
    int mask_ndim = 0;
    if(mask_obj != NULL)
        mask_ndim = (int) PyArray_NDIM(mask_arr);

    if(theta_ndim != 1 || phi_ndim != 1 || t_ndim != 1
            || (mask_obj != NULL && mask_ndim != 1))
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be 1-D.");
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
        Py_DECREF(t_arr);
        if(mask_obj != NULL)
            Py_DECREF(mask_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(theta_arr, 0);
    int Np = (int)PyArray_DIM(phi_arr, 0);
    int Nt = (int)PyArray_DIM(t_arr, 0);
    int Nmask = 0;
    if(mask_obj != NULL)
        Nmask = (int)PyArray_DIM(mask_arr, 0);

    if(N!=Np || N!=Nt || Np!=Nt)
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be same size.");
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
        Py_DECREF(t_arr);
        if(mask_obj != NULL)
            Py_DECREF(mask_arr);
        return NULL;
    }
    if(mask_obj != NULL && Nmask%9 != 0)
    {
        PyErr_SetString(PyExc_RuntimeError, 
                            "Mask length must be multiple of 9.");
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
        Py_DECREF(t_arr);
        if(mask_obj != NULL)
            Py_DECREF(mask_arr);
        return NULL;
    }

    double *theta = (double *)PyArray_DATA(theta_arr);
    double *phi = (double *)PyArray_DATA(phi_arr);
    double *t = (double *)PyArray_DATA(t_arr);
    double *mask = NULL;
    if(mask_obj != NULL)
        mask = (double *)PyArray_DATA(mask_arr);
    int masklen = Nmask/9;

    //Allocate output array

    npy_intp dims[1] = {N};
    PyObject *te_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *R_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *u_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *thj_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    if(te_obj == NULL || R_obj == NULL || u_obj == NULL || thj_obj == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not make arrays.");
        Py_DECREF(theta_arr);
        Py_DECREF(phi_arr);
        Py_DECREF(t_arr);
        if(mask_obj != NULL)
            Py_DECREF(mask_arr);
        return NULL;
    }
    double *te = PyArray_DATA((PyArrayObject *) te_obj);
    double *R = PyArray_DATA((PyArrayObject *) R_obj);
    double *u = PyArray_DATA((PyArrayObject *) u_obj);
    double *thj = PyArray_DATA((PyArrayObject *) thj_obj);

    // Set Up The Parameters
    double ta = t[0];
    double tb = t[0];
    int i;
    for(i=0; i<N; i++)
    {
        if(t[i] < ta)
            ta = t[i];
        else if(t[i] > tb)
            tb = t[i];
    }

    struct fluxParams fp;
    setup_fluxParams(&fp, d_L, theta_obs, E_iso_core, theta_h_core,
                        theta_h_wing, b,
                        L0, q, ts,
                        n_0, p, epsilon_E, epsilon_B, ksi_N, g0, 
                        E_core_global, theta_h_core_global, ta, tb,
                        tRes, latRes, int_type,
                        rtol_struct, rtol_phi, rtol_theta,
                        nmax_phi, nmax_theta,
                        spec_type, mask, masklen,
                        spread, counterjet, gamma_type);

    // Calculate the intensity!
    calc_shockVals(jet_type, theta, phi, t, te, R, u, thj, N, &fp);
   
    if(fp.error)
    {
        PyErr_SetString(PyExc_RuntimeError, fp.error_msg);
        free_fluxParams(&fp);
        return NULL;
    }
    
    //Free the parameters
    free_fluxParams(&fp);

    // Clean up!
    Py_DECREF(theta_arr);
    Py_DECREF(phi_arr);
    Py_DECREF(t_arr);
    if(mask_obj != NULL)
        Py_DECREF(mask_arr);

    //Build output
    PyObject *ret = Py_BuildValue("NNNN", te_obj, R_obj, u_obj, thj_obj);
    
    return ret;
}

static PyObject *jet_shock(PyObject *self, PyObject *args)
{
    double Rt0, Rt1, E0, n0, thetah, L0, q, ts;
    int tRes;
    int spread = 1;

    //Parse Arguments
    if(!PyArg_ParseTuple(args, "ddidddddd", &Rt0, &Rt1, &tRes, &E0, &n0,
                &thetah, &L0, &q, &ts))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    double ta = Rt0;
    double tb = Rt1;

    struct fluxParams pars;
    pars.ta = ta;
    pars.tb = tb;
    pars.n_0 = n0;
    pars.L0 = L0;
    pars.q = q;
    pars.ts = ts;
    pars.tRes = tRes;
    pars.E_tot = -1.0;
    pars.t_table = NULL;
    pars.R_table = NULL;
    pars.u_table = NULL;
    pars.th_table = NULL;
    pars.mu_table = NULL;
    pars.spread = spread;

    set_jet_params(&pars, E0, thetah);
   
    if(pars.error)
    {
        PyErr_SetString(PyExc_RuntimeError, pars.error_msg);
        free_fluxParams(&pars);
        return NULL;
    }
    pars.Rt0 = Rt0;
    pars.Rt1 = Rt1;
    make_R_table(&pars);
    if(pars.error)
    {
        PyErr_SetString(PyExc_RuntimeError, pars.error_msg);
        free_fluxParams(&pars);
        return NULL;
    }

    //Allocate output arrays
    int N = pars.table_entries;
    npy_intp dims[1] = {N};
    PyObject *t_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *R_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *u_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *th_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    if(t_obj == NULL || R_obj == NULL || u_obj == NULL || th_obj == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not make output arrays.");
        Py_XDECREF(t_obj);
        Py_XDECREF(R_obj);
        Py_XDECREF(u_obj);
        Py_XDECREF(th_obj);
        return NULL;
    }

    double *t = PyArray_DATA((PyArrayObject *) t_obj);
    double *R = PyArray_DATA((PyArrayObject *) R_obj);
    double *u = PyArray_DATA((PyArrayObject *) u_obj);
    double *th = PyArray_DATA((PyArrayObject *) th_obj);

    int i;
    for(i=0; i<N; i++)
    {
        t[i] = pars.t_table[i];
        R[i] = pars.R_table[i];
        u[i] = pars.u_table[i];
        th[i] = pars.th_table[i];
    }

    PyObject *ret = Py_BuildValue("NNNN", t_obj, R_obj, u_obj, th_obj);

    free_fluxParams(&pars);
    
    return ret;
}

static PyObject *jet_shockObs(PyObject *self, PyObject *args)
{
    double ta, tb, E0, n0, thetah, L0, q, ts;
    int tRes;
    int spread = 1;

    //Parse Arguments
    if(!PyArg_ParseTuple(args, "ddidddddd", &ta, &tb, &tRes, &E0, &n0, &thetah,
                                            &L0, &q, &ts))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }

    struct fluxParams pars;
    pars.ta = ta;
    pars.tb = tb;
    pars.n_0 = n0;
    pars.L0 = L0;
    pars.q = q;
    pars.ts = ts;
    pars.tRes = tRes;
    pars.E_tot = -1.0;
    pars.t_table = NULL;
    pars.R_table = NULL;
    pars.u_table = NULL;
    pars.th_table = NULL;
    pars.mu_table = NULL;
    pars.table_entries = 0;
    pars.t_table_inner = NULL;
    pars.R_table_inner = NULL;
    pars.u_table_inner = NULL;
    pars.th_table_inner = NULL;
    pars.mu_table_inner = NULL;
    pars.table_entries_inner = 0;
    pars.spread = spread;

    printf("set_jet_params\n");
    set_jet_params(&pars, E0, thetah);
    if(pars.error)
    {
        PyErr_SetString(PyExc_RuntimeError, pars.error_msg);
        free_fluxParams(&pars);
        return NULL;
    }
    printf("done\n");

    //Allocate output arrays
    int N = pars.table_entries;
    npy_intp dims[1] = {N};
    PyObject *t_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *R_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *u_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject *th_obj = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    if(t_obj == NULL || R_obj == NULL || u_obj == NULL || th_obj == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not make output arrays.");
        Py_XDECREF(t_obj);
        Py_XDECREF(R_obj);
        Py_XDECREF(u_obj);
        Py_XDECREF(th_obj);
        return NULL;
    }

    double *t = PyArray_DATA((PyArrayObject *) t_obj);
    double *R = PyArray_DATA((PyArrayObject *) R_obj);
    double *u = PyArray_DATA((PyArrayObject *) u_obj);
    double *th = PyArray_DATA((PyArrayObject *) th_obj);

    int i;
    for(i=0; i<N; i++)
    {
        t[i] = pars.t_table[i];
        R[i] = pars.R_table[i];
        u[i] = pars.u_table[i];
        th[i] = pars.th_table[i];
    }

    PyObject *ret = Py_BuildValue("NNNN", t_obj, R_obj, u_obj, th_obj);

    free_fluxParams(&pars);
    
    return ret;
}

static PyObject *jet_find_jet_edge(PyObject *self, PyObject *args)
{
    double tobs, phi, theta_obs, theta_0;
    PyObject *t_obj = NULL;
    PyObject *R_obj = NULL;
    PyObject *thj_obj = NULL;
    

    //Parse Arguments
    if(!PyArg_ParseTuple(args, "OOOdddd", &t_obj, &R_obj, &thj_obj, &tobs,
                         &phi, &theta_obs, &theta_0))
    {
        //PyErr_SetString(PyExc_RuntimeError, "Could not parse arguments.");
        return NULL;
    }
    
    PyArrayObject *t_arr;
    PyArrayObject *R_arr;
    PyArrayObject *thj_arr;
    t_arr = (PyArrayObject *) PyArray_FROM_OTF(t_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    R_arr = (PyArrayObject *) PyArray_FROM_OTF(R_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    thj_arr = (PyArrayObject *) PyArray_FROM_OTF(thj_obj, NPY_DOUBLE,
                                                NPY_ARRAY_IN_ARRAY);
    
    if(t_arr == NULL || R_arr == NULL || thj_arr == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Could not read input arrays.");
        Py_XDECREF(t_arr);
        Py_XDECREF(R_arr);
        Py_XDECREF(thj_arr);
        return NULL;
    }

    int t_ndim = (int) PyArray_NDIM(t_arr);
    int R_ndim = (int) PyArray_NDIM(R_arr);
    int thj_ndim = (int) PyArray_NDIM(thj_arr);

    if(t_ndim != 1 || R_ndim != 1 || thj_ndim != 1)
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be 1-D.");
        Py_DECREF(t_arr);
        Py_DECREF(R_arr);
        Py_DECREF(thj_arr);
        return NULL;
    }

    int N = (int)PyArray_DIM(t_arr, 0);
    int NR = (int)PyArray_DIM(R_arr, 0);
    int Nthj = (int)PyArray_DIM(thj_arr, 0);

    if(NR != N || Nthj != N)
    {
        PyErr_SetString(PyExc_RuntimeError, "Arrays must be same size.");
        Py_DECREF(t_arr);
        Py_DECREF(R_arr);
        Py_DECREF(thj_arr);
        return NULL;
    }

    double *t = (double *)PyArray_DATA(t_arr);
    double *R = (double *)PyArray_DATA(R_arr);
    double *thj = (double *)PyArray_DATA(thj_arr);

    double *mu = (double *)malloc(N * sizeof(double));
    int i;
    for(i=0; i<N; i++)
        mu[i] = (t[i] - tobs) * v_light / R[i];

    double th = find_jet_edge(phi, cos(theta_obs), sin(theta_obs), theta_0,
                              mu, thj, N);

    free(mu);

    PyObject *ret = Py_BuildValue("d", th);
    
    return ret;
}
