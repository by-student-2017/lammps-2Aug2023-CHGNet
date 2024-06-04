/*
 * Copyright (C) 2023 AdvanceSoft Corporation
 *
 * This source code is licensed under the GNU General Public License Version 2
 * found in the LICENSE file in the root directory of this source tree.
 */

/*
 * The rewrite about Magmom was done by "By Student".
 */

#include "pair_chgnet.h"

using namespace LAMMPS_NS;

PairCHGNet::PairCHGNet(LAMMPS *lmp) : Pair(lmp)
{
    single_enable           = 0;
    restartinfo             = 0;
    one_coeff               = 1;
    manybody_flag           = 1;
    no_virial_fdotr_compute = 1;
    centroidstressflag      = CENTROID_NOTAVAIL;

    this->atomNumMap        = nullptr;
    this->maxinum           = 10;
    this->initializedPython = 0;
    this->cutoff            = 0.0;
    this->npythonPath       = 0;
    this->pythonPaths       = nullptr;
    this->pyModule          = nullptr;
    this->pyFunc            = nullptr;

    if (!PairCHGNet::finalized)
    {
        PairCHGNet::finalized = 1;
        std::atexit(PairCHGNet::finalize);
    }
}

PairCHGNet::~PairCHGNet()
{
    if (copymode)
    {
        return;
    }

    if (this->atomNumMap != nullptr)
    {
        delete[] this->atomNumMap;
    }

    if (allocated)
    {
        memory->destroy(cutsq);
        memory->destroy(setflag);
        memory->destroy(this->cell);
        memory->destroy(this->atomNums);
        memory->destroy(this->positions);
        memory->destroy(this->forces);
        memory->destroy(this->stress);
        memory->destroy(this->magmoms);
        memory->destroy(this->charges);
    }

    if (this->pythonPaths != nullptr)
    {
        for (int i = 0; i < this->npythonPath; ++i)
        {
            delete[] this->pythonPaths[i];
        }

        delete[] this->pythonPaths;
    }

    if (this->initializedPython)
    {
        this->finalizePython();
    }
}

int PairCHGNet::finalized = 0;

void PairCHGNet::finalize()
{
    if(Py_IsInitialized()) Py_Finalize();
}

void PairCHGNet::allocate()
{
    allocated = 1;

    const int ntypes = atom->ntypes;

    memory->create(cutsq,   ntypes + 1, ntypes + 1,   "pair:cutsq");
    memory->create(setflag, ntypes + 1, ntypes + 1,   "pair:setflag");

    memory->create(this->cell,      3, 3,             "pair:cell");
    memory->create(this->atomNums,  this->maxinum,    "pair:atomNums");
    memory->create(this->positions, this->maxinum, 3, "pair:positions");
    memory->create(this->forces,    this->maxinum, 3, "pair:forces");
    memory->create(this->stress,    6,                "pair:stress");
    memory->create(this->magmoms,   this->maxinum,    "pair:magmoms");
	memory->create(this->charges,   this->maxinum,    "pair:charges");
}

void PairCHGNet::compute(int eflag, int vflag)
{
    ev_init(eflag, vflag);

    if (eflag_atom)
    {
        error->all(FLERR, "Pair style CHGNet does not support atomic energy");
    }

    if (vflag_atom)
    {
        error->all(FLERR, "Pair style CHGNet does not support atomic virial pressure");
    }

    this->prepareGNN();

    this->performGNN();
}

void PairCHGNet::prepareGNN()
{
    int i;
    int iatom;

    int*  type = atom->type;
    double** x = atom->x;

    int  inum  = list->inum;
    int* ilist = list->ilist;

    double* boxlo = domain->boxlo;

    double *q = atom->q; // charge

    // grow with inum and nneigh
    if (inum > this->maxinum)
    {
        this->maxinum = inum + this->maxinum / 2;

        memory->grow(this->atomNums,  this->maxinum,    "pair:atomNums");
        memory->grow(this->positions, this->maxinum, 3, "pair:positions");
        memory->grow(this->forces,    this->maxinum, 3, "pair:forces");
        memory->grow(this->magmoms,   this->maxinum,    "pair:magmoms");
        memory->grow(this->charges,   this->maxinum,    "pair:charges");
    }

    // set cell
    this->cell[0][0] = domain->h[0]; // xx
    this->cell[1][1] = domain->h[1]; // yy
    this->cell[2][2] = domain->h[2]; // zz
    this->cell[2][1] = domain->h[3]; // yz
    this->cell[2][0] = domain->h[4]; // xz
    this->cell[1][0] = domain->h[5]; // xy
    this->cell[0][1] = 0.0;
    this->cell[0][2] = 0.0;
    this->cell[1][2] = 0.0;

    // set atomNums and positions
    #pragma omp parallel for private(iatom, i)
    for (iatom = 0; iatom < inum; ++iatom)
    {
        i = ilist[iatom];

        this->atomNums[iatom] = this->atomNumMap[type[i]];

        this->positions[iatom][0] = x[i][0] - boxlo[0];
        this->positions[iatom][1] = x[i][1] - boxlo[1];
        this->positions[iatom][2] = x[i][2] - boxlo[2];
    }

    // set atomic charges
    for (iatom = 0; iatom < inum; ++iatom)
    {
        i = ilist[iatom];

        //this->magmoms[iatom] = q[i];
        this->charges[iatom] = q[i];
    }
}

void PairCHGNet::performGNN()
{
    int i;
    int iatom;

    double** f = atom->f;

    int  inum  = list->inum;
    int* ilist = list->ilist;

    double volume;
    double evdwl = 0.0;

    double *q = atom->q; // charge

    // perform Graph Neural Network Potential of CHGNet
    evdwl = this->calculatePython();

    // set total energy
    if (eflag_global)
    {
        eng_vdwl += evdwl;
    }

    // set atomic forces
    for (iatom = 0; iatom < inum; ++iatom)
    {
        i = ilist[iatom];

        f[i][0] += this->forces[iatom][0];
        f[i][1] += this->forces[iatom][1];
        f[i][2] += this->forces[iatom][2];
    }

    // set virial pressure
    if (vflag_global)
    {
        volume = domain->xprd * domain->yprd * domain->zprd;

        virial[0] -= volume * this->stress[0]; // xx
        virial[1] -= volume * this->stress[1]; // yy
        virial[2] -= volume * this->stress[2]; // zz
        virial[3] -= volume * this->stress[3]; // yz
        virial[4] -= volume * this->stress[4]; // xz
        virial[5] -= volume * this->stress[5]; // xy
    }

    // set atomic charges
    for (iatom = 0; iatom < inum; ++iatom)
    {
        i = ilist[iatom];

        //q[i] = this->magmoms[iatom];
        q[i] = this->charges[iatom];
    }
}

void PairCHGNet::settings(int narg, char **arg)
{
    if (comm->nprocs > 1)
    {
        error->all(FLERR, "Pair style CHGNet does not support MPI parallelization");
    }

    if (narg < 1)
    {
        return;
    }

    this->npythonPath = narg;
    this->pythonPaths = new char*[this->npythonPath];

    for (int i = 0; i < this->npythonPath; ++i)
    {
        this->pythonPaths[i] = new char[512];
        strcpy(this->pythonPaths[i], arg[i]);
    }
}

void PairCHGNet::coeff(int narg, char **arg)
{
    int i, j;
    int iarg;
    int count;

    int ntypes = atom->ntypes;
    int ntypesEff;

    int as_path = 0;
    int dftd3   = withDFTD3();
    int gpu     = withGPU();

    if (narg != (3 + ntypes) && narg != (4 + ntypes))
    {
        error->all(FLERR, "Incorrect number of arguments for pair_coeff.");
    }

    if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    {
        error->all(FLERR, "Only wildcard asterisk is allowed in place of atom types for pair_coeff.");
    }

    if (strcmp(arg[2], "path") == 0)
    {
        iarg    = 4;
        as_path = 1;
    }
    else
    {
        iarg    = 3;
        as_path = 0;
    }

    if (this->atomNumMap != nullptr)
    {
        delete this->atomNumMap;
    }

    this->atomNumMap = new int[ntypes + 1];

    ntypesEff = 0;
    for (i = 0; i < ntypes; ++i)
    {
        if (strcmp(arg[i + iarg], "NULL") == 0)
        {
            this->atomNumMap[i + 1] = 0;
        }
        else
        {
            this->atomNumMap[i + 1] = this->elementToAtomNum(arg[i + iarg]);
            ntypesEff++;
        }
    }

    if (ntypesEff < 1)
    {
        error->all(FLERR, "There are no elements for pair_coeff of CHGNet.");
    }

    if (!allocated)
    {
        allocate();
    }

    if (this->initializedPython)
    {
        this->finalizePython();
    }

    this->cutoff = this->initializePython(arg[iarg - 1], as_path, dftd3, gpu);

    if (this->cutoff <= 0.0)
    {
        error->all(FLERR, "Cutoff is not positive for pair_coeff of CHGNet.");
    }

    count = 0;

    for (i = 1; i <= ntypes; ++i)
    {
        for (j = i; j <= ntypes; ++j)
        {
            if (this->atomNumMap[i] > 0 && this->atomNumMap[j] > 0)
            {
                setflag[i][j] = 1;
                count++;
            }
            else
            {
                setflag[i][j] = 0;
            }
        }
    }

    if (count == 0)
    {
        error->all(FLERR, "Incorrect args for pair coefficients");
    }
}

double PairCHGNet::init_one(int i, int j)
{
    if (setflag[i][j] == 0)
    {
        error->all(FLERR, "All pair coeffs are not set");
    }

    double r, rr;

    r = this->cutoff;
    rr = r * r;

    cutsq[i][j] = rr;
    cutsq[j][i] = rr;

    return r;
}

void PairCHGNet::init_style()
{
    if (strcmp(update->unit_style, "metal") != 0)
    {
        error->all(FLERR, "Pair style CHGNet requires 'units metal'");
    }

    int* periodicity = domain->periodicity;

    if (!(periodicity[0] && periodicity[1] && periodicity[2]))
    {
        error->all(FLERR, "Pair style CHGNet requires periodic boundary condition");
    }

    neighbor->add_request(this, NeighConst::REQ_FULL);
}

int PairCHGNet::withDFTD3()
{
    return 0;
}

int PairCHGNet::withGPU()
{
    return 0;
}

void PairCHGNet::finalizePython()
{
    if (this->initializedPython == 0)
    {
        return;
    }

    this->initializedPython = 0;

    Py_XDECREF(this->pyFunc);
    Py_XDECREF(this->pyModule);

    // Py_Finalize() does not work well, because of the following bug:
    // https://python.readthedocs.io/fr/latest/c-api/init.html#c.Py_FinalizeEx
    //if(Py_IsInitialized()) Py_Finalize();
}

double PairCHGNet::initializePython(const char *name, int as_path, int dftd3, int gpu)
{
    if (this->initializedPython != 0)
    {
        return this->cutoff;
    }

    double cutoff = -1.0;

    PyObject* pySys    = nullptr;
    PyObject* pyPath   = nullptr;
    PyObject* pyName   = nullptr;
    PyObject* pyModule = nullptr;
    PyObject* pyFunc   = nullptr;
    PyObject* pyArgs   = nullptr;
    PyObject* pyArg1   = nullptr;
    PyObject* pyArg2   = nullptr;
    PyObject* pyArg3   = nullptr;
    PyObject* pyArg4   = nullptr;
    PyObject* pyValue  = nullptr;

    if (!Py_IsInitialized()) Py_Initialize();

    pySys  = PyImport_ImportModule("sys");
    pyPath = PyObject_GetAttrString(pySys, "path");

    pyName = PyUnicode_DecodeFSDefault(".");
    if (pyName != nullptr)
    {
        PyList_Append(pyPath, pyName);
        Py_DECREF(pyName);
    }

    if (this->pythonPaths != nullptr)
    {
        for (int i = 0; i < this->npythonPath; ++i)
        {
            pyName = PyUnicode_DecodeFSDefault(this->pythonPaths[i]);
            if (pyName != nullptr)
            {
                PyList_Append(pyPath, pyName);
                Py_DECREF(pyName);
            }
        }
    }

    pyName = PyUnicode_DecodeFSDefault("chgnet_driver");
    if (pyName != nullptr)
    {
        pyModule = PyImport_Import(pyName);
        Py_DECREF(pyName);
    }

    if (pyModule != nullptr)
    {
        pyFunc = PyObject_GetAttrString(pyModule, "chgnet_initialize");

        if (pyFunc != nullptr && PyCallable_Check(pyFunc))
        {
            pyArg1 = PyUnicode_FromString(name);
            pyArg2 = PyBool_FromLong(as_path);
            pyArg3 = PyBool_FromLong(dftd3);
            pyArg4 = PyBool_FromLong(gpu);

            pyArgs = PyTuple_New(4);
            PyTuple_SetItem(pyArgs, 0, pyArg1);
            PyTuple_SetItem(pyArgs, 1, pyArg2);
            PyTuple_SetItem(pyArgs, 2, pyArg3);
            PyTuple_SetItem(pyArgs, 3, pyArg4);

            pyValue = PyObject_CallObject(pyFunc, pyArgs);

            Py_DECREF(pyArgs);

            if (pyValue != nullptr && PyFloat_Check(pyValue))
            {
                this->initializedPython = 1;
                cutoff = PyFloat_AsDouble(pyValue);
            }
            else
            {
                if (PyErr_Occurred()) PyErr_Print();
            }

            Py_XDECREF(pyValue);
        }

        else
        {
            if (PyErr_Occurred()) PyErr_Print();
        }

        Py_XDECREF(pyFunc);

        //pyFunc = PyObject_GetAttrString(pyModule, "chgnet_get_energy_forces_stress");
        pyFunc = PyObject_GetAttrString(pyModule, "chgnet_get_energy_forces_stress_magmoms");

        if (pyFunc != nullptr && PyCallable_Check(pyFunc))
        {
            // NOP
        }
        else
        {
            this->initializedPython = 0;
            if (PyErr_Occurred()) PyErr_Print();
        }

        //Py_XDECREF(pyFunc);
        //Py_DECREF(pyModule);
    }

    else
    {
        if (PyErr_Occurred()) PyErr_Print();
    }

    if (this->initializedPython == 0)
    {
        Py_XDECREF(pyFunc);
        Py_XDECREF(pyModule);

        if(Py_IsInitialized()) Py_Finalize();

        error->all(FLERR, "Cannot initialize python for pair_coeff of CHGNet.");
    }

    this->pyModule = pyModule;
    this->pyFunc   = pyFunc;

    return cutoff;
}

static const double ALL_ELEMENTS_ION[] = { 0.0,
     1.0,  0.0,  1.0,  2.0,  3.0,  0.0, -3.0, -2.0, -1.0,  0.0,  1.0,  2.0,  3.0,  4.0, -3.0, -2.0,
    -1.0,  0.0,  1.0,  2.0,  3.0,  3.5,  4.0,  2.5,  3.0,  2.5,  2.5,  2.5,  1.5,  2.0,  3.0,  0.0,
    -3.0, -2.0, -1.0,  0.0,  1.0,  2.0,  3.0,  4.0,  4.0,  6.0,  7.0,  3.5,  3.0,  3.0,  1.0,  2.0,
     3.0,  3.0,  4.0, -2.0, -1.0,  0.0,  1.0,  2.0,  3.0,  3.0,  3.0,  3.0,  3.0,  2.5,  2.5,  3.0,
     3.0,  3.0,  3.0,  3.0,  3.0,  3.0,  3.0,  4.0,  5.0,  6.0,  7.0,  4.0,  4.0,  3.0,  2.0,  2.0,
     2.0,  3.0,  3.0,  3.0, -1.0,  0.0,  1.0,  2.0,  3.0,  4.0,  4.5,  5.0,  5.0,  4.0,  3.5,  3.0,
     3.0,  3.0,  3.0,  3.0,  2.5,  2.5,  3.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  0.0,  0.0,  0.0,  0.0
};

double PairCHGNet::calculatePython()
{
    int i;
    int iatom;
    int natom = list->inum;

    double energy = 0.0;
    int hasEnergy = 0;
    int hasForces = 0;
    int hasStress = 0;
    int hasMagmoms= 0;

    PyObject* pyFunc  = this->pyFunc;
    PyObject* pyArgs  = nullptr; // call function
    PyObject* pyArg1  = nullptr; // set cell -> pyArg1
    PyObject* pyArg2  = nullptr; // set atomNums -> pyArg2
    PyObject* pyArg3  = nullptr; // set positions -> pyArg3
    PyObject* pyArg4  = nullptr; // set charge(=magmoms) -> pyArg4
    PyObject* pyAsub  = nullptr; // 3 component for lattice and position
    PyObject* pyValue = nullptr;
    PyObject* pyVal1  = nullptr; // get energy <- pyValue
    PyObject* pyVal2  = nullptr; // get forces <- pyValue
    PyObject* pyVal3  = nullptr; // get stress <- pyValue
    PyObject* pyVal4  = nullptr; // get magmoms <- pyValue
    PyObject* pyVsub  = nullptr; // force 3x3
    PyObject* pyVobj  = nullptr; // values from GNN

    // set cell -> pyArg1
    pyArg1 = PyList_New(3);

    for (i = 0; i < 3; ++i)
    {
        pyAsub = PyList_New(3);
        PyList_SetItem(pyAsub, 0, PyFloat_FromDouble(this->cell[i][0]));
        PyList_SetItem(pyAsub, 1, PyFloat_FromDouble(this->cell[i][1]));
        PyList_SetItem(pyAsub, 2, PyFloat_FromDouble(this->cell[i][2]));
        PyList_SetItem(pyArg1, i, pyAsub);
    }

    // set atomNums -> pyArg2
    pyArg2 = PyList_New(natom);

    for (iatom = 0; iatom < natom; ++iatom)
    {
        PyList_SetItem(pyArg2, iatom, PyLong_FromLong(this->atomNums[iatom]));
    }

    // set positions -> pyArg3
    pyArg3 = PyList_New(natom);

    for (iatom = 0; iatom < natom; ++iatom)
    {
        pyAsub = PyList_New(3);
        PyList_SetItem(pyAsub, 0, PyFloat_FromDouble(this->positions[iatom][0]));
        PyList_SetItem(pyAsub, 1, PyFloat_FromDouble(this->positions[iatom][1]));
        PyList_SetItem(pyAsub, 2, PyFloat_FromDouble(this->positions[iatom][2]));
        PyList_SetItem(pyArg3, iatom, pyAsub);
    }

    // set atomNums -> pyArg4
    pyArg4 = PyList_New(natom);
    
    for (iatom = 0; iatom < natom; ++iatom)
    {
        PyList_SetItem(pyArg4, iatom, PyLong_FromLong(this->magmoms[iatom]));
    }

    // call function
    pyArgs = PyTuple_New(4);
    PyTuple_SetItem(pyArgs, 0, pyArg1);
    PyTuple_SetItem(pyArgs, 1, pyArg2);
    PyTuple_SetItem(pyArgs, 2, pyArg3);
    PyTuple_SetItem(pyArgs, 3, pyArg4);

    pyValue = PyObject_CallObject(pyFunc, pyArgs); // get values from GNN

    Py_DECREF(pyArgs);

    if (pyValue != nullptr && PyTuple_Check(pyValue) && PyTuple_Size(pyValue) >= 4)
    {
        // get energy <- pyValue
        pyVal1 = PyTuple_GetItem(pyValue, 0);
        if (pyVal1 != nullptr && PyFloat_Check(pyVal1))
        {
            hasEnergy = 1;
            energy = PyFloat_AsDouble(pyVal1);
        }
        else
        {
            if (PyErr_Occurred()) PyErr_Print();
        }

        // get forces <- pyValue
        pyVal2 = PyTuple_GetItem(pyValue, 1);
        if (pyVal2 != nullptr && PyList_Check(pyVal2) && PyList_Size(pyVal2) >= natom)
        {
            hasForces = 1;

            for (iatom = 0; iatom < natom; ++iatom)
            {
                pyVsub = PyList_GetItem(pyVal2, iatom);
                if (pyVsub != nullptr && PyList_Check(pyVsub) && PyList_Size(pyVsub) >= 3)
                {
                    for (i = 0; i < 3; ++i)
                    {
                        pyVobj = PyList_GetItem(pyVsub, i);
                        if (pyVobj != nullptr && PyFloat_Check(pyVobj))
                        {
                            this->forces[iatom][i] = PyFloat_AsDouble(pyVobj);
                        }
                        else
                        {
                            if (PyErr_Occurred()) PyErr_Print();
                            hasForces = 0;
                            break;
                        }
                    }
                }
                else
                {
                    if (PyErr_Occurred()) PyErr_Print();
                    hasForces = 0;
                    break;
                }

                if (hasForces == 0)
                {
                    break;
                }
            }
        }
        else
        {
            if (PyErr_Occurred()) PyErr_Print();
        }

        // get stress <- pyValue
        pyVal3 = PyTuple_GetItem(pyValue, 2);
        if (pyVal3 != nullptr && PyList_Check(pyVal3) && PyList_Size(pyVal3) >= 6)
        {
            hasStress = 1;

            for (i = 0; i < 6; ++i)
            {
                pyVobj = PyList_GetItem(pyVal3, i);
                if (pyVobj != nullptr && PyFloat_Check(pyVobj))
                {
                    this->stress[i] = PyFloat_AsDouble(pyVobj);
                }
                else
                {
                    if (PyErr_Occurred()) PyErr_Print();
                    hasStress = 0;
                    break;
                }
            }
        }
        else
        {
            if (PyErr_Occurred()) PyErr_Print();
        }

        // get magmoms <- pyValue and convert magmoms to charges
        pyVal4 = PyTuple_GetItem(pyValue, 3);
        if (pyVal4 != nullptr && PyList_Check(pyVal4) && PyList_Size(pyVal4) >= natom)
        {
            hasMagmoms = 1;

            for (iatom = 0; iatom < natom; ++iatom)
            {
                pyVobj = PyList_GetItem(pyVal4, iatom);
                if (pyVobj != nullptr && PyFloat_Check(pyVobj))
                {
                    this->magmoms[iatom] = PyFloat_AsDouble(pyVobj);
                    // An approximate formula for calculating charge from atomic number and MAGMOM can be inserted here.
                    // Currently a bad alternative.
                    if (ALL_ELEMENTS_ION[this->atomNums[iatom]] > 0){
                        this->charges[iatom] = ALL_ELEMENTS_ION[this->atomNums[iatom]] - this->magmoms[iatom]*2.0;
                    }else{
                        this->charges[iatom] = ALL_ELEMENTS_ION[this->atomNums[iatom]] + this->magmoms[iatom]*2.0;
                    }
                    if (this->atomNums[iatom] == 22){ // Ti
                        this->charges[iatom] = 4.0 - this->magmoms[iatom]*2.0;
                    }else if (this->atomNums[iatom] == 23 || this->atomNums[iatom] == 41){ // V, Nb
                        this->charges[iatom] = 5.0 - this->magmoms[iatom]*2.0;
                    }else if (this->atomNums[iatom] == 24 || 
                        this->atomNums[iatom] == 26 || this->atomNums[iatom] == 27 || this->atomNums[iatom] == 28 ||
                        this->atomNums[iatom] == 79){ // Cr, Fe, Co, Ni, Au
                        this->charges[iatom] = this->magmoms[iatom]+1;
                        if (this->magmoms[iatom]+1 > 3.0){
                            this->charges[iatom] = 3.0 - (this->magmoms[iatom]+1) + 2.8;
                        }
                    }else if (this->atomNums[iatom] == 25 || this->atomNums[iatom] == 78 ||
                        this->atomNums[iatom] == 44 || this->atomNums[iatom] == 46){ // Mn, Pt, Ru, Pd
                        this->charges[iatom] = this->magmoms[iatom]+1;
                        if (this->magmoms[iatom]+1 > 4.0){
                            this->charges[iatom] = 4.0 - (this->magmoms[iatom]+1) + 3.8;
                        }
                    }else if (this->atomNums[iatom] == 29){ // Cu
                        this->charges[iatom] = this->magmoms[iatom]+1;
                        if (this->magmoms[iatom]+1 > 2.0){
                            this->charges[iatom] = 2.0 - (this->magmoms[iatom]+1) + 1.8;
                        }
                    }
                    // 
                }
                else
                {
                    if (PyErr_Occurred()) PyErr_Print();
                    hasMagmoms = 0;
                    break;
                }
            }
        }
        else
        {
            if (PyErr_Occurred()) PyErr_Print();
        }
    }

    else
    {
        if (PyErr_Occurred()) PyErr_Print();
    }

    Py_XDECREF(pyValue);

    if (hasEnergy == 0 || hasForces == 0 || hasStress == 0 || hasMagmoms == 0)
    {
        error->all(FLERR, "Cannot calculate energy, forces, stress and magmoms by python of CHGNet.");
    }

    return energy;
}

static const int NUM_ELEMENTS = 118;

static const char* ALL_ELEMENTS[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
    "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};

//static const double ALL_ELEMENTS_ION[] = { 0.0,
//     1.0,  0.0,  1.0,  2.0,  3.0,  0.0, -3.0, -2.0, -1.0,  0.0,  1.0,  2.0,  3.0,  4.0, -3.0, -2.0,
//    -1.0,  0.0,  1.0,  2.0,  3.0,  3.5,  4.0,  2.5,  3.0,  2.5,  2.5,  2.5,  1.5,  2.0,  3.0,  0.0,
//    -3.0, -2.0, -1.0,  0.0,  1.0,  2.0,  3.0,  4.0,  4.0,  6.0,  7.0,  3.5,  3.0,  3.0,  1.0,  2.0,
//     3.0,  3.0,  4.0, -2.0, -1.0,  0.0,  1.0,  2.0,  3.0,  3.0,  3.0,  3.0,  3.0,  2.5,  2.5,  3.0,
//     3.0,  3.0,  3.0,  3.0,  3.0,  3.0,  3.0,  4.0,  5.0,  6.0,  7.0,  4.0,  4.0,  3.0,  2.0,  2.0,
//     2.0,  3.0,  3.0,  3.0, -1.0,  0.0,  1.0,  2.0,  3.0,  4.0,  4.5,  5.0,  5.0,  4.0,  3.5,  3.0,
//     3.0,  3.0,  3.0,  3.0,  2.5,  2.5,  3.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
//     0.0,  0.0,  0.0,  0.0,  0.0,  0.0
//};

int PairCHGNet::elementToAtomNum(const char *elem)
{
    char elem1[16];

    strcpy(elem1, elem);

    this->toRealElement(elem1);

    if (strlen(elem1) > 0)
    {
        for (int i = 0; i < NUM_ELEMENTS; ++i)
        {
            if (strcasecmp(elem1, ALL_ELEMENTS[i]) == 0)
            {
                return (i + 1);
            }
        }
    }

    char estr[256];
    sprintf(estr, "Incorrect name of element: %s", elem);
    error->all(FLERR, estr);

    return 0;
}

void PairCHGNet::toRealElement(char *elem)
{
    int n = strlen(elem);
    n = n > 2 ? 2 : n;

    int m = n;

    for (int i = 0; i < n; ++i)
    {
        char c = elem[i];
        if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
            c == '5' || c == '6' || c == '7' || c == '8' || c == '9' || c == ' ' ||
            c == '_' || c == '-' || c == '+' || c == '*' || c == '~' || c == ':' || c == '#')
        {
            m = i;
            break;
        }

        elem[i] = c;
    }

    elem[m] = '\0';
}

