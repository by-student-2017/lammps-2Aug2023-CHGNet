/*
 * Copyright (C) 2020 AdvanceSoft Corporation
 *
 * This software is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */

#ifndef NNP_NNARCH_H_
#define NNP_NNARCH_H_

#include "nnp_common.h"
#include "nnp_property.h"
#include "nnp_symm_func.h"
#include "nnp_symm_func_manybody.h"
#ifdef _NNP_GPU
#include "nnp_symm_func_gpu_behler.h"
#include "nnp_symm_func_gpu_chebyshev.h"
#else
#include "nnp_symm_func_behler.h"
#include "nnp_symm_func_chebyshev.h"
#endif
#include "nnp_nnlayer.h"
#include "nnp_reax_pot.h"
#include "memory.h"

class NNArch
{
public:
    NNArch(int numElems, const Property* property, LAMMPS_NS::Memory* memory);
    virtual ~NNArch();

    int getAtomNum(int ielem) const
    {
        return this->atomNum[ielem];
    }

    void restoreNN(FILE* fp, char** elemNames, bool zeroEatom, int rank, MPI_Comm world);

    void initGeometry(int numAtoms, int* elements,
                      int* numNeighbor, int** elemNeighbor, nnpreal*** posNeighbor);

    SymmFunc* getSymmFunc();

    void calculateSymmFuncs();

    void initLayers();

    void goForwardOnEnergy();

    void goBackwardOnForce();

    void goForwardOnCharge();

    void obtainEnergies(nnpreal* energies) const;

    void obtainForces(nnpreal*** forces) const;

    void obtainCharges(nnpreal* charges) const;

    const nnpreal* getLJLikeA1() const
    {
        return this->ljlikeA1;
    }

    const nnpreal* getLJLikeA2() const
    {
        return this->ljlikeA2;
    }

    const nnpreal* getLJLikeA3() const
    {
        return this->ljlikeA3;
    }

    const nnpreal* getLJLikeA4() const
    {
        return this->ljlikeA4;
    }

    ReaxPot* getReaxPot() const
    {
        return this->reaxPot;
    }

private:
    int mode;
    int numElems;
    int numAtoms;

    const Property* property;

    LAMMPS_NS::Memory* memory;

    int* atomNum;

    int*       elements;
    int*       numNeighbor;
    int*       idxNeighbor;
    int**      elemNeighbor;
    nnpreal*** posNeighbor;

    int  sizeNumAtom;
    int  sizeTotNeigh;
    int* sizeNbatch;

    int* nbatch;
    int* ibatch;

    nnpreal***  energyData;
    nnpreal***  energyGrad;

    nnpreal*    forceData;

    nnpreal***  chargeData;

    nnpreal*    symmData;
    nnpreal*    symmDiff;
    nnpreal*    symmAve;
    nnpreal*    symmDev;
    SymmFunc*   symmFunc;

    NNLayer**** interLayersEnergy;
    NNLayer***  lastLayersEnergy;

    NNLayer**** interLayersCharge;
    NNLayer***  lastLayersCharge;

    nnpreal* ljlikeA1;
    nnpreal* ljlikeA2;
    nnpreal* ljlikeA3;
    nnpreal* ljlikeA4;

    ReaxPot* reaxPot;

    bool isEnergyMode() const
    {
        return this->mode == NNARCH_MODE_BOTH || this->mode == NNARCH_MODE_ENERGY;
    }

    bool isChargeMode() const
    {
        return this->mode == NNARCH_MODE_BOTH || this->mode == NNARCH_MODE_CHARGE;
    }
};

#endif /* NNP_NNARCH_H_ */
