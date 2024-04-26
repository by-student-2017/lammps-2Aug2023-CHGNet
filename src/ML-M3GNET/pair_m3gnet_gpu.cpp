/*
 * Copyright (C) 2024 AdvanceSoft Corporation
 *
 * This source code is licensed under the GNU General Public License Version 2
 * found in the LICENSE file in the root directory of this source tree.
 */

#include "pair_m3gnet_gpu.h"

using namespace LAMMPS_NS;

PairM3GNetGPU::PairM3GNetGPU(LAMMPS *lmp) : PairM3GNet(lmp)
{
    if (copymode)
    {
        return;
    }

    // NOP
}

PairM3GNetGPU::~PairM3GNetGPU()
{
    // NOP
}

int PairM3GNetGPU::withDFTD3()
{
    return 0;
}

int PairM3GNetGPU::withGPU()
{
    return 1;
}

