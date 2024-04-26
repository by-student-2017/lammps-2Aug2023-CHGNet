/*
 * Copyright (C) 2024 AdvanceSoft Corporation
 *
 * This source code is licensed under the GNU General Public License Version 2
 * found in the LICENSE file in the root directory of this source tree.
 */

#include "pair_m3gnet_d3_gpu.h"

using namespace LAMMPS_NS;

PairM3GNetD3GPU::PairM3GNetD3GPU(LAMMPS *lmp) : PairM3GNet(lmp)
{
    if (copymode)
    {
        return;
    }

    // NOP
}

PairM3GNetD3GPU::~PairM3GNetD3GPU()
{
    // NOP
}

int PairM3GNetD3GPU::withDFTD3()
{
    return 1;
}

int PairM3GNetD3GPU::withGPU()
{
    return 1;
}

