/*
 * Copyright (C) 2024 AdvanceSoft Corporation
 *
 * This source code is licensed under the GNU General Public License Version 2
 * found in the LICENSE file in the root directory of this source tree.
 */

#ifdef PAIR_CLASS

PairStyle(m3gnet/gpu, PairM3GNetGPU)

#else

#ifndef LMP_PAIR_M3GNET_GPU_H_
#define LMP_PAIR_M3GNET_GPU_H_

#include "pair_m3gnet.h"

namespace LAMMPS_NS
{

class PairM3GNetGPU: public PairM3GNet
{
public:
    PairM3GNetGPU(class LAMMPS*);

    virtual ~PairM3GNetGPU() override;

protected:
    int withDFTD3() override;

    int withGPU() override;
};

}  // namespace LAMMPS_NS

#endif /* LMP_PAIR_M3GNET_GPU_H_ */
#endif
