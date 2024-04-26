/*
 * Copyright (C) 2022 AdvanceSoft Corporation
 *
 * This source code is licensed under the GNU General Public License Version 2
 * found in the LICENSE file in the root directory of this source tree.
 */

#ifdef PAIR_CLASS

PairStyle(m3gnet/d3/gpu, PairM3GNetD3GPU)

#else

#ifndef LMP_PAIR_M3GNET_D3_GPU_H_
#define LMP_PAIR_M3GNET_D3_GPU_H_

#include "pair_m3gnet.h"

namespace LAMMPS_NS
{

class PairM3GNetD3GPU: public PairM3GNet
{
public:
    PairM3GNetD3GPU(class LAMMPS*);

    virtual ~PairM3GNetD3GPU() override;

protected:
    int withDFTD3() override;

    int withGPU() override;
};

}  // namespace LAMMPS_NS

#endif /* LMP_PAIR_M3GNET_D3_GPU_H_ */
#endif
