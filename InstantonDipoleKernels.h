// -*- C++ -*-
#ifndef Herwig_InstantonDipoleKernels_H
#define Herwig_InstantonDipoleKernels_H

#include "Herwig/Shower/Dipole/Kernels/FFMqx2qgxDipoleKernel.h"
#include "Herwig/Shower/Dipole/Kernels/FIMqx2qgxDipoleKernel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Preserve outgoing quark ParticleData in final-final radiation.
 *
 * The campaign uses massless canonical quark data for five-flavour ISR and
 * private physical-mass data for hard outgoing quarks. Herwig's stock kernel
 * replaces the emitter with its configured canonical flavour. This generic
 * adapter keeps the actual emitter data and also covers a massless secondary
 * quark recoiling against a massive zero mode. It otherwise inherits the
 * stock massive kernel.
 */
class InstantonFFMqx2qgxDipoleKernel : public FFMqx2qgxDipoleKernel {
public:
  InstantonFFMqx2qgxDipoleKernel() = default;
  InstantonFFMqx2qgxDipoleKernel(
    const InstantonFFMqx2qgxDipoleKernel &) = default;

  bool canHandle(const DipoleIndex &) const override;
  tcPDPtr emitter(const DipoleIndex &) const override;
  static void Init();

protected:
  IBPtr clone() const override;
  IBPtr fullclone() const override;
};

/**
 * Preserve a physical outgoing zero-mode quark in final-initial radiation.
 */
class InstantonFIMqx2qgxDipoleKernel : public FIMqx2qgxDipoleKernel {
public:
  InstantonFIMqx2qgxDipoleKernel() = default;
  InstantonFIMqx2qgxDipoleKernel(
    const InstantonFIMqx2qgxDipoleKernel &) = default;

  bool canHandle(const DipoleIndex &) const override;
  tcPDPtr emitter(const DipoleIndex &) const override;
  static void Init();

protected:
  IBPtr clone() const override;
  IBPtr fullclone() const override;
};

} // namespace Herwig

#endif // Herwig_InstantonDipoleKernels_H
