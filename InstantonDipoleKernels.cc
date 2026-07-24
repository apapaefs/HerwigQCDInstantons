// -*- C++ -*-

#include "InstantonDipoleKernels.h"

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include <cstdlib>

using namespace Herwig;

namespace {

bool isFinalLightQuark(tcPDPtr particle, bool initial) {
  const long id = std::abs(particle->id());
  return !initial && id >= 1 && id <= 5;
}

} // namespace

bool InstantonFFMqx2qgxDipoleKernel::canHandle(
  const DipoleIndex & index) const {
  return useThisKernel()
      && isFinalLightQuark(
           index.emitterData(), index.initialStateEmitter())
      && !(index.emitterData()->mass() == ZERO
           && index.spectatorData()->mass() == ZERO)
      && !index.initialStateSpectator()
      && !index.incomingDecayEmitter()
      && !index.incomingDecaySpectator();
}

tcPDPtr InstantonFFMqx2qgxDipoleKernel::emitter(
  const DipoleIndex & index) const {
  return index.emitterData();
}

bool InstantonFIMqx2qgxDipoleKernel::canHandle(
  const DipoleIndex & index) const {
  return useThisKernel()
      && isFinalLightQuark(
           index.emitterData(), index.initialStateEmitter())
      && index.emitterData()->mass() > ZERO
      && index.initialStateSpectator()
      && index.spectatorData()->mass() == ZERO;
}

tcPDPtr InstantonFIMqx2qgxDipoleKernel::emitter(
  const DipoleIndex & index) const {
  return index.emitterData();
}

#define HERWIG_IMPLEMENT_INSTANTON_DIPOLE_ADAPTER(Class, Base)             \
  IBPtr Class::clone() const {                                              \
    return new_ptr(*this);                                                  \
  }                                                                         \
                                                                            \
  IBPtr Class::fullclone() const {                                          \
    return new_ptr(*this);                                                  \
  }                                                                         \
                                                                            \
  DescribeClass<Class, Base> describe##Class(                               \
    "Herwig::" #Class, "CampaignInstantons.so");                            \
                                                                            \
  void Class::Init() {                                                      \
    static ClassDocumentation<Class> documentation(                         \
      "Preserves physical instanton zero-mode ParticleData while using "    \
      "the corresponding stock Herwig dipole splitting kernel.");          \
  }

HERWIG_IMPLEMENT_INSTANTON_DIPOLE_ADAPTER(
  InstantonFFMqx2qgxDipoleKernel, FFMqx2qgxDipoleKernel)
HERWIG_IMPLEMENT_INSTANTON_DIPOLE_ADAPTER(
  InstantonFIMqx2qgxDipoleKernel, FIMqx2qgxDipoleKernel)

#undef HERWIG_IMPLEMENT_INSTANTON_DIPOLE_ADAPTER
