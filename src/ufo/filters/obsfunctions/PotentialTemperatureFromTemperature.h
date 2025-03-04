/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_POTENTIALTEMPERATUREFROMTEMPERATURE_H_
#define UFO_FILTERS_OBSFUNCTIONS_POTENTIALTEMPERATUREFROMTEMPERATURE_H_

#include <string>
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Calculate model potential temperature near a pressure level (Pa) or at surface.
///  Here is an exampe to derive potential temperature at surface:
///  - filter: Variable Assignment
///    assignments:
///    - name: DerivedMetaData/PotentialTemperatureSurface
///      type: float
///      function:
///        name: PotentialTemperatureFromTemperature@ObsFunction
///        options:
///          use surface pressure: true
///  Here is another exampe to derive potential temperature at model levels near 70000 Pa:
///  - filter: Variable Assignment
///    assignments:
///    - name: DerivedMetaData/PotentialTemperatureSurface
///      type: float
///      function:
///        name: PotentialTemperatureFromTemperature@ObsFunction
///        options:
///          use surface pressure: false
///          pressure to evaluate potential temperature: 70000
///
class PotentialTemperatureFromTemperatureParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(PotentialTemperatureFromTemperatureParameters, Parameters);

 public:
  oops::RequiredParameter<bool> useSurfacePressure{"use surface pressure", this};
  oops::OptionalParameter<float> geovalsPressure{"pressure to evaluate potential temperature",
                                 this};
};

class PotentialTemperatureFromTemperature : public ObsFunctionBase<float> {
 public:
  explicit PotentialTemperatureFromTemperature(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~PotentialTemperatureFromTemperature();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  PotentialTemperatureFromTemperatureParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_POTENTIALTEMPERATUREFROMTEMPERATURE_H_
