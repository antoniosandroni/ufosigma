/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORPARAMETERSBASE_H_
#define UFO_ERRORS_OBSERRORPARAMETERSBASE_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace ufo {

/// \brief Base obs errors parameters class
class ObsErrorParametersBase : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorParametersBase, Parameters)
 public:
  /// \brief Name of the covariance model.
  oops::Parameter<std::string> model{"covariance model", "diagonal", this};
};

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORPARAMETERSBASE_H_
