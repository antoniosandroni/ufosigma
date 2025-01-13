/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_ACTIONS_ROOBSERRINFLATION_H_
#define UFO_FILTERS_ACTIONS_ROOBSERRINFLATION_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "ufo/filters/actions/FilterActionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
}

namespace ufo {

class ObsFilterData;

// -----------------------------------------------------------------------------

class ROobserrInflationParameters : public FilterActionParametersBase {
  OOPS_CONCRETE_PARAMETERS(ROobserrInflationParameters, FilterActionParametersBase);

  // No extra parameters needed
};

// -----------------------------------------------------------------------------

class ROobserrInflation : public FilterActionBase {
 public:
  /// The type of parameters accepted by the constructor of this action.
  /// This typedef is used by the FilterActionFactory.
  typedef ROobserrInflationParameters Parameters_;

  explicit ROobserrInflation(const Parameters_ &);
  ~ROobserrInflation() {}

  void apply(const Variables &, const std::vector<std::vector<bool>> &,
             ObsFilterData &, int,
             ioda::ObsDataVector<int> &, ioda::ObsDataVector<float> &) const override;
  const ufo::Variables & requiredVariables() const override {return allvars_;}
  bool modifiesQCFlags() const override { return false; }

 private:
  Variables allvars_;
  Parameters_ parameters_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_ACTIONS_ROOBSERRINFLATION_H_
