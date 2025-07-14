/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ERRORS_OBSERRORDIAGONAL_H_
#define UFO_ERRORS_OBSERRORDIAGONAL_H_

#include <memory>
#include <string>

#include "ioda/ObsVector.h"

#include "oops/interface/ObsErrorBase.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/ObsTraits.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

/// \brief Parameters for diagonal obs errors
class ObsErrorDiagonalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorDiagonalParameters, Parameters)
 public:
  /// \brief Name of the covariance model.
  oops::Parameter<std::string> model{"covariance model", "diagonal", this};
  /// perturbation amplitude multiplier
  oops::Parameter<double> pert{"random amplitude", 1.0, this};
  oops::Parameter<bool> zeroMeanPerturbations{"zero-mean perturbations", false, this};
  /// 1-based ensemble member index.
  /// Used (and required) only if `zero-mean perturbations` is set to true.
  oops::OptionalParameter<int> member{"member", this, {oops::minConstraint(1)}};
  /// Number of ensemble members.
  /// Used (and required) only if `zero-mean perturbations` is set to true.
  oops::OptionalParameter<int> numberOfMembers{"number of members", this, {oops::minConstraint(1)}};
};

// -----------------------------------------------------------------------------
/// \brief Diagonal observation error covariance matrix.
class ObsErrorDiagonal : public oops::interface::ObsErrorBase<ObsTraits> {
 public:
  /// The type of parameters for this class.
  typedef ObsErrorDiagonalParameters Parameters_;

  static const std::string classname() {return "ufo::ObsErrorDiagonal";}

  ObsErrorDiagonal(const eckit::Configuration &, ioda::ObsSpace &,
                   const eckit::mpi::Comm &timeComm);

/// Update after obs errors potentially changed
  void update(const ioda::ObsVector &) override;

/// Multiply a Departure by \f$R\f$
  void multiply(ioda::ObsVector &) const override;

/// Multiply a Departure by \f$R^{-1}\f$
  void inverseMultiply(ioda::ObsVector &) const override;

/// Generate random perturbation
  void randomize(ioda::ObsVector &) const override;

/// Save obs errors
  void save(const std::string &) const override;

/// Get mean error for Jo table
  double getRMSE() const override {return stddev_.rms();}

/// Get obs errors std deviation
  std::unique_ptr<ioda::ObsVector> getObsErrors() const override;

/// Return inverseVariance
  std::unique_ptr<ioda::ObsVector> getInverseVariance() const override;

 private:
  void print(std::ostream &) const override;
  ioda::ObsVector stddev_;
  ioda::ObsVector inverseVariance_;
  Parameters_ options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ERRORS_OBSERRORDIAGONAL_H_
