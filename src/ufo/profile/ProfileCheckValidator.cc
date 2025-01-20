/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <map>
#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ufo/profile/ProfileCheckValidator.h"
#include "ufo/profile/ProfileVariableNames.h"

#include "ufo/utils/StringUtils.h"

namespace ufo {
  ProfileCheckValidator::ProfileCheckValidator
  (const ConventionalProfileProcessingParameters &options)
    : options_(options)
  {
    // Set offsets due to C++ and Fortran array index starting values
    comparison_offsets_[ufo::ProfileVariableNames::StdLev] = 1;
    comparison_offsets_[ufo::ProfileVariableNames::SigBelow] = 1;
    comparison_offsets_[ufo::ProfileVariableNames::SigAbove] = 1;
    comparison_offsets_[ufo::ProfileVariableNames::IndStd] = 1;
    comparison_offsets_[ufo::ProfileVariableNames::LevErrors] = 1;
    comparison_offsets_[ufo::ProfileVariableNames::Indx] = 1;

    // List of checks performed
    std::vector <std::string> checks = options_.Checks.value();

    // Loop over each check and populate lists of integer and float values to compare
    for (const auto& check : checks) {
      if (check == "Basic") {
      } else if (check == "SamePDiffT") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumAnyErrors,
              ufo::ProfileVariableNames::counter_NumSamePErrObs
              });
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_final_reject_t,
            ufo::ProfileVariableNames::diagflags_interpolation_t});
      } else if (check == "Sign") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumAnyErrors,
              ufo::ProfileVariableNames::counter_NumSignChange});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_final_reject_t,
            ufo::ProfileVariableNames::diagflags_data_correct_t});
      } else if (check == "UnstableLayer") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumAnyErrors,
            ufo::ProfileVariableNames::counter_NumSuperadiabat});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::PBottom});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_superadiabat_t});
      } else if (check == "Interpolation") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumAnyErrors,
            ufo::ProfileVariableNames::counter_NumInterpErrors,
            ufo::ProfileVariableNames::counter_NumInterpErrObs,
            ufo::ProfileVariableNames::NumStd,
            ufo::ProfileVariableNames::NumSig,
            ufo::ProfileVariableNames::StdLev,
            ufo::ProfileVariableNames::SigBelow,
            ufo::ProfileVariableNames::SigAbove,
            ufo::ProfileVariableNames::IndStd,
            ufo::ProfileVariableNames::LevErrors});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::tInterp,
            ufo::ProfileVariableNames::LogP});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_interpolation_t});
      } else if (check == "Hydrostatic") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumAnyErrors,
            ufo::ProfileVariableNames::counter_Num925Miss,
            ufo::ProfileVariableNames::counter_Num100Miss,
            ufo::ProfileVariableNames::counter_NumStdMiss,
            ufo::ProfileVariableNames::counter_NumHydErrObs,
            ufo::ProfileVariableNames::counter_NumIntHydErrors,
            ufo::ProfileVariableNames::HydError});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::DC,
            ufo::ProfileVariableNames::ETol});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_interpolation_t,
            ufo::ProfileVariableNames::diagflags_hydro_t,
            ufo::ProfileVariableNames::diagflags_hydro_z,
            ufo::ProfileVariableNames::diagflags_data_correct_z});
      } else if (check == "UInterp" || check == "UInterpAlternative") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumSamePErrObs,
            ufo::ProfileVariableNames::counter_NumInterpErrObs,
            ufo::ProfileVariableNames::NumStd,
            ufo::ProfileVariableNames::NumSig,
            ufo::ProfileVariableNames::StdLev,
            ufo::ProfileVariableNames::SigBelow,
            ufo::ProfileVariableNames::SigAbove,
            ufo::ProfileVariableNames::LevErrors});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::uInterp,
            ufo::ProfileVariableNames::vInterp,
            ufo::ProfileVariableNames::LogP});
      } else if (check == "RH") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_TotCProfs,
            ufo::ProfileVariableNames::counter_TotHProfs,
            ufo::ProfileVariableNames::counter_TotCFlags,
            ufo::ProfileVariableNames::counter_TotHFlags,
            ufo::ProfileVariableNames::counter_TotLFlags,
            ufo::ProfileVariableNames::FlagH,
            ufo::ProfileVariableNames::Indx});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::Press,
            ufo::ProfileVariableNames::Temp,
            ufo::ProfileVariableNames::rh,
            ufo::ProfileVariableNames::td,
            ufo::ProfileVariableNames::tbk,
            ufo::ProfileVariableNames::rhbk});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_interpolation_rh,
            ufo::ProfileVariableNames::diagflags_final_reject_rh});
      } else if (check == "Time") {
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_perm_reject_u,
            ufo::ProfileVariableNames::diagflags_perm_reject_v});
      } else if (check == "PermanentReject") {
        // todo : need to fix the check code first
        valuesToCompare_bool_.insert({});
      } else if (check.find("Background") != std::string::npos) {
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::pge_air_temperature,
            ufo::ProfileVariableNames::pge_relative_humidity,
            ufo::ProfileVariableNames::pge_eastward_wind,
            ufo::ProfileVariableNames::pge_northward_wind,
            ufo::ProfileVariableNames::pge_geopotential_height});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_back_perf_t,
            ufo::ProfileVariableNames::diagflags_final_reject_t,
            ufo::ProfileVariableNames::diagflags_final_reject_t,
            ufo::ProfileVariableNames::diagflags_back_perf_rh,
            ufo::ProfileVariableNames::diagflags_final_reject_rh,
            ufo::ProfileVariableNames::diagflags_final_reject_rh,
            ufo::ProfileVariableNames::diagflags_back_perf_u,
            ufo::ProfileVariableNames::diagflags_final_reject_u,
            ufo::ProfileVariableNames::diagflags_final_reject_u,
            ufo::ProfileVariableNames::diagflags_back_perf_v,
            ufo::ProfileVariableNames::diagflags_final_reject_v,
            ufo::ProfileVariableNames::diagflags_final_reject_v,
            ufo::ProfileVariableNames::diagflags_back_perf_z,
            ufo::ProfileVariableNames::diagflags_final_reject_z,
            ufo::ProfileVariableNames::diagflags_final_reject_z});
      } else if (check == "Pressure") {
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::obs_air_pressure});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_no_pressure_report
          });
      } else if (check == "AveragePressure") {
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::LogP_derived,
            ufo::ProfileVariableNames::bigPgaps_derived,
            ufo::ProfileVariableNames::modellevels_logP_derived,
            ufo::ProfileVariableNames::modellevels_ExnerP_derived,
            ufo::ProfileVariableNames::modellevels_logP_rho_derived,
            ufo::ProfileVariableNames::modellevels_ExnerP_rho_derived});
      } else if (check == "AverageTemperature") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumGapsT});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::modellevels_air_temperature_derived,
            ufo::ProfileVariableNames::air_temperature_derived});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_final_reject_t,
            ufo::ProfileVariableNames::diagflags_final_reject_rh});
      } else if (check == "AverageWindSpeed") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumGapsU,
            ufo::ProfileVariableNames::counter_NumGapsUWP});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::eastward_wind_derived,
            ufo::ProfileVariableNames::northward_wind_derived});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_final_reject_u,
            ufo::ProfileVariableNames::diagflags_final_reject_v          });
      } else if (check == "AverageRelativeHumidity") {
        valuesToCompare_int_.insert({
            ufo::ProfileVariableNames::counter_NumGapsRH});
        valuesToCompare_float_.insert({
            ufo::ProfileVariableNames::relative_humidity_derived});
        valuesToCompare_bool_.insert({
            ufo::ProfileVariableNames::diagflags_final_reject_rh});
      }
    }
  }

  /// Comparison of single values
  template <typename T>
  void ProfileCheckValidator::compareOutput(const std::string &desc,
                                            const T val1,
                                            const T val2,
                                            const int offset,
                                            const float tol,
                                            int &n)
  {
    if (!differenceWithinTol(val1, val2 + offset, tol)) {
        oops::Log::debug() << "   Mismatch for " << desc << " (OPS, this code): "
                           << val1 << ", " << val2 + offset << std::endl;
        n++;
    }
  }

  /// Comparison of vectors of values
  template <typename T>
  void ProfileCheckValidator::compareOutput(const std::string &desc,
                                            const std::vector <T> &vec1,
                                            const std::vector <T> &vec2,
                                            const int offset,
                                            const float tol,
                                            int &n)
  {
    // Do not compare vectors if at least one is empty
    if (oops::anyVectorEmpty(vec1, vec2))
      {
        if (vec1.empty())
          oops::Log::debug() << "Vector of " << desc << " in OPS output is empty" << std::endl;
        if (vec2.empty())
          oops::Log::debug() << "Vector of " << desc << " in this code is empty" << std::endl;
        return;
      }
    // Compare vector elements up to the smaller of the two sizes.
    const size_t vecsize = std::min(vec1.size(), vec2.size());
    for (size_t jvec = 0; jvec < vecsize; ++jvec) {
      if (!differenceWithinTol(vec1[jvec], vec2[jvec] + offset, tol)) {
        oops::Log::debug() << "   Mismatch for " << desc << "[" << jvec << "] "
                           << "(OPS, this code): " << vec1[jvec] << ", "
                           << vec2[jvec] + offset << std::endl;
        n++;
      }
    }
  }

  void ProfileCheckValidator::validate(ProfileDataHandler &profileDataHandler,
                                       size_t commSize)
  {
    // Reset number of mismatches for this profile
    nMismatches_ = 0;

    float tol = options_.Comparison_Tol.value();  // Comparison tolerance

    // Compare integer values obtained in this code and OPS
    for (const auto& valueToCompare_int : valuesToCompare_int_) {
      oops::Log::debug() << "  " << valueToCompare_int << std::endl;

      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(valueToCompare_int, varname, groupname);
      std::string varname_OPS = groupname + std::string("/OPS_") + varname;
      if (groupname == "Counters") {
        /// Special case: OPS counters have one value per profile level,
        /// and are in the MetaData rather than the Counters group.
        /// This avoids the (default) treatment which assumes
        /// that variables in the Counters group have one value per profile.
        varname_OPS = "MetaData/OPS_" + varname;
      }

      // Obtain values for comparison
      const std::vector <int> &values_thiscode =
        profileDataHandler.get<int>(valueToCompare_int);
      const std::vector <int> &values_OPS =
        profileDataHandler.get<int>(varname_OPS);

      // Account for potential offset between values in this code and OPS
      int offset = 0;

      // Offsets due to C++ and Fortran array indices
      auto comparison_offsets_it = comparison_offsets_.find(valueToCompare_int);
      if (comparison_offsets_it != comparison_offsets_.end())
        offset = comparison_offsets_it->second;

      // Offsets due to particular counters being accumulated over profiles in OPS
      // (and not in this code). NumAnyErrors is not included.
      if (groupname == "Counters" && varname != "NumAnyErrors")
        offset = cumulativeCounters_[valueToCompare_int];

      // Only the first element of each counter is compared;
      // in all other cases tne entire vectors are compared.
      if (groupname == "Counters") {
        // The counter comparison is only performed if there is one processor.
        // With some refactoring it would be possible to use eckit::allGatherv
        // to sum the results over multiple processors but, at present,
        // if two processors have a different number of profiles then the allGatherv
        // routine on the processor with more profiles will hang indefinitely.
        if (commSize == 1) {
          if (!oops::anyVectorEmpty(values_OPS, values_thiscode))
            compareOutput(valueToCompare_int, values_OPS[0], values_thiscode[0],
                          offset, tol, nMismatches_);
        } else {
            oops::Log::debug() << "The counter comparison was not performed "
                               << "because multiple processors are in use." << std::endl;
        }
      } else {
        compareOutput(valueToCompare_int, values_OPS, values_thiscode,
                      offset, tol, nMismatches_);
      }

      // Increment cumulative counters. NumAnyErrors is not included.
      if (groupname == "Counters" && varname != "NumAnyErrors")
        cumulativeCounters_[valueToCompare_int] += values_thiscode[0];
    }

    // Compare float values obtained in this code and OPS
    for (const auto& valueToCompare_float : valuesToCompare_float_) {
      oops::Log::debug() << "  " << valueToCompare_float << std::endl;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(valueToCompare_float, varname, groupname);
      const std::string varname_OPS = groupname + std::string("/OPS_") + varname;
      const std::vector <float> &values_thiscode =
        profileDataHandler.get<float>(valueToCompare_float);
      const std::vector <float> &values_OPS =
        profileDataHandler.get<float>(varname_OPS);
      compareOutput(valueToCompare_float, values_OPS, values_thiscode,
                    0, tol, nMismatches_);
    }

    // Compare bool values obtained in this code and OPS
    for (const auto& valueToCompare_bool : valuesToCompare_bool_) {
      oops::Log::debug() << "  " << valueToCompare_bool << std::endl;
      std::string varname;
      std::string groupname;
      ufo::splitVarGroup(valueToCompare_bool, varname, groupname);
      const std::string varname_OPS = groupname + std::string("/OPS_") + varname;
      const std::vector <bool> &values_thiscode =
        profileDataHandler.get<bool>(valueToCompare_bool);
      const std::vector <bool> &values_OPS =
        profileDataHandler.get<bool>(varname_OPS);
      compareOutput(valueToCompare_bool, values_OPS, values_thiscode,
                    0, tol, nMismatches_);
    }

    oops::Log::debug() << " ... all comparisons done ("
                       << nMismatches_ << " mismatches)" << std::endl;
  }
}  // namespace ufo


