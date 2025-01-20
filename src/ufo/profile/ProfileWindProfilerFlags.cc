/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileWindProfilerFlags.h"
#include "ufo/profile/ProfileVariableNames.h"

namespace ufo {

  static ProfileCheckMaker<ProfileWindProfilerFlags> makerProfileWindProfilerFlags_("WinProFlags");

  ProfileWindProfilerFlags::ProfileWindProfilerFlags
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileWindProfilerFlags::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Wind profiler flag check" << std::endl;

    const int numProfileLevels = profileDataHandler.getNumProfileLevels();

    std::vector <bool> &diagFlagsUFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_u);
    std::vector <bool> &diagFlagsVFinalReject =
      profileDataHandler.get<bool>(ufo::ProfileVariableNames::diagflags_final_reject_v);
    const std::vector <int> &WinProQualInf =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::qualinf_wind_profiler);
    const std::vector <int> &ObsType =
      profileDataHandler.get<int>(ufo::ProfileVariableNames::ObsType);

    if (!oops::allVectorsSameNonZeroSize(diagFlagsUFinalReject, diagFlagsVFinalReject,
                                         WinProQualInf, ObsType)) {
      oops::Log::debug() << "At least one vector is the wrong size. "
                         << "Check will not be performed." << std::endl;
      oops::Log::debug() << "Vector sizes: "
                         << oops::listOfVectorSizes(diagFlagsUFinalReject, diagFlagsVFinalReject,
                                                    WinProQualInf, ObsType)
                         << std::endl;
      return;
    }

    if (ObsType[0] == ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf) {
      for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
        if (WinProQualInf[jlev] > 0) {
          diagFlagsUFinalReject[jlev] = true;
          diagFlagsVFinalReject[jlev] = true;
        }
      }
    }
  }
}  // namespace ufo
