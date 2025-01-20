/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/profile/ProfileVariableNames.h"

// Observation values

constexpr const char* const ufo::ProfileVariableNames::obs_air_pressure;
constexpr const char* const ufo::ProfileVariableNames::obs_air_temperature;
constexpr const char* const ufo::ProfileVariableNames::obs_relative_humidity;
constexpr const char* const ufo::ProfileVariableNames::obs_eastward_wind;
constexpr const char* const ufo::ProfileVariableNames::obs_northward_wind;
constexpr const char* const ufo::ProfileVariableNames::obs_geopotential_height;
constexpr const char* const ufo::ProfileVariableNames::obs_dew_point_temperature;

// Derived observation values

constexpr const char* const ufo::ProfileVariableNames::obs_derived_air_pressure;

// Observation errors

constexpr const char* const ufo::ProfileVariableNames::obserr_air_temperature;
constexpr const char* const ufo::ProfileVariableNames::obserr_relative_humidity;
constexpr const char* const ufo::ProfileVariableNames::obserr_eastward_wind;
constexpr const char* const ufo::ProfileVariableNames::obserr_northward_wind;
constexpr const char* const ufo::ProfileVariableNames::obserr_geopotential_height;
constexpr const char* const ufo::ProfileVariableNames::obserr_dew_point_temperature;

// HofX

constexpr const char* const ufo::ProfileVariableNames::hofx_air_temperature;
constexpr const char* const ufo::ProfileVariableNames::hofx_geopotential_height;
constexpr const char* const ufo::ProfileVariableNames::hofx_relative_humidity;
constexpr const char* const ufo::ProfileVariableNames::hofx_eastward_wind;
constexpr const char* const ufo::ProfileVariableNames::hofx_northward_wind;
constexpr const char* const ufo::ProfileVariableNames::hofx_dew_point_temperature;

// Probability of gross error

constexpr const char* const ufo::ProfileVariableNames::pge_air_temperature;
constexpr const char* const ufo::ProfileVariableNames::pge_relative_humidity;
constexpr const char* const ufo::ProfileVariableNames::pge_eastward_wind;
constexpr const char* const ufo::ProfileVariableNames::pge_northward_wind;
constexpr const char* const ufo::ProfileVariableNames::pge_geopotential_height;

// MetaData

constexpr const char* const ufo::ProfileVariableNames::station_ID;
constexpr const char* const ufo::ProfileVariableNames::ObsType;
constexpr const char* const ufo::ProfileVariableNames::Latitude;
constexpr const char* const ufo::ProfileVariableNames::Longitude;
constexpr const char* const ufo::ProfileVariableNames::Zstation;
constexpr const char* const ufo::ProfileVariableNames::LevelType;
constexpr const char* const ufo::ProfileVariableNames::InstrType;
constexpr const char* const ufo::ProfileVariableNames::extended_obs_space;

// Diagnostic flags for airTemperature

constexpr const char* const ufo::ProfileVariableNames::diagflags_back_perf_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_back_reject_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_data_correct_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_final_reject_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_hydro_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_interpolation_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_partial_layer_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_perm_reject_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_standard_level_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_sig_temp_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_superadiabat_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_surface_level_t;
constexpr const char* const ufo::ProfileVariableNames::diagflags_tropo_t;

// Diagnostic flags for eastwardWind

constexpr const char* const ufo::ProfileVariableNames::diagflags_back_perf_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_back_reject_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_final_reject_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_interpolation_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_max_wind_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_partial_layer_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_perm_reject_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_sig_wind_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_standard_level_u;
constexpr const char* const ufo::ProfileVariableNames::diagflags_surface_level_u;

// Diagnostic flags for northwardWind

constexpr const char* const ufo::ProfileVariableNames::diagflags_back_perf_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_back_reject_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_final_reject_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_interpolation_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_max_wind_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_partial_layer_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_perm_reject_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_sig_wind_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_standard_level_v;
constexpr const char* const ufo::ProfileVariableNames::diagflags_surface_level_v;

// Diagnostic flags for relativeHumidity

constexpr const char* const ufo::ProfileVariableNames::diagflags_back_perf_rh;
constexpr const char* const ufo::ProfileVariableNames::diagflags_back_reject_rh;
constexpr const char* const ufo::ProfileVariableNames::diagflags_final_reject_rh;
constexpr const char* const ufo::ProfileVariableNames::diagflags_interpolation_rh;
constexpr const char* const ufo::ProfileVariableNames::diagflags_partial_layer_rh;
constexpr const char* const ufo::ProfileVariableNames::diagflags_perm_reject_rh;
constexpr const char* const ufo::ProfileVariableNames::diagflags_standard_level_rh;
constexpr const char* const ufo::ProfileVariableNames::diagflags_surface_level_rh;

// Diagnostic flags for height

constexpr const char* const ufo::ProfileVariableNames::diagflags_back_perf_z;
constexpr const char* const ufo::ProfileVariableNames::diagflags_back_reject_z;
constexpr const char* const ufo::ProfileVariableNames::diagflags_data_correct_z;
constexpr const char* const ufo::ProfileVariableNames::diagflags_final_reject_z;
constexpr const char* const ufo::ProfileVariableNames::diagflags_hydro_z;
constexpr const char* const ufo::ProfileVariableNames::diagflags_interpolation_z;
constexpr const char* const ufo::ProfileVariableNames::diagflags_perm_reject_z;

// Diagnostic flags for observationReport

constexpr const char* const ufo::ProfileVariableNames::diagflags_final_reject_report;
constexpr const char* const ufo::ProfileVariableNames::diagflags_interpolation_report;
constexpr const char* const ufo::ProfileVariableNames::diagflags_no_pressure_report;
constexpr const char* const ufo::ProfileVariableNames::diagflags_out_of_area_report;
constexpr const char* const ufo::ProfileVariableNames::diagflags_perm_reject_report;
constexpr const char* const ufo::ProfileVariableNames::diagflags_surface_level_report;
constexpr const char* const ufo::ProfileVariableNames::diagflags_surplus_report;
constexpr const char* const ufo::ProfileVariableNames::diagflags_track_reject_report;

// Counters

constexpr const char* const ufo::ProfileVariableNames::counter_NumAnyErrors;
constexpr const char* const ufo::ProfileVariableNames::counter_NumSamePErrObs;
constexpr const char* const ufo::ProfileVariableNames::counter_NumSuperadiabat;
constexpr const char* const ufo::ProfileVariableNames::counter_Num925Miss;
constexpr const char* const ufo::ProfileVariableNames::counter_Num100Miss;
constexpr const char* const ufo::ProfileVariableNames::counter_NumStdMiss;
constexpr const char* const ufo::ProfileVariableNames::counter_NumHydErrObs;
constexpr const char* const ufo::ProfileVariableNames::counter_NumIntHydErrors;
constexpr const char* const ufo::ProfileVariableNames::counter_NumInterpErrors;
constexpr const char* const ufo::ProfileVariableNames::counter_NumInterpErrObs;
constexpr const char* const ufo::ProfileVariableNames::counter_NumSignChange;
constexpr const char* const ufo::ProfileVariableNames::counter_TotCProfs;
constexpr const char* const ufo::ProfileVariableNames::counter_TotHProfs;
constexpr const char* const ufo::ProfileVariableNames::counter_TotCFlags;
constexpr const char* const ufo::ProfileVariableNames::counter_TotHFlags;
constexpr const char* const ufo::ProfileVariableNames::counter_TotLFlags;
constexpr const char* const ufo::ProfileVariableNames::counter_NumGapsT;
constexpr const char* const ufo::ProfileVariableNames::counter_NumGapsU;
constexpr const char* const ufo::ProfileVariableNames::counter_NumGapsUWP;
constexpr const char* const ufo::ProfileVariableNames::counter_NumGapsRH;

// Corrections

constexpr const char* const ufo::ProfileVariableNames::obscorrection_air_temperature;
constexpr const char* const ufo::ProfileVariableNames::obscorrection_geopotential_height;

// Intermediate values

constexpr const char* const ufo::ProfileVariableNames::DC;
constexpr const char* const ufo::ProfileVariableNames::ETol;
constexpr const char* const ufo::ProfileVariableNames::D;
constexpr const char* const ufo::ProfileVariableNames::E;
constexpr const char* const ufo::ProfileVariableNames::HydError;
constexpr const char* const ufo::ProfileVariableNames::PBottom;
constexpr const char* const ufo::ProfileVariableNames::StdLev;
constexpr const char* const ufo::ProfileVariableNames::SigAbove;
constexpr const char* const ufo::ProfileVariableNames::SigBelow;
constexpr const char* const ufo::ProfileVariableNames::IndStd;
constexpr const char* const ufo::ProfileVariableNames::LevErrors;
constexpr const char* const ufo::ProfileVariableNames::tInterp;
constexpr const char* const ufo::ProfileVariableNames::uInterp;
constexpr const char* const ufo::ProfileVariableNames::vInterp;
constexpr const char* const ufo::ProfileVariableNames::LogP;
constexpr const char* const ufo::ProfileVariableNames::NumStd;
constexpr const char* const ufo::ProfileVariableNames::NumSig;
constexpr const char* const ufo::ProfileVariableNames::Press;
constexpr const char* const ufo::ProfileVariableNames::Temp;
constexpr const char* const ufo::ProfileVariableNames::rh;
constexpr const char* const ufo::ProfileVariableNames::td;
constexpr const char* const ufo::ProfileVariableNames::tbk;
constexpr const char* const ufo::ProfileVariableNames::rhbk;
constexpr const char* const ufo::ProfileVariableNames::FlagH;
constexpr const char* const ufo::ProfileVariableNames::Indx;

// GeoVaLs

constexpr const char* const ufo::ProfileVariableNames::geovals_orog;
constexpr const char* const ufo::ProfileVariableNames::geovals_pressure;
constexpr const char* const ufo::ProfileVariableNames::geovals_pressure_rho;
constexpr const char* const ufo::ProfileVariableNames::geovals_pressure_rho_minus_one;
constexpr const char* const ufo::ProfileVariableNames::geovals_height;
constexpr const char* const ufo::ProfileVariableNames::geovals_height_rho;
constexpr const char* const ufo::ProfileVariableNames::geovals_height_rho_minus_one;
constexpr const char* const ufo::ProfileVariableNames::geovals_air_potential_temperature;
constexpr const char* const ufo::ProfileVariableNames::geovals_air_temperature;
constexpr const char* const ufo::ProfileVariableNames::geovals_air_pressure_at_surface;
constexpr const char* const ufo::ProfileVariableNames::geovals_relative_humidity;

// GeoVaLs used in validation

constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_logP;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_ExnerP;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_logP_rho;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_ExnerP_rho;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_air_temperature;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_eastward_wind;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_northward_wind;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_relative_humidity;
constexpr const char* const
ufo::ProfileVariableNames::geovals_testreference_air_temperature_qcflags;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_eastward_wind_qcflags;
constexpr const char* const ufo::ProfileVariableNames::geovals_testreference_northward_wind_qcflags;
constexpr const char* const
ufo::ProfileVariableNames::geovals_testreference_relative_humidity_qcflags;

// Averaged values on model levels

constexpr const char* const ufo::ProfileVariableNames::air_temperature_derived;
constexpr const char* const ufo::ProfileVariableNames::eastward_wind_derived;
constexpr const char* const ufo::ProfileVariableNames::northward_wind_derived;
constexpr const char* const ufo::ProfileVariableNames::relative_humidity_derived;

// Derived observation values (used in averaging)

constexpr const char* const ufo::ProfileVariableNames::LogP_derived;
constexpr const char* const ufo::ProfileVariableNames::bigPgaps_derived;

// Derived model values (used in averaging)

constexpr const char* const ufo::ProfileVariableNames::modellevels_logP_derived;
constexpr const char* const ufo::ProfileVariableNames::modellevels_ExnerP_derived;
constexpr const char* const ufo::ProfileVariableNames::modellevels_air_temperature_derived;

// Derived model values on rho levels (used in averaging)

constexpr const char* const ufo::ProfileVariableNames::modellevels_logP_rho_derived;
constexpr const char* const ufo::ProfileVariableNames::modellevels_logPWB_rho_derived;
constexpr const char* const ufo::ProfileVariableNames::modellevels_ExnerP_rho_derived;
