/*
 * (C) Copyright 2025- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "./RunPreProcess.h"
#include "oops/runs/Run.h"
#include "ufo/instantiateObsFilterFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory();
  ufo::RunPreProcess preprocess;
  return run.execute(preprocess);
}
