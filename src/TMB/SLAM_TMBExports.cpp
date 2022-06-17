// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_SLAM_TMBExports
#include <TMB.hpp>
#include "optF.hpp"
#include "SLAM.hpp"
#include "include/functions.h"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "optF") {
    return optF(this);
  } else if(model == "SLAM") {
    return SLAM(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
