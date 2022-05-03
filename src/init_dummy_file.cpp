// Dummy file required so that useDynLib(SLAM, .registration=TRUE) doesn't fail on empty 'src'
// Delete 'src/init_dummy_file.cpp' if C++ files are later added to 'src/'.

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

void attribute_visible R_init_SLAM(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


