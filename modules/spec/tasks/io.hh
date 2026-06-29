#ifndef HARD_COMMON_TASKS_IO_HH
#define HARD_COMMON_TASKS_IO_HH

// Main I/O dispatcher header
// Includes appropriate I/O backend implementations based on configuration

#include "io/csv.hh" // CSV output (always available)
#include "io/vti.hh" // VTK ImageData XML output (always available)

#ifdef HARD_ENABLE_HDF5
#include "io/xdmf.hh" // XDMF+HDF5 output (requires parallel HDF5)
#endif

#endif // HARD_COMMON_TASKS_IO_HH
