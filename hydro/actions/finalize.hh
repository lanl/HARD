#ifndef HARD_HYDRO_FINALIZE_HH
#define HARD_HYDRO_FINALIZE_HH

#include "state.hh"

namespace hard {

template<std::size_t D>
void finalize(control_policy<state, D> &);

} // namespace hard

#endif // HARD_HYDRO_FINALIZE_HH
