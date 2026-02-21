#ifndef HARD_HYDRO_FINALIZE_HH
#define HARD_HYDRO_FINALIZE_HH

#include "state.hh"

#include <../modules/spec/io.hh>
#include <flecsi/flog.hh>

namespace hard::action {

template<std::size_t D>
void finalize(control_policy<state, D> &);

} // namespace hard::action

#endif // HARD_HYDRO_FINALIZE_HH
