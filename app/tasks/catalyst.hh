#ifndef FLASTRO_TASK_EXTERNAL_HH
#define FLASTRO_TASK_EXTERNAL_HH

#include <flecsi/flog.hh>

#include "../catalyst/adaptor.hh"
#include "../catalyst/types.hh"
#include "../state.hh"
#include "../types.hh"

namespace flastro::tasks::external {

/*--------------------------------------------------------------------------*
  init tasks

  Call the 'init' methods of our catalust_attributes instance and simple_cubic
  from within the 'init' tasks.
 *--------------------------------------------------------------------------*/
template<std::size_t D>
inline void
get_lattice_size(typename mesh<D>::template accessor<ro> m) {
  if constexpr(D == 3) {
    number_vertices[0] = m.template size<ax::x, dm::global>() +
                         1; // 1 more vertex than cell in simple cubic lattice
    number_vertices[1] = m.template size<ax::y, dm::global>() + 1;
    number_vertices[2] = m.template size<ax::z, dm::global>() + 1;
  }
  else {
    /* Do nothing, catalyst/paraview visualization only for 3D */
  }
}

template<std::size_t D>
inline void
get_color_data(typename mesh<D>::template accessor<ro> m) {
  if constexpr(D == 3) {
    // get number of color in each direction
    auto ccolors = m.axis_colors();
    number_colors[0] = ccolors[ax::x];
    number_colors[1] = ccolors[ax::y];
    number_colors[2] = ccolors[ax::z];

    // get coordinate of color block in each direction
    auto ccoords = m.color_indeces();
    color_coordinate[0] = ccoords[ax::x];
    color_coordinate[1] = ccoords[ax::y];
    color_coordinate[2] = ccoords[ax::z];

    // get global id of first and last cell in
    // local color block (in each direction)
    auto cell_it_x = m.template cells<ax::x, dm::quantities>();
    auto cell_it_y = m.template cells<ax::y, dm::quantities>();
    auto cell_it_z = m.template cells<ax::z, dm::quantities>();

    auto ax_x = m.template axis<ax::x>();
    auto ax_y = m.template axis<ax::y>();
    auto ax_z = m.template axis<ax::z>();

    color_cell_id_begin[0] = ax_x.global_id(*cell_it_x.begin());
    color_cell_id_end[0] = ax_x.global_id(*cell_it_x.end() - 1);
    color_cell_id_begin[1] = ax_y.global_id(*cell_it_y.begin());
    color_cell_id_end[1] = ax_y.global_id(*cell_it_y.end() - 1);
    color_cell_id_begin[2] = ax_z.global_id(*cell_it_z.begin());
    color_cell_id_end[2] = ax_z.global_id(*cell_it_z.end() - 1);

    my_rank = flecsi::process();
  }
  else {
    /* Do nothing, catalyst/paraview visualization only for 3D */
  }
}

inline void
init_attributes(single<catalyst_attributes>::accessor<wo> c_a,
  std::size_t pt,
  std::size_t cl) {
  /* Need to use the '->' operator to "dereference" 'single' accessor types. */
  c_a->initialize(pt, cl);
}

/*--------------------------------------------------------------------------*
  update task

  Call the 'update' method of our catalyst_attributes instance from within the
  'update' task.
 *--------------------------------------------------------------------------*/
template<std::size_t D>
inline void
update_attributes(single<catalyst_attributes>::accessor<rw> c_a,
  single<double>::accessor<ro> time_a,
  typename mesh<D>::template accessor<ro> m,
  field<double>::accessor<ro, ro> r_a,
  typename field<vec<D>>::template accessor<ro, ro> u_a,
  field<double>::accessor<ro, ro> p_a,
  field<double>::accessor<ro, ro> rE_a,
  field<double>::accessor<ro, ro>
    Esf_a) { // <<< add variables here for catalyst
  double time = *time_a;
  if constexpr(D == 3) {
    auto r = m.template mdcolex<is::cells>(
      r_a); // density (here: = mass (fixed volume))
    auto u = m.template mdcolex<is::cells>(u_a); // velocity
    auto p = m.template mdcolex<is::cells>(p_a); // pressure
    auto rE = m.template mdcolex<is::cells>(rE_a);
    auto Esf = m.template mdcolex<is::cells>(
      Esf_a); // <<< add variables here for catalyst
    flog(info) << "Update catalyst fields at time: " << time << std::endl;
    flecsi::flog::flush();

    /*
    idea: make a vals vector, manually fill it with values, like in
    tasks::hydro.cc then pass to update method of catalyst_attributes
    */
    std::vector<double> density_vals;
    std::vector<double> pressure_vals;
    std::vector<double> velocity_vals;
    std::vector<double> tot_energy_vals;
    std::vector<double> rad_energy_vals; // <<< add variables here for catalyst

    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          density_vals.push_back(r(i, j, k));
          pressure_vals.push_back(p(i, j, k));
          tot_energy_vals.push_back(rE(i, j, k));
          rad_energy_vals.push_back(
            Esf(i, j, k)); // <<< add variables here for catalyst
          // vector data for catalyst/paraview stored non-interleaved: first all
          // x, then all y, etc.
          velocity_vals.push_back(u(i, j, k).x);
        } // for
      } // for
    } // for
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          velocity_vals.push_back(
            u(i, j, k).y); // run loop again to add y components
        } // for
      } // for
    } // for
    for(auto k : m.template cells<ax::z, dm::quantities>()) {
      for(auto j : m.template cells<ax::y, dm::quantities>()) {
        for(auto i : m.template cells<ax::x, dm::quantities>()) {
          velocity_vals.push_back(
            u(i, j, k).z); // run loop again to add z components
        } // for
      } // for
    } // for

    c_a->update_fields(velocity_vals.data(),
      density_vals.data(),
      pressure_vals.data(),
      tot_energy_vals.data(),
      rad_energy_vals.data()); // <<< add variables here for catalyst
  }
  else {
    /* Do nothing, catalyst/paraview visualization only for 3D */
  }
}

/*--------------------------------------------------------------------------*
  execute task

  Call the 'execute' method of our catalyst_attributes instance from within the
  'execute' task.
 *--------------------------------------------------------------------------*/
inline void
execute_catalyst(single<catalyst_attributes>::accessor<rw> c_a,
  size_t step,
  single<double>::accessor<ro> time,
  simple_cubic & lattice) {
  c_a->execute(step, *time, lattice);
}

/*--------------------------------------------------------------------------*
  initialize catalyst task

  Call the a 'initialize' method of our catalyst adaptor instance from within
  the 'initialize' task.
 *--------------------------------------------------------------------------*/
inline void
initialize() {
  catalyst_adaptor::initialize();
  flog(info) << "Catalyst Initialize, ok" << std::endl;
}
/*--------------------------------------------------------------------------*
  finalize task

  Call the a 'finalize' method of our catalyst_attributes instance from within
  the 'finalize' task.
 *--------------------------------------------------------------------------*/
inline void
finalize(single<catalyst_attributes>::accessor<rw> c_a) {
  catalyst_adaptor::finalize();
  c_a->finalize();
  flog(info) << "Catalyst Finalize, ok" << std::endl;
}

} // namespace flastro::tasks::external

#endif // FLASTRO_TASKS_EXTERNAL_HH
