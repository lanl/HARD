#ifndef SPEC_BOUNDARIES_HH
#define SPEC_BOUNDARIES_HH

#include "mesh.hh"
#include "types.hh"
#include "utils.hh"

#include <flecsi/flog.hh>

#include <algorithm>
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace spec {

// Boundary condition on one side for one field.
struct bc_side {
  bd::boundary_type type = bd::flow;
  std::array<double, 3> value{}; // dirichlet value; components for vectors
};

// Boundary conditions of one field on every side: [axis][low/high].
template<std::size_t D>
using bc_spec = std::array<std::array<bc_side, 2>, D>;

// Parsed boundary configuration.
template<std::size_t D>
struct boundaries {
  bc_spec<D> defaults;
  typename mesh<D>::periodic_axes periodic;
  std::map<std::string, bc_spec<D>> fields; // resolved per field

  // Spec of one field: fatal if the field is not a boundary field.
  bc_spec<D> const & spec(std::string const & name) const {
    auto it = fields.find(name);
    if(it == fields.end())
      flog_fatal("boundaries: no boundary conditions for field '"
                 << name << "' (is it tagged boundary?)");
    return it->second;
  }

  // Spec with the same condition on every side (e.g. flow for level sets).
  static bc_spec<D> uniform(bd::boundary_type t) {
    bc_spec<D> sp;
    for(auto & axis : sp)
      for(auto & side : axis)
        side = bc_side{t, {}};
    return sp;
  }
};

// Parse config["boundaries"]. A side is a type name, a (type, value) tuple,
// or a list [default, (field, type), (field, type, value), ...]
template<std::size_t D, typename Fields>
boundaries<D>
parse_boundaries(config_py const & config, Fields const & fields) {
  constexpr const char * axes[] = {"x", "y", "z"};
  constexpr const char * sides[] = {"low", "high"};

  // Fields allowed in the input, with their number of components
  std::map<std::string, std::size_t> comps;
  std::size_t max_comps = 1;
  spec::for_each_field(fields, tag::boundary, [&](auto const & f) {
    comps[f.name] = f.components;
    max_comps = std::max(max_comps, f.components);
  });
  // Keep list of all valid fields
  std::string valid;
  for(auto const & [n, c] : comps)
    valid += " " + n;

  // Side default: "flow" or ("dirichlet", value). A number is broadcast to
  // every component, a list needs an entry per component of the widest field.
  auto parse_default = [&](config_py const & v, std::string const & key) {
    bc_side bs;
    if(v.is<py::str>()) {
      bs.type = utils::mesh_boundary<D>(v.cast<std::string>());
      return bs;
    }
    if(!(v.is<py::tuple>() && v.size() == 2))
      flog_fatal("boundaries: side " << key
                                     << " must be a type name, a (type, value) "
                                        "tuple, or a list starting with one");
    bs.type = utils::mesh_boundary<D>(v[0].cast<std::string>());
    if(bs.type != bd::dirichlet)
      flog_fatal(
        "boundaries: side " << key << ", only dirichlet takes a value");
    auto const val = v[1];
    if(val.is<py::float_>() || val.is<py::int_>()) {
      bs.value.fill(val.cast<double>());
    }
    else {
      auto vv = val.cast<std::vector<double>>();
      if(vv.size() < max_comps)
        flog_fatal("boundaries: default value on side "
                   << key << " needs at least " << max_comps
                   << " components, got " << vv.size());
      std::copy_n(
        vv.begin(), std::min<std::size_t>(vv.size(), 3), bs.value.begin());
    }
    return bs;
  };

  boundaries<D> r;
  struct entry {
    std::size_t axis, sd;
    std::string name;
    bc_side bc;
  };
  std::vector<entry> entries;

  for(std::size_t a = 0; a < D; ++a) {
    for(std::size_t sd = 0; sd < 2; ++sd) {
      const std::string key = std::string(axes[a]) + sides[sd];
      if(!config["boundaries"].contains(key))
        flog_fatal("boundaries: missing side " << key);
      auto side = config["boundaries"][key];

      if(side.is<py::str>() || side.is<py::tuple>()) {
        r.defaults[a][sd] = parse_default(side, key);
      }
      else if(side.is<py::list>()) {
        if(side.size() == 0)
          flog_fatal("boundaries: empty list for side " << key);
        r.defaults[a][sd] = parse_default(side[0], key);

        for(std::size_t i = 1; i < side.size(); ++i) {
          auto e = side[i];
          if(!(e.is<py::tuple>() || e.is<py::list>()) || e.size() < 2 ||
             e.size() > 3)
            flog_fatal("boundaries: entry " << i << " of side " << key
                                            << " must be (field, type) or "
                                               "(field, type, value)");
          entry en{a, sd, e[0].cast<std::string>(), {}};
          en.bc.type = utils::mesh_boundary<D>(e[1].cast<std::string>());

          if(!comps.count(en.name))
            flog_fatal("boundaries: unknown field '"
                       << en.name << "' on side " << key
                       << ". Valid fields:" << valid);
          if(en.bc.type == bd::periodic)
            flog_fatal("boundaries: periodic is only allowed as a side "
                       "default (side "
                       << key << ", field " << en.name << ")");
          if(r.defaults[a][sd].type == bd::periodic)
            flog_fatal("boundaries: side " << key
                                           << " is periodic, it cannot "
                                              "have per-field entries");

          if(en.bc.type == bd::dirichlet) {
            if(e.size() == 3) {
              auto v = e[2];
              const std::size_t n = comps[en.name];
              if(v.is<py::float_>() || v.is<py::int_>()) {
                if(n != 1)
                  flog_fatal("boundaries: field " << en.name << " needs " << n
                                                  << " components on side "
                                                  << key);
                en.bc.value[0] = v.cast<double>();
              }
              else if(v.is<py::list>() || v.is<py::tuple>()) {
                auto vv = v.cast<std::vector<double>>();
                if(vv.size() < n)
                  flog_fatal("boundaries: field "
                             << en.name << " needs at least " << n
                             << " components on side " << key << ", got "
                             << vv.size());
                std::copy_n(vv.begin(), n, en.bc.value.begin());
              }
              else {
                flog_fatal("boundaries: unsupported dirichlet value for field "
                           << en.name << " on side " << key
                           << " (number or list; callables not supported yet)");
              }
            }
          }
          else if(e.size() == 3) {
            flog_fatal("boundaries: a value is only allowed with dirichlet "
                       "(side "
                       << key << ", field " << en.name << ")");
          }
          entries.push_back(std::move(en));
        }
      }
      else {
        flog_fatal("boundaries: side "
                   << key
                   << " must be a type name, a (type, value) tuple, "
                      "or a list");
      }
    }

    const bool lo = r.defaults[a][0].type == bd::periodic;
    const bool hi = r.defaults[a][1].type == bd::periodic;
    if(lo != hi)
      flog_fatal(
        "boundaries: axis " << axes[a] << " must be periodic on both sides");
    r.periodic[a] = lo;
  }

  // Every boundary field starts from the side defaults
  for(auto const & [name, c] : comps) {
    bc_spec<D> sp;
    for(std::size_t a = 0; a < D; ++a)
      for(std::size_t sd = 0; sd < 2; ++sd)
        sp[a][sd] = r.defaults[a][sd];
    r.fields[name] = sp;
  }
  for(auto const & en : entries)
    r.fields[en.name][en.axis][en.sd] = en.bc;

  return r;
}

} // namespace spec

#endif // SPEC_BOUNDARIES_HH
