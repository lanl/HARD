#include "flecsi/util/unit.hh"

#include "../../../modules/spec/boundaries.hh"

using namespace flecsi;
using namespace spec;

// Minimal state: two scalars, one vector, and one field that is not tagged
// boundary.
template<std::size_t D>
struct state {
  static inline const tfield<double, mesh<D>> rho{"mass_density",
    tag::boundary};
  static inline const tfield<double, mesh<D>> c{"sound_speed", tag::boundary};
  static inline const tfield<double, mesh<D>> p{"pressure", tag::boundary};
  static inline const tfield<vec<D>, mesh<D>> u{"velocity", tag::boundary};
  static inline const tfield<double, mesh<D>> T{"temperature", tag::none};
  static auto fields() {
    return std::tie(rho, c, p, u, T);
  }
};

template<std::size_t D>
boundaries<D>
parse(std::string const & py_boundaries) {
  py::gil_scoped_acquire _;
  py::object cfg = py::eval("dict(boundaries=" + py_boundaries + ")");
  config_py config(cfg);
  return parse_boundaries<D>(config, state<D>::fields());
}

int
test_driver(scheduler &) {
  UNIT() {
    py::scoped_interpreter guard{};

    // Plain side defaults: every tagged field gets them.
    {
      auto b = parse<1>(R"({"xlow": "reflecting", "xhigh": "flow"})");
      EXPECT_EQ(b.defaults[ax::x][bd::low].type, bd::reflecting);
      EXPECT_EQ(b.defaults[ax::x][bd::high].type, bd::flow);
      EXPECT_FALSE(b.periodic[ax::x]);
      EXPECT_EQ(b.fields.size(), std::size_t{4}); // temperature not tagged
      EXPECT_EQ(b.spec("pressure")[ax::x][bd::low].type, bd::reflecting);
      EXPECT_EQ(b.spec("velocity")[ax::x][bd::high].type, bd::flow);
    }

    // Per-field overrides with dirichlet values; 3-entry list in 1D.
    {
      auto b = parse<1>(R"({
          "xlow": ["flow",
                   ("sound_speed", "dirichlet"),
                   ("pressure", "dirichlet", 2.5),
                   ("velocity", "dirichlet", [1.0, 0.0, 0.0])],
          "xhigh": ["flow", ("velocity", "reflecting")]})");
      EXPECT_EQ(b.defaults[ax::x][bd::low].type, bd::flow);
      EXPECT_EQ(b.spec("mass_density")[ax::x][bd::low].type, bd::flow);
      auto const & c = b.spec("sound_speed")[ax::x][bd::low];
      EXPECT_EQ(c.type, bd::dirichlet);
      EXPECT_EQ(c.value[0], 0.0);
      auto const & p = b.spec("pressure")[ax::x][bd::low];
      EXPECT_EQ(p.type, bd::dirichlet);
      EXPECT_EQ(p.value[0], 2.5);
      auto const & u = b.spec("velocity")[ax::x][bd::low];
      EXPECT_EQ(u.type, bd::dirichlet);
      EXPECT_EQ(u.value[0], 1.0);
      EXPECT_EQ(b.spec("velocity")[ax::x][bd::high].type, bd::reflecting);
      EXPECT_EQ(b.spec("pressure")[ax::x][bd::high].type, bd::flow);
    }

    // Dirichlet as a side default: (type, value), number broadcast, list
    // used component-wise, per-field override still wins.
    {
      auto b = parse<2>(R"({"xlow": "reflecting", "xhigh": ("dirichlet", 0.0),
                              "ylow": [("dirichlet", 4.0),
                                       ("pressure", "reflecting")],
                              "yhigh": ("dirichlet", [1.0, 2.0, 3.0])})");
      EXPECT_EQ(b.defaults[ax::x][bd::high].type, bd::dirichlet);
      EXPECT_EQ(b.spec("velocity")[ax::x][bd::high].value[1], 0.0);
      EXPECT_EQ(b.spec("mass_density")[ax::y][bd::low].value[0], 4.0);
      EXPECT_EQ(b.spec("pressure")[ax::y][bd::low].type, bd::reflecting);
      EXPECT_EQ(b.spec("velocity")[ax::y][bd::high].value[1], 2.0);
      EXPECT_EQ(b.spec("mass_density")[ax::y][bd::high].value[0], 1.0);
    }

    // Periodic axes are flagged per axis and reach every field.
    {
      auto b = parse<2>(R"({"xlow": "flow", "xhigh": "flow",
                              "ylow": "periodic", "yhigh": "periodic"})");
      EXPECT_FALSE(b.periodic[ax::x]);
      EXPECT_TRUE(b.periodic[ax::y]);
      EXPECT_EQ(b.spec("velocity")[ax::y][bd::low].type, bd::periodic);
    }
    {
      auto b = parse<3>(R"({"xlow": "flow", "xhigh": "flow",
                              "ylow": "flow", "yhigh": "flow",
                              "zlow": "periodic", "zhigh": "periodic"})");
      EXPECT_FALSE(b.periodic[ax::x]);
      EXPECT_FALSE(b.periodic[ax::y]);
      EXPECT_TRUE(b.periodic[ax::z]);
    }
    {
      auto b = parse<3>(R"({"xlow": "flow", "xhigh": "flow",
                              "ylow": "periodic", "yhigh": "periodic",
                              "zlow": "flow", "zhigh": "flow"})");
      EXPECT_TRUE(b.periodic[ax::y]);
      EXPECT_FALSE(b.periodic[ax::z]);
    }
  };
}

util::unit::driver<test_driver> driver;
