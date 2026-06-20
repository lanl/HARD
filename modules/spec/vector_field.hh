#ifndef SPEC_VECTOR_FIELD_HH
#define SPEC_VECTOR_FIELD_HH

#include <array>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "mesh.hh"
#include "types.hh"

namespace spec {

// ════════════════════════════════════════════════════════════════════
// SoA MDColex Wrapper - Enables vec-like access to SoA data
// ════════════════════════════════════════════════════════════════════

template<typename MDColexType, std::size_t D>
class soa_mdcolex_accessor {
private:
  std::array<MDColexType, D> components_;

public:
  // ────────────────────────────────────────────────────────────────
  // Element proxy - Returned by operator()(i, j, k)
  // ────────────────────────────────────────────────────────────────
  template<typename... Indices>
  class vector_reference {
  private:
    std::array<MDColexType, D> & components_;
    std::tuple<Indices...> indices_;

    // Helper to retrieve the N-th component at the stored indices
    template<std::size_t ComponentIdx>
    FLECSI_INLINE_TARGET auto & get_component() const {
      return std::apply(
        [this](auto... idx) -> auto & {
          return this->components_[ComponentIdx](idx...);
        },
        indices_);
    }

  public:
    FLECSI_INLINE_TARGET explicit vector_reference(
      std::array<MDColexType, D> & comps,
      Indices... idx)
      : components_(comps), indices_(std::make_tuple(idx...)) {}

    FLECSI_INLINE_TARGET auto & x() const {
      static_assert(D >= 1, "x() requires dimension >= 1");
      return get_component<0>();
    }

    FLECSI_INLINE_TARGET auto & y() const {
      static_assert(D >= 2, "y() requires dimension >= 2");
      return get_component<1>();
    }

    FLECSI_INLINE_TARGET auto & z() const {
      static_assert(D >= 3, "z() requires dimension >= 3");
      return get_component<2>();
    }

    FLECSI_INLINE_TARGET operator vec<D>() const {
      if constexpr(D == 1) {
        return vec<1>(x());
      }
      else if constexpr(D == 2) {
        return vec<2>(x(), y());
      }
      else if constexpr(D == 3) {
        return vec<3>(x(), y(), z());
      }
    }

    FLECSI_INLINE_TARGET auto & operator[](std::size_t component_idx) const {
      return std::apply(
        [this, component_idx](auto... idx) -> auto & {
          return this->components_[component_idx](idx...);
        },
        indices_);
    }

    FLECSI_INLINE_TARGET double norm_squared() const {
      double result = 0;
      for(int i = 0; i < D; i++)
        result += (*this)[i] * (*this)[i];
      return result;
    }

    FLECSI_INLINE_TARGET double norm() const {
      return std::sqrt(norm_squared());
    }

    FLECSI_INLINE_TARGET auto abs() const {
      vec<D> result;
      for(int i = 0; i < D; i++)
        result[i] = std::abs((*this)[i]);
      return result;
    }

    FLECSI_INLINE_TARGET vector_reference & operator=(const vec<D> & v) {
      for(int i = 0; i < D; i++)
        (*this)[i] = v[i];
      return *this;
    }

    FLECSI_INLINE_TARGET vector_reference & operator=(
      const vector_reference & other) {
      for(int i = 0; i < D; i++)
        (*this)[i] = other[i];
      return *this;
    }

    FLECSI_INLINE_TARGET vector_reference & operator+=(const vec<D> & v) {
      for(int i = 0; i < D; i++)
        (*this)[i] += v[i];
      return *this;
    }

    FLECSI_INLINE_TARGET vector_reference & operator+=(
      const vector_reference & other) {
      for(int i = 0; i < D; i++)
        (*this)[i] += other[i];
      return *this;
    }

    template<typename ScalarType>
    FLECSI_INLINE_TARGET friend vec<D> operator*(const ScalarType & scalar,
      const vector_reference & proxy) {
      vec<D> result;
      for(int i = 0; i < D; i++)
        result[i] = scalar * proxy[i];
      return result;
    }

    FLECSI_INLINE_TARGET friend vec<D> operator+(const vector_reference & a,
      const vector_reference & b) {
      vec<D> result;
      for(int i = 0; i < D; i++)
        result[i] = a[i] + b[i];
      return result;
    }

    FLECSI_INLINE_TARGET friend vec<D> operator-(const vector_reference & a,
      const vector_reference & b) {
      vec<D> result;
      for(int i = 0; i < D; i++)
        result[i] = a[i] - b[i];
      return result;
    }

    template<typename ScalarType>
    FLECSI_INLINE_TARGET friend vec<D> operator/(const vector_reference & proxy,
      const ScalarType & scalar) {
      vec<D> result;
      for(int i = 0; i < D; i++)
        result[i] = proxy[i] / scalar;
      return result;
    }
  };

  explicit soa_mdcolex_accessor(std::array<MDColexType, D> && comps)
    : components_(std::move(comps)) {}

  // Operator() - Returns a proxy for the element at the given indices.
  template<typename... Indices>
  FLECSI_INLINE_TARGET vector_reference<Indices...> operator()(
    Indices... idx) const {
    return vector_reference<Indices...>(
      const_cast<std::array<MDColexType, D> &>(components_), idx...);
  }

  FLECSI_INLINE_TARGET MDColexType & operator[](std::size_t d) {
    return components_[d];
  }

  FLECSI_INLINE_TARGET const MDColexType & operator[](std::size_t d) const {
    return components_[d];
  }
};

// ════════════════════════════════════════════════════════════════════
// SoA Structure
// ════════════════════════════════════════════════════════════════════

template<typename TYPE, std::size_t DIM, is::index_space INDEX_SPACE>
class vector_field {
public:
  // ────────────────────────────────────────────────────────────────
  // Type aliases
  // ────────────────────────────────────────────────────────────────
  using value_type = TYPE;
  static constexpr is::index_space index_space = INDEX_SPACE;

  using mesh_type = mesh<DIM>;
  using topology_type = typename mesh<DIM>::topology;

  using field_type =
    typename flecsi::field<TYPE>::template definition<mesh_type, INDEX_SPACE>;
  using storage_type = std::array<field_type, DIM>;

  using field_reference =
    typename flecsi::field<TYPE>::template Reference<mesh_type, INDEX_SPACE>;

  template<flecsi::partition_privilege_t privilege,
    flecsi::partition_privilege_t ghost = flecsi::na>
  using field_accessor =
    typename flecsi::field<TYPE>::template accessor<privilege, ghost>;
  // Accessor for tasks side
  template<flecsi::partition_privilege_t privilege,
    flecsi::partition_privilege_t ghost = flecsi::na>
  using accessor = std::vector<field_accessor<privilege, ghost>>;
  // Accessor for the action side
  using accessor_proxy_type = std::vector<field_reference>;

  // Multi accessor for tasks side
  template<flecsi::partition_privilege_t privilege,
    flecsi::partition_privilege_t ghost = flecsi::na>
  using multi_accessor = std::vector<multi<field_accessor<privilege, ghost>>>;

  // Multi accessor for the action side
  using multi_accessor_proxy_type = std::vector<flecsi::data::
      multi_reference<TYPE, flecsi::data::dense, mesh_type, INDEX_SPACE>>;

private:
  storage_type data_;

  // Helper to call operator() on all components
  template<std::size_t... IS>
  accessor_proxy_type get_accessor_impl(topology_type & mesh,
    std::index_sequence<IS...>) const {
    return accessor_proxy_type{data_[IS](mesh)...};
  }

  template<std::size_t... IS>
  multi_accessor_proxy_type get_multi_accessor_impl(
    flecsi::data::launch::mapping<mesh<DIM>> & mesh,
    std::index_sequence<IS...>) const {
    return multi_accessor_proxy_type{data_[IS](mesh)...};
  }

public:
  vector_field() = default;

  explicit vector_field(storage_type data) : data_(std::move(data)) {}

  // Constructor from variadic field definitions
  template<typename... FieldDefs,
    typename = std::enable_if_t<sizeof...(FieldDefs) == DIM>>
  explicit vector_field(FieldDefs &&... data)
    : data_{std::forward<FieldDefs>(data)...} {}

  field_type & operator[](std::size_t dimension) {
    return data_[dimension];
  }

  const field_type & operator[](std::size_t dimension) const {
    return data_[dimension];
  }

  // Call operator() on all field components at once
  accessor_proxy_type operator()(topology_type & mesh) const {
    return get_accessor_impl(mesh, std::make_index_sequence<DIM>{});
  }

  multi_accessor_proxy_type operator()(
    flecsi::data::launch::mapping<mesh<DIM>> & mesh) const {
    return get_multi_accessor_impl(mesh, std::make_index_sequence<DIM>{});
  }

  storage_type & data() noexcept {
    return data_;
  }

  const storage_type & data() const noexcept {
    return data_;
  }

  field_type & x() noexcept {
    static_assert(DIM >= 1, "x() requires at least 1 dimension");
    return data_[0];
  }

  const field_type & x() const noexcept {
    static_assert(DIM >= 1, "x() requires at least 1 dimension");
    return data_[0];
  }

  field_type & y() noexcept {
    static_assert(DIM >= 2, "y() requires at least 2 dimensions");
    return data_[1];
  }

  const field_type & y() const noexcept {
    static_assert(DIM >= 2, "y() requires at least 2 dimensions");
    return data_[1];
  }

  field_type & z() noexcept {
    static_assert(DIM >= 3, "z() requires at least 3 dimensions");
    return data_[2];
  }

  const field_type & z() const noexcept {
    static_assert(DIM >= 3, "z() requires at least 3 dimensions");
    return data_[2];
  }

  auto begin() noexcept {
    return data_.begin();
  }
  auto end() noexcept {
    return data_.end();
  }
  auto begin() const noexcept {
    return data_.begin();
  }
  auto end() const noexcept {
    return data_.end();
  }
  auto cbegin() const noexcept {
    return data_.cbegin();
  }
  auto cend() const noexcept {
    return data_.cend();
  }

private:
  // Helper to build the mdcolex array
  template<class MeshAcc, class VectorOfRefs, std::size_t... Is>
  static auto make_mdcolex_array_impl(MeshAcc && m,
    VectorOfRefs & vec_refs,
    std::index_sequence<Is...>) {
    using mdcolex_type = decltype(m.template mdcolex<INDEX_SPACE>(vec_refs[0]));

    return std::array<mdcolex_type, DIM>{
      m.template mdcolex<INDEX_SPACE>(vec_refs[Is])...};
  }

public:
  // MDColex - Static function to create the accessor wrapper
  template<class MeshAcc, class VectorOfRefs>
  static auto mdcolex(MeshAcc && m, VectorOfRefs & vec_refs) {
    using mdcolex_type = decltype(m.template mdcolex<INDEX_SPACE>(vec_refs[0]));

    auto mdcolex_array = make_mdcolex_array_impl(
      std::forward<MeshAcc>(m), vec_refs, std::make_index_sequence<DIM>{});

    return soa_mdcolex_accessor<mdcolex_type, DIM>(std::move(mdcolex_array));
  }

  static constexpr std::size_t size() noexcept {
    return DIM;
  }

  template<typename VectorType>
  static auto to_accelerator(const VectorType & acc) {
    std::array<typename VectorType::value_type, DIM> result;
    std::copy(acc.begin(), acc.begin() + DIM, result.begin());
    return result;
  }
};

// ════════════════════════════════════════════════════════════════════
// Alias
// ════════════════════════════════════════════════════════════════════

template<std::size_t D>
using vector_field_double = spec::vector_field<double, D, spec::is::cells>;

} // namespace spec

namespace hard {
using spec::vector_field_double;
} // namespace hard

#endif // SPEC_VECTOR_FIELD_HH