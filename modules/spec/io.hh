#ifndef SPEC_IO_HH
#define SPEC_IO_HH

#include <sstream>

namespace spec::io {

struct name {
  name(std::string const & base) {
    ss_ << base;
  }
  template<typename T>
  name & operator<<(T const & value) {
    ss_ << value;
    return *this;
  }

  name & operator<<(
    ::std::ostream & (*basic_manipulator)(::std::ostream & stream)) {
    ss_ << basic_manipulator;
    return *this;
  } // operator<<

  std::string const str() const {
    return ss_.str();
  }

private:
  std::stringstream ss_;
}; // struct name

} // namespace spec::io

#endif // SPEC_IO_HH
