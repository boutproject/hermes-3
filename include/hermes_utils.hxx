#ifndef HERMES_UTILS_H
#define HERMES_UTILS_H

#include <string>

#include <field.hxx>

template <
    typename T, typename Function,
    typename = decltype(std::declval<Function&>()(std::declval<typename T::ind_type&>()))>
inline T filledFrom(const T& f, Function func, std::string region_string = "RGN_ALL") {
  static_assert(bout::utils::is_Field<T>::value, "filledFrom only works on Fields");
  T result{emptyFrom(f)};
  BOUT_FOR(i, result.getRegion(region_string)) {
    result[i] = func(i);
  }
  return result;
}

#endif // HERMES_UTILS_H
