#include <constfilt/constfilt.hpp>
static constexpr constfilt::Butterworth<double, 6u, constfilt::ZOH> kFilter(
    100.0, 1000.0);
