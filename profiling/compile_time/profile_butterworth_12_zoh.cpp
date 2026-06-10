#include <constfilt/constfilt.hpp>
static constexpr constfilt::Butterworth<double, 12u, constfilt::ZOH> kFilter(
    100.0, 1000.0);
