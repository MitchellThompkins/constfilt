#include <constfilt/constfilt.hpp>
static constexpr constfilt::Butterworth<double, 1u, constfilt::Tustin> kFilter(
    100.0, 1000.0);
