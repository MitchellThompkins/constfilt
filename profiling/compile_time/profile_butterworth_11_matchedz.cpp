#include <constfilt/constfilt.hpp>
static constexpr constfilt::Butterworth<double, 11u, constfilt::MatchedZ> kFilter(
    100.0, 1000.0);
