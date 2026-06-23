#include <constfilt/constfilt.hpp>
static constexpr constfilt::Butterworth<double, 9u, constfilt::TustinPW>
    kFilter(100.0, 1000.0);
