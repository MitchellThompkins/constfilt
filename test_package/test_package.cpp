#include <constfilt/constfilt.hpp>

int main()
{
    static constexpr auto filt =
        constfilt::Butterworth<double, 2>(100.0, 1000.0);
    (void)filt;
    return 0;
}
