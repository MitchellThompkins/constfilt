// bench_accuracy_constfilt.cpp: constfilt accuracy vs Octave reference (C++17)
//
// Outputs CSV rows to stdout (no header, run_profiling.sh writes it).
// Human table to stderr.
//
// The filter-type rows live in separate TUs (_bw.cpp, _el.cpp) so that each
// translation unit instantiates only one filter class and stays within GCC's
// memory budget at N=25.

void run_bw_accuracy();
void run_el_accuracy();

int main()
{
    run_bw_accuracy();
    run_el_accuracy();
    return 0;
}
