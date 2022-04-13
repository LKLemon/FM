#include "mono_block.h"

void mono_process(std::vector<float> &mono_processed_data, const std::vector<float> &demod_data, const parameter &params, 
                    const std::vector<float> &mono_lpf, std::vector<float> &mono_state) {
    // filter and resample data
    // std::cerr << "begin resample\n";
    convolve_ud(mono_processed_data, demod_data, mono_lpf, params.U, params.D1, mono_state);
}
