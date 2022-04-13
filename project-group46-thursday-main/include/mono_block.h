#ifndef DY4_MONO_BLOCK_H
#define DY4_MONO_BLOCK_H

#include "audiofunc.h"
#include "dy4.h"
#include "mode.h"
#include "filter.h"

#define MONO_FC 16e3

void mono_process(const std::vector<float> &demod_data, const parameter &params, std::vector<float> &mono_processed_data);
// void upsample(const std::vector<float> &demod_data, std::vector<float> &us_data, int us_factor);

void resample_convolve(const std::vector<float> &h, const std::vector<float> &x_data, const parameter &params, 
                        std::vector<float> &y, std::vector<float> &state);

#endif
