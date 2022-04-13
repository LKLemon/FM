#ifndef DY4_RF_FRONT_END_H
#define DY4_RF_FRONT_END_H

#include "audiofunc.h"
#include "dy4.h"
#include "filter.h"
#include "mode.h"

#include "logfunc.h"

#define FC 100e3

void RF_front_end(std::vector<float> &demod_data, const std::vector<float> &raw_data, const parameter &params,
				const std::vector<float> &front_end_lpf, std::vector<float> &state_i_lpf_100k, std::vector<float> &state_q_lpf_100k, float &state_I, float &state_Q);
void split_iq_data(const std::vector<float> &iq_data, std::vector<float> &i_data, std::vector<float> &q_data);
void demodulator(std::vector<float> &fm_demod, const std::vector<float> &I,const std::vector<float> &Q,
	float &prev_I, float &prev_Q);

#endif
