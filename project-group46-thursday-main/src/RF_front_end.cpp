#include "RF_front_end.h"

// perform demodulation of i/q data
void RF_front_end(std::vector<float> &demod_data, const std::vector<float> &raw_data, const parameter &params,
				const std::vector<float> &front_end_lpf, std::vector<float> &state_i_lpf_100k, std::vector<float> &state_q_lpf_100k,
				float &state_I, float &state_Q) {
	// split I Q data
	std::vector<float> i_data, q_data;
	// std::cout << "begin iq split\n";
	split_iq_data(raw_data, i_data, q_data);
	// filter out the FM channel
	std::vector<float> i_filt, q_filt;
	convolve_ud(i_filt, i_data, front_end_lpf, 1, params.D0, state_i_lpf_100k);
	convolve_ud(q_filt, q_data, front_end_lpf, 1, params.D0, state_q_lpf_100k);

	// std::cout << "begin demodulator\n";
	demodulator(demod_data, i_filt, q_filt, state_I, state_Q);
}

void demodulator(std::vector<float> &fm_demod, const std::vector<float> &I,const std::vector<float> &Q,
	float &prev_I, float &prev_Q) {
	fm_demod.resize(I.size(),0);
	uint k;
	for (k = 0; k < I.size(); k++) {
		float denom = std::pow(I[k],2) + std::pow(Q[k],2);
		if (k == 0){
			if (denom == 0) {
				fm_demod[k] = 0;
			} else {
				fm_demod[k] = (I[k]*(Q[k]-prev_Q) - Q[k]*(I[k]-prev_I)) / denom;
			}
		} else{
			if (denom == 0) {
				fm_demod[k] = fm_demod[k-1];
			} else {
				fm_demod[k] = (I[k]*(Q[k]-Q[k-1]) - Q[k]*(I[k]-I[k-1])) / denom;
			}
		}
	}
	prev_I = I[k-1];
	prev_Q = Q[k-1];
}

void split_iq_data(const std::vector<float> &iq_data, std::vector<float> &i_data, std::vector<float> &q_data){
	for (uint i=0; i<iq_data.size(); i++) {
		if (i%2==0) {
			i_data.push_back(iq_data[i]);
		} else {
			q_data.push_back(iq_data[i]);
		}
	}
}
