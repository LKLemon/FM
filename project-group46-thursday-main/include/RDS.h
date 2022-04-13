#ifndef DY4_RDS_H
#define DY4_RDS_H

#include "dy4.h"
#include <iostream>
#include <algorithm> 
#include <cmath>
#include <functional>
#include "filter.h"
#include "fourier.h"
#include "RDSutility.h"

#define RDS_channel_Fb  54e3
#define RDS_channel_Fe  60e3
#define RDS_recover_Fb  113.5e3
#define RDS_recover_Fe  114.5e3
#define pll_RDS_f       114e3
#define RDS_Fs          2.4e5
#define RDS_taps        151
#define RDS_demod_Fc    3e3

void RDS_processing(const std::vector <float> &fm_demod, const parameter &params,
                    state_pll_q &state_RDS_pll, std::vector<float> &state_RDS, std::vector <float> &state_allpass,
                    std::vector <float> &state_recover, std::vector <float> &state_rational, std::vector <float> &state_rational_q,
                    std::vector <float> &state_RRC, std::vector <float> &state_RRC_q, std::vector <int> &state_CDR, bool &diff_state,std::vector <bool> &frame_state,
                    std::vector<bool> &group_type, int group_state, int sync_flag, std::vector<bool> &DI_seg, std::vector<char> &char_buf,
                    std::vector<float> &RDS_channel_coeff, std::vector<float> &RDS_recover_coeff, std::vector<float> &rational_coeff,
                    std::vector<float> &RRC_coeff, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs);


#endif

