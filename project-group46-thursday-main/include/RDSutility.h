#ifndef DY4_RDSUTIL_H
#define DY4_RDSUTIL_H


#include <bitset>
#include <iostream>
#include <string>
#include "dy4.h"
#include <vector>
#include <algorithm>
#include <functional>
#include "filter.h"
#include "fourier.h"
#include "RF_front_end.h"

typedef struct state_pll_q {
  float integrator;
  float phaseEst;
  float feedbackI;
  float feedbackQ;
  float	ncoOut_0;
  float	ncoOut_q_0;
  float trigOffset;
} state_pll_q;

void getParity_check_mat(std::vector<bool> &parity_check_mat);

void getSyndrome_vecs(std::vector<bool> &syndrome_vecs);

void sync_start(int &idx, std::string &offset_type, std::vector <bool> &stream, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs);

void check_syndrom(std::string &offset_type, std::vector<bool>::iterator b, std::vector<bool>::iterator e, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs);

void parse_msg_A(std::vector <bool> &msg_block);

void parse_msg_B(std::vector <bool> &msg_block, std::vector <bool> &group_type, std::string &PTY, std::vector <bool> &DI_seg);

void parse_msg_D(std::vector <bool> &msg_block, unsigned long &num1, unsigned long &num2);

void fmPll_q(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &ncoOut_q, float freq,
          float Fs, state_pll_q &state, float ncoScale, float phaseAdjust, float normBandwidth);

int find_sampling_pos(std::vector <float> data_block, int SPS);


void CDR(std::vector <float> preCDR, int SPS, int &start, std::vector <int> &state, std::vector <bool> &output, int &resync_flag);

void diff_decode(std::vector <bool> &encoded_data, bool diff_state, std::vector <bool> &decoded_data);

void vec2str(std::vector<bool> &vec, std::string &str);

void frame_sync(std::vector<bool> &decoded_data, std::vector<bool> &frame_state, int &group_state, int &sync_flag, std::vector<char> &char_buf, std::vector<bool> &DI_seg, std::vector<bool> &group_type, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs);

void fmPll_q(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &ncoOut_q, float freq,
          float Fs, state_pll_q &state, float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01);

#endif
