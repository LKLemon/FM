#ifndef DY4_RF_STEREO_PROCESSING_H
#define DY4_RF_STEREO_PROCESSING_H

#include "dy4.h"
#include <algorithm>
#include <stdio.h>      /* printf */
#include <math.h>       /* atan2 */
#include <functional>
#include "filter.h"
#include "fourier.h"
#include "RF_front_end.h"


#define PLL_F       19e3
#define P_FB        18.5e3
#define P_FE        19.5e3
#define STEREO_FB   22e3
#define STEREO_FE   54e3
#define PILOT_FB    18.5e3
#define PILOT_FE    19.5e3


typedef struct state_pll {
  float integrator;
  float phaseEst;
  float feedbackI;
  float feedbackQ;
  float	ncoOut_0;
  float trigOffset;
} state_pll;

void fmPll(std::vector<float> &pllIn, std::vector<float> &ncoOut, float freq,
          float Fs, state_pll &state, float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01);
void stereo_processing(std::vector<float> &left_channel, std::vector<float> &right_channel, const std::vector<float> &mono_data, const std::vector<float> &demod_data, \
                       std::vector<float> &state_stereo, std::vector<float> &state_pilot, std::vector<float> &state_stereo_filt, \
                       const std::vector<float> &pilot_coeff, const std::vector<float> &stereo_coeff, \
		       		const std::vector<float> &stereo_filt_coeff, state_pll &state, const parameter &params);

#endif
