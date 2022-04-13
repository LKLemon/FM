#include "stereo_Processing.h"

void stereo_processing(std::vector<float> &left_channel, std::vector<float> &right_channel, const std::vector<float> &mono_data, const std::vector<float> &demod_data,
                       std::vector<float> &state_stereo, std::vector<float> &state_pilot, std::vector<float> &state_stereo_filt,
                       const std::vector<float> &pilot_coeff, const std::vector<float> &stereo_coeff,
		       		const std::vector<float> &stereo_filt_coeff, state_pll &state, const parameter &params){

    std::vector<float> stereo_carrier;
    std::vector<float> stereo_channel_data;
    std::vector<float> pilot;

	// static int debug_count;
	static std::vector<float> carrier(40,0);

    // STEREO PATH
    // Carrier Recovery
    convolve_ds(pilot, demod_data, pilot_coeff, 1, state_pilot);
    fmPll(pilot, stereo_carrier, PLL_F, params.IF, state, 2.0);

    // Stereo Channel Extraction
    convolve_ds(stereo_channel_data, demod_data, stereo_coeff, 1, state_stereo);

    // Stereo Processing
    // Mixer
    std::vector<float> mixed_data(stereo_channel_data.size(), 0);
    std::transform(stereo_channel_data.begin(), stereo_channel_data.end(),stereo_carrier.begin(), mixed_data.begin(), std::multiplies<float>()); // assumes v1,v2 of same size > 1,
	for (uint i=0; i < mixed_data.size(); i++) {
		mixed_data[i] *= 2;
	}

    // Filtering and rate conversion
    std::vector<float> stereo_filt;
    convolve_ud(stereo_filt, mixed_data, stereo_filt_coeff, params.U, params.D1, state_stereo_filt);

	left_channel.resize(stereo_filt.size(), 0);
	right_channel.resize(stereo_filt.size(), 0);
	std::transform( mono_data.begin(), mono_data.end(), stereo_filt.begin(), left_channel.begin(), std::plus<float>() );
	std::transform( mono_data.begin(), mono_data.end(), stereo_filt.begin(), right_channel.begin(), std::minus<float>() );

  
	// if (debug_count == 5) {
	// 	std::copy(stereo_carrier.end()-10, stereo_carrier.end(), carrier.begin());
	// }
	// if (debug_count == 6) {
	// 	std::copy(stereo_carrier.begin()+1, stereo_carrier.begin()+11, carrier.begin()+10);
	// 	std::vector<float> vector_index;
	// 	genIndexVector(vector_index, carrier.size());
	// 	logVector("carrier_time", vector_index, pilot);
	// }
	// debug_count++;
}

void fmPll(std::vector<float> &pllIn, std::vector<float> &ncoOut, float freq, float Fs, state_pll &state, float ncoScale, float phaseAdjust, float normBandwidth)
{


	/*
	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)
	freq 			float
					reference frequency to which the PLL locks
	Fs  			float
					sampling rate for the input/output signals
	ncoScale		float
					frequency scale factor for the NCO output
	phaseAdjust		float
					phase adjust to be added to the NCO output only
	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)
	state 			to be added
	*/

  	float Cp = 2.666;
	float Ci = 3.555;


	float Kp = (normBandwidth)*Cp;

	float Ki = (normBandwidth*normBandwidth)*Ci;

	//output array for the NCO
	ncoOut.resize(pllIn.size()+1,0);

	// initialize internal state
  	float integrator = state.integrator;
	float phaseEst = state.phaseEst;
	float feedbackI = state.feedbackI;
	float feedbackQ = state.feedbackQ;
  	ncoOut[0] = state.ncoOut_0;
	float trigOffset = state.trigOffset;
	// note: state saving will be needed for block processing

	for(uint k = 0; k < pllIn.size(); k++){

		// phase detector
		float errorI = pllIn[k] * (+feedbackI);  // complex conjugate of the
		float errorQ = pllIn[k] * (-feedbackQ);  // feedback complex exponential

		// four-quadrant arctangent discriminator for phase error detection
		float errorD = std::atan2(errorQ, errorI);

		// loop filter
		integrator = integrator + Ki*errorD;

		// update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

		// internal oscillator
		trigOffset += 1;
		float trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);
	}

  	state.integrator = integrator;
	state.phaseEst = phaseEst;
	state.feedbackI = feedbackI;
	state.feedbackQ = feedbackQ;
  	state.ncoOut_0 = ncoOut[ncoOut.size()-1];
	state.trigOffset = trigOffset;

	// for stereo only the in-phase NCO component should be returned
	// for block processing you should also return the state


	// for RDS add also the quadrature NCO component to the output
}
