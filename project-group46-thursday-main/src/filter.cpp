/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the impulse response
  h.clear(); h.resize(num_taps, 0.0);
// the rest of the code in this function is to be completed by you
// based on your understanding and the Python code from the first lab
  float normcutoff = Fc/(Fs/2);
  for (int i = 0; i<num_taps; ++i){
    float h_i = 0;
    if (i==(num_taps-1)/2){
      h_i = normcutoff;
    } else {
      h_i= normcutoff * std::sin(PI*normcutoff*(i-(num_taps-1)/2))/(PI*normcutoff*(i-(num_taps-1)/2));
    }
    h[i] = h_i * (std::pow(std::sin(i*PI/num_taps), 2));
  }
}

// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF_amp(float Fs, float Fc, unsigned short int num_taps, int amplify_factor, std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the impulse response
  h.clear(); h.resize(num_taps, 0.0);
// the rest of the code in this function is to be completed by you
// based on your understanding and the Python code from the first lab
  float normcutoff = Fc/(Fs/2);
  for (int i = 0; i<num_taps; ++i){
    float h_i = 0;
    if (i==(num_taps-1)/2){
      h_i = normcutoff;
    } else {
      h_i= normcutoff * std::sin(PI*normcutoff*(i-(num_taps-1)/2))/(PI*normcutoff*(i-(num_taps-1)/2));

    }
    h[i] = amplify_factor * (h_i * (std::pow(std::sin(i*PI/num_taps), 2)));
  }
}

void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h) {
	h.clear(); h.resize(num_taps, 0.0);
  float norm_center = ((Fe + Fb)/2)/(Fs/2);
	float norm_pass = (Fe - Fb)/(Fs/2);

	for(int i = 0; i < num_taps; i++){
		if(i == (num_taps - 1)/2){
			h[i] = norm_pass;
		}
		else{
			h[i] = norm_pass*(std::sin(PI*(norm_pass/2)*(i-(num_taps-1)/2))/(PI*(norm_pass/2)*(i-(num_taps-1)/2)));
		}
		h[i] = h[i]*std::cos(i*PI*norm_center);
		h[i] = h[i]*std::sin(i*PI/num_taps)*std::sin(i*PI/num_taps);
	}

}

void APF(std::vector<float> &output_block, const std::vector<float> &input_block, std::vector<float> &state_block) {
  std::copy(state_block.begin(), state_block.end(), output_block.begin());
  std::copy(input_block.begin(), input_block.end()-state_block.size(), output_block.begin()+state_block.size());
  std::copy(input_block.end()-state_block.size(), input_block.end(), state_block.begin());
}

void RRC(float Fs,int N_taps,std::vector<float> &impulseResponseRRC ){
  float T_symbol = 1/2375.0;
  float beta = 0.9;
  impulseResponseRRC.resize(N_taps,0);

  // std::vector impulseResponseRRC = resize(N_taps,0);
  for (int k=0; k < N_taps; k++){
    float t;
    t = float((k - N_taps/2))/Fs;
    if(t == 0.0)
      impulseResponseRRC[k] = 1.0 + beta *((4/PI)-1);
    else if ((t == -T_symbol/(4*beta)) || (t == T_symbol/(4*beta))) {
      impulseResponseRRC[k] = (beta/std::sqrt(2))*(((1+2/PI)*
					(std::sin(PI/(4*beta)))) + ((1-2/PI)*(std::cos(PI/(4*beta)))));

    } else {
      impulseResponseRRC[k] = (std::sin(PI*t*(1-beta)/T_symbol) + 
  					4*beta*(t/T_symbol)*std::cos(PI*t*(1+beta)/T_symbol))/
  					(PI*t*(1-(4*beta*t/T_symbol)*(4*beta*t/T_symbol))/T_symbol);
    }
  }
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// bring your own functionality
  y.clear(); y.resize(x.size()+h.size()-1, 0.0);

// the rest of the code in this function is to be completed by you
// based on your understanding and the Python code from the first lab
  for (unsigned int i=0; i < y.size(); i++) {
    y[i] = 0;
    for (unsigned int j=0; j < h.size(); j++) {
      if (i-j >= 0 && i-j < x.size()) {
        y[i] += h[j] * x[i-j];
      }
    }
  }
}

void convolve_ds(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,
						const int D, std::vector<float> &state)
{
  y.resize(x.size()/D);
  // std::cout << y.size() << std::endl;
	for (uint i=0; i < y.size(); i++) {
    // std::cout << "calc y" << i << std::endl;
		y[i] = 0;
    int idx;
		for (unsigned int j=0; j < h.size(); j++) {
      idx = (int) (D*i-j);
			if (idx >= 0) {
        // std::cout << "term " << D*i-j << std::endl;
				y[i] += h[j] * x[idx];
			} else { // account for state
				y[i] += h[j] * state[state.size()+idx];
			}
		}
	}
	// save last h.size()-1 elements in current block
	for (unsigned int i=0; i < state.size(); i++) {
		state[i] = x[x.size()-state.size()+i];
	}
}

void convolve_ud(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,
                const int U, const int D, std::vector<float> &state) {
    // prepare output vector
    y.resize(x.size() * U / D);
    // std::cerr << x.size() << " " << y.size() << "\n";
    // account for downsample
    for (uint i=0; i < y.size(); i++) {
        y[i] = 0;
        // account for upsample
        int phase = (D * i) % U;
        int x_idx;
        for (uint j=phase; j < h.size(); j += U) {
            x_idx = (int)(i*D - j) / U;
            if (x_idx >= 0) {
                y[i] += h[j] * x[x_idx];
            } else {
                // account for state
                y[i] += h[j] * state[state.size()+x_idx];
            }
        }
    }
    // store the state
    for (uint i = 0; i < state.size(); i++) {
        state[i] = x[(int)(x.size()-state.size()+i)];
    }
}
