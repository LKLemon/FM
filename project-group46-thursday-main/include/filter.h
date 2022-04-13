/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <cmath>
#include <vector>

// declaration of a function prototypes
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void impulseResponseLPF_amp(float Fs, float Fc, unsigned short int num_taps, int amplify_factor, std::vector<float> &h);
void impulseResponseBPF(float Fs, float Fb, float Fe, unsigned short int num_taps, std::vector<float> &h);
void APF(std::vector<float> &output_block, const std::vector<float> &input_block, std::vector<float> &state_block);
void RRC(float Fs,int N_taps,std::vector<float> &impulseResponseRRC);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void convolve_ud(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, 
                const int U, const int D, std::vector<float> &state);
void convolve_ds(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h,
        const int D, std::vector<float> &state);
#endif // DY4_FILTER_H
