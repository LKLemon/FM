#ifndef DY4_AUDIOFUNC_H
#define DY4_AUDIOFUNC_H

#include <fstream>
#include <iostream>
#include <vector>

// #define BLOCK_SIZE 102400

void read_stdin_block(uint num_samples, uint block_id, std::vector<float> &block_data);
void read_audio_data(const std::string in_fname, std::vector<uint8_t> &audio_data);
void split_audio_into_channels(const std::vector<float> &audio_data, std::vector<float> &audio_left, std::vector<float> &audio_right);
void write_audio_data(const std::string out_fname, const std::vector<float> &audio_left, const std::vector<float> &audio_right);

#endif
