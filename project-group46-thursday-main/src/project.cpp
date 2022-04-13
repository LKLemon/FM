/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "audiofunc.h"
#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include "mode.h"
#include "mono_block.h"
#include "RF_front_end.h"
#include "stereo_Processing.h"
#include "RDS.h"

#include <atomic>
#include <thread>

int mode = 0;
int channel = 1;

parameter params = set_parameter(mode);

std::vector<bool> parity_check_mat;
std::vector<bool> syndrome_vecs;

int block_size;
unsigned short mono_taps;
unsigned short stereo_filt_taps;
int rf_ele_size;
int queue_element = 8;

void RF_thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& audio_read_offset, std::atomic<int>& RDS_read_offset);
void audio_thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset);
void RDS_thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset);


int main(int argc, char *argv[])
{
	// validate input and set mode

	if (argc < 2) {
		std::cerr << "Operating in default mode 0" << std::endl;
	} else if (argc == 2) {
		mode = atoi(argv[1]);
		if (mode > 3) {
			std::cerr << "Wrong mode" << mode << std::endl;
			exit(1);
		}
	} else if (argc == 3) {
		channel = atoi(argv[2]);
	} else {
		std::cerr << "Usage: " << argv[0] << std::endl;
		std::cerr << "or " << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode> " << std::endl;
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;

		exit(1);
	}
	std::cerr << "Operating in mode " << mode << std::endl;

	// matrix for RDS frame sync
	getParity_check_mat(parity_check_mat);
	getSyndrome_vecs(syndrome_vecs);

	// parameter params;
	// params = set_parameter(mode);
	block_size = 1024 * params.D0 * params.D1 * 2;
	// block_size = 1024 * params.D_RDS * 50 * 2;
	mono_taps = TAPS * params.U;
	// stereo process facilities
	stereo_filt_taps = TAPS * params.U;

	rf_ele_size = block_size/2/params.D0;

	std::vector<float> rf_queue(rf_ele_size*queue_element);
	std::atomic<int> RF_write_offset(0);
	std::atomic<int> audio_read_offset(0);
	std::atomic<int> RDS_read_offset(0);
//-----------------------front end----------------------------//
	std::thread RF_producer = std::thread(RF_thread, std::ref(rf_queue),std::ref(RF_write_offset),std::ref(audio_read_offset), std::ref(RDS_read_offset));
//-----------------------audio Path---------------------------//
	std::thread audio_consumer = std::thread(audio_thread, std::ref(rf_queue),std::ref(RF_write_offset),std::ref(audio_read_offset));
//-----------------------RDS path-----------------------------//
	std::thread RDS_consumer = std::thread(RDS_thread, std::ref(rf_queue),std::ref(RF_write_offset),std::ref(RDS_read_offset));


	RF_producer.join();
	audio_consumer.join();
	RDS_consumer.join();

	return 0;
}


/*----------------------------------thread func() ------------------------------------ */

void RF_thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& audio_read_offset, std::atomic<int>& RDS_read_offset){
  std::vector<float> front_end_lpf;
  impulseResponseLPF(params.RF, FC, TAPS, front_end_lpf);
  std::vector<float> state_i_lpf_100k(TAPS-1, 0.0);
  std::vector<float> state_q_lpf_100k(TAPS-1, 0.0);
  std::vector<float> input_demod;
  float state_I = 0;
  float state_Q = 0;
  for (unsigned int block_id=0; ; block_id++){
	std::vector<float> block_data(block_size);
	read_stdin_block(block_size, block_id, block_data);
	// std::cerr << "read block" << block_id << "\n";
	if((std::cin.rdstate()) != 0) {
		std::cerr << "End of input stream reached" << std::endl;
		exit(1);
	}
    RF_front_end(input_demod, block_data, params, front_end_lpf, state_i_lpf_100k, state_q_lpf_100k, state_I, state_Q);
    //copy to queue
	// wait until not "full"
    while(write_offset.load() >= (audio_read_offset.load()+queue_element) || write_offset.load() >= (RDS_read_offset.load()+queue_element)) {
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }    
    std::vector<float>::difference_type address_offset = (write_offset.load() % queue_element)*rf_ele_size;
    std::copy_n(input_demod.begin(), input_demod.size(), rf_queue.begin()+address_offset);
    write_offset.fetch_add(1);
  }
}

void audio_thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset){
  // mono process facilities
  // unsigned short mono_taps = TAPS * params.U;
  std::vector<float> mono_lpf(mono_taps);
  std::vector<float> state_mono_lpf(mono_taps-1, 0.0);
  std::vector<float> mono_processed_data;
  std::vector<float> state_demod((TAPS-1)/2, 0);

  // stereo process facilities
  std::vector<float> pilot_coeff;
  impulseResponseBPF(params.IF, PILOT_FB, PILOT_FE, TAPS, pilot_coeff);
  std::vector<float> stereo_coeff;
  std::vector<float> stereo_filt_coeff;
  unsigned short stereo_filt_taps = TAPS * params.U;

  std::vector<float> state_pilot(TAPS-1);
  std::vector<float> state_stereo(TAPS-1);
  std::vector<float> state_stereo_filt(stereo_filt_taps-1);
  state_pll pll_state;
  pll_state.integrator = 0.0;
  pll_state.phaseEst = 0.0;
  pll_state.feedbackI = 1.0;
  pll_state.feedbackQ = 0.0;
  pll_state.ncoOut_0 = 1.0;
  pll_state.trigOffset = 0.0;

  std::vector<float> left_channel;
  std::vector<float> right_channel;
  std::vector<float> input_demod(rf_ele_size);

  std::vector<float> delayed_demod(input_demod.size(), 0);
  impulseResponseLPF_amp(params.IF * params.U, MONO_FC, mono_taps, params.U, mono_lpf);

  impulseResponseBPF(params.IF, STEREO_FB, STEREO_FE, TAPS, stereo_coeff);

  impulseResponseLPF_amp(params.IF * params.U, MONO_FC, stereo_filt_taps, params.U, stereo_filt_coeff);
  while(1){

    while(write_offset.load() <= read_offset.load()) std::this_thread::sleep_for(std::chrono::milliseconds(10));
    std::vector<float>::difference_type address_offset = (read_offset.load() % queue_element)*rf_ele_size;
    std::copy_n(rf_queue.begin()+address_offset, input_demod.size(),input_demod.begin());
    read_offset.fetch_add(1);

    APF(delayed_demod, input_demod, state_demod);
    convolve_ud(mono_processed_data, delayed_demod, mono_lpf, params.U, params.D1, state_mono_lpf);

    // stereo path
    stereo_processing(left_channel, right_channel, mono_processed_data, input_demod,
              state_stereo, state_pilot, state_stereo_filt,
              pilot_coeff, stereo_coeff, stereo_filt_coeff, pll_state, params);

    if (channel == 0) {
      std::vector<short int> mono_audio_data(mono_processed_data.size());
      for (uint k=0; k < mono_processed_data.size(); k++) {
        if (std::isnan(mono_processed_data[k])) mono_audio_data[k] = 0;
        else mono_audio_data[k] = static_cast<short int> (mono_processed_data[k] * 16384);
      }
      fwrite(&mono_audio_data[0], sizeof(short int), mono_audio_data.size(), stdout);
    } else if (channel == 1) {
      std::vector<short int> stereo_audio_data(left_channel.size()*2);
      for (uint k=0; k < stereo_audio_data.size(); k++) {
        if (k%2 == 0) {
          // left channel
          if (std::isnan(left_channel[k/2])) stereo_audio_data[k] = 0;
          else stereo_audio_data[k] = static_cast<short int> (left_channel[k/2] * 16384);
        } else {
          // right channel
          if (std::isnan(right_channel[(k-1)/2])) stereo_audio_data[k] = 0;
          else stereo_audio_data[k] = static_cast<short int> (right_channel[(k-1)/2] * 16384);
        }

      }

      fwrite(&stereo_audio_data[0], sizeof(short int), stereo_audio_data.size(), stdout);
    }
  }
}

void RDS_thread(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset) {
	// state for RDS
	state_pll_q state_RDS_pll;
	state_RDS_pll.integrator = 0.0;
	state_RDS_pll.phaseEst = 0.0;
	state_RDS_pll.feedbackI = 1.0;
	state_RDS_pll.feedbackQ = 0.0;
	state_RDS_pll.ncoOut_0 = 1.0;
	state_RDS_pll.ncoOut_q_0 = 0.0;
	state_RDS_pll.trigOffset = 0.0;

	std::vector<float> state_RDS(RDS_taps-1, 0);
	std::vector <float> state_allpass((RDS_taps-1)/2, 0);
    std::vector <float> state_recover(RDS_taps-1, 0);
	std::vector <float> state_rational(params.U_RDS*RDS_taps-1, 0);
	std::vector <float> state_rational_q(params.U_RDS*RDS_taps-1, 0);
    std::vector <float> state_RRC(RDS_taps-1, 0);
	std::vector <float> state_RRC_q(RDS_taps-1, 0);
	std::vector <int> state_CDR(1,0);
	bool diff_state = 0;
	std::vector <bool> frame_state;
    std::vector<bool> group_type(5,0);
	int group_state = 0;
	int sync_flag = 1;
	std::vector<bool> DI_seg(2,0);
	std::vector<char> char_buf(8,-1);

	// coefficients of RDS
	// std::cerr << "Prepare coeffs\n"; 
    std::vector<float> RDS_channel_coeff;
	impulseResponseBPF(params.IF, RDS_channel_Fb, RDS_channel_Fe, RDS_taps, RDS_channel_coeff);
	std::vector<float> RDS_recover_coeff;
	impulseResponseBPF(params.IF, RDS_recover_Fb, RDS_recover_Fe, RDS_taps, RDS_recover_coeff);
	std::vector<float> rational_coeff;
	impulseResponseLPF_amp(params.IF*params.U_RDS, RDS_demod_Fc, RDS_taps*params.U_RDS, params.U_RDS, rational_coeff);
	std::vector<float> RRC_coeff;
	RRC(params.SPS*2375, RDS_taps, RRC_coeff);

	std::vector<float> input_demod(rf_ele_size);

	while(1){

		while(write_offset.load() <= read_offset.load()) std::this_thread::sleep_for(std::chrono::milliseconds(10));
		std::vector<float>::difference_type address_offset = (read_offset.load() % queue_element)*rf_ele_size;
		std::copy_n(rf_queue.begin()+address_offset, input_demod.size(),input_demod.begin());
		read_offset.fetch_add(1);
		// std::vector<float> vector_index;
		// vector_index.resize(input_demod.size());
		// genIndexVector(vector_index, input_demod.size());
		// logVector("input_demod_data", vector_index, input_demod);		
		// std::cerr << "Begin rds process\n";
		RDS_processing(input_demod, params,
						state_RDS_pll, state_RDS, state_allpass, state_recover, state_rational, state_rational_q, state_RRC, state_RRC_q,
						state_CDR, diff_state, frame_state, group_type, group_state, sync_flag, DI_seg, char_buf,
						RDS_channel_coeff, RDS_recover_coeff, rational_coeff, RRC_coeff, parity_check_mat, syndrome_vecs);
	}

}
