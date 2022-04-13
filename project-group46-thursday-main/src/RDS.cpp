#include "RDS.h"

void RDS_processing(const std::vector <float> &fm_demod, const parameter &params,
                    state_pll_q &state_RDS_pll, std::vector<float> &state_RDS, std::vector <float> &state_allpass,
                    std::vector <float> &state_recover, std::vector <float> &state_rational, std::vector <float> &state_rational_q,
                    std::vector <float> &state_RRC, std::vector <float> &state_RRC_q, std::vector <int> &state_CDR, bool &diff_state,std::vector <bool> &frame_state,
                    std::vector<bool> &group_type, int group_state, int sync_flag, std::vector<bool> &DI_seg, std::vector<char> &char_buf,
                    std::vector<float> &RDS_channel_coeff, std::vector<float> &RDS_recover_coeff, std::vector<float> &rational_coeff,
                    std::vector<float> &RRC_coeff, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs) {   

   //RDS Channel Extraction

   std::vector <float> RDS_SN_data;
   std::vector <float> RDS_channel_data;
   std::vector <float> RDS_carrier;
   std::vector <float> RDS_carrier_q;
   std::vector <float> RDS_allpass_data;
   std::vector <float> RDS_recover_data;
   std::vector <float> RDS_mixed_data;
   std::vector <float> RDS_mixed_data_q;
   std::vector <float> preCDR;
   std::vector <float> preCDR_q;
   std::vector <bool> CDR_data;
   std::vector <bool> decoded_data;
   std::vector <float> preRRC_q;
   std::vector <float> preRRC;
   std::vector<float> vector_index;


   
   convolve_ds(RDS_channel_data, fm_demod, RDS_channel_coeff, 1.0, state_RDS);
   //All pass filter to create delay
   RDS_allpass_data.resize(RDS_channel_data.size());
   APF(RDS_allpass_data, RDS_channel_data, state_allpass);

   //Carrier Recovery
   RDS_SN_data.resize(RDS_channel_data.size());
   std::transform(RDS_channel_data.begin(), RDS_channel_data.end(), RDS_channel_data.begin(), RDS_SN_data.begin(), std::multiplies<float>());
   convolve_ds(RDS_recover_data, RDS_SN_data, RDS_recover_coeff, 1.0, state_recover);
   fmPll_q(RDS_recover_data, RDS_carrier, RDS_carrier_q, pll_RDS_f, RDS_Fs, state_RDS_pll, 0.5, 1.18, 0.001);
   //RDS Demodulatiom
   //Mixer
   RDS_carrier.pop_back();
   RDS_mixed_data.resize(RDS_carrier.size());
   std::transform(RDS_carrier.begin(), RDS_carrier.end(), RDS_allpass_data.begin(), RDS_mixed_data.begin(), std::multiplies<float>() );
   for(int i = 0; i < RDS_mixed_data.size(); i++){
     RDS_mixed_data[i]*=2;
   }
   // Rational resampler
    convolve_ud(preRRC, RDS_mixed_data, rational_coeff, params.U_RDS, params.D_RDS, state_rational);

   //Root Raised Cosine Filter
   convolve_ds(preCDR, preRRC, RRC_coeff, 1, state_RRC);

   int start = find_sampling_pos(preCDR, params.SPS);
   CDR(preCDR, params.SPS, start, state_CDR, CDR_data, sync_flag);

   // RDS Data Processing
   diff_decode(CDR_data, diff_state, decoded_data);

   // frame_sync
    frame_sync(decoded_data, frame_state, group_state, sync_flag, char_buf, DI_seg, group_type, parity_check_mat, syndrome_vecs);
}
