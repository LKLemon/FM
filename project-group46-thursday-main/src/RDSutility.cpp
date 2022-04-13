#include "RDSutility.h"

void getParity_check_mat(std::vector<bool> &parity_check_mat){
    std::string bits = "10000000000100000000001000000000010000000000100000000001000000000010000000000100000000001000000000011011011100010110111000101101111010000111111001111111000100111101010101110111011001101110111000000001111101110001111011100011110111101010011111100011111100011011";

// parity_check_mat.reserve(bits.size());
// for(auto a : bits){
//     parity_check_mat.push_back(a == '1');
// }
    for(int i = 0; i < bits.size(); i++){
        parity_check_mat.push_back(bits[i] == '1');
    }
}

void getSyndrome_vecs(std::vector<bool> &syndrome_vecs){
    std::string bits = "11110110001111010100100101110011110011001001011000";
    /*parity_check_mat.reserve(bits.size());
    for(auto a : bits){
    parity_check_mat.push_back(a == '1');
}*/
    for(int i = 0; i < bits.size(); i++){
        syndrome_vecs.push_back(bits[i] - '0');
    }
}

std::string program_type[32] = {"No programme type or undefined",
    "News",
    "Information",
    "Sports",
    "Talk",
    "Rock",
    "Classic rock",
    "Adult hits",
    "Soft rock",
    "Top 40",
    "Country",
    "Oldies",
    "Soft music",
    "Nostalgia",
    "Jazz",
    "Classical",
    "Rhythm and blues",
    "Soft rhythm and blues",
    "Language",
    "Religious music",
    "Religious talk",
    "Personality",
    "Public",
    "College",
    "Spanish Talk",
    "Spanish Music",
    "Hip hop",
    "Unassigned",
    "Unassigned",
    "Weather",
    "Emergency test",
    "Emergency"};




void sync_start(int &idx, std::string &offset_type, std::vector <bool> &stream, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs){
    std::cerr << "first check\n";
    check_syndrom(offset_type, stream.begin()+idx, stream.begin()+idx+26, parity_check_mat, syndrome_vecs);
    
    while(offset_type != "A" && (idx+26) < stream.size()){
        idx += 1;
        check_syndrom(offset_type, stream.begin()+idx, stream.begin()+idx+26, parity_check_mat, syndrome_vecs);
        // std::cerr << offset_type << "\n";
    }
    
}

void check_syndrom(std::string &offset_type, std::vector<bool>::iterator b, std::vector<bool>::iterator e, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs){
    std::vector <bool> syndrome_vec;
    std::vector <bool> parity_vec;
    std::vector <bool> prod;
    prod.resize(26,0);
    syndrome_vec.resize(10);
    for(int col = 0; col < 10; col++){
        parity_vec.resize(26,0);
        for(int row = 0; row < 26; row++){
            parity_vec[row] = parity_check_mat[row*10+col];
        }
        // element-wise multiplication (AND)
        std::transform(b,e,parity_vec.begin(), prod.begin(), std::multiplies<bool>());
        syndrome_vec[col] = (count(prod.begin(), prod.end(), 1) % 2);  // sum (XOR) of products
    }
    size_t found = std::search(syndrome_vecs.begin(), syndrome_vecs.end(), syndrome_vec.begin(), syndrome_vec.end()) - syndrome_vecs.begin();
    // std::cerr << found << "\n";
    if(found <= 40){
        size_t idx = found;
        if(idx == 0){
            offset_type =  "A";
        }
        else if(idx == 10){
            offset_type = "B";
        }
        else if(idx == 20){
            offset_type = "Cp";
        }
        else if(idx == 30){
            offset_type = "C";
        }
        else if(idx == 40){
            offset_type = "D";
        }
    }
    offset_type = "None";
}

void parse_msg_A(std::vector <bool> &msg_block){
    std::vector <bool> PI_code;
    PI_code.resize(16,0);
    std::string output;

    std::copy(msg_block.begin(), msg_block.begin()+16, PI_code.begin());
    std::cerr << "PI code: " << std::endl;
    vec2str(PI_code, output);
    std::cerr << std::hex << std::stoi(output.c_str(), nullptr, 2) << std::dec << std::endl;
}

void parse_msg_B(std::vector <bool> &msg_block, std::vector <bool> &group_type, std::string &PTY, std::vector <bool> &DI_seg){
    // group_type = msg_block[:4]
    std::copy(msg_block.begin(), msg_block.begin()+4, group_type.begin());
    // PTY = msg_block[6:11]
    std::vector<bool> temp(5,0);
    std::string str;
    std::copy(msg_block.begin() + 6, msg_block.begin()+11, temp.begin());
    vec2str(temp, str);
    PTY = program_type[std::stoi(str.c_str(), nullptr, 2)];
    // DI_seg = msg_block[14:16]
    std::copy(msg_block.begin() + 14, msg_block.begin() + 16, DI_seg.begin());
}

void parse_msg_D(std::vector <bool> &msg_block, char &num1, char &num2){

    std::vector <bool> temp1;
    std::vector <bool> temp2;

    std::string str1;
    std::string str2;
    
    std::copy(msg_block.begin(), msg_block.begin() + 8, temp1.begin());
    std::copy(msg_block.begin() + 8, msg_block.begin() + 16, temp2.begin());

    vec2str(temp1, str1);
    vec2str(temp2, str2);
    num1 = std::stoi(str1.c_str(), nullptr, 2);
    num2 = std::stoi(str2.c_str(), nullptr, 2);
}

int find_sampling_pos(std::vector <float> data_block, int SPS){
    int offset = SPS/2;
    float mid = 0.0;
    int i, j, pos;
    std::vector <float> chunk;

    std::vector<float> array;
    array.resize(SPS,0);

    for(int i = 0; i < SPS; i++){
        float error = 0.0;
        j = i;
        chunk.resize(SPS+1, 0);
        while(j+SPS+1 < data_block.size()){
            std::copy(data_block.begin()+j, data_block.begin()+j+SPS+1, chunk.begin());
            mid = chunk[offset];
            if(chunk[0] * chunk.back() >= 0){
                // HH/LL mid point should close to the neiboring samples
                error += abs(mid - (chunk[0] + chunk.back())/2);
            }
            else{
                // HL/LH mid point should close to the 0
                error += abs(mid - 0);
            }

            j += SPS;
        }
        array[i] = error;
    }

    pos = std::min_element(array.begin(),array.end()) - array.begin();

    return pos;

}

void CDR(std::vector <float> preCDR, int SPS, int &start, std::vector <int> &state, std::vector <bool> &output, int &resync_flag){
    int i, j;
    std::vector <int> indicator;
    std::vector <int> temp;
    output.clear();
    indicator.resize(preCDR.size()+state.size(),0);
    std::copy(state.begin(), state.end(), indicator.begin());
    state.clear();
    // std::cerr << "indicator init\n";
    for (i = start; i < preCDR.size(); i += SPS){
        if (preCDR[i] > 0 ) indicator.push_back(1);
        else indicator.push_back(-1);
    }
    //first block or sync error happened, need to detect HH or LL
    if (resync_flag){
        // std::cerr << "sync first block\n";
        for (i = 0; i < (indicator.size()-1);i += 2){
            //have pair of HH or LL, move back 1 symbol
            if(indicator[i] == indicator[i+1]){
                output.clear();
                for (j = 1; j < (indicator.size()-1); j += 2){
                    if (indicator[j] == 1) output.push_back(1);
                    else output.push_back(0);
                }
                //if have left symbol, pair it with next block
                if ((indicator.size()-1)%2 != 0) state.push_back(indicator[indicator.size()-1]);
                // return output, state, start
            }
            //normal pair
            if (indicator[i] == 1) output.push_back(1);
            else output.push_back(0);
        }
        //if have left symbol, pair it with next block
        if (indicator.size()%2 != 0) state.push_back(indicator[indicator.size()-1]);
    }
    //normal block, pair directly
    else{
        for (i = 0; i < (indicator.size()-1); i += 2){
            if (indicator[i] == 1) output.push_back(1);
            else output.push_back(0);
        }
        if (indicator.size()%2 != 0) state.push_back(indicator[indicator.size()-1]);
    }
    // return output, state, start
}

void diff_decode(std::vector <bool> &encoded_data, bool diff_state, std::vector <bool> &decoded_data){
    decoded_data.resize(encoded_data.size(),0);

    for(int i = 0; i < encoded_data.size(); i++){
        if(i == 0){
            decoded_data[i] = diff_state ^ encoded_data[i];
        }
        else{
            decoded_data[i] = encoded_data[i] ^ encoded_data[i-1];
        }
    }
    diff_state = encoded_data.back();
}

void vec2str(std::vector<bool> &vec, std::string &str) {
    for(int i = 0; i < vec.size(); i++){
        vec[i] ? str.push_back('1') : str.push_back('0');
    }
}

void frame_sync(std::vector<bool> &decoded_data, std::vector<bool> &frame_state, int &group_state, int &sync_flag, std::vector<char> &char_buf, std::vector<bool> &DI_seg, std::vector<bool> &group_type, std::vector<bool> &parity_check_mat, std::vector<bool> &syndrome_vecs){

    int idx = 0; 
    int initial_flag = 0;
    std::string offset_type;
    // intermeidate message info
    std::vector <bool> msg_block;
    msg_block.resize(26,0);
    std::string PTY;

   // brute-force find the first valid message block, normalized to type A
    if (sync_flag) {
        // std::cerr << "*** Processing block {block_count}" << std::endl;
        sync_start(idx, offset_type, decoded_data, parity_check_mat, syndrome_vecs);
        std::copy(decoded_data.begin()+idx, decoded_data.begin()+idx+26, msg_block.begin());
        if(group_state == 0 && offset_type == "A"){
            group_state = 1;
            parse_msg_A(msg_block);
            idx += 26;
            sync_flag = 0;
            initial_flag = 1;
        }
        else{
            std::cerr << "frame error" << std::endl;
            std::cerr << offset_type << "\n";
            sync_flag = 1;
            return;
        }
    }
    else{
        // no need to re-sync
        std::cerr <<  "*** Processing block {block_count}" << std::endl;
        idx = 26 - frame_state.size();

        std::copy(frame_state.begin(), frame_state.end(), msg_block.begin());
        std::copy(decoded_data.begin(), decoded_data.begin()+idx, msg_block.begin()+frame_state.size());

        check_syndrom(offset_type, msg_block.begin(), msg_block.end(), parity_check_mat, syndrome_vecs);
       if(group_state == 0 && offset_type == "A"){
           group_state = 1;}
       else if (group_state == 1 && offset_type == "B"){
           group_state = 2;
           parse_msg_B(msg_block, group_type, PTY, DI_seg);
           if(initial_flag){
               std::cerr << "Program type: " << PTY << std::endl;
               initial_flag = 0;
               return;
           }
       }
       else if (group_state == 2 && (offset_type == "C" || offset_type == "Cp")){
           group_state = 3;
       }
       else if(group_state == 3 && offset_type == "D"){
           group_state = 0;
            std::string group_type_str;
            vec2str(group_type, group_type_str);

           if (std::stoi(group_type_str.c_str(), nullptr, 2) == 0) {
                std::copy(frame_state.begin(), frame_state.end(), msg_block.begin());
                std::copy(decoded_data.begin(), decoded_data.begin() + idx, msg_block.begin()+frame_state.size());
                std::string DI;
                vec2str(DI_seg, DI);
                int char_addr = std::stoi(DI, nullptr, 2);

                parse_msg_D(msg_block,char_buf[char_addr*2],char_buf[char_addr*2+1]);

                int flag = 0;
                for (int i = 0; i < char_buf.size(); i++) {
                    if (char_buf[i] == -1) {
                        flag = 1;
                        break;
                    }
                }

                if (!flag) {
                    std::cerr <<  "Program service:" << std::endl;
                    for(int k = 0; k < char_buf.size(); k++){
                        std::cerr <<  char_buf[k] << std::endl;
                        char_buf[k] = -1;
                    }
                    std::cout<<"\n"<<std::endl;
                }

           }
           
       }
       else{
        //    std::cerr << "frame error" << std::endl;
            std::cerr <<  "Program service:" << std::endl;
            for(int k = 0; k < char_buf.size(); k++){
                std::cerr <<  char_buf[k] << std::endl;
                char_buf[k] = -1;
            }
            std::cout<<"\n"<<std::endl;
            sync_flag = 1;
       }
   }
   // rest of message blocks
   while (idx+26 < decoded_data.size()){
       std::copy(decoded_data.begin()+idx, decoded_data.begin()+idx+26, msg_block.begin());
       check_syndrom(offset_type, msg_block.begin(), msg_block.end(), parity_check_mat, syndrome_vecs);
       if(group_state == 0 && offset_type == "A"){
           group_state = 1;
       }
       else if(group_state == 1 && offset_type == "B"){
           group_state = 2;
           parse_msg_B(msg_block, group_type, PTY, DI_seg);
           if(initial_flag){
               std::cerr << "Program type: " << PTY << std::endl;
               initial_flag = 0;
           }
       }
       else if(group_state == 2 && (offset_type == "C" || offset_type == "Cp"))
           group_state = 3;

       else if(group_state == 3 && offset_type == "D") {
           group_state = 0;
           std::string group_type_str;
            vec2str(group_type, group_type_str);

           if (std::stoi(group_type_str.c_str(), nullptr, 2) == 0) {
                std::copy(frame_state.begin(), frame_state.end(), msg_block.begin());
                std::copy(decoded_data.begin(), decoded_data.begin() + idx, msg_block.begin()+frame_state.size());
                std::string DI;
                vec2str(DI_seg, DI);
                int char_addr = std::stoi(DI, nullptr, 2);

                parse_msg_D(msg_block,char_buf[char_addr*2],char_buf[char_addr*2+1]);

                int flag = 0;
                for (int i = 0; i < char_buf.size(); i++) {
                    if (char_buf[i] == -1) {
                        flag = 1;
                        break;
                    }
                }

                if (!flag) {
                    std::cerr <<  "Program service:" << std::endl;
                    for(int k = 0; k < char_buf.size(); k++){
                        std::cerr <<  char_buf[k] << std::endl;
                        char_buf[k] = -1;
                    }
                    std::cout<<"\n"<<std::endl;
                }
           }
       }
       else{
        //   std::cerr << "frame error"<< std::endl;
          std::cerr <<  "Program service:" << std::endl;
            for(int k = 0; k < char_buf.size(); k++){
                std::cerr <<  char_buf[k] << std::endl;
                char_buf[k] = -1;
            }
            std::cout<<"\n"<<std::endl;
            sync_flag = 1;
      }

      idx += 26;
    }

    frame_state.resize(decoded_data.size()-idx, 0);
    std::copy(decoded_data.begin() + idx, decoded_data.end(), frame_state.begin());
}

void fmPll_q(std::vector<float> &pllIn, std::vector<float> &ncoOut, std::vector<float> &ncoOut_q, float freq, float Fs, state_pll_q &state, float ncoScale, float phaseAdjust, float normBandwidth)
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
    ncoOut_q.resize(pllIn.size()+1,0);

	// initialize internal state
  	float integrator = state.integrator;
	float phaseEst = state.phaseEst;
	float feedbackI = state.feedbackI;
	float feedbackQ = state.feedbackQ;
  	ncoOut[0] = state.ncoOut_0;
    ncoOut_q[0] = state.ncoOut_q_0;
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
        ncoOut_q[k+1] = std::sin(trigArg*ncoScale + phaseAdjust);
	}

  	state.integrator = integrator;
	state.phaseEst = phaseEst;
	state.feedbackI = feedbackI;
	state.feedbackQ = feedbackQ;
  	state.ncoOut_0 = ncoOut[ncoOut.size()-1];
    state.ncoOut_q_0 = ncoOut_q[ncoOut_q.size()-1];
	state.trigOffset = trigOffset;

	// for stereo only the in-phase NCO component should be returned
	// for block processing you should also return the state


	// for RDS add also the quadrature NCO component to the output
}
