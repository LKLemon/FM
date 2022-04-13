from bitstring import BitStream
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import sys

from fmPll import fmPll
from fmSupportLib import *
from fmRRC import impulseResponseRootRaisedCosine
from RDS_utility import sync_start, check_syndrom, parse_msg

# front end settings
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

# RDS settings
RDS_channel_Fb = 54e3
RDS_channel_Fe = 60e3
RDS_recover_Fb = 113.5e3
RDS_recover_Fe = 114.5e3
pll_RDS_f = 114e3
RDS_Fs =  rf_Fs / rf_decim
RDS_taps = 151
RDS_demod_Fc = 3e3

SPS = 19 #default mode 0
FsIN = 240e3
FsOUT = 2375 * SPS
if SPS == 19: #Mode 0
    GCD = 125 # gdc of 240K & (19 * 2375 = 45125)
    U = int(FsOUT/GCD) # 361
    D = int(FsIN/GCD) # 1920
elif SPS == 42: #Mode 2
    GCD = 750 # gdc of 240K & (42 * 2375 = 99750)
    U = int(FsOUT/GCD) # 133
    D = int(FsIN/GCD) # 32

# Bandpass
def bandpass(Fb, Fe, Fs, Ntaps):
    h = np.zeros(Ntaps)
    norm_center = ((Fe + Fb)/2)/(Fs/2)
    norm_pass = (Fe - Fb)/(Fs/2)
    for i in range(0,Ntaps):
        if i == ((Ntaps - 1)/2):
            h[i] = norm_pass
        else:
            h[i] = norm_pass*(math.sin(math.pi*(norm_pass/2)*(i-(Ntaps-1)/2))/(math.pi*(norm_pass/2)*(i-(Ntaps-1)/2)))
        h[i] = h[i]*math.cos(i*math.pi*norm_center)
        h[i] = h[i]*math.pow(math.sin(i*math.pi/Ntaps),2)

    return h

#Allpass
def allPass(input_block, state_block):
    output_block = np.concatenate((state_block, input_block[:-len(state_block)]))
    state_block = input_block[-len(state_block):]
    return output_block, state_block

#Rational Resampler
def rationalResampler(demod_RDS, U, D, state_rational, coeff):
    #Upsample
    # print("upsample")
    output = np.zeros(int(len(demod_RDS) * U / D))
    for i in range(len(output)):
        phase = (D*i) % U
        for j in range(phase, len(state_rational), U):
            x_idx = int((i*D - j) / U)
            if (x_idx >= 0):
                output[i] += coeff[j] * demod_RDS[x_idx]
            else:
                output[i] += coeff[i] * state_rational[x_idx]
    state_rational = demod_RDS[-len(state_rational):]
    return output, state_rational

# Assume strong signal
def CDR(preCDR, SPS, start, state, resync_flag,block_count):
    indicator = []
    output = []

    indicator += state
    state.clear()

    for i in range(start, len(preCDR), SPS):
        if preCDR[i] > 0:
            indicator.append(1)
        else:
            indicator.append(-1)
    # first block or sync error happened, need to detect HH or LL
    if (block_count == 0 or resync_flag[0]):
        for i in range(0, (len(indicator)-1), 2):
            # have pair of HH or LL, move back 1 symbol
            if(indicator[i] == indicator[i+1]):
                output.clear()
                for j in range(1, (len(indicator)-1), 2):
                    if indicator[j] == 1:
                        output.append(1)
                    else:
                        output.append(0)
                # if have left symbol, pair it with next block 
                if (len(indicator)-1)%2 != 0:
                    state.append(indicator[len(indicator)-1])
                return output, state, start
            # normal pair
            if indicator[i] == 1:
                output.append(1)
            else:
                output.append(0)
        # if have left symbol, pair it with next block 
        if len(indicator)%2 != 0:
            state.append(indicator[len(indicator)-1])

    # normal block, pair directly
    else:
        for i in range(0, (len(indicator)-1), 2):
            if indicator[i] == 1:
                output.append(1)
            else:
                output.append(0)

        if len(indicator)%2 != 0:
            state.append(indicator[len(indicator)-1])

    return output, state, start

def diff_decode(encoded_data, diff_state):
    decoded_data = BitStream(len(encoded_data))
    for i in range(len(encoded_data)):
        if i == 0:
            decoded_data[i] = diff_state ^ encoded_data[i]
        else:
            decoded_data[i] = encoded_data[i] ^ encoded_data[i-1]
    diff_state = encoded_data[-1]
    # print("diff_decode: ",len(decoded_data))
    return decoded_data, diff_state

def frame_sync(decoded_data, frame_state, group_state, sync_flag, char_buf, DI_seg, group_type):
    # index of start of the next block to be checked
    idx = 0
    # whether the first group after sync
    initial_flag = 0
    # print("frame: ",len(decoded_data))
    # brute-force find the first valid message block, normalized to type A
    if sync_flag[0]:
        # print(f"*** Parsing block {block_count} ***")
        idx, offset_type = sync_start(decoded_data)
        if group_state == 0 and offset_type == "A":
            group_state = 1
            # print(" decode_data idx: ", idx)
            print("syncing to block A")
            parse_msg(decoded_data[idx:idx+26], "A")
            idx += 26
            sync_flag[0] = 0
            initial_flag = 1
        else:
            # print(f"FRAME ERROR, got {offset_type}")
            sync_flag[0] = 1
            return [], 0, char_buf, DI_seg, group_type
    else:
        # no need to re-sync
        # print(f"*** Parsing block {block_count} ***")
        idx = 26 - len(frame_state)
        # print(idx)
        offset_type = check_syndrom(frame_state + decoded_data[:idx])
        # print(offset_type)
        if group_state == 0 and offset_type == "A":
            # print("block type:", offset_type)
            group_state = 1
        elif group_state == 1 and offset_type == "B":
            # print("block type:", offset_type)
            group_state = 2
            group_type, PTY, DI_seg = parse_msg(frame_state + decoded_data[:idx], "B")
            # print("group type: ",group_type.uint)
            if initial_flag:
                print(f"Program type: {PTY}")
                initial_flag = 0
        elif group_state == 2 and offset_type in ["C", "Cp"]:
            # print("block type:", offset_type)
            group_state = 3
        elif group_state == 3 and offset_type == "D":
            # print("block type:", offset_type)
            group_state = 0
            if group_type.uint == 0:
                char_buf[DI_seg.uint*2],char_buf[DI_seg.uint*2+1] = parse_msg(frame_state + decoded_data[:idx], "D")
                if (-1 not in char_buf):
                    print("Program service:", end="")
                    for ch in char_buf:
                        print(ch, end="")
                    char_buf = [-1]*8
                    print()
        else:
            # print(f"FRAME ERROR, got {offset_type}")
            print("Program service:", end="")
            for ch in char_buf:
                print(ch, end="")
            print()
            char_buf = [-1]*8
            sync_flag[0] = 1
            return [], 0, char_buf, DI_seg, group_type
    # rest of message blocks
    while (idx+26 < len(decoded_data)):
        offset_type = check_syndrom(decoded_data[idx:idx+26])
        if group_state == 0 and offset_type == "A":
            # print("block type:", offset_type)
            group_state = 1
        elif group_state == 1 and offset_type == "B":
            # print("block type:", offset_type)
            group_state = 2
            group_type, PTY, DI_seg = parse_msg(decoded_data[idx:idx+26], "B")
            # print("group type: ", group_type.uint)
            if initial_flag:
                print(f"Program type: {PTY}")
                initial_flag = 0
        elif group_state == 2 and offset_type in ["C", "Cp"]:
            # print("block type:", offset_type)
            group_state = 3
        elif group_state == 3 and offset_type == "D":
            # print("block type:", offset_type)
            group_state = 0
            if group_type.uint == 0:
                char_buf[DI_seg.uint*2],char_buf[DI_seg.uint*2+1] = parse_msg(decoded_data[idx:idx+26], "D")
                if (-1 not in char_buf):
                    print("Program service:", end="")
                    for ch in char_buf:
                        print(ch, end="")
                    char_buf = [-1]*8
                    print()
        else:
            # print(offset_type)
            # print(f"FRAME ERROR, got {offset_type}")
            print("Program service:", end="")
            for ch in char_buf:
                print(ch, end="")
            print()
            char_buf = [-1]*8
            sync_flag[0] = 1
            return [], 0, char_buf, DI_seg, group_type
        idx += 26

    frame_state = decoded_data[idx:]

    return frame_state, group_state, char_buf, DI_seg, group_type

def find_sampling_pos(data_block, SPS):
    offset = SPS//2
    mid = 0.0
    i = 0
    j = 0
    array = np.zeros(SPS)
    for i in range(SPS):
        error = 0.0
        j = i
        while j+SPS+1 < len(data_block):
            chunk = data_block[j:j+SPS+1]
            mid = chunk[offset]
            if(chunk[0] * chunk[-1] >= 0):
                # HL/LH mid point should close to the neiboring samples
                error += abs(mid - (chunk[0] + chunk[-1])/2)
            else:
                # HH/LL mid point should close to the 0
                error += abs(mid - 0)

            j += SPS
        array[i] = error

    pos = np.argmin(array)

    return pos

if __name__ == "__main__":

    # read the raw IQ data from the recorded file
    # IQ data is assumed to be in 8-bits unsigned (and interleaved)
    # in_fname = "../data/samples3.raw"

    in_fname = "../data/samples8.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
    # IQ data is normalized between -1 and +1 in 32-bit float format
    iq_data = (np.float32(raw_data) - 128.0)/128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

    # select a block_size that is a multiple of KB
    # and a multiple of decimation factors
    block_size = 1024 * rf_decim * 50 * 2
    block_size = int(block_size)
    # print(block_size)
    block_count = 0
    # print(len(iq_data))

    # global facilities
    # FRONT-END
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_phase = 0
    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

    #RDS ******************************************************************
    state_RDS_pll = {
        "integrator" : 0.0,
        "phaseEst" : 0.0,
        "feedbackI" : 1.0,
        "feedbackQ" : 0.0,
        "ncoOut_0" : 1.0,
        "ncoOut_q_0" : 0.0,
        "trigOffset" : 0
    }
    state_RDS = np.zeros(RDS_taps-1)
    state_allpass = np.zeros((RDS_taps-1)//2)
    state_recover = np.zeros(RDS_taps-1)
    state_rational = np.zeros(int(U*RDS_taps)-1)
    state_rational_q = np.zeros(int(U*RDS_taps)-1)
    state_RRC = np.zeros(RDS_taps-1)
    state_RRC_q = np.zeros(RDS_taps-1)

    RDS_channel_coeff = bandpass(RDS_channel_Fb, RDS_channel_Fe,RDS_Fs, RDS_taps)
    RDS_recover_coeff = bandpass(RDS_recover_Fb, RDS_recover_Fe,RDS_Fs, RDS_taps)
    # coeff for RDS demod
    rational_coeff = U * signal.firwin(int(RDS_taps*U), RDS_demod_Fc/(RDS_Fs*U/2), window=('hann'))
    RRC_coeff = impulseResponseRootRaisedCosine(FsOUT, RDS_taps)

    # data processing 
    state_CDR = []
    diff_state = 0
    frame_state = []
    group_state = 0
    start = 0
    sync_flag = [1]
    # intermeidate message info
    group_type = BitStream(4)
    char_buf = [-1]*8
    DI_seg = BitStream(2)
    # print(len(iq_data))
    #**********************************************************************

    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw IQ file
    while (block_count+1)*block_size < len(iq_data):
        # if block_count == 50:
        #     break
        # if you wish to have shorter runtimes while troubleshooting
        # you can control the above loop exit condition as you see fit
        # print('Processing block ' + str(block_count))

        # Front-end Processing
        # filter to extract the FM channel (I samples are even, Q samples are odd)
        # print("i_filt")
        i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
                zi=state_i_lpf_100k)
        # print("q_filt")
        q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
                zi=state_q_lpf_100k)
        # downsample the I/Q data from the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]
        # FM demodulator
        # print("rf demod")
        fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
        # print("//////////////////////////////////")
        #RDS**********************************************************************
        #RDS Channel Extraction
        # print("channel extraction")
        RDS_channel_data, state_RDS = signal.lfilter(RDS_channel_coeff, 1.0, fm_demod, zi = state_RDS)

        #All pass filter to create delay
        # print("all pass")
        RDS_allpass_data, state_allpass = allPass(RDS_channel_data, state_allpass)

        #Carrier Recovery
        # print("carrier recover")
        RDS_SN_data = RDS_channel_data * RDS_channel_data #Squaring Nonlinearity
        RDS_recover_data, state_recover = signal.lfilter(RDS_recover_coeff, 1.0, RDS_SN_data, zi = state_recover) #114KHz filtering
        RDS_carrier, RDS_carrier_q, state_RDS_pll = fmPll(RDS_recover_data, pll_RDS_f, RDS_Fs, ncoScale=0.5,phaseAdjust = 1.18, normBandwidth = 0.001, state=state_RDS_pll)
        # print("after pll")
        # print(f"i_sample: {len(RDS_carrier)}")
        # print(f"q_sample: {len(RDS_carrier_q)}")
        # if (block_count == 0):
        #     plt.plot(RDS_carrier[:20])
        #     plt.plot(RDS_carrier_q[:20])
        #     plt.show()
        #     sys.exit()

        #RDS Demodulatiom
        #Mixer
        RDS_carrier = RDS_carrier[:-1]
        RDS_mixed_data = 2 * RDS_carrier * RDS_allpass_data # pointwise multiplication
        #Rational sampler
        # print("resampler")
        preRRC, state_rational = rationalResampler(RDS_mixed_data, U, D, state_rational, rational_coeff)
        #Root Raised Cosine Filter
        # print("RRC")
        preCDR, state_RRC = signal.lfilter(RRC_coeff, 1.0, preRRC, zi = state_RRC)

        #########################################################################################
        # Mixer
        RDS_carrier_q = np.array(RDS_carrier_q[:-1])
        RDS_mixed_data_q = 2 * RDS_carrier_q * RDS_allpass_data # pointwise multiplication
        # Rational sampler
        preRRC_q, state_rational_q = rationalResampler(RDS_mixed_data_q, U, D, state_rational_q, rational_coeff)
        #Root Raised Cosine Filter
        preCDR_q, state_RRC_q = signal.lfilter(RRC_coeff, 1.0, preRRC_q, zi = state_RRC_q)

        ########################################################################################
        # if block_count == 2:
        #     pos = find_sampling_pos(preCDR, SPS)
        #     recovered_data = [[],[]]
        #     for i in range(pos, len(preCDR), SPS):
        #         recovered_data[0].append(preCDR[i])
        #         recovered_data[1].append(preCDR_q[i])
        #     plt.scatter(recovered_data[0],recovered_data[1],s=10)
        #     plt.xlim(-1.0,1.0)
        #     plt.ylim(-1.0,1.0)
        #     plt.show()
        #     sys.exit()

        #CDR
        # print("CDR")
        start = find_sampling_pos(preCDR,SPS)
        CDR_data, state_CDR, start = CDR(preCDR, SPS, start, state_CDR, sync_flag, block_count)

        # RDS Data Processing
        decoded_data, diff_state = diff_decode(CDR_data, diff_state)
        
        #frame_sync
        frame_state, group_state, char_buf, DI_seg, group_type = frame_sync(decoded_data, frame_state, group_state, block_count, sync_flag, char_buf, DI_seg, group_type)
        # # # if (block_count == 10):
        # #     sys.exit()
        block_count += 1
    # print('Finished processing all the blocks from the recorded I/Q samples')
