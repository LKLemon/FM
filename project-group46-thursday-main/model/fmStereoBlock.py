import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import sys

from fmPll import fmPll
from fmSupportLib import *

# front end settings
rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

# mono path settings
mono_taps = 151
mono_Fc = 16e3
mono_decim = 5
mono_Fs =  rf_Fs / rf_decim
audio_Fs = 48e3

# stereo path settings
stereo_Fs =  rf_Fs / rf_decim
stereo_taps = 151
pll_f = 19e3
pilot_Fb = 18.5e3
pilot_Fe = 19.5e3
stereo_Fb = 22e3
stereo_Fe = 54e3

# Bandpass
def bandpass(Fb, Fe, Fs, Ntaps):
    h = np.empty(Ntaps)
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
# def write_audio_data()
# if (audio_left.size() != audio_right.size()) {
if __name__ == "__main__":

    # read the raw IQ data from the recorded file
    # IQ data is assumed to be in 8-bits unsigned (and interleaved)
    in_fname = "../data/stereo_l0_r9.raw"
    # in_fname = "/home/pi/samples/samples8.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
    # IQ data is normalized between -1 and +1 in 32-bit float format
    iq_data = (np.float32(raw_data) - 128.0)/128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

    # set up the subfigures for plotting
    # subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
    # plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
    # fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
    # fig.subplots_adjust(hspace = .6)

    # select a block_size that is a multiple of KB
    # and a multiple of decimation factors
    block_size = 1024 * rf_decim * mono_decim * 2
    block_count = 0

    # global facilities
    # FRONT-END
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_phase = 0
    # state_I_Q = [0.0, 0.0]
    # coefficients for the front-end low-pass filter
    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

    # MONO CHANNEL
    mono_data = np.array([])    # output of mono path
    # state for mono channel filter
    state_mono = np.zeros(mono_taps-1)
    # mono lowpass filter
    mono_coeff = signal.firwin(mono_taps, mono_Fc/(mono_Fs/2), window=('hann'))

    # STEREO CHANNEL
    stereo_data = np.array([]) # all the processed stereo data
    # state for stereo channel
    state_pll = {
        "integrator" : 0.0,
        "phaseEst" : 0.0,
        "feedbackI" : 1.0,
        "feedbackQ" : 0.0,
        "ncoOut_0" : 1.0,
        "trigOffset" : 0
    }
    state_stereo = np.zeros(stereo_taps-1)
    state_pilot = np.zeros(stereo_taps-1)
    state_stereo_filt = np.zeros(stereo_taps-1)
    # bandpass filter for stereo channel
    # pilot_coeff = signal.firwin(stereo_taps, [pilot_Fb/(stereo_Fs/2),pilot_Fe/(stereo_Fs/2)] , window=('hann'), pass_zero=False)
    # stereo_coeff = signal.firwin(stereo_taps, [stereo_Fb/(stereo_Fs/2),stereo_Fe/(stereo_Fs/2)] , window=('hann'), pass_zero=False)
    pilot_coeff = bandpass(pilot_Fb, pilot_Fe,stereo_Fs, stereo_taps)
    stereo_coeff = bandpass(stereo_Fb, stereo_Fe,stereo_Fs, stereo_taps)
    stereo_filt_coeff = signal.firwin(stereo_taps, mono_Fc/(stereo_Fs/2), window=('hann'))

    # COMBINED AUDIO
    # audio buffer that stores all the audio blocks
    left_channel = np.array([])
    right_channel = np.array([])

    test_carrier = []

    # if the number of samples in the last block is less than the block size
    # it is fine to ignore the last few samples from the raw IQ file
    while (block_count+1)*block_size < len(iq_data):
        # if block_count == 50:
        #     break
        # if you wish to have shorter runtimes while troubleshooting
        # you can control the above loop exit condition as you see fit
        print('Processing block ' + str(block_count))


        # Front-end Processing
        # filter to extract the FM channel (I samples are even, Q samples are odd)
        i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
                zi=state_i_lpf_100k)
        q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
                iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
                zi=state_q_lpf_100k)
        # downsample the I/Q data from the FM channel
        i_ds = i_filt[::rf_decim]
        q_ds = q_filt[::rf_decim]
        # FM demodulator
        fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)


        # MONO PATH
        # extract the mono audio data through filtering
        mono_filt, state_mono = signal.lfilter(mono_coeff, 1.0, fm_demod, zi = state_mono)
        # downsample audio data
        mono_block = mono_filt[::mono_decim]
        # concatenate the most recently processed audio_block
        # to the previous blocks stored already in audio_data
        mono_data = np.concatenate((mono_data, mono_block))


        # STEREO PATH
        # Carrier Recovery
        pilot, state_pilot= conv(pilot_coeff, fm_demod, state_pilot)
        stereo_carrier, state_pll = fmPll(pilot, pll_f, stereo_Fs, ncoScale=2.0, state=state_pll)
        
        # Stereo Channel Extraction
        stereo_channel_data, state_stereo = conv(stereo_coeff, fm_demod, state_stereo)
        # Stereo Processing
        # Mixer
        stereo_carrier = np.array(stereo_carrier[:-1])
        # if (block_count in [3,4,5,6] ):
        #     test_carrier.extend(stereo_carrier)
        # if block_count == 6:
        #     plt.plot(test_carrier)
        #     plt.show()
        #     sys.exit()
        stereo_channel_data = np.array(stereo_channel_data)
        mixed_data = 2 * stereo_carrier * stereo_channel_data # pointwise multiplication
        # Filtering and rate conversion
        stereo_filt, state_stereo_filt = conv(stereo_filt_coeff, mixed_data, state_stereo_filt)
        stereo_block = stereo_filt[::mono_decim]
        stereo_data = np.concatenate((stereo_data, stereo_block))

        # to save runtime select the range of blocks to log data
        # this includes both saving binary files as well plotting PSD
        # below we assume we want to plot for graphs for blocks 10 and 11
        # if block_count >= 10 and block_count < 12:

        #     # plot PSD of selected block after FM demodulation
        #     ax0.clear()
        #     fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
        #             'Demodulated FM (block ' + str(block_count) + ')')
        #     # output binary file name (where samples are written from Python)
        #     fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
        #     # create binary file where each sample is a 32-bit float
        #     fm_demod.astype('float32').tofile(fm_demod_fname)

        #     fmPlotPSD(ax1, mono_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Mono')
        #     # plot PSD of selected block after extracting mono audio
        #     # ... change as needed
        #     #
        #     # plot PSD of selected block after downsampling mono audio
        #     # ... change as needed
        #     fmPlotPSD(ax2, mono_data, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')

        #     # save figure to file
        #     fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

        block_count += 1

    print('Finished processing all the blocks from the recorded I/Q samples')

    # write mono data to file
    out_fname = "../data/fmMonoBlock.wav"
    # out_fname = "/home/pi/samples/fmMonoBlock.wav"
    wavfile.write(out_fname, int(audio_Fs), np.int16((mono_data/2)*32767))
    print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

    # write combined data to file
    # refer to lab2 block processing and the fact that L = M + S, R = M - S
    left_channel = mono_data + stereo_data
    right_channel = mono_data - stereo_data
    out_fname1 = "../data/fmStereoBlock.wav"
    # out_fname1 = "/home/pi/samples/fmMonoStereoBlock.wav"

    if (len(left_channel)!= len(right_channel)):
        print("audio size messed up with size. ")
    else:
        stereo_data = np.vstack((left_channel, right_channel))
        stereo_data = stereo_data.transpose()
        # Produce an audio file that contains stereo sound
        wavfile.write(out_fname1,int(audio_Fs),np.int16((stereo_data/2)*32767))
        print("Written audio samples to \"" + out_fname1 + "\" in signed 16-bit format")

    # uncomment assuming you wish to show some plots
    # estimatePSD(mono_data,len(mono_data),audio_Fs)
    # plt.show()
