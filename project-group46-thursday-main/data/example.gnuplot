# example.gnuplut : configuration for plotting (change as needed)

reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 2,1 scale 1.0,1.0  # set two plots for this figure

# time domain
set ylabel 'I Sample value'               # set y-axis label
set xlabel 'Sample #'                   # set x-axis label
set yrange [-1:1]                       # set y plot range
set xrange [0:50]                      # set x plot range
plot '../data/i_samples.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

set ylabel 'Q Sample value'               # set y-axis label
set xlabel 'Sample #'                   # set x-axis label
set yrange [-0.5:0.5]                       # set y plot range
set xrange [0:50]  
plot '../data/q_samples.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

# set ylabel 'input_demod Sample value'               # set y-axis label
# set xlabel 'Sample #'                   # set x-axis label
# set yrange [-0.5:0.5]                       # set y plot range
# set xrange [0:0.5]  
# plot '../data/input_demod_data.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle


# # freq domain (Fourier)
# set ylabel 'Spectrum (Mag)'              # set y-axis label
# set xlabel 'Frequency bin'               # set x-axis label
# set yrange [0.0:0.25]                    # set y plot range
# set xrange [0:255]                       # set x plot range
# plot '../data/pilot_freq.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

# # time domain
# set ylabel 'Sample value'               # set y-axis label
# set xlabel 'Sample #'                   # set x-axis label
# set yrange [-2:2]                       # set y plot range
# set xrange [0:511]                      # set x plot range
# plot '../data/stereo_channel_time.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

# # freq domain (Fourier)
# set ylabel 'Spectrum (Mag)'              # set y-axis label
# set xlabel 'Frequency bin'               # set x-axis label
# set yrange [0.0:0.25]                    # set y plot range
# set xrange [0:255]                       # set x plot range
# plot '../data/stereo_channel_freq.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#008800' notitle

# # freq domain (PSD)
# set ylabel 'Spectrum (dB/Hz)'            # set y-axis label
# set xlabel 'Frequency (KHz)'             # set x-axis label
# set yrange [-80:0]                       # set y plot range
# set xrange [0:119]                       # set x plot range
# # add your own .dat file for PSD as part of the take-home
# plot '../data/demod_psd.dat' using 1:2 with lines lt 1 lw 3 lc rgb '#880000' notitle

# # freq domain (PSD)
# set ylabel 'Spectrum (dB/Hz)'            # set y-axis label
# set xlabel 'Frequency (KHz)'             # set x-axis label
# set yrange [-80:0]                       # set y plot range
# set xrange [0:19]                       # set x plot range
# # add your own .dat file for PSD as part of the take-home
# plot '../data/mono_psd.dat' using 1:2 with lines lt 1 lw 3 lc rgb '#880000' notitle

# # time domain
# set ylabel 'value'               # set y-axis label
# set xlabel 'term'                   # set x-axis label
# set yrange [-0.02:0.08]                       # set y plot range
# set xrange [0:151]                      # set x plot range
# plot '../data/coeff.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#000088' notitle

unset multiplot
