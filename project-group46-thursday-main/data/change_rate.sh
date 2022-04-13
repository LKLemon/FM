#! /bin/bash

# Usage: change directory to data and run "./change_rate.sh"
for f in *.raw
do
regex="samples[0-9].raw"
if [[ $f =~ $regex ]]; then
python3 ../model/fmRateChange.py $f 6 0 # mode 1
python3 ../model/fmRateChange.py $f 4 0 # mode 3
fi
done
