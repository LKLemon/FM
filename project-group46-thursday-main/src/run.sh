#! /bin/bash
if [[ $# -ne 2 ]]; then
    echo "Usages:   ./run.sh sample mode"
    echo "Example:  ./run.sh ../data/samples0.raw 0"
else
    if [[ $2 == 0 || $2 == 1 ]]; then
        cat $1 | ./project $2 | aplay -c 1 -f S16_LE -r 48000
    else
        cat $1 | ./project $2 | aplay -c 1 -f S16_LE -r 44100
    fi
fi
