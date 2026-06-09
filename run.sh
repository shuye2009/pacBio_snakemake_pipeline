#!/usr/bin/env bash

# Override target at runtime
case $1 in
    cgi)
        TARGET=cgi
        ;;
    centromere)
        TARGET=centromere
        ;;
    repeat)
        TARGET=repeat
        ;;
    enhancer)
        TARGET=enhancer
        ;;
    *)
        echo "Usage: $0 {cgi|centromere|repeat|enhancer}"
        exit 1
        ;;
esac

export TARGET
sbatch --export=TARGET run_pacbio.sh