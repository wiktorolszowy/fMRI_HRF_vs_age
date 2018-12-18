#!/bin/bash

#! Wiktor Olszowy

set -o errexit

#-Folder where output and errors from running 'sbatch' will be saved
mkdir out_err

export part; sbatch --array=1-641 slurm_submit.array.hphi
