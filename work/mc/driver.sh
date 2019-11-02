#!/bin/bash

n=10
no_exec='--no_exec'
batch=0

usage() {
  cat <<END
This program generates input files for cmsRun to produce a Monte Carlo calculation
Usage: ${0##*/} [options] generator
Allowed options:
  -n <number>: number of events (default is $n)
  -e or --exec: execute cmsRun on produced files
  -b or --batch: execute cmsRun after all files are produced
END
}

cmsrun() {
  echo "cmsRun $1"
  cmsRun "$1" 2>&1 | tee ${1%.*}.log
  return ${PIPESTATUS[0]}
}

argv=(`getopt -o hn:eb -l help,exec,batch -- "$@"`) || exit $?
eval "set -- ${argv[@]}"
while [[ $1 != -- ]]; do
  case $1 in
    -h|--help) usage; exit;;
    -n) n=$2; shift 2;;
    -e|--exec) no_exec=''; shift;;
    -b|--batch) batch=1; shift;;
  esac
done
shift

if (($# < 1)); then
  echo "${0##*/}: no generator provided" 2>&1
  exit 1
fi

generator=$1

if (($batch)) && [[ -z $no_exec ]]; then
  echo "${0##*/}: --batch and --exec options are exclusive" >&2
  exit 1
fi

set -xe

cmsDriver.py "${generator}_cff.py" \
             --fileout file:"${generator}_GEN-SIM.root" \
             -s GEN,SIM \
             --mc --datatier GEN-SIM \
             --conditions auto:phase1_2017_realistic \
             --eventcontent RAWSIM \
             --era Run2_2017 \
             --python_filename "${generator}_GEN-SIM.py" \
             -n $n \
             $no_exec

cmsDriver.py -s DIGI,L1,DIGI2RAW,HLT \
             --datatier GEN-SIM-RAW \
             --conditions auto:phase1_2017_realistic \
             --eventcontent RAWSIM \
             --era Run2_2017 \
             --filein file:"${generator}_GEN-SIM.root" \
             --fileout file:"${generator}_GEN-SIM-RAW.root" \
             --python_filename "${generator}_GEN-SIM-RAW.py" \
             -n -1 \
             --mc \
             $no_exec

cmsDriver.py -s RAW2DIGI,L1Reco,RECO \
             --datatier RECO \
             --conditions auto:phase1_2017_realistic \
             --eventcontent AODSIM \
             --era Run2_2017 \
             --filein file:"${generator}_GEN-SIM-RAW.root" \
             --fileout file:"${generator}_GEN-SIM-RAW-RECO.root" \
             --python_filename "${generator}_GEN-SIM-RAW-RECO.py" \
             -n -1 \
             --mc \
             $no_exec

set +x

if (($batch)); then
  for f in GEN-SIM{,-RAW{,-RECO}}; do
    cmsrun "${generator}_$f.py" || break
  done
fi
