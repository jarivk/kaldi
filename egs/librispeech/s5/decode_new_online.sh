#!/bin/bash

# Copyright 2014  Johns Hopkins University (Author: Daniel Povey)
# Apache 2.0

# Begin configuration section.
stage=0
nj=1
#cmd=run.pl
max_active=1000
#max_active=2000
beam=13.0
lattice_beam=6.0
#acwt=0.083333 # note: only really affects adaptation and pruning (scoring is on
              # lattices).
#acwt=0.058823
#acwt=0.0588
acwt=0.142
per_utt=false
do_endpointing=true
#do_endpointing=false
do_speex_compressing=false
scoring_opts=
skip_scoring=false
# End configuration section.

echo "$0 $@"  # Print the command line for logging

[ -f ./path.sh ] && . ./path.sh; # source the path.
#. parse_options.sh || exit 1;

if [ $# != 3 ]; then
   echo "Usage: $0 [options] <graph-dir> <data-dir> <option>"
   exit 1;
fi


graphdir=$1
data=$2
option=$3
srcdir=$graphdir; # The model directory is one level up from decoding directory.
sdata=$data


for f in $graphdir/HCLG.fst $graphdir/words.txt $data/wav.scp; do
  if [ ! -f $f ]; then
    echo "$0: no such file $f"
    exit 1;
  fi
done

spk2utt_rspecifier="ark:$sdata/spk2utt"

wav_rspecifier="scp:$sdata/wav.scp"
 

if [ $option -eq 1 ]; then
    #--capture-device="plughw:CARD=M,DEV=0"
    #--capture-device="hw:CARD=4,DEV=1"
    ../../../src/onlinedecoderlivecapture/online2-wav-gmm-latgen-faster-live --capture-device="plughw:CARD=M,DEV=0" --do-endpointing=$do_endpointing --endpoint.rule2.max-relative-cost=0.5 --endpoint.rule2.min-trailing-silence=0.5 --endpoint.rule2.min-utterance-length=0.2 --endpoint.rule2.must-contain-nonsilence=true --config=$srcdir/conf/online_decoding.conf --max-active=$max_active --beam=$beam --lattice-beam=$lattice_beam --acoustic-scale=$acwt --word-symbol-table=$graphdir/words.txt $graphdir/HCLG.fst $spk2utt_rspecifier "$wav_rspecifier" "ark,t:one-best.tra"
elif [ $option -eq 2 ];then
    ../../../src/onlinedecoder/online2-wav-gmm-latgen-faster-novad --capture-device="plughw:CARD=M,DEV=0" --do-endpointing=$do_endpointing --endpoint.rule2.max-relative-cost=1.0 --endpoint.rule2.min-trailing-silence=0.5 --endpoint.rule2.min-utterance-length=0.5 --endpoint.rule2.must-contain-nonsilence=true --config=$srcdir/conf/online_decoding.conf --max-active=$max_active --beam=$beam --lattice-beam=$lattice_beam --acoustic-scale=$acwt --word-symbol-table=$graphdir/words.txt $graphdir/HCLG.fst $spk2utt_rspecifier "$wav_rspecifier" "ark,t:one-best.tra"
elif [ $option -eq 3 ];then
    ../../../src/onlinedecoderfile/online2-wav-gmm-latgen-faster-file --do-endpointing=$do_endpointing --endpoint.rule2.max-relative-cost=0.5 --endpoint.rule2.min-trailing-silence=0.5 --endpoint.rule2.min-utterance-length=0.5 --endpoint.rule2.must-contain-nonsilence=true --config=$srcdir/conf/online_decoding.conf --max-active=$max_active --beam=$beam --lattice-beam=$lattice_beam --acoustic-scale=$acwt --word-symbol-table=$graphdir/words.txt $graphdir/HCLG.fst $spk2utt_rspecifier "$wav_rspecifier" "ark,t:one-best.tra"
fi

exit 0;
