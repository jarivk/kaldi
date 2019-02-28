#!/bin/bash
# Copyright 2016 Api.ai (Author: Ilya Platonov)
# Apache 2.0

# This script demonstrates kaldi decoding using pretrained model. It will decode list of wav files.
#
# IMPORTANT: wav files must be in 16kHz, 16 bit little-endian format.
#
# This script tries to follow with what other scripts are doing in terms of directory structures and data handling.
#
# Use ./download-model.sh script to download asr model
# See https://github.com/api-ai/api-ai-english-asr-model for details about a model and how to use it.


if [ $# -lt 2 ];then
        echo " <online model>" "<1.live with vad 2.live with no vad 3.file input> <incase of file input provide wav file>"
        exit 1
fi
. ./path.sh
DATA_DIR="test-corpus"
ONLINE_DIR=$1
option=$2
file=${3:-1}


#IMPORTANT: wav files must be in 16kHz, 16 bit little-endian format.
if [ $option -eq 3 ]; then 
	if [ "$file" == "1" ];then
		echo "provide a wav file"
		exit 1;
	fi
	./create-corpus.sh $DATA_DIR $file || exit 1;
fi



	
#online2-wav-gmm-latgen-faster --capture-device="plughw:CARD=M,DEV=0" --do-endpointing=$do_endpointing --endpoint.rule2.min-utterance-length=0.5 --endpoint.rule2.max-relative-cost=1.5 --config=$srcdir/conf/online_decoding.conf --max-active=$max_active --beam=$beam --lattice-beam=$lattice_beam --acoustic-scale=$acwt --word-symbol-table=$graphdir/words.txt $graphdir/HCLG.fst $spk2utt_rspecifier "$wav_rspecifier" "ark,t:one-best.tra"

 sh ./decode_new_online.sh $ONLINE_DIR $DATA_DIR $option

#./rescore.sh
if [ $option -eq 3 ];then
	#echo -n "keywords: "

	##for number in $(cat one-best.tra | sed -r 's/[^0-9 ]+//g'); do word+="$(cat $ONLINE_DIR/words.txt | grep " $number$" | sed -r 's/[0-9]+//g;s/<UNK>//g;s/[ ]+$//g' | tr '\n' ' ' | tr -s " ")"; word+=" "; done; echo
	for number in $(cat one-best.tra | sed -r 's/[^0-9 ]+//g'); do word+="$(cat $ONLINE_DIR/words.txt | grep " $number$" | sed -r 's/[0-9]+//g;s///g;s/[ ]+$//g' | tr '\n' ' ' | tr -s " ")"; word+=" "; done; echo
	echo $word >>resultfile.txt
	rm one-best.tra
fi

