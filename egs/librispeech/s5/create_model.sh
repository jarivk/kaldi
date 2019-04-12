#!/bin/sh
if [ $# -ne 8 ];then
        echo "./create_model.sh <list_of_command> <dictionary_dir> <grammar_order> <output_grammar_in_.gz> <lang_dir_name> <model_dir> <train_data_dir> <Output_model_dir>"
        exit 1
fi
command_list=$1
dictionary_dir=$2
order=$3
grammar=$4
lang_dir=$5
model_dir=$6
train_data=$7
output_model=$8

command=$1
command_list="$(echo $command_list | cut -d'.' -f 1)"
#create a directory contains lexicon and words.txt of command
if [ ! -d lexicon_$command_list ];then
	mkdir lexicon_$command_list
fi

if [ -f lexicon_command ];then
	rm -r lexicon_$command_list/words.txt
fi

cat $command | grep -o -E '\w+' | sort -u -f > lexicon_$command_list/words.txt

if [ -f lexicon_$command_list/lexicon.txt ];then
	rm -r lexicon_$command_list/lexicon.txt
fi
python generate_lxicon.py lexicon_$command_list/words.txt
mv lexicon.txt lexicon_$command_list/

#creating dictionery
if [ ! -d $dictionary_dir ];then
	mkdir $dictionary_dir
fi

for file in dictionary/*;do
	cp $file $dictionary_dir
done
cp lexicon_$command_list/lexicon.txt $dictionary_dir

echo '<UNK> SPN' | cat - $dictionary_dir/lexicon.txt > $dictionary_dir/temp && mv $dictionary_dir/temp $dictionary_dir/lexicon.txt
echo '<SPOKEN_NOISE> SPN' | cat - $dictionary_dir/lexicon.txt > $dictionary_dir/temp && mv $dictionary_dir/temp $dictionary_dir/lexicon.txt
echo '!SIL SIL' | cat - $dictionary_dir/lexicon.txt > $dictionary_dir/temp && mv $dictionary_dir/temp $dictionary_dir/lexicon.txt

#create arpaLM.. set grammer order
../../../tools/srilm/bin/i686-m64/ngram-count -text $command -order $order -limit-vocab -vocab lexicon_$command_list/words.txt -unk   -map-unk "<UNK>" -wbdiscount -interpolate -lm lexicon_$command_list/srilm_$grammar

rm $dictionary_dir/lexiconp.txt
#prepare language
utils/prepare_lang.sh $dictionary_dir "<UNK>" lang_tmp $lang_dir

#format LM
utils/format_lm.sh $lang_dir lexicon_$command_list/srilm_$grammar $dictionary_dir/lexicon.txt $lang_dir

#generate decoding graph
utils/mkgraph.sh $lang_dir $model_dir $model_dir/graph_$lang_dir

#prepare for online decoding
steps/online/prepare_online_decoding.sh $train_data $lang_dir $model_dir online_$lang_dir

#creating final model
mkdir $output_model
cp -r online_$lang_dir/* $output_model
cp $model_dir/graph_$lang_dir/HCLG.fst $output_model
cp $model_dir/graph_$lang_dir/words.txt $output_model
