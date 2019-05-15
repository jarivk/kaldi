#!/bin/bash
data=data

# base url for downloads.
#dataset librispeech
data_url=www.openslr.org/resources/12
#language files 
lm_url=www.openslr.org/resources/11

. ./cmd.sh
. ./path.sh
. parse_options.sh

# download the data.  Note: we're using the 100 hour setup 
#dev-clean test-clean dev-other test-other is testing data set
#download train-clean-100
#If the data is available locally it will not download again make sure path is same 
local/download_and_untar.sh $data $data_url train-clean-100

# download the LM resources
local/download_lm.sh $lm_url data/local/lm
          
# format the data as Kaldi data directories
# use underscore-separated names in data directories.
local/data_prep.sh $data/LibriSpeech/train-clean-100 $data/train_clean_100
    


#local/prepare_dict.sh --nj 30 --cmd "$train_cmd" data/local/lm data/local/lm data/local/dict_nosp

utils/prepare_lang.sh data/local/dict_nosp "<UNK>" data/local/lang_tmp_nosp data/lang_nosp

local/format_lms.sh --src-dir data/lang_nosp data/local/lm
    



#MFCC feature
./make_mfcc.sh --cmd "$train_cmd" --nj 40 data/train_clean_100 exp/make_mfcc/train-clean-100 mfccdir
#CMVN feature
./compute_cmvn_stats.sh data/train_clean_100 exp/make_mfcc/train_clean_100 mfccdir




# For the monophone stages we select the shortest utterances, which should make it
# easier to align the data from a flat start.

utils/subset_data_dir.sh --shortest data/train_clean_100 2000 data/train_2kshort
utils/subset_data_dir.sh data/train_clean_100 5000 data/train_5k
utils/subset_data_dir.sh data/train_clean_100 10000 data/train_10k
    



# train a monophone system
steps/train_mono.sh --boost-silence 1.25 --nj 20 --cmd "$train_cmd" data/train_2kshort data/lang_nosp exp/mono
    


steps/align_si.sh --boost-silence 1.25 --nj 10 --cmd "$train_cmd" data/train_5k data/lang_nosp exp/mono exp/mono_ali_5k

# train a first delta + delta-delta triphone system on a subset of 5000 utterances
steps/train_deltas.sh --boost-silence 1.25 --cmd "$train_cmd" 2000 10000 data/train_5k data/lang_nosp exp/mono_ali_5k exp/tri1
    


steps/align_si.sh --nj 10 --cmd "$train_cmd" data/train_10k data/lang_nosp exp/tri1 exp/tri1_ali_10k

# train an LDA+MLLT system.
steps/train_lda_mllt.sh --cmd "$train_cmd" --splice-opts "--left-context=3 --right-context=3" 2500 15000 \
                                                                   ##        data/train_10k data/lang_nosp exp/tri1_ali_10k exp/tri2b
    



# align the entire train_clean_100 subset using the tri2b model
steps/align_si.sh --nj 20 --cmd "$train_cmd" \
    data/train_clean_100 data/lang_nosp \
    exp/tri2b exp/tri2b_ali_clean_100

# train another LDA+MLLT+SAT system on the entire 100 hour subset
  steps/train_lda_mllt.sh  --cmd "$train_cmd" 4200 40000 \
                      data/train_clean_100 data/lang_nosp \
                      exp/tri2b_ali_clean_100 exp/tri3b
    



local/download_and_untar.sh $data $data_url train-clean-360

# now add the "clean-360" subset to the mix ...
local/data_prep.sh $data/LibriSpeech/train-clean-360 $data/train_clean_360

#extract features
./make_mfcc.sh --cmd "$train_cmd" --nj 40 data/train_clean_360 exp/make_mfcc/train_clean_360 mfccdir
./compute_cmvn_stats.sh data/train_clean_360 exp/make_mfcc/train_clean_360 mfccdir

# ... and then combine the two sets into a 460 hrs dataset
utils/combine_data.sh data/train_clean_460 data/train_clean_100 data/train_clean_360
      



# align the new, combined set, using the tri3b model
steps/align_si.sh --nj 40 --cmd "$train_cmd" \
                       data/train_clean_460 data/lang_nosp exp/tri3b exp/tri3b_ali_clean_460

# create a larger model, trained on the 460 hours of data.
steps/train_lda_mllt.sh  --cmd "$train_cmd" 5000 100000 \
                      data/train_clean_460 data/lang_nosp exp/tri3b_ali_clean_460 exp/tri4b
    



# prepare the remaining 500 hours of data
local/download_and_untar.sh $data $data_url train-other-500

# prepare the 500 hrs subset.
local/data_prep.sh $data/LibriSpeech/train-other-500 data/train_other_500

#feature extraction
steps/make_mfcc.sh --cmd "$train_cmd" --nj 40 data/train_other_500 exp/make_mfcc/train_other_500 mfccdir
steps/compute_cmvn_stats.sh data/train_other_500 exp/make_mfcc/train_other_500 mfccdir

# combine all the data
utils/combine_data.sh data/train_960 data/train_clean_460 data/train_other_500
    



steps/align_si.sh --nj 40 --cmd "$train_cmd" data/train_960 data/lang_nosp exp/tri4b exp/tri4b_ali_960

# train a LDA+MLTT model on the 960 hour mixed data
steps/train_lda_mllt.sh --cmd "$train_cmd" 7000 150000 data/train_960 data/lang_nosp exp/tri4b_ali_960 exp/tri960
    
