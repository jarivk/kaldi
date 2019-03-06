// online2bin/online2-wav-gmm-latgen-faster.cc

// Copyright 2014  Johns Hopkins University (author: Daniel Povey)

// See ../../COPYING for clarification regarding multiple authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//  http://www.apache.org/licenses/LICENSE-2.0
//
// THIS CODE IS PROVIDED *AS IS* BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, EITHER EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION ANY IMPLIED
// WARRANTIES OR CONDITIONS OF TITLE, FITNESS FOR A PARTICULAR PURPOSE,
// MERCHANTABLITY OR NON-INFRINGEMENT.
// See the Apache 2 License for the specific language governing permissions and
// limitations under the License.

#include "feat/wave-reader.h"
#include "online2/online-feature-pipeline.h"
#include "online2/online-gmm-decoding.h"
#include "online2/onlinebin-util.h"
#include "online2/online-timing.h"
#include "online2/online-endpoint.h"
#include "fstext/fstext-lib.h"
#include "lat/lattice-functions.h"

#include "fstext/fstext-lib.h"
#include "fstext/kaldi-fst-io.h"

#include<fstream>

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "lat/kaldi-lattice.h"

#include  <alsa/asoundlib.h>

#include<regex>

#include "../webrtcvad/webrtc/common_audio/vad/include/webrtc_vad.h"


namespace kaldi {

//search KW
int searchKW(string text)
{
  //change the regex according to your need
  string pattern ="(\\b\\s*MENU\\s*\\b|\\b\\s*SHOW MAP\\s*\\b|\\b\\s*BECOME DRIVER\\s*\\b|\\b\\s*BOOST\\s*\\b|\\b\\s*BUILD\\s*\\b|\\b\\s*BUILD MEDIVAC\\s*\\b|\\b\\s*BUILD MINE\\s*\\b|\\b\\s*BUILD THREE DROPSHIPS\\s*\\b|\\b\\s*CHANGE CAMERA\\s*\\b|\\b\\s*CHANGE MODE\\s*\\b|\\b\\s*CROUCH\\s*\\b|\\b\\s*DRIVE MOTORBIKE\\s*\\b|\\b\\s*DRIVE VEHICLE\\s*\\b|\\b\\s*ENTER CAR\\s*\\b|\\b\\s*EQUIP MELEE\\s*\\b|\\b\\s*EXPLORE\\s*\\b|\\b\\s*EXTRA SUPPLIES\\s*\\b|\\b\\s*FIRST AID KIT\\s*\\b|\\b\\s*FOUND CITY\\s*\\b|\\b\\s*GET ARMY\\s*\\b|\\b\\s*GET BATTLE UNITS\\s*\\b|\\b\\s*GET IDLE WORKER\\s*\\b|\\b\\s*GRENADE\\s*\\b|\\b\\s*JUMP\\s*\\b|\\b\\s*MAKE BARRACKS\\s*\\b|\\b\\s*MOVE TO PASSENGER\\s*\\b|\\b\\s*NEXT TURN\\s*\\b|\\b\\s*NITRO\\s*\\b|\\b\\s*OPEN CELLPHONE\\s*\\b|\\b\\s*OPEN CHAT\\s*\\b|\\b\\s*OPEN MENU\\s*\\b|\\b\\s*OPTIONS\\s*\\b|\\b\\s*PRIMARY WEAPON\\s*\\b|\\b\\s*PRONE\\s*\\b|\\b\\s*QUICK SAVE\\s*\\b|\\b\\s*RECORD\\s*\\b|\\b\\s*RELOAD\\s*\\b|\\b\\s*REMOVE FOREST\\s*\\b|\\b\\s*SCREENSHOT\\s*\\b|\\b\\s*SET UP ARTILLERY\\s*\\b|\\b\\s*SHOW DIPLOMACY\\s*\\b|\\b\\s*SHOW ECONOMY\\s*\\b|\\b\\s*SHOW HUD\\s*\\b|\\b\\s*SHOW INVENTORY\\s*\\b|\\b\\s*SHOW MILITARY\\s*\\b|\\b\\s*STEALTH MODE\\s*\\b|\\b\\s*SWITCH TO HANDGUN\\s*\\b|\\b\\s*SWITCH TO WEAPON ONE\\s*\\b|\\b\\s*TAKE JEEP\\s*\\b|\\b\\s*TAKE PICTURE\\s*\\b|\\b\\s*TECH TREE\\s*\\b|\\b\\s*TIMEOUT\\s*\\b|\\b\\s*WEAPON TWO\\s*\\b)";
  std::smatch m;
  std::regex r(pattern);
  int count=0;
  while (std::regex_search (text,m,r)) {
    std::cout << m[0] << " ";
    
    text = m.suffix().str();
    count++;
  }
  std::cout<<count<<std::endl;
  if(count>0)
    return 1;
  else
    return 0;
}


//add lattice scale
void LatScale(const std::string rclat,const std::string wclat)
{
  try{
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    bool write_compact = true;
    BaseFloat acoustic_scale = 1.0;
    BaseFloat inv_acoustic_scale = 7.0;
    BaseFloat lm_scale = 1.0;
    BaseFloat acoustic2lm_scale = 0.0;
    BaseFloat lm2acoustic_scale = 0.0;
    std::string lats_rspecifier = rclat,
        lats_wspecifier = wclat;
    int32 n_done = 0;
    KALDI_ASSERT(acoustic_scale == 1.0 || inv_acoustic_scale == 1.0);
    if (inv_acoustic_scale != 1.0)
      acoustic_scale = 1.0 / inv_acoustic_scale;

    std::vector<std::vector<double> > scale(2);
    scale[0].resize(2);
    scale[1].resize(2);
    scale[0][0] = lm_scale;
    scale[0][1] = acoustic2lm_scale;
    scale[1][0] = lm2acoustic_scale;
    scale[1][1] = acoustic_scale;
    if (write_compact) {
      SequentialCompactLatticeReader compact_lattice_reader(lats_rspecifier);

      // Write as compact lattice.
      CompactLatticeWriter compact_lattice_writer(lats_wspecifier);

      for (; !compact_lattice_reader.Done(); compact_lattice_reader.Next()) {
        CompactLattice lat = compact_lattice_reader.Value();
        ScaleLattice(scale, &lat);
        compact_lattice_writer.Write(compact_lattice_reader.Key(), lat);
        n_done++;
      }
    } else {
      SequentialLatticeReader lattice_reader(lats_rspecifier);

      // Write as regular lattice.
      LatticeWriter lattice_writer(lats_wspecifier);

      for (; !lattice_reader.Done(); lattice_reader.Next()) {
        Lattice lat = lattice_reader.Value();
        ScaleLattice(scale, &lat);
        lattice_writer.Write(lattice_reader.Key(), lat);
        n_done++;
      }
    }

    //KALDI_LOG << "Done " << n_done << " lattices.";
  }
  catch(const std::exception &e)
  {
    std::cerr << e.what();
  }
}

//add penalty
void LatPenalty(const std::string rclat, const std::string wclat)
{
  using namespace kaldi;
  typedef kaldi::int64 int64;
  try {


    BaseFloat word_ins_penalty = 0.5;

    std::string lats_rspecifier = rclat,
        lats_wspecifier = wclat;

    SequentialCompactLatticeReader clat_reader(lats_rspecifier);
    CompactLatticeWriter clat_writer(lats_wspecifier); // write as compact.

    int64 n_done = 0;

    for (; !clat_reader.Done(); clat_reader.Next()) {
      CompactLattice clat(clat_reader.Value());
      AddWordInsPenToCompactLattice(word_ins_penalty, &clat);
      clat_writer.Write(clat_reader.Key(), clat);
      n_done++;
    }
    //KALDI_LOG << "Done adding word insertion penalty to " << n_done << " lattices.";
  } catch(const std::exception &e) {
    std::cerr << e.what();
  }

}

//find best path
string LatBestPath(const std::string rclat, const std::string wclat, std::string word_syms_filename, const std::string wali)
{

  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    std::string transcript="";

    BaseFloat acoustic_scale = 0.142;
    BaseFloat lm_scale = 1.0;

    std::string lats_rspecifier = rclat,
        transcriptions_wspecifier = wclat,
        alignments_wspecifier = wali;

    SequentialCompactLatticeReader clat_reader(lats_rspecifier);

    Int32VectorWriter transcriptions_writer(transcriptions_wspecifier);

    Int32VectorWriter alignments_writer(alignments_wspecifier);

    fst::SymbolTable *word_syms = NULL;
    if (word_syms_filename != "")
      if (!(word_syms = fst::SymbolTable::ReadText(word_syms_filename)))
        KALDI_ERR << "Could not read symbol table from file "
                   << word_syms_filename;


    int32 n_done = 0, n_fail = 0;
    int64 n_frame = 0;
    LatticeWeight tot_weight = LatticeWeight::One();

    for (; !clat_reader.Done(); clat_reader.Next()) {
      std::string key = clat_reader.Key();
      CompactLattice clat = clat_reader.Value();
      clat_reader.FreeCurrent();
      fst::ScaleLattice(fst::LatticeScale(lm_scale, acoustic_scale), &clat);
      CompactLattice clat_best_path;
      CompactLatticeShortestPath(clat, &clat_best_path);  // A specialized
      // implementation of shortest-path for CompactLattice.
      Lattice best_path;
      ConvertLattice(clat_best_path, &best_path);
      if (best_path.Start() == fst::kNoStateId) {
        KALDI_WARN << "Best-path failed for key " << key;
        n_fail++;
      } else {
        std::vector<int32> alignment;
        std::vector<int32> words;
        LatticeWeight weight;
        GetLinearSymbolSequence(best_path, &alignment, &words, &weight);
        

        if (transcriptions_wspecifier != "")
          transcriptions_writer.Write(key, words);
        if (alignments_wspecifier != "")
          alignments_writer.Write(key, alignment);
        if (word_syms != NULL) {
          //std::cerr << key << ' ====>';
          for (size_t i = 0; i < words.size(); i++) {
            std::string s = word_syms->Find(words[i]);
            if (s == "")
              KALDI_ERR << "Word-id " << words[i] <<" not in symbol table.";
            transcript += s;
            transcript += " ";
            //std::cerr << s << ' ';

          }
          //std::cerr << '\n';
        }
        n_done++;
        n_frame += alignment.size();
        tot_weight = Times(tot_weight, weight);
      }
    }

    BaseFloat tot_weight_float = tot_weight.Value1() + tot_weight.Value2();
    return transcript;
    delete word_syms;
  } catch(const std::exception &e) {
    std::cerr << e.what();
    return "";
  }
}
}

void LatBest_1_Path(const std::string rclat, const std::string wclat)
{
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::VectorFst;
    using fst::StdArc;

    //Compute best path through lattices and write out as FSTs
    //Note: differs from lattice-nbest with --n=1 because we won't
    //append -1 to the utterance-ids.  Differs from lattice-best-path
    //because output is FST
    //e.g.: lattice-1best --acoustic-scale=0.1 ark:1.lats ark:1best.lats
 
    BaseFloat acoustic_scale = 0.142;
    BaseFloat lm_scale = 1.0;
    BaseFloat word_ins_penalty = 0.0;
    

    std::string lats_rspecifier = rclat,
        lats_wspecifier = wclat;

    SequentialCompactLatticeReader clat_reader(lats_rspecifier);
    
    // Write as compact lattice.
    CompactLatticeWriter compact_1best_writer(lats_wspecifier); 

    int32 n_done = 0, n_err = 0;

    /*if (acoustic_scale == 0.0 || lm_scale == 0.0)
      KALDI_ERR << "Do not use exactly zero acoustic or LM scale (cannot be inverted)";*/
    for (; !clat_reader.Done(); clat_reader.Next()) {
      std::string key = clat_reader.Key();
      CompactLattice clat = clat_reader.Value();
      clat_reader.FreeCurrent();
      fst::ScaleLattice(fst::LatticeScale(lm_scale, acoustic_scale), &clat);
      if (word_ins_penalty > 0.0) {
        AddWordInsPenToCompactLattice(word_ins_penalty, &clat);
      }

      CompactLattice best_path;
      CompactLatticeShortestPath(clat, &best_path);
      
      if (best_path.Start() == fst::kNoStateId) {
        KALDI_WARN << "Possibly empty lattice for utterance-id " << key
                   << "(no output)";
        n_err++;
      } else {
        fst::ScaleLattice(fst::LatticeScale(1.0 / lm_scale, 1.0/acoustic_scale),
                          &best_path);
        if (word_ins_penalty > 0.0) {
          AddWordInsPenToCompactLattice(word_ins_penalty, &clat);
        }
        compact_1best_writer.Write(key, best_path);
        n_done++;
      }
    }
    /*KALDI_LOG << "Done converting " << n_done << " to best path, "
              << n_err << " had errors.";*/
  } catch(const std::exception &e) {
    std::cerr << e.what();
  }
}


void lattice_to_fst(const std::string rclat, const std::string wclat)
{
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    using std::vector;
    BaseFloat acoustic_scale = 0.142;
    BaseFloat lm_scale = 1.0;
    bool rm_eps = true;

    vector<vector<double> > scale = fst::LatticeScale(lm_scale, acoustic_scale);
    
    std::string lats_rspecifier = rclat,
        fsts_wspecifier = wclat;
    
    SequentialCompactLatticeReader lattice_reader(lats_rspecifier);
    TableWriter<fst::VectorFstHolder> fst_writer(fsts_wspecifier);
    
    int32 n_done = 0; // there is no failure mode, barring a crash.
    
    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      std::string key = lattice_reader.Key();
      CompactLattice clat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
      ScaleLattice(scale, &clat); // typically scales to zero.
      RemoveAlignmentsFromCompactLattice(&clat); // remove the alignments...
      fst::VectorFst<StdArc> fst;
      {
        Lattice lat;
        ConvertLattice(clat, &lat); // convert to non-compact form.. won't introduce
        // extra states because already removed alignments.
        ConvertLattice(lat, &fst); // this adds up the (lm,acoustic) costs to get
        // the normal (tropical) costs.
        Project(&fst, fst::PROJECT_OUTPUT); // Because in the standard Lattice format,
        // the words are on the output, and we want the word labels.
      }
      if (rm_eps) RemoveEpsLocal(&fst);
      
      fst_writer.Write(key, fst);
      n_done++;
    }
  } catch(const std::exception &e) {
    std::cerr << e.what();
  }

}




bool OpenAudioDevice(const char *devicename, snd_pcm_t **handle, snd_pcm_stream_t streamType, int channels, int rate, int flags)
{
    int rc;
    static snd_output_t *log;

    printf("Open device : %s %s channels %d rate %d\n", devicename,
        streamType == SND_PCM_STREAM_CAPTURE ? "CAPTURE" : "PLAYBACK",
        channels, rate);

    if ((rc = snd_pcm_open(handle, devicename, streamType, flags)) < 0) {
        printf("unable to open pcm device for recording: %s\n",snd_strerror(rc));
        return false;
    }

    if ((rc = snd_output_stdio_attach(&log, stderr, 0)) < 0) {
        printf("unable to attach log output: %s\n",snd_strerror(rc));
        return false;
    }

    if ((rc = snd_pcm_set_params(*handle, SND_PCM_FORMAT_S16_LE,  SND_PCM_ACCESS_RW_INTERLEAVED, channels, rate, 1, 1000000)) < 0) {
        printf("snd_pcm_set_params error: %s\n", snd_strerror(rc));
        return false;
    }

    if (streamType == SND_PCM_STREAM_CAPTURE) {
        snd_pcm_sw_params_t *sw_params = NULL;
        if((rc = snd_pcm_sw_params_malloc(&sw_params)) < 0) {
            printf("snd_pcm_sw_params_malloc error: %s\n", snd_strerror(rc));
            return false;
        }

        if((rc = snd_pcm_sw_params_current(*handle, sw_params)) < 0) {
            printf("snd_pcm_sw_params_current error: %s\n", snd_strerror(rc));
            return false;
        }

        if ((rc = snd_pcm_sw_params_set_start_threshold(*handle, sw_params, 1)) < 0) {
            printf("snd_pcm_sw_params_set_start_threshold failed: %s\n",snd_strerror(rc));
            return false;
        }

        if ((rc = snd_pcm_sw_params(*handle, sw_params)) < 0) {
            printf("snd_pcm_sw_params failed: %s\n",snd_strerror(rc));
            return false;
        }

        snd_pcm_sw_params_free(sw_params);
    }

    snd_pcm_dump(*handle, log);
    snd_output_close(log);
    return true;
}

bool StopAudioDevice(snd_pcm_t **handle)
{
    if (*handle != NULL) {
        snd_pcm_drain(*handle);
        snd_pcm_close(*handle);
    }
    return true;
}

int main(int argc, char *argv[]) 
{
  try 
  {

    using namespace kaldi;
    using namespace fst;

    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;

    const char *usage =
        "Reads in wav file(s) and simulates online decoding, including\n"
        "basis-fMLLR adaptation and endpointing.  Writes lattices.\n"
        "Models are specified via options.\n"
        "\n"
        "Usage: online2-wav-gmm-latgen-faster [options] <fst-in> "
        "<spk2utt-rspecifier> <wav-rspecifier> <lattice-wspecifier>\n"
        "Run egs/rm/s5/local/run_online_decoding.sh for example\n";

    ParseOptions po(usage);

    std::string word_syms_rxfilename;
    std::string capture_device;
    snd_pcm_t *inhandle;

    OnlineEndpointConfig endpoint_config;
    OnlineFeaturePipelineCommandLineConfig feature_cmdline_config;
    OnlineGmmDecodingConfig decode_config;

    BaseFloat chunk_length_secs = 0.02;
    bool do_endpointing = false;
    std::string use_gpu = "no";

    po.Register("chunk-length", &chunk_length_secs,
                "Length of chunk size in seconds, that we process.");
    po.Register("word-symbol-table", &word_syms_rxfilename,
                "Symbol table for words [for debug output]");
    po.Register("do-endpointing", &do_endpointing,
                "If true, apply endpoint detection");
    po.Register("capture-device", &capture_device,
                "ALSA capture device");

    feature_cmdline_config.Register(&po);
    decode_config.Register(&po);
    endpoint_config.Register(&po);

    po.Read(argc, argv);

    if (po.NumArgs() != 4) 
    {
      po.PrintUsage();
      return 1;
    }

    std::string fst_rxfilename = po.GetArg(1),
        spk2utt_rspecifier = po.GetArg(2),
        wav_rspecifier = po.GetArg(3),
        clat_wspecifier = po.GetArg(4);
    //CompactLattice clat;
    OnlineFeaturePipelineConfig feature_config(feature_cmdline_config);
    OnlineFeaturePipeline pipeline_prototype(feature_config);
    // The following object initializes the models we use in decoding.
    OnlineGmmDecodingModels gmm_models(decode_config);


    fst::Fst<fst::StdArc> *decode_fst = ReadFstKaldiGeneric(fst_rxfilename);

    fst::SymbolTable *word_syms = NULL;
    if (word_syms_rxfilename != "")
      if (!(word_syms = fst::SymbolTable::ReadText(word_syms_rxfilename)))
        KALDI_ERR << "Could not read symbol table from file "
                  << word_syms_rxfilename;

    int32 num_done = 0, num_err = 0;
    double tot_like = 0.0;
    int64 num_frames = 0;

    SequentialTokenVectorReader spk2utt_reader(spk2utt_rspecifier);
    RandomAccessTableReader<WaveHolder> wav_reader(wav_rspecifier);

    if (capture_device != "") 
    {
        if (OpenAudioDevice(capture_device.c_str(), &inhandle, SND_PCM_STREAM_CAPTURE, 1, 16000, 0) == false) 
        {
            KALDI_ERR << "Could not open capture PCM audio "
                  << capture_device;
            return 1;
        }
    }

    //OnlineTimingStats timing_stats;

    int count=0;
    
    OnlineGmmAdaptationState adaptation_state;
    VadInst *vad=WebRtcVad_Create();
    WebRtcVad_Init(vad);
    WebRtcVad_set_mode(vad, 1);
    
    //const int16_t * temp = sample.data();
    while(capture_device != "") 
    {
      string utt="uttr1";
      cout<<"==============================wait for next utterance "<<++count<<"==============================\n";
      CompactLatticeWriter clat_writer("ark:lat_in");
      CompactLatticeWriter clat_writer_partial("ark:lat_partial");

      SingleUtteranceGmmDecoder decoder(decode_config,
                                        gmm_models,
                                        pipeline_prototype,
                                        *decode_fst,
                                        adaptation_state);

      //OnlineTimer decoding_timer(utt);

      //BaseFloat samp_freq = wave_data.SampFreq();
      float32 samp_freq = 16000;
      int32 chunk_length = int32(samp_freq * chunk_length_secs);
      if (chunk_length == 0) chunk_length = 1;
      int cn=0;
      int flag1 = 0,flag0 = 0;
      int32 j=0;
      int32 buffer_length = chunk_length*10;
      Vector<BaseFloat> buffer(buffer_length, kSetZero);
      int isValid = 1;
      while(1)
      {
        //cout<<"====================speak start "<<++cn<<"====================\n";
        //cout<<"===========================>> "<<j<<endl;
        int16_t samples[chunk_length];
        Vector<BaseFloat> wave_part(chunk_length, kSetZero);
        int rc;
        if ((rc = snd_pcm_readi(inhandle, samples, chunk_length)) < 0) {
          printf("read failed (%s) ======>>>\n", snd_strerror (rc));
          if ((rc = snd_pcm_recover(inhandle, rc, 0)) < 0) {
            printf("recover failed (%s)\n", snd_strerror (rc));
            break;
          }
        }
        //check voice activity
        int isActive = WebRtcVad_Process(vad, 16000, samples, 320);
        //cout<<isActive<<endl;
        //check contiguous 10 frames for nonsilence, if silence reset everything
        if(isActive && flag1<10)
        {
          //cout<<"here1"<<endl;
          flag0 = 0;
          flag1++;
          //store everything untill silence
          for (int32 i= 0; i < chunk_length; ++i,++j) 
          {
            buffer(j) = samples[i];
          }
          //cout<<j<<endl;
          continue;
        }
        else if(!isActive && flag1<10) //
        {
          //cout<<"here2"<<endl;
          flag0++;
          //cout<<"reset1"<<endl;
          flag1 = 0;
          j = 0;
          buffer.SetZero();
          continue;
        }
        else if(!isActive)
        {
          flag0++;
          //cout<<flag0<<endl;
        }
        //else if(isActive)
          //flag0 = 0;
        if(flag1==10)//not silence send everything from buffer for decoding
        { 
          for(int32 i=0;i<buffer_length;i+=chunk_length)
          {
            for(int32 k=0;k<chunk_length;k++)
            {
              wave_part(k)=buffer(i+k);
            }
            decoder.FeaturePipeline().AcceptWaveform(samp_freq, wave_part);
            decoder.AdvanceDecoding();
            if (do_endpointing && decoder.EndpointDetected(endpoint_config)) 
            {
              cout<<"===================detected end of speech small====================\n";
              break;
            }
          }
          buffer.SetZero();
          j = 0;
          flag1++;
          flag0 = 0;
          continue;
        }
        //afterthat normal decoding
        for (uint32 i = 0; i < chunk_length; ++i) 
        {
          wave_part(i) = samples[i];
        }
        decoder.FeaturePipeline().AcceptWaveform(samp_freq, wave_part);
        decoder.AdvanceDecoding();

        //cout<<"f0: "<<flag0<<"  flag1: "<<flag1<<endl;
        
        if (do_endpointing && decoder.EndpointDetected(endpoint_config)) 
        {
          //decoder.FeaturePipeline().InputFinished();
          cout<<"f0: "<<flag0<<"flag1: "<<flag1<<endl;
          //if((flag1*4)<flag0)
            //isValid = 0;
          cout<<"===================detected end of speech====================\n";
          break;
        }
        else if(flag0>100 && (2*flag1)<flag0)
        {
          //isValid = 0;
          cout<<"f00: "<<flag0<<"flag11: "<<flag1<<endl;
          cout<<"exiting"<<endl;
          break;
        }
        /*else
        {
          //decoder.PruneActiveTokens(1.5);
          CompactLattice clat_partial;
          decoder.GetLattice(true, true, &clat_partial);
          clat_writer_partial.Write(utt, clat_partial);
          //decoder.PruneActiveTokens(1.5);
          //BestPathEnd(bool true,
                 //              BaseFloat *final_cost = NULL)
          LatScale("ark:lat_partial","ark:lat_scale_partial");
          LatPenalty("ark:lat_scale_partial","ark:lat_pen_partial");
          //LatBest_1_Path("ark:lat_pen_partial","ark:one-best.lat");
          string transcript_partial=LatBestPath("ark:lat_pen_partial","ark,t:test.tra", word_syms_rxfilename,"ark:partial.ali");
          cout<<"partial: "<<transcript_partial<<endl;
          //Lattice partial_best;
          //decoder.GetBestPath(true,&partial_best);
          //lattice_to_fst("ark:one_best","ark,t:one_best.tra");
          //system("cat one_best.tra");

        }*/
      }
      decoder.FeaturePipeline().InputFinished();
      decoder.FinalizeDecoding();

      bool end_of_utterance = true;
      decoder.EstimateFmllr(end_of_utterance);
      bool rescore_if_needed = true;
      CompactLattice clat;
      decoder.GetLattice(rescore_if_needed, end_of_utterance, &clat);

      // we want to output the lattice with un-scaled acoustics.
      if (decode_config.acoustic_scale != 0.0) {
        BaseFloat inv_acoustic_scale = 1.0 / decode_config.acoustic_scale;
        ScaleLattice(AcousticLatticeScale(inv_acoustic_scale), &clat);
      }
      clat_writer.Write(utt, clat);

      LatScale("ark:lat_in","ark:lat_scale");

      LatPenalty("ark:lat_scale","ark:lat_pen");
      if(!isValid)
      {
        cout<<"last"<<endl;
        continue;
      }
      //LatBest_1_Path("ark:lat_pen","ark,t:one-best.tra");
      string transcript=LatBestPath("ark:lat_pen",clat_wspecifier, word_syms_rxfilename,"ark:best.ali");
      cout<<"Transcript: "<< transcript<< endl;
      cout<<"key words:  ";
      if(!searchKW(transcript))
        cout<<"No keyword found"<<endl;
      //system("pwd");
      //system("../../../src/bin/ali-to-phones exp/custom_model_460_online/final.mdl ark:best.ali ark,t:best");
      //system("less best | utils/int2sym.pl -f 2- data/lang_3g_command/phones.txt");
      //remove("lat_in");
      //remove("lat_scale");
      //remove("lat_pen");
      //remove("one-best.tra");
    }
    if (capture_device != "")
        StopAudioDevice(&inhandle);
    //delete decode_fst;
    //delete word_syms; // will delete if non-NULL.
    return (num_done != 0 ? 0 : 1);
  } 
  catch(const std::exception& e) {
    std::cerr << e.what();
    std::cout<<"bye\n";
    return -1;
  }
} // main()
