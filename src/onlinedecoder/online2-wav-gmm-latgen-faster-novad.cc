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

namespace kaldi {

//search KW
int searchKW(string text)
{
  string pattern ="(\\b\\s*MENU\\s*\\b|\\b\\s*SHOW MAP\\s*\\b|\\b\\s*BECOME DRIVER\\s*\\b|\\b\\s*BOOST\\s*\\b|\\b\\s*BUILD\\s*\\b|\\b\\s*BUILD MEDIVAC\\s*\\b|\\b\\s*BUILD MINE\\s*\\b|\\b\\s*BUILD THREE DROPSHIPS\\s*\\b|\\b\\s*CHANGE CAMERA\\s*\\b|\\b\\s*CHANGE MODE\\s*\\b|\\b\\s*CROUCH\\s*\\b|\\b\\s*DRIVE MOTORBIKE\\s*\\b|\\b\\s*DRIVE VEHICLE\\s*\\b|\\b\\s*ENTER CAR\\s*\\b|\\b\\s*EQUIP MELEE\\s*\\b|\\b\\s*EXPLORE\\s*\\b|\\b\\s*EXTRA SUPPLIES\\s*\\b|\\b\\s*FIRST AID KIT\\s*\\b|\\b\\s*FOUND CITY\\s*\\b|\\b\\s*GET ARMY\\s*\\b|\\b\\s*GET BATTLE UNITS\\s*\\b|\\b\\s*GET IDLE WORKER\\s*\\b|\\b\\s*GRENADE\\s*\\b|\\b\\s*JUMP\\s*\\b|\\b\\s*MAKE BARRACKS\\s*\\b|\\b\\s*MOVE TO PASSENGER\\s*\\b|\\b\\s*NEXT TURN\\s*\\b|\\b\\s*NITRO\\s*\\b|\\b\\s*OPEN CELLPHONE\\s*\\b|\\b\\s*OPEN CHAT\\s*\\b|\\b\\s*OPEN MENU\\s*\\b|\\b\\s*OPTIONS\\s*\\b|\\b\\s*PRIMARY WEAPON\\s*\\b|\\b\\s*PRONE\\s*\\b|\\b\\s*QUICK SAVE\\s*\\b|\\b\\s*RECORD\\s*\\b|\\b\\s*RELOAD\\s*\\b|\\b\\s*REMOVE FOREST\\s*\\b|\\b\\s*SCREENSHOT\\s*\\b|\\b\\s*SET UP ARTILLERY\\s*\\b|\\b\\s*SHOW DIPLOMACY\\s*\\b|\\b\\s*SHOW ECONOMY\\s*\\b|\\b\\s*SHOW HUD\\s*\\b|\\b\\s*SHOW INVENTORY\\s*\\b|\\b\\s*SHOW MILITARY\\s*\\b|\\b\\s*STEALTH MODE\\s*\\b|\\b\\s*SWITCH TO HANDGUN\\s*\\b|\\b\\s*SWITCH TO WEAPON ONE\\s*\\b|\\b\\s*TAKE JEEP\\s*\\b|\\b\\s*TAKE PICTURE\\s*\\b|\\b\\s*TECH TREE\\s*\\b|\\b\\s*TIMEOUT\\s*\\b|\\b\\s*WEAPON TWO\\s*\\b)";
  std::smatch m;
  std::regex r(pattern);
  int count=0;
  while (std::regex_search (text,m,r)) {
    std::cout << m[0] << " ";
    std::cout << std::endl;
    text = m.suffix().str();
    count++;
  }
  std::cout<<count<<std::endl;
  if(count>0)
    return 1;
  else
    return 0;
}

//fst-project
/*void fstproject(const std::string rclat, const std::string wclat) {
  try {
    using namespace kaldi;
    using namespace fst;
    typedef kaldi::int32 int32;
    typedef kaldi::uint64 uint64;
    bool project_output = false;

    std::string fsts_rspecifier = rclat,
        fsts_wspecifier = wclat;


    SequentialTableReader<VectorFstHolder> fst_reader(fsts_rspecifier);
    TableWriter<VectorFstHolder> fst_writer(fsts_wspecifier);

    int32 n_done = 0;
    for (; !fst_reader.Done(); fst_reader.Next()) {
      std::string key = fst_reader.Key();
      VectorFst<StdArc> fst(fst_reader.Value());

      Project(&fst, project_output ? PROJECT_OUTPUT : PROJECT_INPUT);

      fst_writer.Write(key, fst);
      n_done++;
    }

    //KALDI_LOG << "Projected " << n_done << " FSTs";
  } catch(const std::exception &e) {
    std::cerr << e.what();
  }
}*/



//Kaldi rescore
/*void rescore(const std::string rclat,const std::string lm, const std::string wclat) {
  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    using fst::ReadFstKaldi;
    BaseFloat lm_scale = 1.0;
    int32 num_states_cache = 50000;

    std::string lats_rspecifier = rclat,
        fst_rxfilename = lm,
        lats_wspecifier = wclat;

    VectorFst<StdArc> *std_lm_fst = ReadFstKaldi(fst_rxfilename);
    if (std_lm_fst->Properties(fst::kILabelSorted, true) == 0) {
      // Make sure LM is sorted on ilabel.
      fst::ILabelCompare<StdArc> ilabel_comp;
      fst::ArcSort(std_lm_fst, ilabel_comp);
    }

    // mapped_fst is the LM fst interpreted using the LatticeWeight semiring,
    // with all the cost on the first member of the pair (since it's a graph
    // weight).
    fst::CacheOptions cache_opts(true, num_states_cache);
    fst::MapFstOptions mapfst_opts(cache_opts);
    fst::StdToLatticeMapper<BaseFloat> mapper;
    fst::MapFst<StdArc, LatticeArc, fst::StdToLatticeMapper<BaseFloat> >
        lm_fst(*std_lm_fst, mapper, mapfst_opts);
    delete std_lm_fst;

    // The next fifteen or so lines are a kind of optimization and
    // can be ignored if you just want to understand what is going on.
    // Change the options for TableCompose to match the input
    // (because it's the arcs of the LM FST we want to do lookup
    // on).
    fst::TableComposeOptions compose_opts(fst::TableMatcherOptions(),
                                          true, fst::SEQUENCE_FILTER,
                                          fst::MATCH_INPUT);

    // The following is an optimization for the TableCompose
    // composition: it stores certain tables that enable fast
    // lookup of arcs during composition.
    fst::TableComposeCache<fst::Fst<LatticeArc> > lm_compose_cache(compose_opts);

    // Read as regular lattice-- this is the form we need it in for efficient
    // composition and determinization.
    SequentialLatticeReader lattice_reader(lats_rspecifier);

    // Write as compact lattice.
    CompactLatticeWriter compact_lattice_writer(lats_wspecifier);

    int32 n_done = 0, n_fail = 0;

    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      std::string key = lattice_reader.Key();
      Lattice lat = lattice_reader.Value();
      lattice_reader.FreeCurrent();
      if (lm_scale != 0.0) {
        // Only need to modify it if LM scale nonzero.
        // Before composing with the LM FST, we scale the lattice weights
        // by the inverse of "lm_scale".  We'll later scale by "lm_scale".
        // We do it this way so we can determinize and it will give the
        // right effect (taking the "best path" through the LM) regardless
        // of the sign of lm_scale.
        fst::ScaleLattice(fst::GraphLatticeScale(1.0 / lm_scale), &lat);
        ArcSort(&lat, fst::OLabelCompare<LatticeArc>());

        Lattice composed_lat;
        // Could just do, more simply: Compose(lat, lm_fst, &composed_lat);
        // and not have lm_compose_cache at all.
        // The command below is faster, though; it's constant not
        // logarithmic in vocab size.
        TableCompose(lat, lm_fst, &composed_lat, &lm_compose_cache);

        Invert(&composed_lat); // make it so word labels are on the input.
        CompactLattice determinized_lat;
        DeterminizeLattice(composed_lat, &determinized_lat);
        fst::ScaleLattice(fst::GraphLatticeScale(lm_scale), &determinized_lat);
        if (determinized_lat.Start() == fst::kNoStateId) {
          KALDI_WARN << "Empty lattice for utterance " << key << " (incompatible LM?)";
          n_fail++;
        } else {
          compact_lattice_writer.Write(key, determinized_lat);
          n_done++;
        }
      } else {
        // zero scale so nothing to do.
        n_done++;
        CompactLattice compact_lat;
        ConvertLattice(lat, &compact_lat);
        compact_lattice_writer.Write(key, compact_lat);
      }
    }

    //KALDI_LOG << "Done " << n_done << " lattices, failed for " << n_fail;
   
  } catch(const std::exception &e) {
    std::cerr << e.what();
  }
}*/


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
    BaseFloat inv_acoustic_scale = 17.0;
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


    BaseFloat word_ins_penalty = 0.5882;

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
string LatBestPath(const std::string rclat, const std::string wclat, std::string word_syms_filename)
{

  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;
    std::string transcript="";

    BaseFloat acoustic_scale = 0.5882;
    BaseFloat lm_scale = 1.0;

    std::string lats_rspecifier = rclat,
        transcriptions_wspecifier = wclat,
        alignments_wspecifier = "";

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

int main(int argc, char *argv[]) {
  try {

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

    BaseFloat chunk_length_secs = 0.05;
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

    if (po.NumArgs() != 4) {
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

    if (capture_device != "") {
        if (OpenAudioDevice(capture_device.c_str(), &inhandle, SND_PCM_STREAM_CAPTURE, 1, 16000, 0) == false) {
            KALDI_ERR << "Could not open capture PCM audio "
                  << capture_device;
            return 1;
        }
    }

    //OnlineTimingStats timing_stats;

    int count=0;
    
    OnlineGmmAdaptationState adaptation_state;
    while(capture_device != "") {
      string utt="uttr1";
      cout<<"==============================wait for next utterance "<<++count<<"==============================\n";
      CompactLatticeWriter clat_writer("ark:lat_in");

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
      while(1){
        //cout<<"====================speak start "<<++cn<<"====================\n";
        //int32 num_samp;
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
        for (uint32 i = 0; i < chunk_length; ++i) {
          wave_part(i) = samples[i];
        }
        decoder.FeaturePipeline().AcceptWaveform(samp_freq, wave_part);
        decoder.AdvanceDecoding();

        if (do_endpointing && decoder.EndpointDetected(endpoint_config)) {
          cout<<"===================detected end of speech====================\n";
          break;
        }
      }
      decoder.FeaturePipeline().InputFinished();
      decoder.FinalizeDecoding();

      bool end_of_utterance = true;
      //decoder.EstimateFmllr(end_of_utterance);
      bool rescore_if_needed = true;
      CompactLattice clat;
      decoder.GetLattice(rescore_if_needed, end_of_utterance, &clat);

      // we want to output the lattice with un-scaled acoustics.
      if (decode_config.acoustic_scale != 0.0) {
        BaseFloat inv_acoustic_scale = 1.0 / decode_config.acoustic_scale;
        ScaleLattice(AcousticLatticeScale(inv_acoustic_scale), &clat);
      }
      clat_writer.Write(utt, clat);


      //string oldlmcommand="fstproject --project_output=true G1.fst |";
      //string newlmcommand="fstproject --project_output=true G2.fst |";

      //cout<<"============="<<oldlmcommand<<endl;

      //rescore("ark:lat_in","fstproject --project_output=true G1.fst |","ark:lat_res");
      //rescore("ark:lat_res","fstproject --project_output=true G2.fst |","ark:lat_res_new");

      LatScale("ark:lat_in","ark:lat_scale");

      LatPenalty("ark:lat_scale","ark:lat_pen");
      string transcript=LatBestPath("ark:lat_pen",clat_wspecifier, word_syms_rxfilename);
      cout<<"Transcript: "<< transcript<< endl;
      cout<<"key words:"<<endl;
      if(!searchKW(transcript))
        cout<<"No keyword found"<<endl;
      //remove("lat_in");
      //remove("lat_scale");
      //remove("lat_pen");
    }
    if (capture_device != "")
        StopAudioDevice(&inhandle);
    delete decode_fst;
    delete word_syms; // will delete if non-NULL.
    return (num_done != 0 ? 0 : 1);
  } 
  catch(const std::exception& e) {
    std::cerr << e.what();
    return -1;
  }
} // main()