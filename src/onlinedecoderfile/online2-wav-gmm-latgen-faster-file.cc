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


#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "lat/kaldi-lattice.h"

#include  <alsa/asoundlib.h>

namespace kaldi {

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


    BaseFloat word_ins_penalty = 1.0;

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
void LatBestPath(const std::string rclat, const std::string wclat, std::string word_syms_filename)
{

  try {
    using namespace kaldi;
    typedef kaldi::int32 int32;
    typedef kaldi::int64 int64;
    using fst::SymbolTable;
    using fst::VectorFst;
    using fst::StdArc;

    BaseFloat acoustic_scale = 1.0;
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
        /*KALDI_LOG << "For utterance " << key << ", best cost "
                  << weight.Value1() << " + " << weight.Value2() << " = "
                  << (weight.Value1() + weight.Value2())
                  << " over " << alignment.size() << " frames.";*/
        if (transcriptions_wspecifier != "")
          transcriptions_writer.Write(key, words);
        if (alignments_wspecifier != "")
          alignments_writer.Write(key, alignment);
        if (word_syms != NULL) {
          std::cerr << key << ' ';
          for (size_t i = 0; i < words.size(); i++) {
            std::string s = word_syms->Find(words[i]);
            if (s == "")
              KALDI_ERR << "Word-id " << words[i] <<" not in symbol table.";
            std::cerr << s << ' ';
          }
          std::cerr << '\n';
        }
        n_done++;
        n_frame += alignment.size();
        tot_weight = Times(tot_weight, weight);
      }
    }

    BaseFloat tot_weight_float = tot_weight.Value1() + tot_weight.Value2();
    /*KALDI_LOG << "Overall cost per frame is " << (tot_weight_float/n_frame)
              << " = " << (tot_weight.Value1()/n_frame) << " [graph]"
              << " + " << (tot_weight.Value2()/n_frame) << " [acoustic]"
              << " over " << n_frame << " frames.";
    KALDI_LOG << "Done " << n_done << " lattices, failed for " << n_fail;*/

    delete word_syms;
  } catch(const std::exception &e) {
    std::cerr << e.what();
  }



}

void GetDiagnosticsAndPrintOutput(const std::string &utt,
                                  const fst::SymbolTable *word_syms,
                                  const CompactLattice &clat,
                                  int64 *tot_num_frames,
                                  double *tot_like) {
  if (clat.NumStates() == 0) {
    KALDI_WARN << "Empty lattice.";
    return;
  }
  CompactLattice best_path_clat;
  CompactLatticeShortestPath(clat, &best_path_clat);

  Lattice best_path_lat;
  ConvertLattice(best_path_clat, &best_path_lat);

  double likelihood;
  LatticeWeight weight;
  int32 num_frames;
  std::vector<int32> alignment;
  std::vector<int32> words;
  GetLinearSymbolSequence(best_path_lat, &alignment, &words, &weight);
  num_frames = alignment.size();
  likelihood = -(weight.Value1() + weight.Value2());
  *tot_num_frames += num_frames;
  *tot_like += likelihood;
  /*KALDI_VLOG(2) << "Likelihood per frame for utterance " << utt << " is "
                << (likelihood / num_frames) << " over " << num_frames
                << " frames.";*/

  if (word_syms != NULL) {
    std::string result;
    std::cerr << utt << ' ';
    for (size_t i = 0; i < words.size(); i++) {
      std::string s = word_syms->Find(words[i]);
      if (s == "")
        KALDI_ERR << "Word-id " << words[i] << " not in symbol table.";
      //std::cerr << s << ' ';
      result.append(s);
      result.append(" ");
    }
    KALDI_LOG << "Result " << result;
    std::cerr << std::endl;
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

    if ((rc = snd_pcm_set_params(*handle, SND_PCM_FORMAT_S16_LE,  SND_PCM_ACCESS_RW_INTERLEAVED, channels, rate, 1, 500000)) < 0) {
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

    CompactLattice clat;
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

    OnlineTimingStats timing_stats;


    for (; !spk2utt_reader.Done(); spk2utt_reader.Next()) {
      std::string spk = spk2utt_reader.Key();
      const std::vector<std::string> &uttlist = spk2utt_reader.Value();
      OnlineGmmAdaptationState adaptation_state;
      for (size_t i = 0; i < uttlist.size() || (capture_device != ""); i++) {
        std::string utt = uttlist[0];
        CompactLatticeWriter clat_writer("ark:lat_in");
        if (!wav_reader.HasKey(utt)) {
          KALDI_WARN << "Did not find audio for utterance " << utt;
          num_err++;
          continue;
        }
        const WaveData &wave_data = wav_reader.Value(utt);
        // get the data for channel zero (if the signal is not mono, we only
        // take the first channel).
        SubVector<BaseFloat> data(wave_data.Data(), 0);

        SingleUtteranceGmmDecoder decoder(decode_config,
                                          gmm_models,
                                          pipeline_prototype,
                                          *decode_fst,
                                          adaptation_state);

        OnlineTimer decoding_timer(utt);

        BaseFloat samp_freq = wave_data.SampFreq();
        int32 chunk_length = int32(samp_freq * chunk_length_secs);
        if (chunk_length == 0) chunk_length = 1;

        int32 samp_offset = 0;
        while ((capture_device != "") || samp_offset < data.Dim()) {
          int32 num_samp;
          if (capture_device == "") {
              int32 samp_remaining = data.Dim() - samp_offset;
              num_samp = chunk_length < samp_remaining ? chunk_length
                                                             : samp_remaining;

              SubVector<BaseFloat> wave_part(data, samp_offset, num_samp);
              decoder.FeaturePipeline().AcceptWaveform(samp_freq, wave_part);
          }
          else {
            int16_t samples[chunk_length];
            Vector<BaseFloat> wave_part(chunk_length, kSetZero);
            int rc;
            if ((rc = snd_pcm_readi(inhandle, samples, chunk_length)) < 0) {
                printf("read failed (%s)\n", snd_strerror (rc));
                if ((rc = snd_pcm_recover(inhandle, rc, 0)) < 0) {
                    printf("recover failed (%s)\n", snd_strerror (rc));
                    break;
                }
            }
            for (uint32 i = 0; i < chunk_length; ++i) {
                wave_part(i) = samples[i];
            }
            decoder.FeaturePipeline().AcceptWaveform(samp_freq, wave_part);
            num_samp = chunk_length;
          }
          samp_offset += num_samp;
          decoding_timer.WaitUntil(samp_offset / samp_freq);
          if (samp_offset == data.Dim() && (capture_device == "")) {
            // no more input. flush out last frames
            decoder.FeaturePipeline().InputFinished();
          }

          decoder.AdvanceDecoding();

          if (do_endpointing && decoder.EndpointDetected(endpoint_config)) {
            KALDI_LOG << "Endpoint detected. Processed samples " << samp_offset;
            break;
          }
        }
        decoder.FinalizeDecoding();

        bool end_of_utterance = true;
        decoder.EstimateFmllr(end_of_utterance);
        bool rescore_if_needed = true;
        decoder.GetLattice(rescore_if_needed, end_of_utterance, &clat);

        //GetDiagnosticsAndPrintOutput(utt, word_syms, clat,
          //                           &num_frames, &tot_like);

        //decoding_timer.OutputStats(&timing_stats);

        // In an application you might avoid updating the adaptation state if
        // you felt the utterance had low confidence.  See lat/confidence.h
        decoder.GetAdaptationState(&adaptation_state);

        // we want to output the lattice with un-scaled acoustics.
        if (decode_config.acoustic_scale != 0.0) {
          BaseFloat inv_acoustic_scale = 1.0 / decode_config.acoustic_scale;
          ScaleLattice(AcousticLatticeScale(inv_acoustic_scale), &clat);
        }
        clat_writer.Write(utt, clat);
        //KALDI_LOG << "Decoded utterance " << utt;

        LatScale("ark:lat_in","ark:lat_scale");
        LatPenalty("ark:lat_scale","ark:lat_pen");
        LatBestPath("ark:lat_pen",clat_wspecifier, word_syms_rxfilename);
        //std::remove("lat_in");
        //std::remove("lat_scale");
        //std::remove("lat_pen");

        num_done++;
      }
    }

    if (capture_device != "")
        StopAudioDevice(&inhandle);

    timing_stats.Print();
    /*KALDI_LOG << "Decoded " << num_done << " utterances, "
              << num_err << " with errors.";
    KALDI_LOG << "Overall likelihood per frame was " << (tot_like / num_frames)
              << " per frame over " << num_frames << " frames.";*/
    delete decode_fst;
    delete word_syms; // will delete if non-NULL.
    return (num_done != 0 ? 0 : 1);
  } catch(const std::exception& e) {
    std::cerr << e.what();
    return -1;
  }
} // main()
