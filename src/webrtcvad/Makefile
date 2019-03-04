EXTRA_CXXFLAGS = -Wno-sign-compare

all:

include ../kaldi.mk

OBJFILES = ./webrtc/common_audio/third_party/spl_sqrt_floor/spl_sqrt_floor.o \
	./webrtc/common_audio/signal_processing/spl_init.o \
	./webrtc/common_audio/signal_processing/downsample_fast.o \
	./webrtc/common_audio/signal_processing/cross_correlation.o \
	./webrtc/common_audio/signal_processing/spl_inl.o \
	./webrtc/common_audio/signal_processing/resample_48khz.o \
	./webrtc/common_audio/signal_processing/complex_fft.o \
	./webrtc/common_audio/signal_processing/get_scaling_square.o \
	./webrtc/common_audio/signal_processing/complex_bit_reverse.o \
	./webrtc/common_audio/signal_processing/min_max_operations.o \
	./webrtc/common_audio/signal_processing/vector_scaling_operations.o \
	./webrtc/common_audio/signal_processing/spl_sqrt.o \
	./webrtc/common_audio/signal_processing/division_operations.o \
	./webrtc/common_audio/signal_processing/resample_fractional.o \
	./webrtc/common_audio/signal_processing/energy.o \
	./webrtc/common_audio/signal_processing/resample_by_2_internal.o \
	./webrtc/common_audio/vad/vad_sp.o \
	./webrtc/common_audio/vad/webrtc_vad.o \
	./webrtc/common_audio/vad/vad_filterbank.o \
	./webrtc/common_audio/vad/vad_gmm.o \
	./webrtc/common_audio/vad/vad_core.o


LDFLAGS = -fPIC

CFLAGS = -I. -fPIC -DNDEBUG -DWEBRTC_POSIX

LIBNAME = webrtcvad

include ../makefiles/default_rules.mk
