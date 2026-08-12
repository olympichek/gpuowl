# Use "make CUDA=1" for a CUDA build, use "make DEBUG=1" for a debug build

# The build artifacts are put in the "build-release" subfolder (or "build-debug" for a debug build).

# On Windows invoke with "make exe" or "make all"

DEBUG = 0
CUDA = 0
STATIC_RUNTIME = 0
STATIC_CUDA = 0

# Uncomment below as desired to set a particular compiler or force a debug build:
# CXX = g++-12
# DEBUG = 1
# or export those into environment, or pass on the command line e.g.
# make all DEBUG=1 CXX=g++-12

HOST_OS = $(shell uname -s)

CXX ?= g++

ifeq ($(CUDA), 1)
 BIN=build-cuda
 CUDASRCS1 = clwrap_cuda.cpp cudawrap.cpp
 CUDAFLAGS = -DCUDA_BACKEND -Isrc/cuda -I/usr/local/cuda/include
 CUDAOBJS = $(CUDASRCS1:%.cpp=$(BIN)/%.o)
 ifeq ($(STATIC_CUDA), 1)
  OPENCL_LIBS = -L/usr/local/cuda/lib64 -Wl,--start-group -lnvrtc_static -lnvrtc-builtins_static -lnvptxcompiler_static -Wl,--end-group -lcuda -lpthread -ldl
 else
  OPENCL_LIBS = -L/usr/local/cuda/lib64 -Wl,-rpath,'$$ORIGIN' -lnvrtc -lcuda -lpthread
 endif
else
 BIN=build-release
 CUDAFLAGS =
 CUDAOBJS =
 ifeq ($(HOST_OS), Darwin)
  OPENCL_LIBS = -framework OpenCL
 else
  OPENCL_LIBS = -lOpenCL -lpthread
 endif
endif

COMMON_FLAGS = -Wall -Wextra $(CUDAFLAGS) -std=c++20

ifeq ($(STATIC_RUNTIME),1)
 LDFLAGS += -static-libstdc++ -static-libgcc

 ifeq ($(findstring MINGW, $(HOST_OS)), MINGW)
# For mingw-64 use this:
  LDFLAGS += -static
 endif
endif

ifeq ($(findstring MINGW, $(HOST_OS)), MINGW)
 CPPFLAGS += -DWINVER=0x0601 -D_WIN32_WINNT=0x0601
 LDFLAGS += -Wl,--subsystem,console:6.01
endif
# -fext-numeric-literals

ifeq ($(DEBUG), 1)

BIN=build-debug
CXXFLAGS = -g -Og $(COMMON_FLAGS)

else

CXXFLAGS = -O3 -flto -DNDEBUG $(COMMON_FLAGS)

endif

SRCS1 = fs.cpp Trig.cpp TuneEntry.cpp Primes.cpp tune.cpp CycleFile.cpp TrigBufCache.cpp Event.cpp Queue.cpp TimeInfo.cpp Profile.cpp bundle.cpp Saver.cpp KernelCompiler.cpp Kernel.cpp gpuid.cpp File.cpp Proof.cpp log.cpp Worktodo.cpp common.cpp main.cpp Gpu.cpp clwrap.cpp Task.cpp timeutil.cpp Args.cpp state.cpp Signal.cpp FFTConfig.cpp AllocTrac.cpp sha3.cpp md5.cpp version.cpp

SRCS2 = test.cpp

RNS31_TEST_OBJS = $(BIN)/RNS31.o $(BIN)/rns31_test.o
CARRY_TRANSFER_TEST_OBJS = $(BIN)/CarryTransfer.o $(BIN)/carry_transfer_test.o
RIESEL_ALGEBRA_TEST_OBJS = $(BIN)/riesel_algebra_test.o

# SRCS=$(addprefix src/, $(SRCS1))

OBJS = $(CUDAOBJS) $(SRCS1:%.cpp=$(BIN)/%.o)
DEPDIR := $(BIN)/.d
$(shell mkdir -p $(DEPDIR) >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c
POSTCOMPILE = @mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

all: prpll

rns31-test: $(BIN)/rns31-test

carry-transfer-test: $(BIN)/carry-transfer-test

folded-syndrome-test: $(BIN)/folded-syndrome-test

gold-fold-decoder-test: $(BIN)/gold-fold-decoder-test

riesel-algebra-test: $(BIN)/riesel-algebra-test

rns31-cuda-bench: build-cuda/rns31-cuda-bench

rns31-ntt-bench: build-cuda/rns31-ntt-bench

rns-multi-ntt-bench: build-cuda/rns-multi-ntt-bench

rns-tensor-bench: build-cuda/rns-tensor-bench

rns-mont-ntt-bench: build-cuda/rns-mont-ntt-bench

rns31-crt-bench: build-cuda/rns31-crt-bench

q24-m61-charged-bench: build-cuda/q24-m61-charged-bench

build-cuda/q24-m61-charged-bench: src/cuda/q24_m61_charged_bench.cu
	/usr/local/cuda/bin/nvcc -ccbin g++ -O3 -std=c++20 -arch=sm_120 -Isrc \
	  -o build-cuda/q24-m61-charged-bench src/cuda/q24_m61_charged_bench.cu

riesel-lazy-tile-bench: build-cuda/riesel-lazy-tile-bench

warp-specialized-edge-bench: build-cuda/warp-specialized-edge-bench

m61-limb-bench: build-cuda/m61-limb-bench

near-m61-radix33-bench: build-cuda/near-m61-radix33-bench

m61-pfa33-edge-bench: build-cuda/m61-pfa33-edge-bench

m61-batch-bench: build-cuda/m61-batch-bench

m61-middle-batch-bench: build-cuda/m61-middle-batch-bench

q24-tensor-bench: build-cuda/q24-tensor-bench

rns31-shape-bench: build-cuda/rns31-shape-bench

m61-resident-tile-bench: build-cuda/m61-resident-tile-bench

m61-warp-tail-bench: build-cuda/m61-warp-tail-bench

m31-m61-composite-bench: build-cuda/m31-m61-composite-bench

m31-m19-composite-bench: build-cuda/m31-m19-composite-bench

m31-prime-power-bench: build-cuda/m31-prime-power-bench

m31-prime-power-tile-bench: build-cuda/m31-prime-power-tile-bench

m31-cubic-power-bench: build-cuda/m31-cubic-power-bench

m31-m61-direct-carry-bench: build-cuda/m31-m61-direct-carry-bench

m31-m61-q3m-tile-bench: build-cuda/m31-m61-q3m-tile-bench

m31-m61-q3m-viable-tile-bench: build-cuda/m31-m61-q3m-viable-tile-bench

m31-m61-q3m-small-tile-bench: build-cuda/m31-m61-q3m-small-tile-bench

m31-m61-m19-3m-tile-bench: build-cuda/m31-m61-m19-3m-tile-bench

m31-m61-q3m-carry-bench: build-cuda/m31-m61-q3m-carry-bench

m31-m61-q3m-viable-carry-bench: build-cuda/m31-m61-q3m-viable-carry-bench

m31-m61-q3m-small-carry-bench: build-cuda/m31-m61-q3m-small-carry-bench

m31-m61-q-radix3-bench: build-cuda/m31-m61-q-radix3-bench

m31-m61-q-viable-radix3-bench: build-cuda/m31-m61-q-viable-radix3-bench

m31-m61-q-small-radix3-bench: build-cuda/m31-m61-q-small-radix3-bench

m31-m61-m19-radix3-bench: build-cuda/m31-m61-m19-radix3-bench

m31-m61-m19-radix7-bench: build-cuda/m31-m61-m19-radix7-bench

m31-m61-radix31-edge-bench: build-cuda/m31-m61-radix31-edge-bench

m31-m61-radix63-edge-bench: build-cuda/m31-m61-radix63-edge-bench

m61-hartley-tile-bench: build-cuda/m61-hartley-tile-bench

m31-m61-m19-3m-carry-bench: build-cuda/m31-m61-m19-3m-carry-bench

m89-limb-bench: build-cuda/m89-limb-bench

m127-limb-bench: build-cuda/m127-limb-bench

m31-gold-crt-bench: build-cuda/m31-gold-crt-bench

fp-compensated-tensor-bench: build-cuda/fp-compensated-tensor-bench

q31-fixed-fft-bench: build-cuda/q31-fixed-fft-bench

q31-fixed-overlap-bench: build-cuda/q31-fixed-overlap-bench

gold-fold-ntt-bench: build-cuda/gold-fold-ntt-bench

gold-pair-shape-bench: build-cuda/gold-pair-shape-bench

prpll: $(BIN)/prpll

amd: $(BIN)/prpll-amd

#$(BIN)/test: $(BIN)/test.o
#	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBPATH)

$(BIN)/prpll: ${OBJS}
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ ${OBJS} $(LIBPATH) $(OPENCL_LIBS)

$(BIN)/rns31-test: $(RNS31_TEST_OBJS)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $(RNS31_TEST_OBJS)

$(BIN)/carry-transfer-test: $(CARRY_TRANSFER_TEST_OBJS)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $(CARRY_TRANSFER_TEST_OBJS)

$(BIN)/folded-syndrome-test: src/folded_syndrome_test.cpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -std=c++20 -o $@ $<

$(BIN)/gold-fold-decoder-test: src/gold_fold_decoder_test.cpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -std=c++20 -o $@ $<

$(BIN)/riesel-algebra-test: $(RIESEL_ALGEBRA_TEST_OBJS)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $(RIESEL_ALGEBRA_TEST_OBJS)

build-cuda/rns31-cuda-bench: src/cuda/rns31_cuda_bench.cu src/RNS31Config.h
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ src/cuda/rns31_cuda_bench.cu

build-cuda/rns31-ntt-bench: src/cuda/rns31_ntt_bench.cu src/cuda/rns31.cuh src/RNS31Config.h
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ src/cuda/rns31_ntt_bench.cu

build-cuda/gold-fold-ntt-bench: src/cuda/gold_fold_ntt_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ $<

build-cuda/gold-pair-shape-bench: src/cuda/gold_pair_shape_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ $<

build-cuda/rns-multi-ntt-bench: src/cuda/rns_multi_ntt_bench.cu src/cuda/rns31.cuh src/RNS31Config.h
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ src/cuda/rns_multi_ntt_bench.cu

build-cuda/rns-tensor-bench: src/cuda/rns_tensor_bench.cu src/cuda/rns31.cuh
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ src/cuda/rns_tensor_bench.cu

build-cuda/rns-mont-ntt-bench: src/cuda/rns_mont_ntt_bench.cu src/cuda/rns31.cuh src/RNS31Config.h
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ src/cuda/rns_mont_ntt_bench.cu

build-cuda/rns31-crt-bench: src/cuda/rns31_crt_bench.cu src/cuda/rns31.cuh src/RNS31Config.h
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ src/cuda/rns31_crt_bench.cu

build-cuda/riesel-lazy-tile-bench: src/cuda/riesel_lazy_tile_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Isrc -o $@ src/cuda/riesel_lazy_tile_bench.cu

build-cuda/warp-specialized-edge-bench: src/cuda/warp_specialized_edge_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=compute_120a -code=sm_120a -Isrc -o $@ src/cuda/warp_specialized_edge_bench.cu

build-cuda/m31-fp-limb-bench: src/cuda/m31_fp_limb_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -lineinfo -Isrc -o $@ src/cuda/m31_fp_limb_bench.cu

m31-fp-limb-bench: build-cuda/m31-fp-limb-bench
	$<

build-cuda/q24-m61-overlap-bench: src/cuda/q24_m61_overlap_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -lineinfo -Isrc -o $@ src/cuda/q24_m61_overlap_bench.cu

q24-m61-overlap-bench: build-cuda/q24-m61-overlap-bench
	$<

build-cuda/m61-tensor-constant-bench: src/cuda/m61_tensor_constant_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

m61-tensor-constant-bench: build-cuda/m61-tensor-constant-bench
	$<

build-cuda/m61-tensor-bit-bench: src/cuda/m61_tensor_bit_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

m61-tensor-bit-bench: build-cuda/m61-tensor-bit-bench
	$<

build-cuda/q53-m61-bench: src/cuda/q53_m61_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

q53-m61-bench: build-cuda/q53-m61-bench
	$<

build-cuda/m61-persistent-limb-bench: src/cuda/m61_persistent_limb_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

m61-persistent-limb-bench: build-cuda/m61-persistent-limb-bench
	$<

build-cuda/m61-limb-bench: src/cuda/m61_limb_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -Isrc -o $@ src/cuda/m61_limb_bench.cu

build-cuda/near-m61-radix33-bench: src/cuda/near_m61_radix33_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m61-pfa33-edge-bench: src/cuda/m61_pfa33_edge_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m61-batch-bench: src/cuda/m61_batch_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m61-middle-batch-bench: src/cuda/m61_middle_batch_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/q24-tensor-bench: src/cuda/q24_tensor_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/rns31-shape-bench: src/cuda/rns31_shape_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m61-resident-tile-bench: src/cuda/m61_resident_tile_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m61-warp-tail-bench: src/cuda/m61_warp_tail_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-composite-bench: src/cuda/m31_m61_composite_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m19-composite-bench: src/cuda/m31_m19_composite_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-prime-power-bench: src/cuda/m31_prime_power_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-prime-power-tile-bench: src/cuda/m31_prime_power_tile_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-cubic-power-bench: src/cuda/m31_cubic_power_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-direct-carry-bench: src/cuda/m31_m61_direct_carry_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-q3m-tile-bench: src/cuda/m31_m61_q3m_tile_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-q3m-viable-tile-bench: src/cuda/m31_m61_q3m_tile_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=2141192191u -DENABLE_Q3_LAZY=0 -o $@ $<

build-cuda/m31-m61-q3m-small-tile-bench: src/cuda/m31_m61_q3m_tile_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=36700159u -DENABLE_Q3_LAZY=0 -o $@ $<

build-cuda/m31-m61-m19-3m-tile-bench: src/cuda/m31_m61_q3m_tile_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=524287u -DENABLE_Q3_LAZY=0 -o $@ $<

build-cuda/m31-m61-q3m-carry-bench: src/cuda/m31_m61_q3m_carry_bench.cu src/cuda/rns31.cuh
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-q3m-viable-carry-bench: src/cuda/m31_m61_q3m_carry_bench.cu src/cuda/rns31.cuh
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=2141192191u -o $@ $<

build-cuda/m31-m61-q3m-small-carry-bench: src/cuda/m31_m61_q3m_carry_bench.cu src/cuda/rns31.cuh
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=36700159u -o $@ $<

build-cuda/m31-m61-q-radix3-bench: src/cuda/m31_m61_q_radix3_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-q-viable-radix3-bench: src/cuda/m31_m61_q_radix3_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=2141192191u -o $@ $<

build-cuda/m31-m61-q-small-radix3-bench: src/cuda/m31_m61_q_radix3_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=36700159u -o $@ $<

build-cuda/m31-m61-m19-radix3-bench: src/cuda/m31_m61_q_radix3_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=524287u -o $@ $<

build-cuda/m31-m61-m19-radix7-bench: src/cuda/m31_m61_m19_radix7_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-radix31-edge-bench: src/cuda/m31_m61_radix31_edge_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-radix63-edge-bench: src/cuda/m31_m61_radix63_edge_bench.cu src/cuda/m31_m61_m19_radix7_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m61-hartley-tile-bench: src/cuda/m61_hartley_tile_bench.cu src/cuda/rns31_cuda_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/m31-m61-m19-3m-carry-bench: src/cuda/m31_m61_q3m_carry_bench.cu src/cuda/rns31.cuh
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -DQ3_MODULUS=524287u -o $@ $<

build-cuda/m89-limb-bench: src/cuda/m89_limb_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -Isrc -o $@ $<

build-cuda/m127-limb-bench: src/cuda/m127_limb_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -Isrc -o $@ $<

build-cuda/m31-gold-crt-bench: src/cuda/m31_gold_crt_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -Isrc -o $@ $<

build-cuda/fp-compensated-tensor-bench: src/cuda/fp_compensated_tensor_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -Isrc -o $@ src/cuda/fp_compensated_tensor_bench.cu

build-cuda/q31-fixed-fft-bench: src/cuda/q31_fixed_fft_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

build-cuda/q31-fixed-overlap-bench: src/cuda/q31_fixed_overlap_bench.cu src/cuda/q31_fixed_fft_bench.cu
	/usr/local/cuda/bin/nvcc -O3 -std=c++20 -arch=sm_120 -Xptxas=-v -lineinfo -Isrc -o $@ $<

# Instead of linking with libOpenCL, link with libamdocl64
$(BIN)/prpll-amd: ${OBJS}
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ ${OBJS} $(LIBPATH) -lamdocl64 -L/opt/rocm/lib

clean:
	rm -rf build-debug build-release build-cuda

$(BIN)/%.o : src/%.cpp $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

$(BIN)/%.o : src/cuda/%.cpp $(DEPDIR)/%.d
	$(COMPILE.cc) $(OUTPUT_OPTION) $<
	$(POSTCOMPILE)

# src/bundle.cpp is just a wrapping of the OpenCL sources (*.cl) as a C string (as well as the CUDA OpenCL translation code)

src/bundle.cpp: genbundle.sh src/cuda/*.cuh src/cl/*.cl
	bash genbundle.sh $^ > src/bundle.cpp

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

src/version.cpp : src/version.inc

src/version.inc: FORCE
	echo \"`basename \`git describe --tags --long --dirty --always --match v/prpll/*\``\" > $(BIN)/version.new
	diff -q -N $(BIN)/version.new $@ >/dev/null || mv $(BIN)/version.new $@
	echo Version: `cat $@`

FORCE:

include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS1))))
# include $(wildcard $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS2))))
