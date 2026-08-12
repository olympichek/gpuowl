// Copyright (C) Mihai Preda

#include "base.cl"
#include "fftwidth.cl"
#include "carryutil.cl"
#include "weight.cl"
#include "middle.cl"
#include "rieselfused.cl"

void spin() {
#if defined(__has_builtin) && __has_builtin(__builtin_amdgcn_s_sleep)
  __builtin_amdgcn_s_sleep(0);
#elif HAS_ASM
  __asm("s_sleep 0");
#else
  // nothing: just spin
  // on Nvidia: see if there's some brief sleep function
#endif
}

// Increasing WMUL to 2 reduces carryShuttle activity.  This led to a 1% speedup on Titan V.  Testing on other GPUs is needed.
#ifndef WMUL
#define WMUL 2
#endif

#if AMDGPU
#define CarryShuttleAccess(me,i)        ((me) * NW + (i))                       // Generates denser global_load_dwordx4 instructions
//#define CarryShuttleAccess(me,i)      ((me) * 4 + (i)%4 + (i)/4 * 4*G_W)      // Also generates global_load_dwordx4 instructions and unit stride when NW=8
#else
#define CarryShuttleAccess(me,i)        ((me) + (i) * G_W)                      // nVidia likes this unit stride better
#endif

// The last WMUL workgroup's carries have been written to global memory.  Now we shuffle WMUL-1 workgroups carries up using local memory.
void OVERLOAD shufl_carries_up(local void *lds2, i64 *carry, u32 me, u32 lowMe) {
  // If WMUL is one, there is no shuffling of carries
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead of the clean looking if statement below we use the uglier #if
  //if (WMUL == 1) return;
#if WMUL > 1

  const u32 lds_i64s = LDS_BYTES / sizeof(i64);                 // Number of i64s in LDS used by shufl for each WMUL line
  local i64 *lds = (local i64 *) lds2;

  // Handle nasty case where we are writing 8-byte quantities but SHUFL_BYTES_W is only 4 bytes
  if (SHUFL_BYTES_W == 4) {
    if (WMUL == 2) {
      // Full barrier needed as we are using the entire LDS buffer.
      bar();
      // Write the carries.  This will use the entire LDS buffer.
      if (me < G_W) for (i32 i = 0; i < NW; ++i) lds[i * G_W + lowMe] = carry[i];
      // Read carries from previous WMUL workgroup
      bar();
      if (me >= G_W) for (i32 i = 0; i < NW; ++i) carry[i] = lds[i * G_W + lowMe];
      // Full barrier needed as one workgroup just read data from two workgroups LDS buffer.  Not compatible with shufl().
      bar();
    }

    // The really nasty case where all the carries will not fit in LDS memory
    else {
      lds += (me / G_W) * lds_i64s + lowMe;                         // This WMUL workgroup's LDS area
      // Write half the carries to next WMUL's workgroup LDS area
      bar();
      if (me < (WMUL-1) * G_W) for (i32 i = 0; i < NW/2; ++i) lds[lds_i64s + i * G_W] = carry[i];
      // Read carries from our WMUL workgroup LDS area
      bar();
      if (me >= G_W) for (i32 i = 0; i < NW/2; ++i) carry[i] = lds[i * G_W];
      // Write the other half of the carries
      bar();
      if (me < (WMUL-1) * G_W) for (i32 i = 0; i < NW/2; ++i) lds[lds_i64s + i * G_W] = carry[i + NW/2];
      // Read carries from our WMUL workgroup LDS area.  Compatible with shufl, no trailing bar() needed.
      bar();
      if (me >= G_W) for (i32 i = 0; i < NW/2; ++i) carry[i + NW/2] = lds[i * G_W];
    }
  }

  // Easy case.  Write carries to local memory (except last WMUL workgroup which was written to global memory).
  else {
    lds += (me / G_W) * lds_i64s + lowMe;                         // This WMUL workgroup's LDS area
    // Full barrier needed as we are moving data to next WMUL workgroup's LDS area
    bar();
    if (me < (WMUL-1) * G_W) for (i32 i = 0; i < NW; ++i) lds[lds_i64s + i * G_W] = carry[i];
    // Full barrier needed as we just moved data from one WMUL workgroup LDS area to the another WMUL workgroup's LDS area
    bar();
    // Read carries from our WMUL workgroup's LDS area.  This is compatible with shufl and no trailing bar() is required.
    if (me >= G_W) for (i32 i = 0; i < NW; ++i) carry[i] = lds[i * G_W];
  }

#endif
}

// The last WMUL workgroup's carries have been written to global memory.  Now we shuffle WMUL-1 workgroup carries up using local memory.
void OVERLOAD shufl_carries_up(local void *lds2, i32 *carry, u32 me, u32 lowMe) {
  // If WMUL is one, there is no shuffling of carries
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead of the clean looking if statement below we use the uglier #if
  //if (WMUL == 1) return;
#if WMUL > 1

  const u32 lds_i32s = LDS_BYTES / sizeof(i32);                 // Number of i32s in LDS used by shufl for each WMUL line
  local i32 *lds = (local i32 *) lds2;
  lds += (me / G_W) * lds_i32s + lowMe;                         // This WMUL workgroup's LDS area

  // Write carries to local memory (except last WMUL workgroup which was written to global memory)
  // Full barrier needed as we are moving data to next WMUL workgroup's LDS area
  bar();
  if (me < (WMUL-1) * G_W) for (i32 i = 0; i < NW; ++i) lds[lds_i32s + i * G_W] = carry[i];
  // Full barrier needed as we just moved data from one WMUL workgroup LDS area to the another WMUL workgroup's LDS area
  bar();
  // Read carries from our WMUL workgroup's LDS area.  This is compatible with shufl and no trailing bar() is required.
  if (me >= G_W) for (i32 i = 0; i < NW; ++i) carry[i] = lds[i * G_W];

#endif
}


#if FFT_TYPE == FFT64

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(T2) out, CP(T2) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, Trig smallTrig,
                              ConstBigTab CONST_THREAD_WEIGHTS, BigTab THREAD_WEIGHTS, P(uint) bufROE) {
  local T2 lds[WMUL * LDS_BYTES / sizeof(T2)];

  T2 u[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

#if HAS_ASM
  __asm("s_setprio 3");
#endif

  dependentLaunchWait();   // Previous kernel was fftMiddleOutFP64

  readCarryFusedLine(in, u, line, lowMe);

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;
  fft_WIDTH1(lds + zerohack, u, smallTrig + zerohack, WMUL, lowMe);

  Word2 wu[NW];
#if !NVIDIAGPU || CUDA_BACKEND
  T2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), TSLOAD(&THREAD_WEIGHTS[G_W + line]));
#else
  T2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), CONST_THREAD_WEIGHTS[line % 64]);
  weights.x = optionalDouble(weights.x);
  weights.y = optionalHalve(weights.y);
  weights = fancyMul(weights, CONST_THREAD_WEIGHTS[64 + line / 64]);
#endif

#if MUL3
  P(i64) carryShuttlePtr = (P(i64)) carryShuttle;
  i64 carry[NW+1];
#else
  P(CFcarry) carryShuttlePtr = (P(CFcarry)) carryShuttle;
  CFcarry carry[NW+1];
#endif

  float roundMax = 0;
  float carryMax = 0;

  // Calculate the most significant 32-bits of FRAC_BPW * the word index.  Also add FRAC_BPW_HI to test first biglit flag.
  u32 word_index = (lowMe * H + line) * 2;
  u32 frac_bits = fracBits(word_index) + FRAC_BPW_HI;
  const u32 frac_bits_bigstep = fracBits(G_W * H * 2);
  u32 starting_frac_bits = frac_bits;

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  T invBase = optionalDouble(weights.x);
  for (u32 i = 0; i < NW; ++i) {
    T invWeight1 = i == 0 ? invBase : optionalDouble(fancyMul(invBase, iweightStep(i)));
    T invWeight2 = optionalDouble(fancyMul(invWeight1, IWEIGHT_STEP));

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(u[i]), U2(invWeight1, invWeight2),
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    // Generate frac_bits for next pair
    frac_bits += frac_bits_bigstep;
  }
  frac_bits = starting_frac_bits;     // Restore starting frac_bits for applying weights after carry propagation

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) {
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, roundMax);
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Calculate inverse weights
  T base = optionalHalve(weights.y);
  for (u32 i = 0; i < NW; ++i) {
    T weight1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)));
    T weight2 = optionalHalve(fancyMul(weight1, WEIGHT_STEP));
    u[i] = U2(weight1, weight2);
  }

  // Shuffle carries up
  shufl_carries_up(lds, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H/WMUL we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words
  for (i32 i = 0; i < NW; ++i) {
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
    u[i] = U2(u[i].x * wu[i].x, u[i].y * wu[i].y);

    // Generate frac_bits for next pair
    frac_bits += frac_bits_bigstep;
  }

  dependentLaunch();   // Next kernel will be fftMiddleInFP64

  fft_WIDTH2(lds, u, smallTrig, WMUL, lowMe);
  writeCarryFusedLine(u, out, line, lowMe);
}


/**************************************************************************/
/*            Similar to above, but for an FFT based on FP32              */
/**************************************************************************/

#elif FFT_TYPE == FFT32

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(F2) out, CP(F2) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, TrigFP32 smallTrig,
                              ConstBigTabFP32 CONST_THREAD_WEIGHTS, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE) {
  local F2 lds[WMUL * LDS_BYTES / sizeof(F2)];

  F2 u[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

#if HAS_ASM
  __asm("s_setprio 3");
#endif

  dependentLaunchWait();   // Previous kernel was fftMiddleOutFP32

  readCarryFusedLine(in, u, line, lowMe);

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;
  fft_WIDTH1(lds + zerohack, u, smallTrig + zerohack, WMUL, lowMe);

  Word2 wu[NW];
  u32 me_frac_bits = fracBits(lowMe * H * 2);
#if !NVIDIAGPU || CUDA_BACKEND
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), TSLOAD(&THREAD_WEIGHTS[G_W + line]));
  u32 line_frac_bits = fracBits(line * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > line_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > line_frac_bits);
#else
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), CONST_THREAD_WEIGHTS[line % 64]);
  u32 partialLine_frac_bits = fracBits((line % 64) * 2);
  u32 base_frac_bits = me_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
  weights = fancyMul(weights, CONST_THREAD_WEIGHTS[64 + line / 64]);
  partialLine_frac_bits = fracBits(((line / 64) * 64) * 2);
  base_frac_bits = base_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
#endif

  P(CFcarry) carryShuttlePtr = (P(CFcarry)) carryShuttle;
  CFcarry carry[NW+1];

  float roundMax = 0;
  float carryMax = 0;

  // Calculate the most significant 32-bits of FRAC_BPW * the word index (it's the same as base_frac_bits).
  u32 word_index = (lowMe * H + line) * 2;
  const u32 frac_bits_bigstep = fracBits(G_W * H * 2 - 1);

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  F invBase = weights.x;
  u32 frac_bits = base_frac_bits;
  for (u32 i = 0; i < NW; ++i) {
    F invWeight1 = i == 0 ? invBase : optionalDouble(fancyMul(invBase, iweightStep(i)), frac_bits > base_frac_bits);
    F invWeight2 = optionalDouble(fancyMul(invWeight1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    frac_bits += FRAC_BPW_HI;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(u[i]), U2(invWeight1, invWeight2),
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    // Generate frac_bits for next pair
    frac_bits += frac_bits_bigstep;
  }

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) { 
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, roundMax);
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Shuffle carries up
  shufl_carries_up(lds, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words
  F base = weights.y;
  frac_bits = base_frac_bits;
  for (i32 i = 0; i < NW; ++i) {
    // Calculate inverse weights
    F weight1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F weight2 = optionalHalve(fancyMul(weight1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    frac_bits += FRAC_BPW_HI;
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
    u[i] = U2(weight1 * wu[i].x, weight2 * wu[i].y);

    // Generate frac_bits for next pair
    frac_bits += frac_bits_bigstep;
  }

  dependentLaunch();   // Next kernel will be fftMiddleInFP32

  fft_WIDTH2(lds, u, smallTrig, WMUL, lowMe);
  writeCarryFusedLine(u, out, line, lowMe);
}


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M31^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT31

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(GF31) out, CP(GF31) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, TrigGF31 smallTrig, P(uint) bufROE) {
  local GF31 lds[WMUL * LDS_BYTES / sizeof(GF31)];

  GF31 u[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

#if HAS_ASM
  __asm("s_setprio 3");
#endif

  dependentLaunchWait();   // Previous kernel was fftMiddleOutGF31

  readCarryFusedLine(in, u, line, lowMe);

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;
  fft_WIDTH1(lds + zerohack, u, smallTrig + zerohack, WMUL, lowMe);

  Word2 wu[NW];

  P(CFcarry) carryShuttlePtr = (P(CFcarry)) carryShuttle;
  CFcarry carry[NW+1];

  u32 roundMax = 0;
  float carryMax = 0;

  u32 word_index = (lowMe * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Weights can be applied with shifts because 2 is the 60th root GF31.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 31;
  u64 starting_combo_counter = combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = weight_shift + log2_NWORDS + 1;
  if (weight_shift > 31) weight_shift -= 31;

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  for (u32 i = 0; i < NW; ++i) {
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(u[i]), weight_shift0, weight_shift1,
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }
  combo_counter = starting_combo_counter;     // Restore starting counter for applying weights after carry propagation

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) { 
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  float fltRoundMax = (float) roundMax / (float) M31;      // For speed, roundoff was computed as 32-bit integer.  Convert to float.
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, fltRoundMax);
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Shuffle carries up
  shufl_carries_up(lds, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words.  Apply weights.
  for (i32 i = 0; i < NW; ++i) {
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;
    // Generate big-word/little-word flag, propagate final carry
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
    u[i] = U2(shl(make_Z31(wu[i].x), weight_shift0), shl(make_Z31(wu[i].y), weight_shift1));
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }

  dependentLaunch();   // Next kernel will be fftMiddleInGF31

  fft_WIDTH2(lds, u, smallTrig, WMUL, lowMe);
  writeCarryFusedLine(u, out, line, lowMe);
}


/**************************************************************************/
/*          Similar to above, but for an NTT based on GF(M61^2)           */
/**************************************************************************/

#elif FFT_TYPE == FFT61

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(GF61) out, CP(GF61) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, TrigGF61 smallTrig, P(uint) bufROE) {
  local GF61 lds[WMUL * LDS_BYTES / sizeof(GF61)];

  GF61 u[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

#if HAS_ASM
  __asm("s_setprio 3");
#endif

  dependentLaunchWait();   // Previous kernel was fftMiddleOutGF61

  readCarryFusedLine(in, u, line, lowMe);

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;
  fft_WIDTH1(lds + zerohack, u, smallTrig + zerohack, WMUL, lowMe);

  Word2 wu[NW];

#if MUL3
  P(i64) carryShuttlePtr = (P(i64)) carryShuttle;
  i64 carry[NW+1];
#else
  P(CFcarry) carryShuttlePtr = (P(CFcarry)) carryShuttle;
  CFcarry carry[NW+1];
#endif

  u32 roundMax = 0;
  float carryMax = 0;

  u32 word_index = (lowMe * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Weights can be applied with shifts because 2 is the 60th root GF61.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 61.
  const u32 log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 61;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 61;
  u64 starting_combo_counter = combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = weight_shift + log2_NWORDS + 1;
  if (weight_shift > 61) weight_shift -= 61;

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  for (u32 i = 0; i < NW; ++i) {
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(u[i]), weight_shift0, weight_shift1,
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 61) weight_shift -= 61;
  }
  combo_counter = starting_combo_counter;     // Restore starting counter for applying weights after carry propagation

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) { 
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  float fltRoundMax = (float) roundMax / (float) (M61 >> 32);      // For speed, roundoff was computed as 32-bit integer.  Convert to float.
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, fltRoundMax);
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Shuffle carries up
  shufl_carries_up(lds, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words.  Apply weights.
  for (i32 i = 0; i < NW; ++i) {
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;
    // Generate big-word/little-word flag, propagate final carry
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
    u[i] = U2(shl(make_Z61(wu[i].x), weight_shift0), shl(make_Z61(wu[i].y), weight_shift1));
    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 61) weight_shift -= 61;
  }

  dependentLaunch();   // Next kernel will be fftMiddleInGF61

  fft_WIDTH2(lds, u, smallTrig, WMUL, lowMe);
  writeCarryFusedLine(u, out, line, lowMe);
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP64 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT6431

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(T2) out, CP(T2) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, Trig smallTrig,
                              ConstBigTab CONST_THREAD_WEIGHTS, BigTab THREAD_WEIGHTS, P(uint) bufROE) {
  local T2 lds[WMUL * LDS_BYTES / sizeof(T2)];
  local GF31 *lds31 = (local GF31 *) lds;

  T2 u[NW];
  GF31 u31[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);
  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);

#if HAS_ASM
  __asm("s_setprio 3");
#endif

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;

  readCarryFusedLine(in, u, line, lowMe);
  fft_WIDTH1(lds + zerohack, u, smallTrig + zerohack, WMUL, lowMe);

  dependentLaunchWait();   // Previous kernel was fftMiddleOutGF31

  readCarryFusedLine(in31, u31, line, lowMe);
  fft_WIDTH1(lds31 + zerohack, u31, smallTrig31 + zerohack, WMUL, lowMe);

  Word2 wu[NW];
#if !NVIDIAGPU || CUDA_BACKEND
  T2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), TSLOAD(&THREAD_WEIGHTS[G_W + line]));
#else
  T2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), CONST_THREAD_WEIGHTS[line % 64]);
  weights.x = optionalDouble(weights.x);
  weights.y = optionalHalve(weights.y);
  weights = fancyMul(weights, CONST_THREAD_WEIGHTS[64 + line / 64]);
#endif

  P(i64) carryShuttlePtr = (P(i64)) carryShuttle;
  i64 carry[NW+1];

  float roundMax = 0;
  float carryMax = 0;

  u32 word_index = (lowMe * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 31;
  u64 starting_combo_counter = combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = weight_shift + log2_NWORDS + 1;
  if (weight_shift > 31) weight_shift -= 31;

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  T invBase = optionalDouble(weights.x);
  for (u32 i = 0; i < NW; ++i) {
    // Generate the FP64 weights and second GF31 weight shift
    T invWeight1 = i == 0 ? invBase : optionalDouble(fancyMul(invBase, iweightStep(i)));
    T invWeight2 = optionalDouble(fancyMul(invWeight1, IWEIGHT_STEP));
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(u[i]), SWAP_XY(u31[i]), invWeight1, invWeight2, weight_shift0, weight_shift1,
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      LL != 0, (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }
  combo_counter = starting_combo_counter;     // Restore starting counter for applying weights after carry propagation

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

  // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) { 
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, roundMax);
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Calculate inverse weights
  T base = optionalHalve(weights.y);
  for (u32 i = 0; i < NW; ++i) {
    T weight1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)));
    T weight2 = optionalHalve(fancyMul(weight1, WEIGHT_STEP));
    u[i] = U2(weight1, weight2);
  }

  // Shuffle carries up
  shufl_carries_up(lds, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H/WMUL we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words.  Apply weights.
  for (i32 i = 0; i < NW; ++i) {
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;
    // Generate big-word/little-word flag, propagate final carry
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
    u[i] = U2(u[i].x * wu[i].x, u[i].y * wu[i].y);
    u31[i] = U2(shl(make_Z31(wu[i].x), weight_shift0), shl(make_Z31(wu[i].y), weight_shift1));

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }

  fft_WIDTH2(lds, u, smallTrig, WMUL, lowMe);
  writeCarryFusedLine(u, out, line, lowMe);

  dependentLaunch();   // Next kernel will be fftMiddleInFP32

  fft_WIDTH2(lds31, u31, smallTrig31, WMUL, lowMe);
  writeCarryFusedLine(u31, out31, line, lowMe);
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M31^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3231

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(T2) out, CP(T2) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, Trig smallTrig,
                              ConstBigTabFP32 CONST_THREAD_WEIGHTS, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE) {
  local F2 ldsF2[WMUL * LDS_BYTES / sizeof(F2)];
  local GF31 *lds31 = (local GF31 *) ldsF2;

  F2 uF2[NW];
  GF31 u31[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

  CP(F2) inF2 = (CP(F2)) in;
  P(F2) outF2 = (P(F2)) out;
  TrigFP32 smallTrigF2 = (TrigFP32) smallTrig;
  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);
  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);

#if HAS_ASM
  __asm("s_setprio 3");
#endif

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;

  readCarryFusedLine(inF2, uF2, line, lowMe);
  fft_WIDTH1(ldsF2 + zerohack, uF2, smallTrigF2 + zerohack, WMUL, lowMe);

  dependentLaunchWait();   // Previous kernel was fftMiddleOutGF31

  readCarryFusedLine(in31, u31, line, lowMe);
  fft_WIDTH1(lds31 + zerohack, u31, smallTrig31 + zerohack, WMUL, lowMe);

  Word2 wu[NW];
  u32 me_frac_bits = fracBits(lowMe * H * 2);
#if !NVIDIAGPU || CUDA_BACKEND
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), TSLOAD(&THREAD_WEIGHTS[G_W + line]));
  u32 line_frac_bits = fracBits(line * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > line_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > line_frac_bits);
#else
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), CONST_THREAD_WEIGHTS[line % 64]);
  u32 partialLine_frac_bits = fracBits((line % 64) * 2);
  u32 base_frac_bits = me_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
  weights = fancyMul(weights, CONST_THREAD_WEIGHTS[64 + line / 64]);
  partialLine_frac_bits = fracBits(((line / 64) * 64) * 2);
  base_frac_bits = base_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
#endif

  P(i32) carryShuttlePtr = (P(i32)) carryShuttle;
  i32 carry[NW+1];

  float roundMax = 0;
  float carryMax = 0;

  u32 word_index = (lowMe * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
  const u32 log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 31;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 30) % 31;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 31;
  u64 starting_combo_counter = combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = weight_shift + log2_NWORDS + 1;
  if (weight_shift > 31) weight_shift -= 31;

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  F invBase = weights.x;
  for (u32 i = 0; i < NW; ++i) {
    // Generate the FP32 weights and second GF31 weight shift
    F invWeight1 = i == 0 ? invBase : optionalDouble(fancyMul(invBase, iweightStep(i)), frac_bits > base_frac_bits);
    F invWeight2 = optionalDouble(fancyMul(invWeight1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(uF2[i]), SWAP_XY(u31[i]), invWeight1, invWeight2, weight_shift0, weight_shift1,
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }
  combo_counter = starting_combo_counter;     // Restore starting counter for applying weights after carry propagation

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) { 
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  updateStats((local u32 *) ldsF2, G_W * WMUL, H / WMUL, bufROE, posROE, roundMax);
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) ldsF2, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Shuffle carries up
  shufl_carries_up(ldsF2, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H/WMUL we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words.  Apply weights.
  F base = weights.y;
  for (i32 i = 0; i < NW; ++i) {
    F weight1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F weight2 = optionalHalve(fancyMul(weight1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 31) weight_shift -= 31;
    u32 weight_shift1 = weight_shift;
    // Generate big-word/little-word flag, propagate final carry
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
    uF2[i] = U2(weight1 * wu[i].x, weight2 * wu[i].y);
    u31[i] = U2(shl(make_Z31(wu[i].x), weight_shift0), shl(make_Z31(wu[i].y), weight_shift1));

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 31) weight_shift -= 31;
  }

  fft_WIDTH2(ldsF2, uF2, smallTrigF2, WMUL, lowMe);
  writeCarryFusedLine(uF2, outF2, line, lowMe);

  dependentLaunch();   // Next kernel will be fftMiddleInFP32

  fft_WIDTH2(lds31, u31, smallTrig31, WMUL, lowMe);
  writeCarryFusedLine(u31, out31, line, lowMe);
}


/**************************************************************************/
/*    Similar to above, but for a hybrid FFT based on FP32 & GF(M61^2)    */
/**************************************************************************/

#elif FFT_TYPE == FFT3261

#if PARITY_PACKED
u32 parityBallot(bool predicate) {
#if CUDA_BACKEND
  return __ballot_sync(0xffffffffu, predicate);
#else
  return predicate ? 1u : 0u;
#endif
}
#endif

#if PARITY_DIRECT
void storeNextExpectedParity(P(u32) parity, u32 sourceWord, u32 bit) {
  u32 const outputPair = sourceWord & (NWORDS / 2 - 1);
  u32 const frac = ((u32) EXP * outputPair) & (NWORDS - 1);
  bool const selectsUpper = frac != 0 && frac <= NWORDS / 2;
  bool const isUpper = sourceWord >= NWORDS / 2;
  if (selectsUpper == isUpper) {
    u32 const x = outputPair / BIG_HEIGHT;
    u32 const line = outputPair - x * BIG_HEIGHT;
    parity[line * WIDTH + x] = bit & 1u;
  }
}
#endif

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(T2) out, CP(T2) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, Trig smallTrig,
                              ConstBigTabFP32 CONST_THREAD_WEIGHTS, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE
#if PARITY_SQUARE
                              , CP(u32) parityIn, P(u32) parityOut
#endif
#if FOLD_SYNDROME
                              , P(GF31) foldOut
#if FOLD_SYNDROME_VALIDATE
                              , P(Word2) foldWords
#endif
#endif
                              ) {
  local GF61 lds61[WMUL * LDS_BYTES / sizeof(GF61)];
  local F2 *ldsF2 = (local F2 *) lds61;

  F2 uF2[NW];
  GF61 u61[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

  CP(F2) inF2 = (CP(F2)) in;
  P(F2) outF2 = (P(F2)) out;
  TrigFP32 smallTrigF2 = (TrigFP32) smallTrig;
  CP(GF61) in61 = (CP(GF61)) (in + DISTGF61);
  P(GF61) out61 = (P(GF61)) (out + DISTGF61);
  TrigGF61 smallTrig61 = (TrigGF61) (smallTrig + DISTWTRIGGF61);

#if HAS_ASM
  __asm("s_setprio 3");
#endif

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;

  readCarryFusedLine(inF2, uF2, line, lowMe);
  fft_WIDTH1(ldsF2 + zerohack, uF2, smallTrigF2 + zerohack, WMUL, lowMe);

  dependentLaunchWait();   // Previous kernel was fftMiddleOutGF61

  readCarryFusedLine(in61, u61, line, lowMe);
  fft_WIDTH1(lds61 + zerohack, u61, smallTrig61 + zerohack, WMUL, lowMe);

  Word2 wu[NW];
  u32 me_frac_bits = fracBits(lowMe * H * 2);
#if !NVIDIAGPU || CUDA_BACKEND
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), TSLOAD(&THREAD_WEIGHTS[G_W + line]));
  u32 line_frac_bits = fracBits(line * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > line_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > line_frac_bits);
#else
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), CONST_THREAD_WEIGHTS[line % 64]);
  u32 partialLine_frac_bits = fracBits((line % 64) * 2);
  u32 base_frac_bits = me_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
  weights = fancyMul(weights, CONST_THREAD_WEIGHTS[64 + line / 64]);
  partialLine_frac_bits = fracBits(((line / 64) * 64) * 2);
  base_frac_bits = base_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
#endif

  P(i64) carryShuttlePtr = (P(i64)) carryShuttle;
  i64 carry[NW+1];

  float roundMax = 0;
  float carryMax = 0;
#if ROE && ROE_COUNT
  u32 riskyPairs = 0;
#endif

  u32 word_index = (lowMe * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 61.
  const u32 log2_root_two = (u32) (((1ULL << 60) / NWORDS) % 61);
  const u32 bigword_weight_shift = (NWORDS - EXP % NWORDS) * log2_root_two % 61;
  const u32 bigword_weight_shift_minus1 = (bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } combo;
#define frac_bits       combo.a[0]
#define weight_shift    combo.a[1]
#define combo_counter   combo.b

  const u64 combo_step = make_u64(bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  combo_counter = comboFracBits(word_index) + make_u64(word_index * bigword_weight_shift_minus1, 0xFFFFFFFF);
  weight_shift = weight_shift % 61;
  u64 starting_combo_counter = combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  weight_shift = weight_shift + log2_NWORDS + 1;
  if (weight_shift > 61) weight_shift -= 61;

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  F invBase = weights.x;
  for (u32 i = 0; i < NW; ++i) {
    // Generate the FP32 weights and second GF61 weight shift
    F invWeight1 = i == 0 ? invBase : optionalDouble(fancyMul(invBase, iweightStep(i)), frac_bits > base_frac_bits);
    F invWeight2 = optionalDouble(fancyMul(invWeight1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);

    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

#if PARITY_SQUARE
#if PARITY_PHYSICAL
    u32 parityIndex = (lowMe + i * G_W) * H + line;
#elif PARITY_PREPARED
    u32 parityIndex = line * WIDTH + lowMe + i * G_W;
#else
    u32 parityIndex = (lowMe + i * G_W) * H + line;
#endif
#endif

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
#if ROE && ROE_COUNT
    float pairRoundMax = 0;
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(uF2[i]), SWAP_XY(u61[i]), invWeight1, invWeight2, weight_shift0, weight_shift1,
#if PARITY_SQUARE
#if PARITY_LAZY
                      parityIn, parityIndex,
#else
                      squareCoefficientParity(parityIn, parityIndex), 0,
#endif
#endif
                      LL != 0, (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &pairRoundMax, &carryMax);
    roundMax = max(roundMax, pairRoundMax);
    riskyPairs += pairRoundMax >= (float) ROE_COUNT * 0.001f;
#else
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(uF2[i]), SWAP_XY(u61[i]), invWeight1, invWeight2, weight_shift0, weight_shift1,
#if PARITY_SQUARE
#if PARITY_LAZY
                      parityIn, parityIndex,
#else
                      squareCoefficientParity(parityIn, parityIndex), 0,
#endif
#endif
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      LL != 0, (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);
#endif

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 61) weight_shift -= 61;
  }
  combo_counter = starting_combo_counter;     // Restore starting counter for applying weights after carry propagation

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) { 
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
#if ROE_COUNT
  updateCountStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, riskyPairs);
#else
  updateStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, roundMax);
#endif
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Shuffle carries up
  shufl_carries_up(lds61, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H/WMUL we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words.  Apply weights.
  F base = weights.y;
#if FOLD_SYNDROME
  // Build the next iteration's full-N M31-weighted input, folded by eight.
  // The production carry geometry already gives this thread exactly the
  // eight aliases separated by NWORDS/8 words, so no atomics are required.
  const u32 fold_log2_root_two = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 fold_bigword_shift = (NWORDS - EXP % NWORDS) * fold_log2_root_two % 31;
  const u32 fold_shift_step = (fold_bigword_shift + 30) % 31;
  const u64 fold_combo_step = make_u64(fold_shift_step, FRAC_BPW_HI);
  const u64 fold_combo_bigstep =
    (comboFracBits(G_W * H * 2 - 1) +
     make_u64((G_W * H * 2 - 1) * fold_shift_step, 0)) % (31ULL << 32);
  u64 fold_combo = comboFracBits(word_index) +
                   make_u64(word_index * fold_shift_step, 0xFFFFFFFF);
  fold_combo = make_u64(hi32(fold_combo) % 31, lo32(fold_combo));
#if FOLD_FACTOR == 4
  GF31 folded0 = U2((Z31)0, (Z31)0);
  GF31 folded1 = U2((Z31)0, (Z31)0);
#else
  GF31 folded = U2((Z31)0, (Z31)0);
#endif
#endif
  for (i32 i = 0; i < NW; ++i) {
    // Calculate inverse weights
    F weight1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F weight2 = optionalHalve(fancyMul(weight1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    // Generate the second weight shift
    u32 weight_shift0 = weight_shift;
    combo_counter += combo_step;
    if (weight_shift > 61) weight_shift -= 61;
    u32 weight_shift1 = weight_shift;
    // Generate big-word/little-word flag, propagate final carry
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
#if FOLD_SYNDROME
    u32 fold_shift0 = hi32(fold_combo);
    fold_combo += fold_combo_step;
    if (hi32(fold_combo) > 31) fold_combo -= 31ULL << 32;
    u32 fold_shift1 = hi32(fold_combo);
    GF31 const foldValue =
      U2(shl(make_Z31(wu[i].x), fold_shift0),
         shl(make_Z31(wu[i].y), fold_shift1));
#if FOLD_FACTOR == 4
    if (i & 1) folded1 = add(folded1, foldValue);
    else       folded0 = add(folded0, foldValue);
#else
    folded = add(folded, foldValue);
#endif
#if FOLD_SYNDROME_VALIDATE
    foldWords[(lowMe + i * G_W) * H + line] = wu[i];
#endif
    fold_combo += fold_combo_bigstep;
    if (hi32(fold_combo) > 31) fold_combo -= 31ULL << 32;
#endif
    uF2[i] = U2(weight1 * wu[i].x, weight2 * wu[i].y);
    u61[i] = U2(shl(make_Z61(wu[i].x), weight_shift0), shl(make_Z61(wu[i].y), weight_shift1));
#if PARITY_SQUARE
#if PARITY_PACKED
    u32 const physicalPair = line * WIDTH + lowMe + i * G_W;
    u32 parityMask = parityBallot(((u32) wu[i].x & 1u) != 0);
    if ((me & 31u) == 0) parityOut[physicalPair >> 5] = parityMask;
    parityMask = parityBallot(((u32) wu[i].y & 1u) != 0);
    if ((me & 31u) == 0) parityOut[ND / 32 + (physicalPair >> 5)] = parityMask;
#else
    u32 const logicalPair = (lowMe + i * G_W) * H + line;
#if PARITY_DIRECT
    storeNextExpectedParity(parityOut, 2 * logicalPair, (u32) wu[i].x);
    storeNextExpectedParity(parityOut, 2 * logicalPair + 1, (u32) wu[i].y);
#elif PARITY_PREPARED || PARITY_PHYSICAL
    parityOut[line * WIDTH + lowMe + i * G_W] =
      ((u32) wu[i].x & 1u) | (((u32) wu[i].y & 1u) << 1);
#else
    parityOut[logicalPair] = ((u32) wu[i].x & 1u) | (((u32) wu[i].y & 1u) << 1);
#endif
#endif
#endif

    // Generate weight shifts and frac_bits for next pair
    combo_counter += combo_bigstep;
    if (weight_shift > 61) weight_shift -= 61;
  }

#if FOLD_SYNDROME
#if FOLD_TRANSFORM
  const u32 fold_pair = lowMe * H + line;
#if FOLD_FACTOR == 4
  // The main carry owns eight aliases spaced N/8 words apart.  A factor-four
  // fold produces two bins: even i and odd i.  Reinterpret both natural
  // indices in the side transform's 512x2x512 geometry and transpose directly
  // into the forward-width input layout.
  const u32 fold_pair1 = fold_pair + NWORDS / 16;
  foldOut[(fold_pair & 1023u) * 512u + (fold_pair >> 10)] = folded0;
  foldOut[(fold_pair1 & 1023u) * 512u + (fold_pair1 >> 10)] = folded1;
#elif FOLD_FACTOR == 16
  // Preserve the natural N/8 partial-fold order.  The side stream combines
  // pairs into N/16 and transposes them after this critical carry has exited.
  foldOut[fold_pair] = folded;
#else
  // The N/8 side transform's forward width kernel consumes contiguous lines,
  // so transpose while all indices are already available.
#if FOLD_SHAPE == 2
  foldOut[(fold_pair & 511u) * 512u + (fold_pair >> 9)] = folded;
#else
  foldOut[(fold_pair & 1023u) * 256u + (fold_pair >> 10)] = folded;
#endif
#endif
#else
  foldOut[lowMe * H + line] = folded;
#endif
#endif

  fft_WIDTH2(ldsF2, uF2, smallTrigF2, WMUL, lowMe);
  writeCarryFusedLine(uF2, outF2, line, lowMe);

  dependentLaunch();   // Next kernel will be fftMiddleInFP32

  fft_WIDTH2(lds61, u61, smallTrig61, WMUL, lowMe);
  writeCarryFusedLine(u61, out61, line, lowMe);
}


/**************************************************************************/
/*    Similar to above, but for an NTT based on GF(M31^2)*GF(M61^2)       */
/**************************************************************************/

#elif FFT_TYPE == FFT3161

#ifndef CARRY_NOWAIT
#define CARRY_NOWAIT 0     // Timing scaffold: skip the carry shuttle (WRONG results)
#endif
#ifndef CARRY_EARLY
#define CARRY_EARLY 0      // BROKEN: early spin livelocks the grid; do not enable (kept as evidence)
#endif
#ifndef CARRY_ACQREL
#define CARRY_ACQREL 0     // Release/acquire flag handshake instead of full device fences (exact)
#endif

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(T2) out, CP(T2) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, Trig smallTrig
#if GOLD_PAIR
                              , BigTab THREAD_WEIGHTS
#endif
                              , P(uint) bufROE
#if DETACHED_M31_EDGE
                              // Exact external M31 width transforms: fftWGF31 writes this
                              // buffer flat before this kernel, and fftWOut31 consumes the
                              // flat pre-width residues this kernel writes back to it.
                              , P(T2) detach
#endif
                              ) {
  local GF61 lds61[WMUL * LDS_BYTES / sizeof(GF61)];
  local GF31 *lds31 = (local GF31 *) lds61;

#if !DETACHED_M31_EDGE || DETACHED_M31_EDGE == 2
  GF31 u31[NW];
#endif
  GF61 u61[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);
  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);
  CP(GF61) in61 = (CP(GF61)) (in + DISTGF61);
  P(GF61) out61 = (P(GF61)) (out + DISTGF61);
  TrigGF61 smallTrig61 = (TrigGF61) (smallTrig + DISTWTRIGGF61);

#if HAS_ASM
  __asm("s_setprio 3");
#endif

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;

#if !DETACHED_M31_EDGE
  readCarryFusedLine(in31, u31, line, lowMe);
  fft_WIDTH1(lds31 + zerohack, u31, smallTrig31 + zerohack, WMUL, lowMe);
#endif

  dependentLaunchWait();   // Previous kernel was fftMiddleOutGF61

  readCarryFusedLine(in61, u61, line, lowMe);
  fft_WIDTH1(lds61 + zerohack, u61, smallTrig61 + zerohack, WMUL, lowMe);

#if DETACHED_M31_EDGE
  // Exact external-width designs supply this already transformed line.  This
  // opt-in timing gate deliberately skips that producer until the actual
  // carry saving is large enough to justify the scheduler integration.
#endif

  Word2 wu[NW];

  P(i64) carryShuttlePtr = (P(i64)) carryShuttle;
  i64 carry[NW+1];

  u32 roundMax = 0;
  float carryMax = 0;

  u32 word_index = (lowMe * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
#if GOOD_THOMAS3
  const u32 m31_log2_root_two = 21;
#else
  const u32 m31_log2_root_two = (u32)(((1ULL << 30) / NWORDS) % 31);
#endif
  const u32 m31_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_bigword_weight_shift_minus1 = (m31_bigword_weight_shift + 30) % 31;
#if !GOLD_PAIR
#if GOOD_THOMAS3
  const u32 m61_log2_root_two = 45;
#else
  const u32 m61_log2_root_two = (u32)(((1ULL << 60) / NWORDS) % 61);
#endif
  const u32 m61_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m61_log2_root_two % 61;
  const u32 m61_bigword_weight_shift_minus1 = (m61_bigword_weight_shift + 60) % 61;
#endif

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } m31_combo;
#define frac_bits           m31_combo.a[0]
#define m31_weight_shift    m31_combo.a[1]
#define m31_combo_counter   m31_combo.b
#if !GOLD_PAIR
  union { uint2 a; u64 b; } m61_combo;
#define m61_weight_shift    m61_combo.a[1]
#define m61_combo_counter   m61_combo.b
#endif

  const u64 m31_combo_step = make_u64(m31_bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 m31_combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * m31_bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  m31_combo_counter = comboFracBits(word_index) + make_u64(word_index * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m31_weight_shift = m31_weight_shift % 31;
  u64 m31_starting_combo_counter = m31_combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation
#if GOLD_PAIR
  u32 goldInvExponent;
  Z61 goldInvWeight = goldStartingWeight(
    THREAD_WEIGHTS, lowMe, line, true, &goldInvExponent);
  goldInvWeight = mul(goldInvWeight, (Z61)GOLD_INV_ND);
#else
  const u64 m61_combo_step = make_u64(m61_bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 m61_combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * m61_bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  m61_combo_counter = comboFracBits(word_index) + make_u64(word_index * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m61_weight_shift = m61_weight_shift % 61;
  u64 m61_starting_combo_counter = m61_combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation
#endif

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift + log2_NWORDS + 1);
#if !GOLD_PAIR
  m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift + log2_NWORDS + 1);
#endif

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  for (u32 i = 0; i < NW; ++i) {
    // Generate the second weight shifts
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
#if GOLD_PAIR
    Z61 const goldInvWeight0 = goldInvWeight;
    u32 goldOddExponent = goldInvExponent;
    Z61 const goldInvWeight1 = advanceGoldInverse(
      goldInvWeight, GOLD_DELTA_ONE, (Z61)GOLD_INV_ONE,
      &goldOddExponent);
#else
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;
#endif

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
#if DETACHED_M31_EDGE
    GF31 const detached31 = ((CP(GF31)) detach)[(line * NW + i) * G_W + lowMe];
#endif
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(
#if DETACHED_M31_EDGE
                      detached31
#else
                      u31[i]
#endif
                      ),
#if GOLD_PAIR
                      u61[i],
#else
                      SWAP_XY(u61[i]),
#endif
                      m31_weight_shift0, m31_weight_shift1,
#if GOLD_PAIR
                      goldInvWeight0, goldInvWeight1,
#else
                      m61_weight_shift0, m61_weight_shift1,
#endif
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      LL != 0, (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    // Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_bigstep;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
#if GOLD_PAIR
    goldInvWeight = advanceGoldInverse(
      goldInvWeight, GOLD_DELTA_X, (Z61)GOLD_INV_X,
      &goldInvExponent);
#else
    m61_combo_counter += m61_combo_bigstep;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
#endif
  }
  m31_combo_counter = m31_starting_combo_counter;     // Restore starting counter for applying weights after carry propagation
#if !GOLD_PAIR
  m61_combo_counter = m61_starting_combo_counter;
#endif

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H && !CARRY_NOWAIT) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W && !CARRY_NOWAIT) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
#if CARRY_ACQREL
    if (lowMe % WAVEFRONT == 0) {
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store_explicit((atomic_uint *) &ready[pos], 1, memory_order_release, memory_scope_device);
    }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) {
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

#if CARRY_EARLY
  // Early shuttle wait+load: overlap the L2 round-trip latency with the
  // weights/stats/shuffle work below.  Same protocol and values as the
  // original late wait; incoming carries land in inCarry and are merged
  // after shufl_carries_up.
  CFcarry inCarry[NW];
  if (me < G_W) {
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) { inCarry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]); }
    }
  }
  // Full-block barrier for the rotated-read group keeps barrier pairing
  // uniform across warps (a conditional bar() here deadlocks the block).
  if (gr >= H / WMUL) {
    bar();
    if (me < G_W) {
      for (i32 i = 0; i < NW; ++i) { inCarry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i)]); }
      if (me == 0) {
        CFcarry t = inCarry[NW-1];
        for (i32 i = NW-1; i; --i) { inCarry[i] = inCarry[i-1]; }
        inCarry[0] = t;
      }
    }
  }
#endif

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  float fltRoundMax = (float) roundMax / (float) 0x1FFFFFFF;      // For speed, roundoff was computed as 32-bit integer.  Convert to float - divide by M61.
  updateStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, fltRoundMax);
#elif SPIN_STATS
  // Deferred until after the readiness wait below.
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Shuffle carries up
  shufl_carries_up(lds61, carry, me, lowMe);

  // Wait until our carries are ready
#if SPIN_STATS && !ROE
  u32 spinCount = 0;
#endif
#if CARRY_EARLY
  if (me < G_W) { for (i32 i = 0; i < NW; ++i) carry[i] = inCarry[i]; }
  if (false) {
#elif CARRY_NOWAIT
  if (me < G_W) { for (i32 i = 0; i < NW; ++i) carry[i] = 0; }
  if (false) {
#else
  if (me < G_W) {
#endif
#if OLD_FENCE
    if (me == 0) {
      do {
        spin();
#if SPIN_STATS && !ROE
        ++spinCount;
#endif
      } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device));
    }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do {
        spin();
#if SPIN_STATS && !ROE
        ++spinCount;
#endif
      } while(atomic_load_explicit((atomic_uint *) &ready[pos],
                 CARRY_ACQREL ? memory_order_acquire : memory_order_relaxed, memory_scope_device) == 0);
    }
    if (!CARRY_ACQREL) mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H/WMUL we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

#if SPIN_STATS && !ROE
  updateStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, (float) spinCount);
#endif

#if GOLD_PAIR
  u32 goldFwdExponent;
  Z61 goldFwdWeight = goldStartingWeight(
    THREAD_WEIGHTS, lowMe, line, false, &goldFwdExponent);
#endif

  // Apply each 32 or 64 bit carry to the 2 words.  Apply weights.
  for (i32 i = 0; i < NW; ++i) {
    // Generate the second weight shifts
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
    // Generate big-word/little-word flag, propagate final carry
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
    GF31 const next31 = U2(shl(make_Z31(wu[i].x), m31_weight_shift0),
                           shl(make_Z31(wu[i].y), m31_weight_shift1));
#if DETACHED_M31_EDGE == 1
    ((P(GF31)) detach)[(line * NW + i) * G_W + lowMe] = next31;
#else
    u31[i] = next31;
#endif
#if GOLD_PAIR
    u32 goldOddExponent = goldFwdExponent;
    Z61 const goldFwdWeight1 = advanceGoldForward(
      goldFwdWeight, GOLD_DELTA_ONE, (Z61)GOLD_FWD_ONE,
      &goldOddExponent);
    u61[i] = U2(mul(make_Z61_word(wu[i].x), goldFwdWeight),
                  mul(make_Z61_word(wu[i].y), goldFwdWeight1));
#else
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;
    u61[i] = U2(shl(make_Z61(wu[i].x), m61_weight_shift0), shl(make_Z61(wu[i].y), m61_weight_shift1));
#endif

    // Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_bigstep;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
#if GOLD_PAIR
    goldFwdWeight = advanceGoldForward(
      goldFwdWeight, GOLD_DELTA_X, (Z61)GOLD_FWD_X,
      &goldFwdExponent);
#else
    m61_combo_counter += m61_combo_bigstep;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
#endif
  }

#if !DETACHED_M31_EDGE || DETACHED_M31_EDGE == 2
  fft_WIDTH2(lds31, u31, smallTrig31, WMUL, lowMe);
  writeCarryFusedLine(u31, out31, line, lowMe);
#endif

  dependentLaunch();   // Next kernel will be fftMiddleInGF31

  fft_WIDTH2(lds61, u61, smallTrig61, WMUL, lowMe);
  writeCarryFusedLine(u61, out61, line, lowMe);
}


/**************************************************************************/
/*      Fused carry for M31 plus two independent 32-bit Riesel fields     */
/**************************************************************************/

#elif FFT_TYPE == FFT31R2

KERNEL(G_W * WMUL) carryFused(P(T2) out, CP(T2) in, u32 posROE,
                              P(i64) carryShuttle, P(u32) ready,
                              Trig smallTrig, BigTab THREAD_WEIGHTS,
                              P(uint) bufROE) {
  local GF31 lds[WMUL * LDS_BYTES / sizeof(GF31)];
  GF31 u[NW];

  u32 const gr = get_group_id(0);
  u32 const me = get_local_id(0);
  u32 const H = BIG_HEIGHT;
#if WMUL == 1
  u32 const lowMe = me;
  u32 line = gr;
#else
  u32 const lowMe = me % G_W;
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

  CP(GF31) in31 = (CP(GF31))(in + DISTGF31);
  CP(GF31) in0 = (CP(GF31))(in + DISTR0);
  CP(GF31) in1 = (CP(GF31))(in + DISTR1);
  P(GF31) out31 = (P(GF31))(out + DISTGF31);
  P(GF31) out0 = (P(GF31))(out + DISTR0);
  P(GF31) out1 = (P(GF31))(out + DISTR1);
  TrigGF31 trig31 = (TrigGF31)(smallTrig + DISTWTRIGGF31);
  TrigGF31 trig0 = (TrigGF31)(smallTrig + DISTWTR0);
  TrigGF31 trig1 = (TrigGF31)(smallTrig + DISTWTR1);
  P(GF31) zeroScratch31 = (P(GF31))(carryShuttle + ND + WIDTH);
  P(GF31) zeroScratch0 = zeroScratch31 + G_W * WMUL * NW;

#if HAS_ASM
  __asm("s_setprio 3");
#endif

  u32 const zerohack = ZEROHACK_W * (u32)get_group_id(0) / 131072;
  // Spill each completed inverse width transform to its own global plane.
  // This keeps only one eight-value vector live instead of three, avoiding
  // the register explosion of the straightforward fused implementation.
  readCarryFusedLine(in31, u, line, lowMe);
  fft_WIDTH1(lds + zerohack, u, trig31 + zerohack, WMUL, lowMe);
  if (gr == 0) {
    for (u32 i = 0; i < NW; ++i) zeroScratch31[me * NW + i] = u[i];
  } else {
    writeCarryFusedLine(u, out31, line, lowMe);
  }

  readCarryFusedLine(in0, u, line, lowMe);
  rfFftWidth(lds + zerohack, u, trig0 + zerohack, WMUL, lowMe,
             RIESEL0_T8, RIESEL_Q0, RIESEL_Q0_NEG_INV);
  if (gr == 0) {
    for (u32 i = 0; i < NW; ++i) zeroScratch0[me * NW + i] = u[i];
  } else {
    writeCarryFusedLine(u, out0, line, lowMe);
  }

  dependentLaunchWait();
  readCarryFusedLine(in1, u, line, lowMe);
  rfFftWidth(lds + zerohack, u, trig1 + zerohack, WMUL, lowMe,
             RIESEL1_T8, RIESEL_Q1, RIESEL_Q1_NEG_INV);
  if (gr != 0) writeCarryFusedLine(u, out1, line, lowMe);
  mem_fence(CLK_GLOBAL_MEM_FENCE);

  Word2 wu[NW];
  P(CFcarry) carryShuttlePtr = (P(CFcarry))carryShuttle;
  CFcarry carry[NW + 1];
  u32 roundMax = 0;
  float carryMax = 0;

  u32 const word_index = (lowMe * H + line) * 2;
  const u32 m31_log2_root_two = (u32)(((1ULL << 30) / NWORDS) % 31);
  const u32 m31_bigword_weight_shift =
    (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_shift_step = (m31_bigword_weight_shift + 30) % 31;

  union { uint2 a; u64 b; } m31_combo;
#define r2_frac_bits         m31_combo.a[0]
#define r2_m31_shift         m31_combo.a[1]
#define r2_m31_counter       m31_combo.b
  const u64 m31_combo_step = make_u64(m31_shift_step, FRAC_BPW_HI);
  const u64 m31_combo_bigstep =
    (comboFracBits(G_W * H * 2 - 1) +
     make_u64((G_W * H * 2 - 1) * m31_shift_step, 0)) % (31ULL << 32);
  r2_m31_counter = comboFracBits(word_index) +
                   make_u64(word_index * m31_shift_step, 0xFFFFFFFF);
  r2_m31_shift %= 31;
  u64 const m31_starting_counter = r2_m31_counter;

  const u32 log2_NWORDS =
    (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
    (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
    (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
  r2_m31_shift = adjust_m31_weight_shift(r2_m31_shift + log2_NWORDS + 1);

  u32 invExponent0, invExponent1;
  Z31 inv0 = rieselCarryStartingInverse(THREAD_WEIGHTS, 0, lowMe, line,
                                         &invExponent0);
  Z31 inv1 = rieselCarryStartingInverse(THREAD_WEIGHTS, 1, lowMe, line,
                                         &invExponent1);
  inv0 = riesel0Mul(inv0, (Z31)RIESEL0_INV_SCALE);
  inv1 = riesel1Mul(inv1, (Z31)RIESEL1_INV_SCALE);

  for (u32 i = 0; i < NW; ++i) {
    u32 const shift0 = r2_m31_shift;
    r2_m31_counter += m31_combo_step;
    r2_m31_shift = adjust_m31_weight_shift(r2_m31_shift);
    u32 const shift1 = r2_m31_shift;

    u32 oddExponent0 = invExponent0;
    u32 oddExponent1 = invExponent1;
    Z31 const oddInv0 = advanceRieselCarryInverse(
      inv0, 0, RIESEL_DELTA_ONE, (Z31)RIESEL0_INV_ONE, &oddExponent0);
    Z31 const oddInv1 = advanceRieselCarryInverse(
      inv1, 1, RIESEL_DELTA_ONE, (Z31)RIESEL1_INV_ONE, &oddExponent1);

    bool const biglit0 = r2_frac_bits <= FRAC_BPW_HI;
    bool const biglit1 = r2_frac_bits >= -FRAC_BPW_HI;
    GF31 const value31 = gr == 0 ? zeroScratch31[me * NW + i] :
                                  rfReadCarryValue(out31, line, lowMe, i);
    GF31 const value0 = gr == 0 ? zeroScratch0[me * NW + i] :
                                 rfReadCarryValue(out0, line, lowMe, i);
    GF31 const value1 = gr == 0 ? u[i] :
                                 rfReadCarryValue(out1, line, lowMe, i);
    wu[i] = weightAndCarryPairSloppy(
      SWAP_XY(value31), SWAP_XY(value0), SWAP_XY(value1),
      shift0, shift1, inv0, oddInv0, inv1, oddInv1,
      LL != 0, (LL & (i == 0) & (line == 0) & (me == 0)) ? -2 : 0,
      biglit0, biglit1, &carry[i], &roundMax, &carryMax);

    r2_m31_counter += m31_combo_bigstep;
    r2_m31_shift = adjust_m31_weight_shift(r2_m31_shift);
    inv0 = advanceRieselCarryInverse(inv0, 0, RIESEL_DELTA_X,
                                     (Z31)RIESEL0_INV_X, &invExponent0);
    inv1 = advanceRieselCarryInverse(inv1, 1, RIESEL_DELTA_X,
                                     (Z31)RIESEL1_INV_X, &invExponent1);
  }
  r2_m31_counter = m31_starting_counter;

#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL - 1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) {
      CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)],
              carry[i]);
    }
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
#if OLD_FENCE
    bar(G_W);
    if (lowMe == 0) atomic_store((atomic_uint *)&ready[gr], 1);
#else
    if (lowMe % WAVEFRONT == 0) {
      u32 const pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *)&ready[pos], 1);
    }
#endif
  }

  if (gr == 0) return;

#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
  updateStats((local u32 *)lds, G_W * WMUL, H / WMUL, bufROE, posROE, 0.0f);
#elif SPIN_STATS
  // Deferred until after the readiness wait.
#elif STATS & (1 << MUL3)
  updateStats((local u32 *)lds, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  shufl_carries_up(lds, carry, me, lowMe);

#if SPIN_STATS && !ROE
  u32 spinCount = 0;
#endif
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) {
      do {
        spin();
#if SPIN_STATS && !ROE
        ++spinCount;
#endif
      } while (!atomic_load_explicit((atomic_uint *)&ready[gr - 1],
                                     memory_order_relaxed, memory_scope_device));
    }
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 const pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do {
        spin();
#if SPIN_STATS && !ROE
        ++spinCount;
#endif
      } while (atomic_load_explicit((atomic_uint *)&ready[pos],
                                    memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (me % WAVEFRONT == 0) ready[pos] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[
          (gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {
#if !OLD_FENCE
      bar();
#endif
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[
          (gr - 1) * WIDTH +
          CarryShuttleAccess((me + G_W - 1) % G_W, i)]);
      }
      if (me == 0) {
        carry[NW] = carry[NW - 1];
        for (i32 i = NW - 1; i; --i) carry[i] = carry[i - 1];
        carry[0] = carry[NW];
      }
    }
  }

#if SPIN_STATS && !ROE
  updateStats((local u32 *)lds, G_W * WMUL, H / WMUL,
              bufROE, posROE, (float)spinCount);
#endif

  // Finalize carried words and produce the M31 forward width transform.
  for (i32 i = 0; i < NW; ++i) {
    u32 const shift0 = r2_m31_shift;
    r2_m31_counter += m31_combo_step;
    r2_m31_shift = adjust_m31_weight_shift(r2_m31_shift);
    u32 const shift1 = r2_m31_shift;

    wu[i] = carryFinal(wu[i], carry[i], r2_frac_bits <= FRAC_BPW_HI);
    u[i] = U2(shl(make_Z31(wu[i].x), shift0),
              shl(make_Z31(wu[i].y), shift1));
    r2_m31_counter += m31_combo_bigstep;
    r2_m31_shift = adjust_m31_weight_shift(r2_m31_shift);
  }
  fft_WIDTH2(lds, u, trig31, WMUL, lowMe);
  writeCarryFusedLine(u, out31, line, lowMe);

  // Produce q0.  Its weight state dies before q1 is constructed, further
  // shortening live ranges in this already register-heavy kernel.
  u32 fwdExponent0;
  Z31 fwd0 = rieselCarryStartingForward(THREAD_WEIGHTS, 0, lowMe, line,
                                         &fwdExponent0);
  for (i32 i = 0; i < NW; ++i) {
    u32 oddExponent0 = fwdExponent0;
    Z31 const oddFwd0 = advanceRieselCarryForward(
      fwd0, 0, RIESEL_DELTA_ONE, (Z31)RIESEL0_FWD_ONE, &oddExponent0);
    u[i] = U2(riesel0Mul(rfMakeWord(wu[i].x, RIESEL_Q0,
                                    RIESEL_Q0_NEG_INV, RIESEL0_R2), fwd0),
              riesel0Mul(rfMakeWord(wu[i].y, RIESEL_Q0,
                                    RIESEL_Q0_NEG_INV, RIESEL0_R2), oddFwd0));
    fwd0 = advanceRieselCarryForward(fwd0, 0, RIESEL_DELTA_X,
                                     (Z31)RIESEL0_FWD_X, &fwdExponent0);
  }
  rfFftWidth(lds, u, trig0, WMUL, lowMe,
             RIESEL0_T8, RIESEL_Q0, RIESEL_Q0_NEG_INV);
  writeCarryFusedLine(u, out0, line, lowMe);

  // Produce q1 using the same register vector and LDS allocation.
  u32 fwdExponent1;
  Z31 fwd1 = rieselCarryStartingForward(THREAD_WEIGHTS, 1, lowMe, line,
                                         &fwdExponent1);
  for (i32 i = 0; i < NW; ++i) {
    u32 oddExponent1 = fwdExponent1;
    Z31 const oddFwd1 = advanceRieselCarryForward(
      fwd1, 1, RIESEL_DELTA_ONE, (Z31)RIESEL1_FWD_ONE, &oddExponent1);
    u[i] = U2(riesel1Mul(rfMakeWord(wu[i].x, RIESEL_Q1,
                                    RIESEL_Q1_NEG_INV, RIESEL1_R2), fwd1),
              riesel1Mul(rfMakeWord(wu[i].y, RIESEL_Q1,
                                    RIESEL_Q1_NEG_INV, RIESEL1_R2), oddFwd1));
    fwd1 = advanceRieselCarryForward(fwd1, 1, RIESEL_DELTA_X,
                                     (Z31)RIESEL1_FWD_X, &fwdExponent1);
  }
  rfFftWidth(lds, u, trig1, WMUL, lowMe,
             RIESEL1_T8, RIESEL_Q1, RIESEL_Q1_NEG_INV);
  writeCarryFusedLine(u, out1, line, lowMe);
  dependentLaunch();

#undef r2_frac_bits
#undef r2_m31_shift
#undef r2_m31_counter
}


/******************************************************************************/
/*  Similar to above, but for a hybrid FFT based on FP32*GF(M31^2)*GF(M61^2)  */
/******************************************************************************/

#elif FFT_TYPE == FFT323161

// The "carryFused" is equivalent to the sequence: fftW, carryA, carryB, fftPremul.
// It uses "stairway forwarding" (forwarding carry data from one workgroup to the next)
KERNEL(G_W * WMUL) carryFused(P(T2) out, CP(T2) in, u32 posROE, P(i64) carryShuttle, P(u32) ready, Trig smallTrig,
                              ConstBigTabFP32 CONST_THREAD_WEIGHTS, BigTabFP32 THREAD_WEIGHTS, P(uint) bufROE
#if PARITY_SQUARE
                              , CP(u32) parityIn, P(u32) parityOut
#endif
                              ) {
  local GF61 lds61[WMUL * LDS_BYTES / sizeof(GF61)];
  local F2 *ldsF2 = (local F2 *) lds61;
  local GF31 *lds31 = (local GF31 *) lds61;

  F2 uF2[NW];
  GF31 u31[NW];
  GF61 u61[NW];

  u32 gr = get_group_id(0);
  u32 me = get_local_id(0);

  u32 H = BIG_HEIGHT;
#if WMUL == 1
  u32 lowMe = me;
  u32 line = gr;
#else
  u32 lowMe = me % G_W;           // lane-id in one of the WMUL sub-workgroups.
  u32 line = gr * WMUL + me / G_W;
#endif
  if (line >= H) line -= H;

  CP(F2) inF2 = (CP(F2)) in;
  P(F2) outF2 = (P(F2)) out;
  TrigFP32 smallTrigF2 = (TrigFP32) smallTrig;
  CP(GF31) in31 = (CP(GF31)) (in + DISTGF31);
  P(GF31) out31 = (P(GF31)) (out + DISTGF31);
  TrigGF31 smallTrig31 = (TrigGF31) (smallTrig + DISTWTRIGGF31);
  CP(GF61) in61 = (CP(GF61)) (in + DISTGF61);
  P(GF61) out61 = (P(GF61)) (out + DISTGF61);
  TrigGF61 smallTrig61 = (TrigGF61) (smallTrig + DISTWTRIGGF61);

#if HAS_ASM
  __asm("s_setprio 3");
#endif

// Try this weird FFT_width call that adds a "hidden zero" when unrolling.  This prevents the compiler from finding
// common sub-expressions to re-use in the second fft_WIDTH call.  Re-using this data requires dozens of VGPRs
// which causes a terrible reduction in occupancy.
  u32 zerohack = ZEROHACK_W * (u32) get_group_id(0) / 131072;

  readCarryFusedLine(inF2, uF2, line, lowMe);
  fft_WIDTH1(ldsF2 + zerohack, uF2, smallTrigF2 + zerohack, WMUL, lowMe);

#if !DETACHED_M31_EDGE
  readCarryFusedLine(in31, u31, line, lowMe);
  fft_WIDTH1(lds31 + zerohack, u31, smallTrig31 + zerohack, WMUL, lowMe);
#endif

  dependentLaunchWait();   // Previous kernel was fftMiddleOutGF61

  readCarryFusedLine(in61, u61, line, lowMe);
  fft_WIDTH1(lds61 + zerohack, u61, smallTrig61 + zerohack, WMUL, lowMe);

#if DETACHED_M31_EDGE
  // Timing gate for a decoupled exact side transform.  The completed design
  // supplies this line from an early standalone inverse-width kernel.  Until
  // that scheduler is attached, the values are deliberately not a correct
  // transform and this opt-in mode is performance evidence only.
  readCarryFusedLine(in31, u31, line, lowMe);
#endif

  Word2 wu[NW];
  u32 me_frac_bits = fracBits(lowMe * H * 2);
#if !NVIDIAGPU || CUDA_BACKEND
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), TSLOAD(&THREAD_WEIGHTS[G_W + line]));
  u32 line_frac_bits = fracBits(line * 2);
  u32 base_frac_bits = me_frac_bits + line_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > line_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > line_frac_bits);
#else
  F2 weights = fancyMul(TFLOAD(&THREAD_WEIGHTS[lowMe]), CONST_THREAD_WEIGHTS[line % 64]);
  u32 partialLine_frac_bits = fracBits((line % 64) * 2);
  u32 base_frac_bits = me_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
  weights = fancyMul(weights, CONST_THREAD_WEIGHTS[64 + line / 64]);
  partialLine_frac_bits = fracBits(((line / 64) * 64) * 2);
  base_frac_bits = base_frac_bits + partialLine_frac_bits;
  weights.x = optionalDouble(weights.x, base_frac_bits > partialLine_frac_bits);
  weights.y = optionalHalve(weights.y, base_frac_bits > partialLine_frac_bits);
#endif

  P(i64) carryShuttlePtr = (P(i64)) carryShuttle;
  i64 carry[NW+1];

  float roundMax = 0;
  float carryMax = 0;
#if ROE && ROE_COUNT
  u32 riskyPairs = 0;
#endif

  u32 word_index = (lowMe * H + line) * 2;

  // Weight is 2^[ceil(qj / n) - qj/n] where j is the word index, q is the Mersenne exponent, and n is the number of words.
  // Let s be the shift amount for word 1.  The shift amount for word x is ceil(x * (s - 1) + num_big_words_less_than_x) % 31.
#if GOOD_THOMAS3
  const u32 m31_log2_root_two = 21;
#else
  const u32 m31_log2_root_two = (u32)(((1ULL << 30) / NWORDS) % 31);
#endif
  const u32 m31_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m31_log2_root_two % 31;
  const u32 m31_bigword_weight_shift_minus1 = (m31_bigword_weight_shift + 30) % 31;
#if GOOD_THOMAS3
  const u32 m61_log2_root_two = 45;
#else
  const u32 m61_log2_root_two = (u32)(((1ULL << 60) / NWORDS) % 61);
#endif
  const u32 m61_bigword_weight_shift = (NWORDS - EXP % NWORDS) * m61_log2_root_two % 61;
  const u32 m61_bigword_weight_shift_minus1 = (m61_bigword_weight_shift + 60) % 61;

  // Derive the big vs. little flags from the fractional number of bits in each word.
  // Create a 64-bit counter that tracks both weight shifts and frac_bits (adding 0xFFFFFFFF to effect the ceil operation required for weight shift).
  union { uint2 a; u64 b; } m31_combo, m61_combo;
#define frac_bits           m31_combo.a[0]
#define m31_weight_shift    m31_combo.a[1]
#define m31_combo_counter   m31_combo.b
#define m61_weight_shift    m61_combo.a[1]
#define m61_combo_counter   m61_combo.b

  const u64 m31_combo_step = make_u64(m31_bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 m31_combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * m31_bigword_weight_shift_minus1, 0)) % (31ULL << 32);
  m31_combo_counter = comboFracBits(word_index) + make_u64(word_index * m31_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m31_weight_shift = m31_weight_shift % 31;
  u64 m31_starting_combo_counter = m31_combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation
  const u64 m61_combo_step = make_u64(m61_bigword_weight_shift_minus1, FRAC_BPW_HI);
  const u64 m61_combo_bigstep = (comboFracBits(G_W * H * 2 - 1) + make_u64((G_W * H * 2 - 1) * m61_bigword_weight_shift_minus1, 0)) % (61ULL << 32);
  m61_combo_counter = comboFracBits(word_index) + make_u64(word_index * m61_bigword_weight_shift_minus1, 0xFFFFFFFF);
  m61_weight_shift = m61_weight_shift % 61;
  u64 m61_starting_combo_counter = m61_combo_counter;     // Save starting counter before adding log2_NWORDS+1 for applying weights after carry propagation

  // We also adjust shift amount for the fact that NTT returns results multiplied by 2*NWORDS.
#if GOOD_THOMAS3
  // Each Good--Thomas exact channel transforms 2^20 scalar words.  The
  // common +1 below accounts for the packed DGT's additional factor of two.
  const u32 log2_NWORDS = 20;
#else
  const u32 log2_NWORDS = (WIDTH == 256 ? 8 : WIDTH == 512 ? 9 : WIDTH == 1024 ? 10 : 12) +
                          (MIDDLE == 1 ? 0 : MIDDLE == 2 ? 1 : MIDDLE == 4 ? 2 : MIDDLE == 8 ? 3 : 4) +
                          (SMALL_HEIGHT == 256 ? 8 : SMALL_HEIGHT == 512 ? 9 : SMALL_HEIGHT == 1024 ? 10 : 12) + 1;
#endif
  m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift + log2_NWORDS + 1);
  m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift + log2_NWORDS + 1);

  // Apply the inverse weights and carry propagate pairs to generate the output carries

  F invBase = weights.x;
  for (u32 i = 0; i < NW; ++i) {
    // Generate the FP32 weights and second GF31 and GF61 weight shift
    F invWeight1 = i == 0 ? invBase : optionalDouble(fancyMul(invBase, iweightStep(i)), frac_bits > base_frac_bits);
    F invWeight2 = optionalDouble(fancyMul(invWeight1, IWEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;

    // Generate big-word/little-word flags
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    bool biglit1 = frac_bits >= -FRAC_BPW_HI;   // Same as frac_bits + FRAC_BPW_HI <= FRAC_BPW_HI;

    // Apply the inverse weights, optionally compute roundoff error, and convert to integer.  Also apply MUL3 here.
    // Then propagate carries through two words (the first carry does not have to be accurately calculated because it will
    // be accurately calculated by carryFinal later on).  The second carry must be accurate for output to the carry shuttle.
#if ROE && ROE_COUNT
    float pairRoundMax = 0;
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(uF2[i]), SWAP_XY(u31[i]), SWAP_XY(u61[i]), invWeight1, invWeight2, m31_weight_shift0, m31_weight_shift1, m61_weight_shift0, m61_weight_shift1,
#if PARITY_SQUARE
                      squareCoefficientParity(parityIn, (lowMe + i * G_W) * H + line), 0,
#endif
                      LL != 0, (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &pairRoundMax, &carryMax);
    roundMax = max(roundMax, pairRoundMax);
    riskyPairs += pairRoundMax >= (float) ROE_COUNT * 0.001f;
#else
    wu[i] = weightAndCarryPairSloppy(SWAP_XY(uF2[i]), SWAP_XY(u31[i]), SWAP_XY(u61[i]), invWeight1, invWeight2, m31_weight_shift0, m31_weight_shift1, m61_weight_shift0, m61_weight_shift1,
#if PARITY_SQUARE
                      squareCoefficientParity(parityIn, (lowMe + i * G_W) * H + line), 0,
#endif
                      // For an LL test, add -2 as the very initial "carry in"
                      // We'd normally use logical &&, but the compiler whines with warning and bitwise fixes it
                      LL != 0, (LL & (i == 0) & (line==0) & (me == 0)) ? -2 : 0, biglit0, biglit1, &carry[i], &roundMax, &carryMax);
#endif

    // Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_bigstep;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    m61_combo_counter += m61_combo_bigstep;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
  }
  m31_combo_counter = m31_starting_combo_counter;     // Restore starting counter for applying weights after carry propagation
  m61_combo_counter = m61_starting_combo_counter;

  // Write out our carries for the last line in this group. Only groups 0 to H/WMUL-1 need to write carries out.
  // Group H/WMUL is a duplicate of group 0 (producing the same results) so we don't care about that group writing out,
  // but it's fine either way.
  // AMD's OpenCL Windows compiler generates warnings about always true if statements for WMUL-1.  So instead an #if is required
#if WMUL == 1
  if (gr < H) {
#else
  if (gr < H / WMUL && me >= (WMUL-1) * G_W) {
#endif
    for (i32 i = 0; i < NW; ++i) { CSSTORE(&carryShuttlePtr[gr * WIDTH + CarryShuttleAccess(lowMe, i)], carry[i]); }

    // Tell next group that its carries are ready
#if OLD_FENCE
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    bar(G_W);
    if (lowMe == 0) { atomic_store((atomic_uint *) &ready[gr], 1); }
#else
    write_mem_fence(CLK_GLOBAL_MEM_FENCE);
    if (lowMe % WAVEFRONT == 0) { 
      u32 pos = gr * (G_W / WAVEFRONT) + lowMe / WAVEFRONT;
      atomic_store((atomic_uint *) &ready[pos], 1);
    }
#endif
  }

  // Group zero will be redone when gr == H / WMUL
  if (gr == 0) { return; }

  // Do some work while our carries may not be ready
#if HAS_ASM
  __asm("s_setprio 0");
#endif

#if ROE
#if ROE_COUNT
  updateCountStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, riskyPairs);
#else
  updateStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, roundMax);
#endif
#elif STATS & (1 << MUL3)
  updateStats((local u32 *) lds61, G_W * WMUL, H / WMUL, bufROE, posROE, carryMax);
#endif

  // Shuffle carries up
  shufl_carries_up(lds61, carry, me, lowMe);

  // Wait until our carries are ready
  if (me < G_W) {
#if OLD_FENCE
    if (me == 0) { do { spin(); } while(!atomic_load_explicit((atomic_uint *) &ready[gr - 1], memory_order_relaxed, memory_scope_device)); }
    // work_group_barrier(CLK_GLOBAL_MEM_FENCE, memory_scope_device);
    bar();
    read_mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me == 0) ready[gr - 1] = 0;
#else
    u32 pos = (gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT;
    if (me % WAVEFRONT == 0) {
      do { spin(); } while(atomic_load_explicit((atomic_uint *) &ready[pos], memory_order_relaxed, memory_scope_device) == 0);
    }
    mem_fence(CLK_GLOBAL_MEM_FENCE);
    // Clear carry ready flag for next iteration
    if (me % WAVEFRONT == 0) ready[(gr - 1) * (G_W / WAVEFRONT) + me / WAVEFRONT] = 0;
#endif
#if HAS_ASM
    __asm("s_setprio 1");
#endif

    // Read from the carryShuttle carries produced by the previous WIDTH group.  Rotate carries from the last WIDTH line.
    // The new carry layout lets the AMD compiler generate global_load_dwordx4 instructions.
    if (gr < H / WMUL) {
      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess(me, i)]);
      }
    } else {

#if !OLD_FENCE
      // For gr==H/WMUL we need the barrier since the carry reading is shifted, thus the per-wavefront trick does not apply.
      bar();
#endif

      for (i32 i = 0; i < NW; ++i) {
        carry[i] = CSLOAD(&carryShuttlePtr[(gr - 1) * WIDTH + CarryShuttleAccess((me + G_W - 1) % G_W, i) /* ((me!=0) + NW - 1 + i) % NW*/]);
      }

      if (me == 0) {
        carry[NW] = carry[NW-1];
        for (i32 i = NW-1; i; --i) { carry[i] = carry[i-1]; }
        carry[0] = carry[NW];
      }
    }
  }

  // Apply each 32 or 64 bit carry to the 2 words.  Apply weights.
  F base = weights.y;
  for (i32 i = 0; i < NW; ++i) {
    // Calculate inverse weights
    F weight1 = i == 0 ? base : optionalHalve(fancyMul(base, fweightStep(i)), frac_bits > base_frac_bits);
    F weight2 = optionalHalve(fancyMul(weight1, WEIGHT_STEP), frac_bits + FRAC_BPW_HI > FRAC_BPW_HI);
    // Generate the second weight shifts
    u32 m31_weight_shift0 = m31_weight_shift;
    m31_combo_counter += m31_combo_step;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    u32 m31_weight_shift1 = m31_weight_shift;
    u32 m61_weight_shift0 = m61_weight_shift;
    m61_combo_counter += m61_combo_step;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
    u32 m61_weight_shift1 = m61_weight_shift;
    // Generate big-word/little-word flag, propagate final carry
    bool biglit0 = frac_bits <= FRAC_BPW_HI;
    wu[i] = carryFinal(wu[i], carry[i], biglit0);
#if PARITY_SQUARE
    u32 logicalPair = (lowMe + i * G_W) * H + line;
    parityOut[logicalPair] = ((u32) wu[i].x & 1u) | (((u32) wu[i].y & 1u) << 1);
#endif
    uF2[i] = U2(weight1 * wu[i].x, weight2 * wu[i].y);
    u31[i] = U2(shl(make_Z31(wu[i].x), m31_weight_shift0), shl(make_Z31(wu[i].y), m31_weight_shift1));
    u61[i] = U2(shl(make_Z61(wu[i].x), m61_weight_shift0), shl(make_Z61(wu[i].y), m61_weight_shift1));

    // Generate weight shifts and frac_bits for next pair
    m31_combo_counter += m31_combo_bigstep;
    m31_weight_shift = adjust_m31_weight_shift(m31_weight_shift);
    m61_combo_counter += m61_combo_bigstep;
    m61_weight_shift = adjust_m61_weight_shift(m61_weight_shift);
  }

  fft_WIDTH2(ldsF2, uF2, smallTrigF2, WMUL, lowMe);
  writeCarryFusedLine(uF2, outF2, line, lowMe);

  dependentLaunch();   // Next kernel will be fftMiddleInFP32

#if !DETACHED_M31_EDGE
  fft_WIDTH2(lds31, u31, smallTrig31, WMUL, lowMe);
#endif
  writeCarryFusedLine(u31, out31, line, lowMe);

  fft_WIDTH2(lds61, u61, smallTrig61, WMUL, lowMe);
  writeCarryFusedLine(u61, out61, line, lowMe);
}


#else
error - missing CarryFused kernel implementation
#endif
