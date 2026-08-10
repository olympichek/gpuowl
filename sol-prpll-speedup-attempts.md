# Architectural speedup strategies for PRPLL

This document records possible architectural changes for increasing the aggregate
throughput of PRPLL on the NVIDIA RTX PRO 6000 Blackwell Max-Q. Single-exponent
latency is useful but is not the primary objective: unrelated exponents may be
tested concurrently whenever that improves total iterations per second or
exponents completed per day. The document intentionally concentrates on
representation, algorithms, execution pipelines, and hardware facilities rather
than launch-parameter or instruction-level micro-optimizations.

The initial reference target is exponent `136279841`, currently using a 4M
`GF(M31^2) * GF(M61^2)` NTT. The requested 1M-iteration run averaged 205.5
microseconds per iteration; shorter repeat runs are normally about 200--206
microseconds. Any estimated improvement below is a research target, not a
measured result.

For throughput experiments, this exponent can be combined with unrelated
exponents having the same or similar transform geometry. Results should be
reported both as per-exponent latency and as aggregate useful iterations per
second.

## Result summary

No experimental architecture has yet beaten the tuned production `M31+M61`
path **end to end** for exponent `136279841`. The first campaign's strongest
standalone direct-prime result was a correct 4M x three-prime Montgomery NTT at
0.708 ms before CRT and carry, versus about 0.204 ms for the complete production
iteration. Running multiple unchanged exponent streams reduced aggregate
throughput by 13--35%, and the batch-shaped exact Tensor Core tile was 42% slower
than sparse SIMT.

A second campaign found one materially different component-level lead: retain
PRPLL's quadratic/conjugate packing and replace the 16-byte `GF(M61^2)` plane
with two 8-byte quadratic Riesel-prime planes.  The first 92.992-bit prime pair
later failed the norm-one transform-root gate.  The valid integrated pair keeps
the residue state at 48 MiB and provides 92.957 bits of exact CRT range.  Its
standalone radix-16 tile was 15.7--15.8% faster than the current M31+M61 tile,
but the complete three-plane PRPLL path measured 302.5 us/iteration versus
205.5 us/iteration for production.  The component lead therefore did not
become an end-to-end speedup.

The later shorter-transform and lazy-field experiments also did not beat
production.  A mathematically valid 2M four-field
`M31+M61+q0+q1` schedule reached 225.0 us/iteration while still using an
optimistically cheap placeholder carry, already 9.5% behind production.  Two
sub-`2^30` Riesel fields then enabled redundant `[0,2q)` Montgomery arithmetic:
this reduced the generic-q production kernels by 8--17% in static instruction
count and passed a full 1M-iteration PRPLL run, but improved the old-prime
M31R2 path by only about 0.6% in paired controls.  Its sustained result was
285.4 us/iteration with the correct `52b03a7cc55e677d` residue.  Thus the
campaign still has no end-to-end speedup over the pre-existing 205.5 us
M31+M61 code.

A production-sized middle/width/carry fusion proxy subsequently appeared to
offer a 19.5% component gain with a compact 33-bit bridge.  Integrating the
actual M31/M61 CRT and carry made the prototype exact, but reversed the result:
the compact fused path measured 311.0 us/iteration against a nearby 230.3 us
production control.  The proxy had omitted the arithmetic and live registers
that determine occupancy in the real carry.  This is now another measured
no-go, not a surviving speedup lead.

Two later clean-room gates also failed.  Keeping a complete 64-KiB M61 tile
resident through middle/height/square was exact and spill-free, but ran 2--3%
slower than three cache-resident kernels.  A novel 3M M61 times M31-prime-power
ring design first passed its scalar arithmetic gate; after lifting real roots
and adding both required fields, its exact radix-16 tile was 23.5--24.1% slower
than the architecture-sized production M31+M61 tile.  It was rejected before
radix 3, wider CRT, or carry.  A standalone program remains fully acceptable,
but it is held to the same complete 180-us/1M-iteration correctness gate.

The first forum-derived 3M `M31+M61+(11*2^26-1)` implementation produced useful
component timings, but a subsequent whole-algorithm gate invalidated that
specific prime: 2 is not a cube in its quadratic field, so the required 3M
Crandall--Fagin weight does not exist.  Its exact tile, radix-3, lazy-arithmetic,
and carry measurements remain lower-bound data, not a viable architecture.
The compatible replacement `q=1021*2^21-1=2141192191` is now being re-gated;
it admits both the radix-3 root and an exact 3M weight.  No result from either
version is yet an end-to-end speedup or an established route below 180 us.
Its tile, radix-3, exact 121-bit carry, and Good--Thomas carry-layout gates now
pass.  In particular, a shared-memory three-channel transpose limits the
otherwise severe CRT-index gather/scatter penalty to about 4.5--5.4 us.  The
production-shaped transform core then measured 144.0 us on the same timeline
method that gives 129.6 us for production.  Core plus the 38.3--38.6-us exact
shared carry already exceeds 182 us before either width edge, so the current
3M integration fails the 180-us lower-bound gate and is no longer a surviving
end-to-end route.

A subsequent clean-room single-field `GF((2^127-1)^2)` proposal retained the
same 48-MiB 3M state and eliminated CRT, but failed its first exact arithmetic
gate.  At architecture populations it was 1.81x M61 for one quadratic multiply
and 3.52x for an arithmetic-dense chain; even the theoretical benefit of
two-limb Karatsuba cannot close that gap.  It was rejected before roots or a
full transform.

The 3M FP32+M31+M61 hybrid has now also been rejected.  It uses one FP32
estimate plane plus the production M31 and M61 exact planes, decomposed as
three independent 1M Good--Thomas channels for the NTT fields.  The initially
attractive `1024x3x512` profile reported an 88.480-us transform core, but its
one-channel exact trig tables were accidentally zero.  After implementing the
exact natural/channel mapping, repairing those tables, and correcting the
packed-DGT scale, the main p136 square recurrence matches production at
iteration 2,000 (`05d6515c416b83e2`).  On representative nonzero residues the
same core is actually 136.160 us median, and the direct split implementation
is about 265.3 us/iteration.  Even an optimistic coalesced three-kernel edge
has a measured-component lower bound near 104.6 us.  Thus the apparent 40-us
lead was a zero-data artifact, not a surviving speedup route; a fused edge and
Gerbicz `tailMul` integration were deliberately not implemented after the
corrected lower-bound failure.

A newer 3M all-integer experiment superseded the generic-q arithmetic forecast:
`M31*M61*M19`, where `M19=2^19-1`.  M19 is prime and its quadratic field has
exactly the norm-one 2-adicity required by each 1M-word Good--Thomas channel;
it also has a base-field cube root.  A Beatty-distribution coefficient proof
gives 0.854 bit of p136 headroom.  Its first exact architecture-sized
radix-16/square/inverse tile is 10.1--10.7% faster than the 4M production-sized
M31+M61 control, and its radix-3 and exact CRT/carry component gates pass.
However, the production-shaped transform core improves by only 2.048 us over
the paired generic-q control: 137.728 us rather than 139.776 us.  Core plus the
38.400-us exact shared-layout carry is already 176.128 us before either width
edge, leaving less than 3.9 us for two mandatory transforms.  The route is
therefore rejected at the complete-budget gate; it is not an end-to-end
speedup.

The standalone permission also prompted the first campaign timing of the
Lucas--Lehmer recurrence.  It does not help: a CUDA timeline gives an
approximately 206-us steady cycle, with the same roughly 129.6-us M31/M61
transform core and a roughly 74.1-us dedicated LL fused edge.  The recurrence
removes PRP-specific checking/proof work but does not remove an NTT, and its
subtract-two carry is slightly slower than the approximately 70.5-us PRP fused
edge.  LL is therefore another measured no-go for the 180-us gate, not an
unexplored clean-room shortcut.

The latest standalone screening closed two more apparent Mersenne-field leads.
`M29=2^29-1` is composite, so it failed the field/root gate before timing.  A
genuinely distinct single-ring `M31^3` design has the same 48-MiB state and one
more nominal modulus bit than production, but its exact arithmetic-dense GPU
gate is 2.24--2.29x the same-population M31/M61 comparator.  It is rejected
before root lifting or a transform tile.  A subsequent production-stream-
priority experiment was also negative:
equal-priority controls measured 198.5 and 199.3 us/iteration, while giving the
critical M61 stream high CUDA priority measured 203.1 and 203.3 us/iteration.
The incumbent equal-priority overlap remains best, and the campaign still has
no speedup over the pre-existing code or exact result at or below 180 us.

An exact FP32-offloaded M31 representation was then tested as a resource-
complementarity architecture.  It preserves the production M31*M61 range but
stores each live M31 scalar in `11+11+9`-bit FP32 limbs, making all six
Karatsuba limb products exactly representable.  The architecture-population
quadratic gate is exact but decisive: it is 3.11--3.14x the integer path for one
multiply and 5.13--5.15x for an eight-multiply chain.  With no spills, this is
arithmetic expansion rather than an occupancy accident.  It cannot fit under
the M61 critical chain and is rejected before a full transform.

A narrower exact-FP32 replacement was also re-gated rather than inferred from
the limb result.  `qC=7*2^21-1` times M61 has just enough empirical range at
p150: an exact coefficient oracle found no value outside its balanced CRT range
over 100,000 iterations.  Its FP32 arithmetic overlaps M61 slightly better, but
an optimistic production-population proxy saves only 4.0--4.3 us at an
eight-multiply density.  That is one sixth of the required 25.5-us gain before
charging generic roots, weights, and a less favorable CRT.  The q24*M61
replacement is therefore rejected before integration.

The work below supplies tested host and CUDA arithmetic, complete transform
baselines, a multi-exponent scheduler, an exact Tensor Core prototype, an exact
FP32-estimated CRT prototype, and a proved carry-transfer reference. Each rejected
branch has a measured or mathematical go/no-go reason. Known-slower experimental
code is kept as a benchmark. The packed Riesel branch remains behind explicit
integration gates and has not replaced the production arithmetic.

## Retrospective: why the initial speedup expectation was wrong

The central forecasting error was to overweight the theoretical advantage of
lower-width INT32 arithmetic while underweighting the mature, fused dataflow of
the existing M31+M61 implementation. The initial 20--40% research target was not
supported by a complete-transform measurement.

Specific mistakes were:

- **Extrapolating from arithmetic-dense microbenchmarks.** Three q31 channels are
  30--35% faster than the current arithmetic proxy in long multiplication chains,
  but a complete NTT also pays for butterflies, twiddle access, transposes,
  synchronization, reductions, and stores. The multiplication advantage was not
  the dominant end-to-end cost.
- **Reducing the comparison to operand width and plane bytes.** The existing path
  benefits from cheap Mersenne reduction, quadratic-field packing, specialized
  width/middle/height geometry, cache-aware scheduling, and fusion of transform,
  CRT, carry, weighting, and premultiplication. Replacing its representation also
  discards those advantages.
- **Using an incomplete memory model.** The 48 MiB residue state fits in L2, but
  state plus root tables and multiple active planes create substantial L2 traffic
  and contention. Naive stage fusion reduced nominal global passes while harming
  coalescing enough to lose performance.
- **Using `prime_count * N * log2(N)` as too much of the cost model.** It omits
  modulus-specific reduction cost, twiddle storage, CRT width, carry width,
  packing, register pressure, occupancy, and synchronization.
- **Assuming narrower integers imply proportional GPU throughput.** Exact q31
  multiplication still needs multiply-high, reduction, and correction work, and
  three residues require three complete transform streams. There is no native
  modular INT32 instruction that turns peak INT32 throughput directly into NTT
  throughput.
- **Misreading coarse utilization telemetry.** Low reported DRAM utilization did
  not imply that data movement was unimportant; much of the relevant traffic was
  through L2, shared memory, registers, and transpose paths. At the same time, the
  existing kernel already reached the 300 W board-power limit.
- **Attributing the fused-region time to carry.** The 71.5 us `carryFused` region
  also contains a transform edge, CRT, reweighting, and premultiplication.
  Instrumentation found only 18.4 short readiness polls on average, showing that
  the stairway carry dependency itself was not a large removable bottleneck.
- **Expecting batching to consume idle hardware.** A single production stream was
  already power-limited. Additional unrelated streams divided clock, cache, and
  memory resources, reducing aggregate throughput by 13--35% instead of filling
  unused capacity.
- **Comparing Tensor Cores by peak arithmetic rate.** NTT butterflies are sparse,
  whereas Tensor Cores reward dense matrix work. Exact limb packing and sixteen
  dense MMA contributions made the Tensor tile 42% slower than the sparse SIMT
  butterfly despite beating a deliberately dense SIMT matrix.
- **Setting the whole-operation gate too late.** A replacement transform needed
  to be roughly 0.13 ms or faster to leave useful room for CRT and carry inside a
  0.204--0.206 ms iteration. The best direct-prime transform was already 0.708 ms
  before either, making further integration unable to recover the gap.

Lessons for future campaigns:

1. Establish a hard end-to-end time budget for every proposed subsystem before
   optimizing its arithmetic primitives.
2. Build the smallest *complete* forward/square/inverse pipeline early; use
   microbenchmarks to explain its result, not to predict it in isolation.
3. Model the incumbent's representation, fusion boundaries, twiddle working set,
   cache traffic, carry range, and power behavior—not only operation counts and
   persistent-state bytes.
4. Report aggregate throughput and energy under sustained thermal conditions;
   distinguish cold short-run results from long-run performance.
5. Apply early abandonment gates. If a correct standalone lower bound already
   exceeds the complete production iteration by several times, preserve it as
   research infrastructure but do not spend effort integrating it.
6. Treat hardware peak rates as useful only after mapping the algorithm's actual
   sparsity, precision, packing, reduction, and synchronization costs to that
   hardware.
7. Validate every proposed modulus as a field and construct the maximum required
   root before treating arithmetic timings as evidence. The early `2^23-1`
   control was mistakenly called M23; it is `47 * 178481`, not a prime.
8. Search for representation-preserving substitutions before replacing the whole
   transform. The useful Riesel lead keeps conjugate packing and total bytes;
   the failed direct-prime design discarded both and therefore answered a much
   less favorable architectural question.
9. Reducing static normalization instructions does not imply proportional NTT
   speed.  The sub-`2^30` experiment removed up to one third of a q kernel's
   `VIMNMX` operations without changing its `IMAD` count or memory traffic; an
   arithmetic-dense tile improved several percent, but complete PRPLL improved
   only about 0.6%.

The main durable lesson is that PRPLL's advantage is not one fast modular
operation. It comes from the interaction of its number representation,
specialized transform structure, cache scheduling, and cross-stage fusion.

## Current constraints and diagnosis

For transform length `N = 4,194,304`, the representation has

```text
bits per word = 136279841 / 4194304 = 32.491646
```

A useful first-order upper bound for an uncarried convolution coefficient is

```text
log2(N) + 2 * bits_per_word = 86.983 bits
```

Weighting, multiplication by three, signed balancing, and a safety margin increase
the required range. The current product

```text
(2^31 - 1) * (2^61 - 1)
```

provides approximately 92 bits. PRPLL reconstructs this into an `i96` in
[`src/cl/carryutil.cl`](src/cl/carryutil.cl). Consequently, neither one 61-bit
prime nor FP32 alone can safely replace the current pair at this transform length.

The current data footprint is approximately:

```text
M31 channel: 16 MiB
M61 channel: 32 MiB
Total:       48 MiB
```

The GPU reports 128 MiB of L2 cache. During the measured workload it reaches
approximately 99% GPU utilization and the 300 W power limit, while reported
memory utilization is only around 9--10%. This suggests that the current path is
primarily arithmetic/power limited, although proper Nsight Compute counters should
be collected before treating the coarse memory-utilization figure as conclusive.

The relevant implementation points are:

- Transform-type selection in [`src/FFTConfig.cpp`](src/FFTConfig.cpp).
- `GF31` and `GF61` representations in [`src/cl/base.cl`](src/cl/base.cl).
- M31 and M61 modular arithmetic in [`src/cl/math.cl`](src/cl/math.cl).
- 92-bit CRT reconstruction in [`src/cl/carryutil.cl`](src/cl/carryutil.cl).
- L2-aware queue scheduling in [`src/Gpu.cpp`](src/Gpu.cpp).
- Fused transform/carry/premultiplication and stairway carry forwarding in
  [`src/cl/carryfused.cl`](src/cl/carryfused.cl).

## Strategy 1: replace M31+M61 with three direct 31-bit NTTs

This was the highest-confidence architectural candidate before measurement. The
implementation log below records why it is now a no-go for this GPU and target.

### Candidate primes

The following primes fit in signed 31-bit values and have ample power-of-two root
orders:

| Prime | Factorization of `q - 1` | 2-adicity |
|---:|---:|---:|
| 1,811,939,329 | `27 * 2^26` | 26 |
| 2,013,265,921 | `15 * 2^27` | 27 |
| 2,113,929,217 | `63 * 2^25` | 25 |

Their product has 92.639 bits of range, about 1.557 times the range of the
current M31*M61 product. Each supports a primitive root of order `2^23`, which is
enough for a 4M transform and the expected negacyclic weighting requirement.

### Proposed representation

- Use three scalar `uint32_t` residue planes in structure-of-arrays form.
- Use a direct base-field NTT over each prime rather than copying the existing
  quadratic extension-field representation.
- Keep the total state at 48 MiB: three 16 MiB planes.
- Use Shoup multiplication for fixed twiddles.
- Use Montgomery or specialized Barrett reduction for arbitrary products and
  pointwise squaring.
- Use lazy reduction inside radix kernels where bounds allow values in `[0, 2q)`.
- Reconstruct a balanced result with three-way Garner CRT into a 96-bit type.

The direct-prime design removes all 64x64-to-128 modular products from the main
transform. It also potentially reduces base multiplications because the current
quadratic-field complex multiplication requires several base-field products.

### Important risks

- M31 and M61 are Mersenne primes, so multiplication by powers of two used in
  Crandall--Fagin weighting is exceptionally cheap. The proposed primes require
  ordinary modular multiplication for those weights.
- Three-way CRT is more complicated than the current two-prime reconstruction.
- Three simultaneous INT32 streams may contend for the same execution pipelines,
  whereas the current M31 and M61 queues may overlap somewhat complementary work.
- The existing transform uses extension-field packing and Hermitian-style pair
  handling. The direct scalar transform and its convolution sign convention must
  be rederived, not inferred from buffer layouts.

### Go/no-go criterion

Before integrating a full transform, benchmark:

1. Two direct 31-bit residue channels against the existing GF61 channel.
2. Forward twiddle multiplication, pointwise square, and inverse multiplication.
3. Weight application and the additional CRT cost.

Proceed only if the result projects at least a 15--20% whole-iteration improvement.
A 20--40% latency reduction was the premeasurement research goal. The complete
measurements below reject this path on the tested GPU.

## Strategy 2: shorter transforms with more 32-bit residue channels

Transform length and RNS width can be traded against each other. Approximate
first-order costs are:

| Design | Raw coefficient bound | Relative butterfly work | Main complication |
|---|---:|---:|---|
| 4M x 3 primes | 87 bits | 1.00 | Three-way CRT |
| 3M x 4 primes | 108 bits | 0.98 | Mixed radix and about 65-bit carries |
| 2M x 5 primes | 151 bits | 0.80 | About 160-bit CRT and about 86-bit carries |
| 8M x 2 primes | 55.5 bits | 1.39 | More computation and state traffic |

The work estimate is proportional to

```text
number_of_primes * N * log2(N)
```

and ignores radix, reduction, CRT, and carry costs.

### Experimental 2M prime set

For a 2M design, five near-32-bit primes can provide 159.61 bits of total range:

```text
3,892,314,113
3,942,645,761
4,076,863,489
4,194,304,001
4,253,024,257
```

They are of the form `k * 2^m + 1` with sufficient root order. This design would
use about 40 MiB of residue state and approximately 20% fewer first-order scalar
butterfly operations than 4M x 3 primes.

Its carry path is the likely blocker:

- Nearly every digit is 65 bits.
- The reconstructed coefficient is roughly 151 bits before safety margin.
- The propagated carry can be roughly 86 bits, rather than fitting in `int64_t`.
- The current `i64` carry shuttle would need a multiword or RNS replacement.
- Unsigned primes close to `2^32` require reductions that correctly preserve the
  extra carry bit in Montgomery or Barrett arithmetic.

The 2M design should therefore be investigated only after a working direct-prime
4M implementation provides reusable modular arithmetic and CRT infrastructure.

The 3M x 4-prime design has nearly identical first-order work and storage to 4M x
3 primes. Its extra CRT lane, mixed-radix transform, and slightly-too-wide carry
make it unlikely to win unless a four-lane layout maps unusually well to Tensor
Cores or vectorized memory operations.

## Strategy 3: use FP32 only as an exact algorithm's estimator

FP32 has 24 significant bits, far below the approximately 87 coefficient bits
needed by the 4M transform. FP16, TF32, FP8, and FP4 have still less directly
usable precision. They must not determine the persistent PRP state without an
exact recovery mechanism.

PRPLL already contains FP32+M31, FP32+M61, and FP32+M31+M61 transform modes. The
current tuner chose M31+M61 for this exponent, so simply enabling an existing
FP32 path is not a new solution.

A more useful role for FP32 is quotient estimation during RNS conversion and
carry:

1. Estimate the signed CRT representative or `x / P` from the residues.
2. Estimate the quotient by the current mixed-radix word base.
3. Recover the exact digit and carry using residue checks.
4. Prove a strict error bound limiting correction to a small number of cases.

This could make three-way or five-way CRT cheaper. The existing
`FFT323161` reconstruction in [`src/cl/carryutil.cl`](src/cl/carryutil.cl) already
uses FP32 to estimate how many multiples of M31*M61 must be added and can serve as
a conceptual starting point.

Possible variants include:

- FP32 estimate plus three exact 31-bit residues.
- FP64 estimate plus two exact residues if FP64 cost is confined to carry.
- An extra redundant small modulus used solely to validate and correct the
  approximate quotient.

Every variant needs a proof that an estimation error cannot silently corrupt the
PRP state.

## Strategy 4: exact Tensor Core NTT using integer limbs

Blackwell Tensor Cores do not natively implement 31-bit modular multiplication.
The relevant exact mode is INT8 input with INT32 accumulation through
`tcgen05.mma`.

A possible exact construction is:

1. Split every centered 31-bit residue into four signed 8-bit limbs.
2. Fuse several butterfly stages into a radix-16 or radix-32 matrix tile.
3. Arrange independent transform tiles as matrix columns.
4. Evaluate limb products using INT8 Tensor Core MMA.
5. Keep every partial dot product within signed INT32 range.
6. Recombine limb columns and reduce modulo the selected prime using CUDA cores.

FP16 with roughly 10-bit integer limbs and FP32 accumulation is another
possibility when exact accumulator bounds can be proven. INT8 is preferable
because it has explicit integer semantics.

### Risks

- An NTT butterfly network is sparse, while Tensor Cores are optimized for dense
  matrix multiplication. A naive dense DFT tile performs excessive arithmetic.
- Four-limb by four-limb multiplication can require up to sixteen MMA
  contributions before recombination.
- Packing, layout conversion, modular reduction, and Tensor Memory management may
  cost more than the MMA saves.
- Research results for Tensor-Core NTTs commonly rely on batching many independent
  homomorphic-encryption operations. That is directly relevant when PRPLL is
  allowed to process multiple unrelated exponents together. Each individual chain
  remains serial, but independent chains can populate the columns of dense MMA
  tiles and amortize packing and scheduling costs.
- The current OpenCL-compatible CUDA layer is not an appropriate abstraction for
  `tcgen05`, Tensor Memory, or specialized shared-memory layouts.

The implementation should be a CUDA-native fast path for one complete radix tile,
not a rewrite of the whole program. It should be tested both with tiles from one
large transform and with matrix columns drawn from several exponents. Continue
only if the tile, including packing, reduction, and stores, improves aggregate
throughput substantially over the corresponding SIMT implementation.

Useful references:

- [NVIDIA PTX ISA](https://docs.nvidia.com/cuda/parallel-thread-execution/contents.html)
- [TensorNTT](https://research.polyu.edu.hk/en/publications/tensorntt-architecture-aware-optimizations-for-number-theoretic-t/)
- [TensorFHE](https://arxiv.org/abs/2212.14191)

## Strategy 5: batch unrelated exponents

This is now a primary strategy because the objective is total throughput rather
than the latency of one exponent.

Each PRP chain is sequential, but different exponents are independent. A batched
engine can advance many chains by one iteration in the same launch:

```text
for every active exponent e in the batch:
    state[e] = square_and_reduce(state[e], 2^p[e] - 1)
```

### Batch exponents by transform geometry

The simplest batches contain exponents using the same:

- Transform length and factorization.
- RNS prime set.
- Radix decomposition and kernel variant.
- Carry implementation and buffer layout.

Crandall--Fagin weights and big-word positions still depend on the exponent, so
they remain per-exponent metadata. The roots of unity and most kernel code can be
shared by every exponent in the bucket.

If the input work contains many exponent sizes, maintain several geometry buckets
and launch the fullest bucket first. Sparse buckets can fall back to the existing
single-exponent path.

### Batched memory layout

For scalar CUDA-core kernels, viable layouts include:

- `prime -> exponent -> coefficient`, which gives contiguous coefficients for
  each transform and preserves the existing access patterns.
- `prime -> coefficient tile -> exponent`, which makes the same tile from several
  exponents contiguous and is better for Tensor Core matrix columns.
- An AoSoA layout that groups a small fixed number of exponents per tile while
  retaining coalesced coefficient accesses.

The best layout may differ between the outer transform stages and Tensor Core
tiles. Layout conversion must therefore be included in all performance results.

The 96 GB GPU has ample capacity for many 48 MiB transform states, even after
allowing for work buffers, proof data, twiddles, and checkpoint state. Capacity is
unlikely to limit practical batch sizes; register pressure, active CTAs, power,
and the desired checkpoint interval will set the useful limit first.

### Throughput mechanisms

Batches can improve throughput by:

- Filling Tensor Core MMA columns with matching tiles from independent exponents.
- Amortizing kernel-launch, graph, and host scheduling costs.
- Sharing read-only twiddle tables across all exponents in a geometry bucket.
- Hiding short or imbalanced tail and carry phases from one exponent with useful
  work from another.
- Allowing persistent kernels to draw the next `(exponent, tile)` item from a
  device-side work queue.
- Using separate streams for different geometry buckets when one bucket alone
  cannot occupy the device efficiently.
- Batching residue checks, proof generation, and checkpoint transfers.

The current single-exponent run already reaches the board power limit and high
reported GPU utilization, so merely launching several unchanged instances may
not help. The larger opportunity is to change the kernel shape: batch work into
dense Tensor Core operations, reduce per-exponent metadata traffic, and remove
scheduler bubbles.

### Batched carry handling

Carry propagation is sequential within one exponent but independent across the
batch. Assign different warps, CTAs, or clusters to different exponents during
the carry phase. This provides useful parallelism even if carry propagation along
one number cannot be made fully parallel.

A persistent batched carry kernel could switch to another exponent while one
chain waits for a carry boundary, reducing the cost of the current stairway
dependency. Exponent identity and carry-region identity should be explicit work
queue dimensions.

### Multi-GPU execution

Splitting one NTT across PCIe-connected GPUs is unattractive, but assigning
independent exponents to different GPUs is embarrassingly parallel and should
scale well. A multi-GPU scheduler should:

- Keep each exponent resident on one GPU between checkpoints.
- Balance work using measured iteration cost rather than exponent count.
- Group compatible exponents locally for batching.
- Avoid migrating live transform state unless a GPU fails or becomes unavailable.

### Measurements

For batch sizes `1, 2, 4, 8, 16, ...`, collect:

- Aggregate useful iterations per second.
- Iterations per second per exponent.
- Exponents or candidate work completed per day.
- Joules per useful iteration.
- Median and tail time for advancing a whole batch once.
- Achieved Tensor Core, integer-pipe, L2, and memory throughput.
- Clock rate and power-limit behavior.

An increase in per-exponent latency is acceptable when aggregate useful throughput
or energy efficiency improves.

## Strategy 6: CUDA-native tile fusion, clusters, DSM, and TMA

Reducing GDDR7 traffic alone is unlikely to produce a major improvement because
the current state fits in L2 and the run appears compute/power limited. Memory
work is valuable when it also removes synchronization, address generation,
launches, or redundant transformations.

Potential changes include:

- Fuse compatible portions of `fftMiddleIn`, `tailSquare`, and `fftMiddleOut` so
  that an intermediate tile stays in registers, shared memory, or distributed
  shared memory.
- Use thread-block clusters when a tail tile requires data produced by multiple
  CTAs.
- Use Tensor Memory Accelerator operations for multidimensional transfers and
  reduce address-generation instructions.
- Use asynchronous double buffering between L2 and shared memory.
- Use structure-of-arrays prime planes and aligned vector loads; avoid a 12-byte
  `uint3` array-of-structures layout.
- Partition the 128 MiB L2 working set deliberately among residue channels and
  twiddle tables.

PRPLL already uses cache-aware replay, multiple queues, dependent launches, and a
fused carry kernel. A new implementation must eliminate actual global round trips
or synchronization; merely replacing launches with another launch mechanism is
unlikely to help.

Relevant hardware documentation:

- [NVIDIA Blackwell Tuning Guide](https://docs.nvidia.com/cuda/blackwell-tuning-guide/)
- [RTX Blackwell architecture](https://www.nvidia.com/content/dam/en-zz/Solutions/design-visualization/quadro-product-literature/NVIDIA-RTX-Blackwell-PRO-GPU-Architecture-v1.0.pdf)

## Strategy 7: redesign carry propagation

The profiled `carryFused` region is the largest named region, although timing is
complicated by overlap between the M31 and M61 command queues. It includes the
ending width transform, CRT, carry propagation, and premultiplication, so its
entire time cannot be attributed to carry forwarding.

The current kernel uses stairway forwarding through a global `carryShuttle` and
readiness markers. Architectural alternatives are:

### Block transfer functions and prefix scan

- Normalize each block assuming a canonical incoming carry.
- Represent the block's effect on a bounded incoming carry as a small transfer
  function.
- Compose transfer functions with a prefix scan.
- Apply incoming carries with a short correction pass.

This is only attractive if the post-local incoming carry has a provably small
state space. The alternating big-word/little-word radix must be included in the
proof.

### Persistent cluster carry

- Assign persistent CTAs or clusters to ordered carry regions.
- Keep shuttle data in distributed shared memory where possible.
- Use decoupled lookback or cluster barriers rather than global spin polling.
- Fuse the correction into the next premultiplication stage.

### Carry in RNS form

For shorter-transform designs with carries wider than 64 bits, retain part of the
carry in RNS form and reconstruct only the quotient information needed to emit a
digit. Approximate FP32 quotient estimation plus an exact redundant modulus may
help, but this is a mathematical research project rather than a straightforward
kernel rewrite.

## Strategy 8: twiddle representation and generation

Large GPU NTTs can be bandwidth-limited by twiddle tables, and published work has
shown benefits from generating some roots on the fly. For this workload, extra
modular multiplication is expensive and the state plus current tables should fit
in L2, so full on-the-fly generation is not the default recommendation.

Possibilities to measure are:

- Store one starting root and one step root per tile, then generate a short chain.
- Store three aligned 32-bit twiddles for the proposed RNS channels in one compact
  record.
- Generate trivial roots and quarter/eighth rotations algebraically.
- Share exponent/index calculations across residue channels while keeping the
  residue values in separate planes.
- Fuse twiddle generation with TMA tile loading to hide some arithmetic latency.

Any experiment must include the effect on clocks and power, not only instruction
count, because the current run is at the board power limit.

Reference: [Accelerating Number Theoretic Transformations for Bootstrappable
Homomorphic Encryption on GPUs](https://arxiv.org/abs/2012.01968).

## Strategy 9: carry-save or delayed normalization

Avoiding a forward/inverse transform by retaining the state in the NTT domain
would be extremely valuable, but normal carry propagation is nonlinear in that
domain. The next modular square depends on the normalized mixed-radix state.
Without normalization, coefficient widths approximately double on every square
and quickly exceed any practical RNS range.

Research directions that could revisit this conclusion include:

- A redundant digit system with a provable multi-iteration growth bound.
- Carry-save normalization that resolves only enough carries to keep the next
  convolution bounded.
- A polynomial-ring representation in which reduction modulo `2^p - 1` becomes a
  cheap linear operation.

These ideas have potentially enormous upside but very low near-term confidence.
Any proposal must account for the fact that `p` is not divisible by the transform
length, which is why Crandall--Fagin weighting and mixed word sizes are necessary.

## Strategies unlikely to help this target

- **Pure FP32/TF32/FP16/FP8 state:** insufficient precision without exact
  residues and correction.
- **Direct Tensor Core use on M61:** no native 61-bit modular MMA operation.
- **RT cores:** fixed-function ray traversal/intersection hardware is not
  programmable for modular arithmetic.
- **NVENC/NVDEC:** unrelated fixed-function media engines.
- **Splitting one exponent across multiple GPUs:** every NTT would require
  communication, and PCIe synchronization is likely to dominate. Assigning whole,
  unrelated exponents to different GPUs remains highly attractive.
- **Pruned transforms:** after a few iterations the state is dense and all output
  coefficients are required.
- **Blindly reducing transform length:** the extra RNS range and carry width must
  be counted; the transform alone is not the complete cost.
- **Only reducing DRAM bytes:** the present transform appears to fit in L2 and is
  already power/compute saturated.

## Original implementation sequence (now evaluated)

### Phase 0: correctness and cost model

1. Derive a rigorous coefficient bound for exponent `136279841`, including
   weights, `MUL3`, signed balancing, and maximum carry.
2. Confirm the exact required root order and cyclic/negacyclic convention.
3. Implement CPU reference arithmetic for the candidate primes and Garner CRT.
4. Model transform, CRT, carry, and memory costs separately.

### Phase 1: modular arithmetic microkernels

1. Implement direct-prime add, subtract, Shoup multiply, Montgomery multiply, and
   square in CUDA.
2. Compare two candidate 31-bit channels against the existing GF61 arithmetic.
3. Measure instructions, integer-pipe utilization, achieved clocks, power, L2 hit
   rate, and elapsed time with Nsight Compute.

### Phase 2: one complete transform tile

1. Implement a CUDA-native radix tile for one prime.
2. Extend it to three residue planes.
3. Benchmark columns from one exponent and from batches of unrelated exponents.
4. Include loads, twiddles, reduction, pointwise squaring, inverse scaling, and
   stores in the benchmark.
5. Validate every output against a CPU reference.

### Phase 3: experimental 4M RNS31x3 transform

1. Add a new transform type without disturbing existing FFT/NTT paths.
2. Add three-plane buffer and twiddle management.
3. Add exact three-way CRT and the current carry behavior.
4. Compare complete iteration residues against M31+M61.
5. Run progressively longer comparisons before a 1M-iteration validation.

### Phase 4: carry and fusion work

1. Profile the new path and separate CRT, local carry, shuttle waiting, final
   transform, and premultiplication costs.
2. Prototype FP32 quotient estimation with exact correction.
3. Prototype cluster tile fusion and a carry prefix/transfer-function design.

### Phase 5: higher-risk experiments

1. Add a persistent multi-exponent scheduler and geometry buckets.
2. Implement a batched INT8 Tensor Core radix tile.
3. Compare single-exponent and multi-exponent layouts at several batch sizes.
4. Prototype the 2M x 5-prime transform with multiword carry.
5. Investigate delayed or redundant carry representations.

## Validation requirements

Performance changes must not weaken correctness. Every new exact path should:

- Match current M31+M61 residues on every tested iteration.
- Test CRT values around `0`, `P/2`, and `P-1` for each combined modulus `P`.
- Test maximum big-word and little-word digits and positive/negative carries.
- Detect coefficient-range exhaustion before ambiguity is possible.
- Preserve save/resume and proof-generation formats, or introduce a versioned
  conversion path.
- Run the existing self-tests and known-residue tests.
- Complete at least one 1M-iteration A/B validation on exponent `136279841`.
- For batch tests, validate every exponent independently against its unbatched
  reference path.
- Report aggregate throughput, per-exponent throughput, median, mean, standard
  deviation, percentiles, power, energy per iteration, and clock rate after
  warm-up; one short timing is insufficient on a power-limited Max-Q GPU.

## Original priority order (all entries now evaluated)

1. Direct 4M, three-prime INT32 RNS NTT.
2. Multi-exponent batching by transform geometry.
3. Exact batched INT8 Tensor Core radix tiles.
4. Three-way CRT/carry optimized with an FP32 estimate and exact correction.
5. CUDA-native cluster and shared-memory fusion of transform stages.
6. Shorter 2M transform with five 32-bit primes and multiword carry.
7. Delayed normalization or a new redundant representation.

## Implementation status and experiment log

This section is the authoritative running record of implementation work. A
candidate is not a production speedup until it is integrated, correctness-tested,
and measured in the relevant full PRPLL workload. A branch may instead be closed
as a no-go when a correct standalone lower bound is already slower than the whole
production iteration or when a mathematical bound excludes it.

The log is chronological. Its `Next work` lists record the gate identified at
that point in the investigation; later entries and the completion matrix record
how each gate was resolved.

### 2026-08-09: production timing distribution

The original requested 1M-iteration run reports 205.5 us/iteration after its
2,000-iteration startup. Fourteen subsequent 100k runs with the identical
production FFT geometry and compile settings give:

```text
sample count:       14
mean:          203.514 us/iteration
median:        203.900 us/iteration
sample SD:       2.301 us/iteration
minimum:        199.500 us/iteration
maximum:        206.400 us/iteration
```

A final 30k confirmation after the GPU had been idle measured 198.5 us/iteration
and produced the expected residue. The sequential warm runs rose from 202.3 to
206.4 us/iteration, so the 198.5 us result is a useful cold/short-run lower bound,
not a sustainable long-run mean. Architectural comparisons therefore use about
0.204--0.206 ms as the production baseline. This thermal variability does not
affect any direct-prime go/no-go conclusion because the best alternative already
takes 0.708 ms before adding CRT or carry.

### 2026-08-09: three-prime reference arithmetic

Implemented:

- Added shared prime configuration in [`src/RNS31Config.h`](src/RNS31Config.h).
- Added reusable host reference arithmetic in [`src/RNS31.h`](src/RNS31.h) and
  [`src/RNS31.cpp`](src/RNS31.cpp).
- Implemented modular add/subtract/reference multiply, primitive power-of-two
  roots, Montgomery conversion/multiplication, Shoup constant multiplication,
  residue splitting, three-way Garner reconstruction, and balanced CRT output.
- Added [`src/rns31_test.cpp`](src/rns31_test.cpp), a standalone test executable
  available as `make CUDA=1 rns31-test` and as CTest `rns31-reference`.

Measured outcome:

```text
RNS31 reference tests passed
combined modulus: 0x18eac000a2d8000162000001
combined modulus bits: 92.639058023
```

The test covers:

- Root orders through `2^23` for every prime.
- 200,000 randomized Montgomery and Shoup products per prime, for 600,000 of
  each multiplication method in total.
- Modular addition and subtraction on the same randomized inputs.
- 200,000 randomized CRT round trips.
- CRT boundary values around zero, `P/2`, and `P-1`.
- Balanced positive and negative reconstruction.

CTest also passes from a clean temporary CMake build.

Outcome: the proposed three-prime set and 92.639-bit CRT range are now backed by
executable reference code. This completes the host arithmetic portion of Phase 0,
but not the rigorous Crandall--Fagin coefficient/carry bound or GPU integration.

### 2026-08-09: native Blackwell arithmetic benchmark

Implemented:

- Added [`src/cuda/rns31_cuda_bench.cu`](src/cuda/rns31_cuda_bench.cu).
- Added the `make rns31-cuda-bench` target, compiled specifically for `sm_120`.
- Implemented and CPU-validated four arithmetic-chain kernels:
  - Current GF61 quadratic-extension multiplication.
  - Two direct 31-bit prime channels, representing a possible GF61 replacement.
  - Current combined GF31+GF61 arithmetic.
  - Three direct 31-bit prime channels.
- The direct kernels use `__umulhi` Shoup multiplication; the comparison kernels
  use Mersenne reductions and `__umul64hi` for M61.
- Implemented a complete three-prime radix-16 cyclic-square tile: forward DIF
  NTT, Montgomery pointwise square, inverse DIT NTT, and inverse scaling. The tile
  uses 16-lane shuffle groups and processes two packed scalar transforms per
  thread group.

For 1,048,576 pairs, 21 timing samples, and increasing chain lengths, median
kernel-time ratios were:

| Multiplies per loaded value | `2q31 / GF61` | `3q31 / current` |
|---:|---:|---:|
| 1 | 1.2128 | 1.1442 |
| 2 | 1.0704 | 1.1627 |
| 4 | 0.8656 | 0.9011 |
| 8 | 0.7648 | 0.7078 |
| 16 | 0.6906 | 0.6956 |
| 32 | 0.6472 | 0.6724 |
| 64 | 0.6296 | 0.6611 |
| 128 | 0.6201 | 0.6560 |

At chain length 32, the representative measurement was:

```text
GF61 extension          0.088 ms
two direct q31          0.057 ms
current GF31+GF61       0.118 ms
three direct q31        0.079 ms
```

Every arithmetic-chain run validated 256 output pairs against CPU arithmetic. The
radix-16 kernel was independently validated against direct cyclic convolution in
all three primes and both packed streams. For 1,048,576 pairs it took 0.044 ms,
or 47.4 billion coefficients per second through a forward-transform, square, and
inverse-transform tile. `ptxas` reports:

```text
current combined arithmetic: 38 registers, no spills
three-prime arithmetic:       38 registers, no spills
three-prime radix-16 square:  36 registers, no spills
```

Interpretation:

- With only one or two modular products per load, the direct-prime path loses due
  to extra lanes and memory/launch overhead.
- From four products per load onward, it wins.
- In arithmetic-dense chains it is approximately 30--35% faster than the current
  combined arithmetic proxy.
- This supports continuing to a full radix tile, but does not yet predict a
  30--35% PRPLL speedup. Real transforms include butterflies, shared-memory
  exchanges, twiddle loads, transpose stages, weighting, CRT, and carry.

Next work:

1. Compose radix tiles into a complete multi-tile direct-prime NTT and establish a
   correct, deliberately simple full-transform baseline.
2. Compare scalar SoA and batched-exponent AoSoA layouts.
3. Replace global radix-2 passes with radix-8/radix-16 stage fusion and transpose
   tiles modeled on PRPLL's width/middle/height decomposition.
4. Extend cyclic-square validation to the exact weighted convolution convention
   needed by PRPLL.
5. Integrate the winning transform behind a new experimental transform type.

### Completion matrix

| Architectural strategy | State | Evidence and final gate |
|---|---|---|
| 4M direct three-prime INT32 NTT | Complete standalone, no-go | Full correct transforms reach 0.707 ms with Montgomery roots before CRT/carry, versus about 0.200 ms for the complete current iteration |
| Transform-length/RNS-width trade | Measured, no-go | Exact 2M x five-prime is 4% slower than 4M x three-prime before wider carry; exact 8M x two-prime takes 0.966 ms before CRT/carry |
| FP32-assisted exact quotient/CRT | Complete standalone, neutral/no-go | Exact estimator matches Garner; both take 0.161 ms for 8.39M reconstructions, so it does not rescue the slower RNS transform |
| Exact Tensor Core NTT | Tile measured, no-go | Exact INT8 limb matrix is 42% slower than sparse SIMT radix-16 with packing/reduction included |
| Multi-exponent batching | Scheduler complete, scalar no-go | `-prps` implemented; 2/3/4 unchanged workers lose 13%/28%/35% aggregate throughput; Tensor batch tile also loses |
| CUDA cluster/DSM/TMA fusion | Assessed, no integration gate | Carry waits only 18.4 short polls; direct CUDA tile is over 3x too slow, leaving no plausible TMA/DSM recovery margin |
| Carry propagation redesign | Reference + instrumented, no-go | Block transfers are exact constant/binary functions, but the current fused stairway already performs local speculative carry and waits only 18.4 short polls |
| Twiddle generation/layout experiments | Measured, no-go overall | Full Montgomery roots give a 24% one-plane gain; 16x generation is 39% slower; three full-root planes still take 0.708 ms |
| Delayed/redundant carry representation | Bounded, no-go | One unnormalized iteration raises the next convolution to about 200 bits and at least seven q31 channels |
| Packed M31 + two Riesel quadratic planes | Complete, no-go in present form | Same 48 MiB state and 92.957-bit CRT; the tile lead did not survive production weighting, normalization, CRT, and scheduling |
| Compact two-phase M31/M61 middle/carry fusion | Exact end to end, no-go | Correct through 10k, but 311.0 us/iteration versus a nearby 230.3 us control; exact halves cost about 210 us and use 98/90 registers with a 512-thread, 96 KiB block |
| Radix-`2^31` limb representation of M61 | Exact arithmetic gate, no-go | Same M61 field and bytes, but explicit three-product limb arithmetic is 55--57% slower than the generated 64-bit path in arithmetic-dense chains |
| Riesel pseudo-Mersenne shift/add reduction | Exact tile, no-go | Four serial folds make the tile about 2x slower; generic Montgomery wins |
| M31 + Goldilocks | Exact algebraic tile, no-go | Fast scalar roots violate conjugate packing; correct unit-norm tile is 11.6% slower than current |
| Proth-prime packed replacement | Mathematical/search no-go | No relevant sub-2^31 candidate gives 2 a 2^22-th root in the base or quadratic field |

### 2026-08-09: complete direct-prime NTT baseline and first fusion attempt

Implemented:

- Added reusable native CUDA modular arithmetic in
  [`src/cuda/rns31.cuh`](src/cuda/rns31.cuh).
- Added [`src/cuda/rns31_ntt_bench.cu`](src/cuda/rns31_ntt_bench.cu) and the
  `make rns31-ntt-bench` target.
- Implemented a complete direct three-prime power-of-two NTT baseline using one
  fused-prime radix-2 CUDA kernel per stage.
- Implemented forward DIF, Montgomery pointwise square, inverse DIT, and inverse
  scaling for sizes through `2^25` where supported by the selected primes.
- Added generated Shoup root tables for every stage and prime.
- Added CUDA Graph measurement of the complete operation.
- Added a first radix-fusion experiment combining up to four consecutive stages
  with 4- or 16-lane shuffle groups.

Correctness results:

- Dense forward/inverse round trips pass for every coefficient and all three
  primes at sizes 1K, 1M, and 4M.
- Sparse cyclic squares pass against independently calculated CPU convolution at
  the same sizes.
- The fused-stage path independently passes both validations.

Median timing with 15 samples:

| Transform size | Radix-2 launches | Radix-2 graph | Fused radix | Fused graph |
|---:|---:|---:|---:|---:|
| 1M | 0.244 ms | 0.236 ms | 0.212 ms | 0.211 ms |
| 4M | 0.949 ms | 0.942 ms | 0.937 ms | 0.935 ms |

At 4M the radix-2 operation uses 46 kernels. The fused version uses only 14, but
does not materially improve elapsed time. Its nominal full-state traffic is much
lower, yet its effective traffic rate falls from approximately 4.58 TiB/s to 1.40
TiB/s.

Outcome:

- The direct-prime arithmetic, roots, full NTT, pointwise square, and inverse are
  now executable and correct at the target transform length.
- The simple implementation is approximately 4.6 times slower than the complete
  current PRPLL iteration and cannot be integrated as a performance path.
- CUDA Graphs save less than 1% at 4M, confirming that launch overhead is not the
  main problem at this size.
- Naively fusing radix stages while accessing widely strided elements sacrifices
  coalescing and L2 efficiency. Fewer global passes do not help when each pass has
  a poor memory layout.

Next work:

1. Replace strided shuffle fusion with coalesced shared-memory transpose tiles.
2. Factor the 4M scalar transform into width/middle/height dimensions modeled on
   the existing `512 x 8 x 512` extension-field implementation.
3. Keep several local stages in registers/shared memory, transpose once, then run
   the next local dimension contiguously.
4. Compare separate-prime kernels with fused-three-prime kernels; the latter may
   lose occupancy or cache efficiency despite sharing address calculations.
5. Integrate only after the standalone direct-prime square approaches or beats
   the current full-iteration time with enough headroom for CRT and carry.

### 2026-08-09: coalesced radix fusion and high-radix limits

Implemented three more complete, independently selectable transform paths in
[`src/cuda/rns31_ntt_bench.cu`](src/cuda/rns31_ntt_bench.cu):

- A coalesced radix-16 shared-memory transpose tile.
- A radix-32 shared-memory tile.
- A radix-256 decomposition, first with shared-memory local transforms and then
  with a register/warp-shuffle implementation.

Every path passes dense forward/inverse round trips and sparse cyclic-square
checks against the CPU reference for all three primes at 1K, 1M, and 4M.  The
4M operation timings (forward transform, pointwise square, inverse transform,
and inverse scale; 21 samples) were:

| Complete 4M three-prime operation | Global kernels | Median time |
|---|---:|---:|
| One radix-2 stage per kernel | 46 | 0.949 ms |
| Naive fused shuffle stages | 14 | 0.937 ms |
| Coalesced radix-16 tiles | 12 | **0.795 ms** |
| Coalesced radix-32 tiles | 10 | 0.844 ms |
| Radix-256 shared-memory decomposition | 8 | 1.226 ms |
| Radix-256 register/warp decomposition | 8 | 0.991 ms |

At 1M, the corresponding best radix-16 result was 0.165 ms, versus 0.244 ms for
the radix-2 baseline. CUDA Graph capture did not improve the radix-16 result and
was somewhat variable. `ptxas` reports approximately 40 registers and no spills
for the tiled kernels; the register-heavy radix-256 forward/inverse kernels use
62/56 registers without spills.

Outcome:

- Restoring coalescing makes stage fusion worthwhile: radix-16 is about 16%
  faster than the simple 4M radix-2 path and about 21% faster than it at 1M.
- Kernel/pass count is not a sufficient optimization target. Radix 32 and 256
  spend progressively more time in modular products, root addressing,
  synchronization, and shared/register exchange than they save in global
  traffic.
- Even the best generic radix-16 path is roughly four times the current tuned
  PRPLL iteration time. A generic full scalar NTT is therefore not an integration
  candidate. Further direct-prime work must use a PRPLL-style specialized
  width/middle/height factorization and fuse weighting/carry work, or it must gain
  throughput by batching independent exponents.

Next work:

1. Measure aggregate throughput with one through four independent exponents on
   the existing CUDA path; this is now a first-class target because latency of an
   individual exponent is not important.
2. Identify whether concurrent workers fill unused execution resources or only
   split the existing 300 W power budget.
3. Use that result to choose between an exponent-interleaved AoSoA transform and
   independent CUDA streams/worker contexts before doing more scalar NTT work.

### 2026-08-09: multi-exponent scheduler and throughput experiment

Implemented:

- Added `-prps <e1,e2,...>` to [`src/Args.cpp`](src/Args.cpp) and a synchronized
  in-memory batch source in [`src/Worktodo.cpp`](src/Worktodo.cpp).
- A command such as `-prps 136279841,136279879 -workers 2` now assigns unrelated
  exponents directly to PRPLL's existing worker threads and CUDA queues. It does
  not require per-worker `worktodo-N.txt` setup.
- Built the CUDA backend and exercised the new path with one through four nearby
  prime exponents, all forced to the tuned `4M 1:512:8:512:202` geometry.

Steady timings exclude the initial 2,000-iteration warm-up.  The one- and
two-worker runs used 100,000 requested iterations per exponent; three workers
used 50,000 and four workers used 30,000.  Each stopped at the requested count
with a valid residue/check value.

| Workers | Tail mode | Per-exponent us/iteration | Aggregate iterations/s | Change from one worker |
|---:|---:|---|---:|---:|
| 1 | 2 | 199.6 | **5,010** | baseline |
| 2 | 2 | 486.5, 486.5 | 4,111 | -17.9% |
| 2 | 3 | 457.1, 457.2 | 4,375 | -12.7% |
| 3 | 3 | 838.0, 835.9, 836.3 | 3,585 | -28.4% |
| 4 | 3 | 1227, 1226, 1228, 1229 | 3,259 | -35.0% |

The documented two-worker `TAIL_KERNELS=3` split is 6.4% better than tail mode 2
for two workers, but it still loses materially to a single exponent. With two
workers, telemetry showed 100% GPU activity, about 52% memory-controller
activity, a sustained 300 W board-power limit, and an SM clock around 1.22 GHz.
Three workers raised memory-controller activity to roughly 70%, but aggregate
useful work fell further.

Outcome:

- Independent CUDA streams running the existing kernels are a no-go on this
  300 W Max-Q GPU. They increase resource activity but exceed the efficient power
  envelope and reduce total useful iterations per second.
- The scheduler remains useful operationally and as the launch mechanism for
  future batch-native kernels.
- Batching is only worth further implementation if it changes the computation:
  for example, exponent-interleaved AoSoA tiles that share metadata, Tensor Core
  matrix columns populated by different exponents, or a persistent queue that
  removes otherwise idle phases. Merely overlapping unchanged exponent streams
  must not be presented as a throughput optimization.

Next work:

1. Keep batch size one as the default on this device.
2. Prototype batch-interleaved radix tiles in the standalone harness and include
   layout conversion in their timing.
3. Revisit batch size only when a tile uses a hardware path unavailable to the
   scalar single-exponent implementation.

### 2026-08-09: existing FP32 precision boundary and exact hybrid cost

Two existing architectures were forced at exponent `136279841` and the same 4M
geometry:

- `FFT3261` (FP32+M61) has a documented maximum exponent of 132,791,664 for this
  variant. At the target it warmed to approximately 171.0 us/iteration, but
  failed at iteration 2,000, reloaded, and produced the identical failure again.
  Both runs reported `ROEmax=0.500` and PRPLL stopped on the consistent error.
- `FFT323161` (FP32+M31+M61) restores ample exact range and produced the same
  checked residues as the baseline. Its steady 30,000-iteration timing was
  258.7 us/iteration, about 29.6% slower than the 199.6 us comparison baseline.

Outcome: lower precision is not merely blocked by a conservative table entry;
the faster existing hybrid crosses a hard rounding boundary at this exponent.
Adding exact residues as a complete parallel transform removes the speedup. An
FP32 contribution is only promising if it replaces work inside CRT/quotient
recovery rather than adding another complete transform channel.

### 2026-08-09: steady kernel profile

A 10,000-iteration `-time` run on the tuned baseline measured:

| Named region | Average kernel time | Steady profile share |
|---|---:|---:|
| `kCarryFused` | 71.5 us | 37.87% |
| `kfftMidOutGF31` | 51.1 us | 28.35% |
| `ktailSquareGF31` | 40.0 us | 22.17% |
| `kfftMidInGF31` | 17.8 us | 9.87% |
| periodic ROE carry | 66.7 us when invoked | 1.64% amortized |

`kCarryFused` includes the final width transform, two-prime CRT, carry
propagation, reweighting, and the first transform of the next iteration. The
profile confirms that carry redesign has a large ceiling, but it does not imply
that 71.5 us is pure serial carry overhead.

### 2026-08-09: exact 2M x five-prime transform

Implemented:

- Added five verified unsigned-32-bit NTT primes and primitive roots to
  [`src/RNS31Config.h`](src/RNS31Config.h).
- Added overflow-aware `addFull32` and `montgomeryReduceFull32` CUDA primitives.
  REDC explicitly retains the carry from the 64-bit `value + multiplier*q` sum.
- Added [`src/cuda/rns_multi_ntt_bench.cu`](src/cuda/rns_multi_ntt_bench.cu) and
  `make CUDA=1 rns-multi-ntt-bench`.
- The same templated coalesced radix-16 pipeline can run either 4M x three
  signed-31-bit primes or 2M x five full-32-bit primes.

Both configurations pass a dense full-array forward/inverse round trip and a
sparse cyclic square for every residue plane. Timings are complete forward,
pointwise square, inverse, and scale operations with 21 samples:

| Architecture in generic harness | CRT range | Median | Mean |
|---|---:|---:|---:|
| 4M x three q31 | 92.639 bits | 0.929 ms | 0.925 ms |
| 2M x five q32 | 159.611 bits | 0.826 ms | 0.823 ms |

The shorter transform is 11.1% faster under identical generic code. The earlier
three-prime-specialized radix-16 path, however, is 0.795 ms, making the five-prime
path about 3.9% slower than the best relevant comparator. `ptxas` reports 40
registers and no spills for all five-prime transform kernels.

Outcome: the five-prime path realizes part of the theoretical 20% butterfly-work
reduction, but not enough. It loses before reconstructing 160-bit coefficients,
emitting 65-bit digits, or propagating approximately 87-bit carries. It is a
no-go for integration unless a future specialized layout improves it by a large
margin and simultaneously supplies a cheap wide-carry representation.

### 2026-08-09: exact INT8 Tensor Core radix-tile experiment

Implemented [`src/cuda/rns_tensor_bench.cu`](src/cuda/rns_tensor_bench.cu) and
`make CUDA=1 rns-tensor-bench`:

- Split each residue into four unsigned base-256 limbs.
- Evaluate a 16x16 modular DFT matrix over 16 columns with sixteen
  `u8 * u8 -> s32` MMA operations.
- Every MMA dot product is bounded by `16 * 255^2 = 1,040,400` and is exact.
- Defer modular reduction of the sixteen limb matrices because their weighted
  sum is below `2^55`.
- Include input limb packing, MMA, recombination, reduction, and output stores.
- Validate every output against both a dense Shoup CUDA kernel and an independent
  sparse four-stage radix-16 NTT kernel.

For 32,768 tiles and 8,388,608 input/output coefficients:

| Exact radix-16 operation | Median time |
|---|---:|
| INT8 Tensor Core dense matrix | 0.069 ms |
| SIMT dense Shoup matrix | 0.075 ms |
| SIMT sparse radix-16 butterflies | **0.049 ms** |

SASS contains `IMMA.16816.U8.U8`; this is genuine Tensor Core execution, not a
compiler scalarization. Both Tensor and dense SIMT kernels use 40 registers with
no spills; sparse radix-16 uses 11.

Outcome: limb MMA is 8% faster than a deliberately dense scalar matrix, proving
the exact construction works, but it is 42% slower than the actual sparse NTT
butterfly network. Supplying sixteen simultaneous columns already captures the
batching advantage that unrelated exponents would provide. Dense Tensor Core
work does too much arithmetic and packing to replace radix-16 SIMT on this GPU.
Larger limb tiles should not be integrated unless they eliminate other work in
addition to the NTT itself.

### 2026-08-09: fused versus separate residue planes

The generic radix-16 harness can also instantiate a one-prime pipeline. One 4M
plane repeatedly transformed in isolation takes 0.294 ms, suggesting 0.882 ms
for three planes if each plane's state and root tables remained cache-hot.
Running all three one-prime pipelines sequentially in the same allocation instead
takes 0.968 ms, while the fused-three-prime generic kernel takes 0.929 ms.

Outcome: the apparent isolated-plane advantage disappears when the three root
table working sets actually compete for L2. Fusing planes shares indexing and
launch work and is about 4% faster than realistic sequential plane scheduling.
The 0.795 ms explicitly specialized fused kernel remains the best direct-prime
result. Splitting primes into separate kernels is therefore not a memory-
architecture solution; future work would need genuinely compressed/generated
twiddles rather than relying on one plane staying resident.

### 2026-08-09: Montgomery state and twiddle compression

Implemented [`src/cuda/rns_mont_ntt_bench.cu`](src/cuda/rns_mont_ntt_bench.cu)
and `make CUDA=1 rns-mont-ntt-bench` with two complete NTT paths:

1. Keep the residue state and all roots in Montgomery form. Each twiddle entry
   is one 32-bit value instead of a normal root plus a 32-bit Shoup factor.
2. Store only every sixteenth Montgomery root plus sixteen small powers per
   stage, generating the exact root with one extra Montgomery product.

Both paths pass dense full-array round trips and sparse cyclic squares for all
three selected primes. Timings for one 4M plane are:

| Montgomery root representation | Forward+inverse root storage | Complete operation |
|---|---:|---:|
| Full root tables | 32 MiB | **0.224 ms** |
| 16x compressed/generated | about 2 MiB | 0.312 ms |

The comparable one-prime Shoup path in the generic harness is 0.294 ms. Thus,
full Montgomery roots improve a complete one-plane operation by about 24%, while
generating each twiddle makes it about 39% slower than full Montgomery tables.

With all three primes and their working sets allocated, three sequential
full-root Montgomery pipelines take 0.707 ms. Three compressed pipelines take
0.931 ms. The full-root result is 11% faster than the best 0.795 ms specialized
fused Shoup path even without fusing address calculations across primes.

Outcome:

- Montgomery-domain persistent state is the winning direct-prime representation.
  Halving twiddle bytes matters more than the arithmetic difference between
  Montgomery and Shoup multiplication on Blackwell.
- Aggressive on-the-fly generation is a no-go: one additional modular product
  per twiddle overwhelms the 30 MiB/plane table saving.
- At 0.707 ms before CRT/carry, even the winning three-plane scalar transform is
  still more than three times a complete tuned PRPLL iteration. It is therefore
  valuable infrastructure and a measured design improvement, but not a path that
  should replace the current extension-field implementation.

Next work:

1. If direct primes are revisited, fuse the three full-Montgomery planes and keep
   only forward or inverse roots live in the active cache phase.
2. Do not pursue finer root compression unless generation is fused with an
   otherwise-required modular product.

### 2026-08-09: bounded carry transfer-function reference

Implemented [`src/CarryTransfer.h`](src/CarryTransfer.h),
[`src/CarryTransfer.cpp`](src/CarryTransfer.cpp), and
[`src/carry_transfer_test.cpp`](src/carry_transfer_test.cpp), available through
`make CUDA=1 carry-transfer-test` and CTest `carry-transfer-reference`.

The reference models PRPLL's exact centered power-of-two recurrence:

```text
digit_i = centered_low_bits(coefficient_i + carry_i, bits_i)
carry_(i+1) = (coefficient_i + carry_i - digit_i) / 2^bits_i
```

The recurrence is monotone in its incoming carry. After at least two 32/33-bit
words, an incoming interval of width at most `2^60` is divided by a product base
of at least `2^64`. Consequently, a block transfer over the conservative
`[-2^59, 2^59]` input interval has either one output or two adjacent outputs. In
the latter case it is represented exactly by a single threshold.

Tests generated 20,000 random blocks of 2--16 approximately 90-bit coefficients,
using the actual big/little-word sequence for exponent `136279841`. Results were:

```text
19,959 constant transfer functions
    41 binary threshold transfer functions
```

Every boundary and randomized input sampled from the carry interval matched
direct word-by-word propagation.

Compiled-kernel inspection of the actual tuned `carryFused` PTX shows 96
registers, 9,088 bytes of shared memory, and no spills for the steady kernel.
The ROE form has one 8-byte spill. Nsight Compute counter collection was attempted
but the driver rejected performance-counter access with `ERR_NVGPUCTRPERM`.

Outcome: the mathematical condition needed for a block-transfer scan is valid,
and almost all realistic randomized blocks are input-independent under this
conservative range. A production design would perform local normalization,
scan the small constant/binary transfer descriptors, and correct only the block
prefix affected by the true incoming carry.

The remaining architectural issue is synchronization and data lifetime. The
current fused kernel retains normalized words plus both transform channels in
registers while using stairway readiness flags. A conventional two-kernel scan
would spill that large intermediate state to global memory, likely losing the
benefit. The next viable implementation should therefore use a CUDA cooperative
grid or persistent clustered kernel with a grid-wide phase boundary; it should
not split the current fused region without measuring the extra state traffic.

### 2026-08-09: direct carry-wait instrumentation and cluster decision

Added an opt-in `SPIN_STATS=1` instrumentation path to
[`src/cl/carryfused.cl`](src/cl/carryfused.cl). Normal kernels are unchanged.
With `STATS=1,SPIN_STATS=1`, each steady fused carry invocation reduces the
readiness-poll counts across all carry groups and records their maximum through
the existing statistics buffer.

Over the 8,598 steady instrumented iterations after warm-up:

```text
maximum observed readiness polls:       29
mean of per-iteration group maxima:  18.389
standard deviation:                   1.980
```

On NVIDIA the `spin()` body has no sleep instruction, so one poll is an atomic
load and short loop. The instrumented run remained correct and took 202.2
us/iteration, close to the uninstrumented 199.6 us comparison.

Outcome: global readiness polling is not a large hidden fraction of the 71.5 us
fused region. A cooperative-grid or DSM cluster rewrite aimed only at replacing
those waits cannot provide an architectural speedup and risks losing occupancy
or forcing intermediate stores. TMA/DSM remain useful only if a future kernel
fusion keeps an otherwise-global transform tile on chip; the current direct-
prime path is more than 3x too slow for transfer machinery alone to close its
gap, so no production cluster/TMA integration is justified by these results.

### 2026-08-09: delayed normalization bound

For a normalized 33-bit signed digit the magnitude is at most `2^32`. A 4M
convolution therefore has a basic centered bound of

```text
2^22 * (2^32)^2 = 2^86
```

before Crandall--Fagin weight and multiply-by-three margins. If normalization is
skipped, the next iteration would square values of roughly 88--89 bits. Its raw
coefficient bound becomes roughly

```text
2^22 * (2^89)^2 = 2^200
```

which requires at least seven 31-bit residue channels, versus three today.
Keeping digit and carry terms separately does not avoid this expansion: the next
square needs digit-square, cross, and carry-square convolutions.

Outcome: delaying normalization across a complete PRP iteration increases
transform width and work by more than the carry pass it removes. Delayed
reductions *within* a radix kernel remain useful range tracking and are already
employed extensively by the current M31/M61 code, but a cross-iteration
carry-save representation is an algorithmic no-go for this exponent.

### 2026-08-09: exact FP32-estimated three-prime CRT

Implemented:

- Added exact FP32-assisted reconstruction to [`src/RNS31.cpp`](src/RNS31.cpp).
  FP32 estimates which multiple of the combined modulus to remove, while the
  96-bit basis products, subtraction, range check, and balanced result remain
  integer-exact.
- Extended [`src/rns31_test.cpp`](src/rns31_test.cpp) with uniform CRT values,
  coefficient-range inputs, boundaries, and exact comparison with Garner CRT.
- Added [`src/cuda/rns31_crt_bench.cu`](src/cuda/rns31_crt_bench.cu) and
  `make CUDA=1 rns31-crt-bench`. The benchmark compares GPU Garner reconstruction
  against the estimator and validates every 96-bit output.

The normalized exact CRT basis sum is in `[0,3)`. A conservative bound on the
three FP32 conversion, reciprocal, multiply, and FMA errors is far below one, so
the estimated floor can differ from the exact CRT multiple by at most one. The
candidate is therefore corrected once against `[0,P)` and is never trusted
without the exact range check.

Correctness results:

```text
host: 0 corrections / 200,000 uniform CRT values
host: 1 correction  / 200,007 coefficient-range and boundary values
GPU:  1 correction  / 8,388,608 target-range values
all estimated reconstructions exactly match Garner
```

With 8,388,608 GPU reconstructions and 31 timing samples:

| Reconstruction | Median | Mean |
|---|---:|---:|
| Integer Garner CRT | 0.161 ms | 0.161--0.162 ms |
| FP32-estimated exact CRT | 0.161 ms | 0.161 ms |

Outcome: the estimator is exact and essentially cost-neutral in this
memory-resident standalone pass, but it is not faster than Garner. It may still
be a convenient quotient primitive in a future RNS carry design, yet it cannot
offset a direct-prime transform that already takes 0.708 ms before CRT/carry.
There is no end-to-end integration gate on this target.

### 2026-08-09: 8M x two-prime range/work trade

Parameterized [`src/cuda/rns_mont_ntt_bench.cu`](src/cuda/rns_mont_ntt_bench.cu)
by transform size and plane count, then ran an exact 8M x two-q31 operation. Both
planes pass dense round-trip and sparse cyclic-square validation.

| 8M x two-prime operation | Median | Mean |
|---|---:|---:|
| Full Montgomery roots | 0.966 ms | 0.965 ms |
| 16x-compressed/generated roots | 1.253 ms | 1.253 ms |

One 8M full-root plane takes 0.444 ms. The two-plane combined modulus has about
61.66 bits, versus the approximately 55.5-bit raw coefficient estimate, but the
operation is already about 4.8 times slower than the complete tuned production
iteration before CRT, carry, weighting, or safety checks. Increasing the
transform to reduce RNS width is therefore a measured no-go.

### 2026-08-09: second architectural campaign

The second campaign started from the first campaign's main lesson: preserve the
incumbent packing, working-set size, and transform geometry, and replace only a
costly arithmetic plane. It also used exact production kernels as an oracle
rather than treating a faster approximate residue as acceptable.

#### FP32 repair and longer-transform controls

An exact `FFT323161` sidecar confirmed that the FP32 quotient-parity mapping was
correct over 8,998 fused iterations: there were zero predicted-versus-exact
parity mismatches, and the 2,000/10,000 iteration residues matched. The failure
was instead information loss in the stored FP32 coefficients:

```text
odd quotient errors:                  31,686.6 coefficient pairs/iteration
wrong +/-1 direction from residual:       87.1 pairs/iteration
```

Using FP64 for the FP32 complex multiply or FMA did not reduce that population;
it only raised the exact triple-path time to about 0.354--0.449 ms. The lost bits
must be retained as state, not recovered by locally more accurate operations.
Residual-window repair was also rejected because wrong-direction cases were
spread across both small and large residuals.

Longer exact transforms did not create a faster low-precision path:

| Architecture | Correctness | Steady time |
|---|---:|---:|
| 8M FP32+M61 | correct | 0.416 ms |
| 8M M61 only | correct | 0.314 ms |
| 8M FP32+M31, best tested geometry | deterministic failure | 0.233 ms |
| 4M FP64+M31 | correct | 0.521 ms |

The 4M FP64+M31 result used `51:512:8:512:202`, produced the expected residues at
2,000 and 10,000 iterations, and had low ROE; its performance rejects it.

#### Coefficient range and the false M23 lead

The exact M31+M61 oracle observed no balanced coefficient outside `[-2^78,2^78)`
over 8,998 iterations. Coefficients at or above `2^77` occurred only about 2.4
pairs per iteration, while values at or above `2^76` occurred about 9,787 times
per iteration. This explains why an 84-bit empirical range can appear safe even
though the conservative all-input convolution bound is about 87 bits before
additional margins.

An arithmetic control initially labeled `M23` was invalid as an NTT proposal:

```text
2^23 - 1 = 8,388,607 = 47 * 178481
```

It is not a prime field and cannot supply the required quadratic-extension
roots. M31 is the smallest actual Mersenne prime whose quadratic extension can
support this construction. The old `2^23-1 + M61` timing remains only as a
throughput control in the CUDA benchmark and must not be cited as a candidate.

#### Representation-preserving three-plane design

The viable substitution is:

```text
current:     GF(M31^2), 8 bytes/pair + GF(M61^2), 16 bytes/pair
replacement: GF(M31^2), 8 bytes/pair + two GF(q^2), 8 bytes/pair each
```

Both use 24 bytes per packed pair, or 48 MiB for the 4M-word transform. Unlike
the failed three-direct-prime transform, the replacement retains the
two-real/conjugate packing and performs the same number of packed transforms.

The first valid pair was:

```text
q0 = 1,769,996,287 = 211 * 2^23 - 1
q1 = 1,887,436,799 = 225 * 2^23 - 1
```

The product with M31 has 92.535 bits. Exact combined arithmetic and complete
radix-16 forward/square/inverse tiles gave:

```text
current M31+M61 tile:     0.053 ms
M31+two-q tile:           0.044--0.045 ms
replacement/current:      0.844--0.850
```

Every plane was validated against direct cyclic convolution. At chain lengths
1 and 2 the replacement was tied with current because traffic dominated. It was
about 19--25% faster from chain length 4 through 64, and the complete tile ratio
remained about 0.842--0.843 across those tests.

#### Riesel-prime search

The two generic primes above are Riesel primes. A broader search showed that
`2^23` divisibility was unnecessarily strong. The production packed transform
has `ND = 2^21`, so a quadratic field with a `2^21` root is sufficient. For an
odd Riesel prime `q = k*2^n-1`, `v2(q^2-1)=n+1`; therefore `n >= 20` is the root
gate. Also, `q == 7 (mod 8)`, so `-1` is nonsquare (preserving the `a+bi` field)
and 2 is a quadratic residue of odd order. Raising to `2^22` is consequently an
automorphism on the subgroup containing 2, so the exact Crandall--Fagin root of
2 exists.

The first high-product sub-`2^31` candidates found were:

| Prime | Riesel form | Pseudo-Mersenne form | `v2(q^2-1)` |
|---:|---:|---:|---:|
| 2,142,240,767 | `2043*2^20-1` | `2^31-(5*2^20+1)` | 21 |
| 2,141,192,191 | `1021*2^21-1` | `2^31-(3*2^21+1)` | 22 |

Deterministic 32-bit primality testing passed. Their product with M31 has
92.9922406 bits, and explicit base-field roots `theta` satisfy
`theta^(2^22) == 2` for both moduli. Generic Montgomery arithmetic has the same
measured cost as the first pair:

```text
M31+two Riesel arithmetic, chain 32:  0.090 ms
current M31+M61 arithmetic:           0.118 ms
Riesel/current:                       0.758

Riesel forward/square/inverse tile:   0.044 ms
current tile:                         0.053 ms
Riesel/current:                       0.843
```

The `2043*2^20-1` candidate was subsequently rejected: it has a root in the
full quadratic multiplicative group but not the required order-`2^21`
norm-one root used by PRPLL's conjugate packing.  The integrated replacement
uses `997*2^21-1` and `1021*2^21-1`; its product with M31 is 92.9572168 bits.

Their low-Hamming pseudo-Mersenne complements suggested multiplication-free
shift/add reduction. The exact four-fold implementation validated, but serial
folding was much worse on Blackwell:

```text
Riesel shift/add arithmetic:  0.300--0.306 ms
Riesel shift/add tile:        0.104--0.105 ms
tile/current:                 1.97--2.00
```

This is an important hardware-specific result: a reduction attractive in custom
hardware is not automatically attractive on a GPU with efficient multiply-high
instructions. Generic Montgomery remains the selected representation.

#### CRT and Crandall--Fagin weighting gate

At 4,194,304 balanced coefficient reconstructions, both paths exactly matched
random signed targets below `2^78`:

| Path | Median |
|---|---:|
| Current M31+M61 CRT | 0.024 ms |
| M31+two-Riesel Garner CRT | 0.026 ms |

A synthetic fused inverse-weight/CRT/forward-weight pass measured 0.130 ms for
current and 0.129 ms for the Riesel replacement. Adding one recurrence multiply
to model generation of each generic inverse and forward weight also measured
0.129 ms. These kernels are traffic-heavy, so equality means only that weight
generation did not move this standalone lower bound; register pressure and
weight recurrence inside the production fused carry still require measurement.
M31/M61 rotations remain cheaper operations in isolation.

No full-size weight table should be used: two tables would add 32 MiB and the
associated traffic. A production candidate should factor starting weights by
thread/group and advance them with two fixed recurrence ratios, or absorb the
scales into otherwise-required transform-edge constants where algebra permits.

#### Other prime shapes

An exhaustive sub-`2^31` Proth-prime search over relevant 2-adicities found no
candidate for which 2 has a `2^22`-th root in either the base field or its
quadratic extension. A Proth direct NTT would also abandon the conjugate packing
and return to the already-rejected direct-prime architecture.

`M31 + Goldilocks`, with `Goldilocks = 2^64-2^32+1`, gives almost exactly 95 bits
and the same 48 MiB data footprint. A scalar-root tile was 16.2% faster, but that
root is fixed by conjugation and does not implement PRPLL's Hermitian packing.
The correct unit-norm root satisfies `conj(root)=root^-1` and requires a general
three-product complex multiply. Its exact results were:

```text
M31+Goldilocks general arithmetic: 0.137 ms (1.161x current)
unit-norm Goldilocks tile:         0.059 ms (1.116x current)
```

Thus Goldilocks is a no-go; the superficially faster scalar-root result is an
example of benchmarking the wrong algebra.

#### Current decision and next gates

The production path remains unchanged because no end-to-end speedup has been
measured. The packed Riesel design is the only surviving architectural lead.
The next gates, in order, are:

1. Construct and validate primitive roots through order `2^21` plus the exact
   Crandall--Fagin weight recurrence for both Riesel fields.
2. Build the smallest complete 4M packed forward/square/inverse pipeline using
   production-like width/middle/height layouts. A generic scalar NTT would again
   benchmark the wrong architecture.
3. Integrate a three-way balanced CRT and generated weights into the fused carry,
   measuring registers, spills, and the complete fused-region time.
4. Require exact 2,000/10,000 iteration residues, then run a sustained A/B long
   enough to compare against the 0.204--0.206 ms thermal baseline.

A reasonable but unproven projection is a 10--15% whole-iteration opportunity:
the exact tile saves about 15.7%, standalone CRT adds about 2 microseconds, and
total residue bytes are unchanged. The projection must not be reported as an
achieved speedup.

## Final verification

The final source tree was checked with:

```text
make -j4 CUDA=1 all rns31-test carry-transfer-test \
  rns31-cuda-bench rns31-ntt-bench rns-multi-ntt-bench \
  rns-tensor-bench rns-mont-ntt-bench rns31-crt-bench
./build-cuda/rns31-test
./build-cuda/carry-transfer-test
fresh CMake Release build followed by ctest --output-on-failure
git diff --check
```

The CUDA application and all six standalone CUDA benchmark targets build for
`sm_120`. Both host reference tests pass; the clean CMake build reports 2/2
CTest tests passed. The final uninstrumented 30k production confirmation for
exponent `136279841` produced residue `9139db3046e846d4` at 198.5 us/iteration,
and the original 1M run produced residue `52b03a7cc55e677d` at 205.5
us/iteration. `git diff --check` reports no whitespace errors.

Supporting workflow changes also remain in place: [`src/main.cpp`](src/main.cpp)
disables stdout buffering with `setvbuf(stdout, nullptr, _IONBF, 0)`, and the
CUDA backend is rebuilt in [`build-cuda/prpll`](build-cuda/prpll).

## Final decision and remaining research boundary

The production configuration still stays on the tuned 4M `M31+M61`
implementation. Integrating any direct RNS, Tensor Core, compressed-root,
multi-stream, Goldilocks, FP32-repair, or delayed-carry prototype would knowingly
reduce throughput or correctness. Those experimental targets remain buildable
as comparison oracles.

The original strategies are exhausted for this GPU under their measured gates.
The packed Riesel substitution is new evidence and is the sole active research
branch. It must pass the complete-transform and fused-carry gates described
above before it is eligible for production integration. Other future work would
require one of:

- Permission to collect Nsight Compute hardware counters. The current driver
  rejects them with `ERR_NVGPUCTRPERM`; counters could identify a previously
  invisible bottleneck inside the existing extension-field kernels.
- A new sparse Tensor Core formulation that performs asymptotically fewer than
  the dense 16x16 limb construction. Merely increasing batch size cannot fix its
  42% tile deficit.
- A different exact convolution formulation whose complete transform is below
  roughly 0.13 ms, leaving useful headroom for CRT/carry under the 0.204 ms
  production budget. The best direct-prime result is currently 0.708 ms; the
  packed Riesel path does not yet have a complete-transform measurement.
- Additional GPUs. Assigning whole unrelated exponents to separate devices is
  embarrassingly parallel and should scale throughput, unlike concurrent streams
  sharing this already power-limited device.

The Riesel branch has now passed its algebraic and end-to-end correctness gates,
but not its performance gate.  The results below supersede the prospective
Riesel plan above.  The highest-throughput demonstrated action therefore remains
one exponent at a time per GPU with the tuned production configuration.

## End-to-end three-prime M31R2 architecture

The packed-Riesel lead was integrated as a proper three-prime architecture,
rather than routing the two Riesel residues through the old packed `M61` data
type.  FFT type `54` (`M31R2`) stores three independent complex `uint32` planes:

- `M31 = 2^31-1`
- `q0 = 997*2^21-1 = 2090860543`
- `q1 = 1021*2^21-1 = 2141192191`

Each plane has its own roots, weights, transforms, and CUDA queue.  The carry
stage reconstructs the signed coefficient with an exact 93-bit three-way CRT.
The complete `54:512:8:512` implementation reproduced the known residues at
iterations 2,000, 10,000, and 1,000,000 for exponent `136279841`; the final
million-iteration residue was `52b03a7cc55e677d`.

With `INPLACE=1,MULTI_Q=1`, the independent-plane long-carry implementation ran
at 302.5 microseconds/iteration over the sustained million-iteration run.  Using
three queues was useful: short measurements improved from about 310.3 to
295.3--295.8 microseconds/iteration.  It nevertheless remained about 24% slower
than the matched current long-carry `M31*M61` result of 243.5--244.0
microseconds/iteration, and much slower than the approximately 205.5
microseconds/iteration tuned production fused path.

### Fused M31R2 carry kernel

A direct analogue of the `M31M61` fused carry was also implemented.  It performs
the three inverse width transforms, exact CRT reconstruction, carry staircase,
forward reweighting, and three forward width transforms in one kernel.  Since
both Riesel fields are simultaneously visible, it uses explicit parameterized
32-bit Montgomery arithmetic rather than pretending the two primes are one
packed `uint64` field.

The first version kept all three transform vectors live.  It was correct, but
register pressure/spilling produced about 425.7 microseconds/iteration.  A
lower-register version reused one vector and staged intermediate residues in
global memory.  That exposed an important correctness hazard: workgroup zero
can overwrite the original first lines before the duplicate wraparound
workgroup reads them.  Dedicated scratch for those first lines removed the race.
The resulting kernel reproduced `6a5c8b8989125413` at iteration 20,000.

Final matched short-run measurements were:

```text
M31R2 split, three queues: 303.0 us/iteration
M31R2 fused carry:         402.6 us/iteration
fused kCarryFused alone:   258.1 us/invocation
```

Thus fusion increased iteration latency by about 33%, or reduced throughput by
about 25%, relative to the split M31R2 path.  It removes intermediate kernel
boundaries but serializes work from three independently schedulable residue
fields.  Keeping all fields live exceeds the attractive register footprint;
reducing that footprint reintroduces global spill/reload traffic.  On this GPU,
those two costs are much larger than the launch and intermediate-traffic savings
that made the two-field `M31M61` fused kernel successful.

The fused implementation remains available for experiments with
`-use RIESEL_FUSED=1`.  M31R2 defaults to the faster split carry.  This result is
also a general lesson for future RNS designs: fusion must be evaluated against
field-level concurrency and live-state size, not only against eliminated global
memory passes and kernel launches.

## Follow-up M31R2 layout, scheduling, and SIMD experiments

The independent-plane implementation was tested across the remaining obvious
layout controls.  All listed runs reproduced the known checkpoint residue.

| Experiment | Result |
|---|---:|
| CUDA graphs enabled | 298.9 us/iteration |
| Matched graphs-disabled control | 300.3 us/iteration |
| `54:512:8:512` | approximately 295--300 us/iteration |
| `54:1024:4:512` | 311.3 us/iteration |
| `54:512:4:1024` | 319.2 us/iteration |
| `54:256:16:512` | 307.9 us/iteration |
| `54:1024:8:256` | 325.9 us/iteration |
| `54:512:16:256` | 317.6 us/iteration |
| `54:256:8:1024` | 330.9 us/iteration |

FFT variants `202`, `112`, `211`, and `210` measured 299.2, 300.1, 301.4,
and 302.2 us/iteration respectively in short runs; the existing `212`/`202`
choices remained within the same noise band.  Long-carry block lengths 4, 8,
16, and 32 similarly gave 299.6, approximately 299.2, 299.2, and 301.4
us/iteration.  These controls do not expose a missing layout win.

### Packed-within-thread q0/q1 SIMD

An opt-in `RIESEL_PACKED=1` path put q0 and q1 into the low and high halves of
one `uint64`, retained independent Montgomery arithmetic for the two moduli,
and reduced the three transform queues to an M31 stream plus a paired-q stream.
It is exact at iterations 2,000 and 10,000, but measured 407.9 us/iteration.
The paired tail alone took approximately 148.5 us.

Blackwell's scalar CUDA integer cores do not execute two independent 32-bit
modular multiplications as a packed `uint64` SIMD operation.  Pairing merely
doubles each thread's dependency chain and register footprint.  This is direct
end-to-end evidence that any useful field-level SIMD must assign fields to
different warps/blocks or GPUs, not pack them into one scalar thread.

## Generated PTX/SASS and the M31-as-Riesel question

M31 is indeed the special Riesel prime `1*2^31-1`, and the production code
already exploits the important consequence.  Crandall--Fagin weights and NTT
constants that are powers of two become cyclic shifts/folds modulo M31.  The
proper Riesel fields still need general Montgomery products for those constants.
Treating M31 as a generic member of a uniform three-prime loop would therefore
discard the most valuable property of the prime.

Disassembly of matched 4M kernels on the RTX PRO 6000 Blackwell Max-Q made the
difference explicit:

| Kernel | Registers | SASS instructions | `IMAD` | `ISETP` | `VIMNMX` |
|---|---:|---:|---:|---:|---:|
| M31 tail | 62 | about 2,450 | 238 | 11 | 0 |
| q0 tail, portable Montgomery C | 61 | about 2,960 | 533 | 258 | 395 |
| q0 tail, inline-PTX Montgomery | 60--63 | about 2,765 | 532 | 265 | 395 |
| M61 tail | 96 | about 3,845 | 614 | 26 | 0 |

NVVM's portable-C Montgomery sequence initially looked inefficient in PTX: it
contained two `mul.wide.u32` operations, a mask, and a 64-bit add.  Final SASS
was better: `ptxas` recognized most instances and formed an `IMAD.WIDE`, a low
`IMAD`, and an `IMAD.HI`.  Inline PTX spelling the low/high carry chain did not
remove a multiply, but it prevented some conservative 64-bit carry bookkeeping
and reduced the static q-tail instruction count by about 6.5%.

The `RIESEL_PTX_MONT=1` path reproduced all known residues, including the
million-iteration residue `52b03a7cc55e677d`.  Short A/B runs initially favored
it by 2--7 us/iteration, but the sustained million-iteration result was 302.1
us/iteration versus 302.5 us/iteration for the existing portable-C run.  That
0.1% difference is not a demonstrated speedup.  The option remains useful for
compiler experiments, but inline PTX does not address the architectural gap.

## Why production M31/M61 is about 100 us faster

The original expectation was that replacing a 64-bit M61 plane with two
32-bit planes would expose more integer throughput while keeping the same
48 MiB residue footprint.  That model omitted four production advantages of
the M31/M61 path:

1. Both moduli are Mersenne primes.  Multiplication, reduction, roots, inverse
   weights, and forward weights use shifts, rotations, adds, and delayed folds.
   The q fields pay Montgomery multiplication for arbitrary roots and for every
   weight recurrence.
2. The mature M31/M61 butterflies track redundant value ranges and postpone
   normalization.  Their M61 tail has only 26 predicate instructions and no
   `VIMNMX`; one q tail has roughly 265 predicates and 395 `VIMNMX` reductions.
3. The two-modulus CRT is unusually cheap.  After decoding M31, it forms
   `u61-n31`, adds a 31-bit shift, performs one M61 fold, and constructs the
   signed 96-bit value.  The three-prime Garner path decodes and inverse-weights
   two generic residues, performs another modular multiply, two wide products,
   and a 96-bit balanced-range comparison.
4. The fused M31/M61 kernel can keep both field vectors live while doing inverse
   width transforms, CRT/carry, reweighting, and forward width transforms.  Its
   generated kernel has 4,640 SASS instructions, 96 registers, and no local
   stack.  The reduced-register M31R2 analogue has 7,982 instructions, 118
   registers, a 272-byte local stack, and additional explicit global staging.

The measured approximately 97 us end-to-end gap separates cleanly into two
parts.  Like-for-like long carry is 302.5 us for M31R2 versus 243.5--244.0 us
for M31/M61, so approximately 59 us comes from generic-prime transforms,
normalization, weights, and three-way CRT.  Production M31/M61 then saves about
another 38 us through its successful fusion and reaches 205.5 us.  M31R2 cannot
copy that fusion mechanically: its measured fused path is 402.6 us because it
loses inter-field concurrency and crosses the register/local-memory threshold.

Replacing M31 with another proper Riesel prime is therefore contraindicated:
it would turn the cheapest plane (238 tail `IMAD`s and almost no normalization)
into another generic plane (about 532 tail `IMAD`s plus hundreds of reductions).
A uniform block-level scheduler may still retain specialized M31 blocks, but
uniform arithmetic is the wrong objective.

## New prime-shape lead: lazy ranges below 2^30

The disassembly suggests a more relevant use for different Riesel primes.  If
`q < 2^30`, a Harvey butterfly can keep residues in `[0,4q)` in one `uint32`.
Montgomery multiplication can consume the lazy input and return a bounded lazy
result, allowing reductions to be deferred across butterflies.  The selected
near-`2^31` primes cannot do this because `4q` overflows 32 bits.

Two verified candidates are:

| Prime | Riesel form | Norm-one root gate | root-of-two gate |
|---:|---:|---:|---:|
| 1,038,090,239 | `495*2^21-1` | pass | `theta^(2^22)=2` |
| 1,031,798,783 | `123*2^23-1` | pass | `theta^(2^22)=2` |

Together with M31 they provide 90.8937992 bits of CRT range.  This is smaller
than production's 92 bits, but balanced 32/33-bit digits give a conservative
target-exponent square bound below `2 * NWORDS * 2^64`; even after multiplication
by three this is below 88.585 bits before the much smaller incoming carry.  Half
the candidate CRT product is 89.894 bits, leaving about 1.3 bits of conservative
headroom for exponent `136279841`.  This bound must be formalized against every
carry/input invariant before integration.

This is an architectural rather than cosmetic prime substitution: the go/no-go
test is whether a range-annotated lazy radix tile materially reduces final SASS
predicates and elapsed time.  Simply substituting the smaller primes while
retaining canonical reduction after every operation is not expected to help.

## Forum and literature survey: fusion boundaries and modulus choice

The external survey found no published GPU implementation of the complete
PRPLL-specific chain

```text
inverse NTT edge -> exact multi-prime CRT -> cross-coefficient carry -> forward NTT edge
```

in one kernel.  The closest exact-integer GPU work performs transforms per
prime and then a separate CRT kernel; its carry discussion is prospective and
its measurements omit digit adjustment and carry ([Emeliyanenko 2009](https://workasm.github.io/pdf/poly_mul_2009.pdf)).
This absence matters because CRT is elementwise, while PRPLL carry has a
cross-coefficient dependency that the usual FHE fusion papers do not face.

Published GPU NTT fusion results support a narrower boundary.  They fuse final
inverse scaling or other elementwise work while transform values are already in
registers, and fuse the first forward stage with adjacent elementwise work.
They do not retain several complete residue-field tiles across CRT.  Jung et
al. eliminate memory passes by applying scaling immediately before store
([TCHES 2021](https://d-nb.info/1246458438/34)); Shivdikar et al. fuse the final
CT stage, Hadamard product, and first GS stage
([SEED 2022](https://arxiv.org/abs/2209.01290)).  Zhai et al. also report that
radix 16 loses to radix 8 once register spilling exceeds the saved traffic
([IPDPS 2022](https://arxiv.org/pdf/2109.14704)).  This is consistent with the
local 118-register, 272-byte-stack M31R2 fused kernel.  The resulting design
rule is:

- Keep field transforms independent.
- Normalize and inverse-weight at the final inverse edge if this does not
  increase its live tile.
- Perform CRT/carry scalar by scalar.
- Apply forward weighting and the first forward edge immediately before
  storing each field again.

### Why Riesel form alone did not make q0/q1 cheap

Harvey's redundant-butterfly paper records earlier use of redundant
representation for a transform over `GF((2^61-1)^2)`
([Harvey 2014](https://arxiv.org/pdf/1205.2926)).  Together with the
hierarchical DGT analysis of Alves, Ortiz, and Aranha
([DGT paper](https://eprint.iacr.org/2020/861.pdf)), this explains the local
result.  A Riesel field `k*2^s-1` supplies the required large two-adic subgroup,
but only `k=1` also supplies Mersenne folds, power-of-two roots/weights, and
wide redundant ranges with almost no normalization.  M31 is already the best
case Riesel prime.  Percival's DWT for a Riesel-shaped *target modulus* is a
different benefit and does not make an internal Riesel NTT field cheap
([Percival 2003](https://www.daemonology.net/papers/fft.pdf)).

This distinction was missing from the initial M31R2 projection.  That
projection counted nominal 32-bit arithmetic and bytes, but did not price the
generic roots, two Montgomery weight recurrences, normalization predicates,
and loss of the production carry fusion.

### Forum-derived shorter-transform architectures

Yves Gallot proposed combining special fields so shorter transforms fit cache,
including `Goldilocks*M61*M31` at roughly 156 bits
([proposal](https://www.mersenneforum.org/node/1110517?p=1110569#post1110569),
[cache rationale](https://www.mersenneforum.org/node/1110517?p=1111489#post1111489)).
He also proposed a 3M `M31*M61*(7*2^29-1)` system
([post](https://www.mersenneforum.org/node/1110517?p=1112086#post1112086))
and listed alternatives such as `11*2^26-1 = 738197503`
([prime list](https://www.mersenneforum.org/node/1110517?p=1112082#post1112082)).
The 3M design is sound but has limited upside: four 32-bit limbs at 3M retain
the current 48 MiB payload and have about `0.98x` the current first-order
limb-stage work.  Prime95 independently estimated about a 2% improvement for
an analogous mode
([analysis](https://www.mersenneforum.org/node/1110517?p=1113714#post1113714)).

The forum also warns that a standalone odd-radix pass increases global-memory
traffic; a radix-3 edge should be embedded in an existing kernel
([radix-3 discussion](https://www.mersenneforum.org/node/1110517?p=1112564#post1112564)).
Prime95 keeps transforms separate by field and combines fields only at carry
([implementation note](https://www.mersenneforum.org/node/1110517?p=1112740#post1112740)).
Corrected Good-Thomas/radix-9 comparisons and measured radix-7 attempts likewise
did not beat the proper comparator
([Good-Thomas result](https://www.mersenneforum.org/node/1110517?p=1120952#post1120952),
[radix-7 result](https://www.mersenneforum.org/node/1110517?p=1112489#post1112489)).

#### Reopened local gate: exact 3M M31+M61+one-q architecture

The 3M proposal above had only a forum/first-order work estimate locally; it
was never subjected to the same exact architecture-sized tile gate as M31R2,
Gold, and the prime-power designs.  It is materially different from the failed
M31R2 path: retain both production fields, add only
`q=11*2^26-1=738197503`, and shorten all populations from 2M to 1.5M quadratic
values.  The product has about 121.46 bits; the conservative p150 coefficient
bound of about `2^118.17` appears to leave over two bits below its balanced CRT
limit.  The q field is below `2^30` and permits Harvey lazy ranges.  Those
properties made it a sensible arithmetic gate, but they do not prove that the
complete weighted Mersenne convolution exists.

[`src/cuda/m31_m61_q3m_tile_bench.cu`](src/cuda/m31_m61_q3m_tile_bench.cu)
and `make m31-m61-q3m-tile-bench` implement a real forward radix-16,
pointwise quadratic square, and inverse radix-16 in all three proposed fields.
The 3M candidate is compared with the current M31+M61 tile at its 4M
population.  An independent direct cyclic convolution validates all three
candidate outputs, including conversion into and out of q Montgomery form.
The candidate uses 40 registers and has no stack or spills.

```text
current 4M M31+M61 tile:   109.088 us
candidate 3M three-field:  103.712 us
candidate/current:           0.9507
```

This passed the arithmetic tile gate, but not the later whole-algorithm root and
weight gate described below.  Its
fixed-radix component saves only 4.9%; the complete transform also replaces 21
power-of-two stages by 19 such stages plus a Good--Thomas radix 3.  Conversely,
the exposed M61 population becomes 25% shorter and the M31/q work may be
scheduled separately.

[`src/cuda/m31_m61_q_radix3_bench.cu`](src/cuda/m31_m61_q_radix3_bench.cu)
implements the previously missing exact radix-3 gate.  It uses
`s=b+c`, `base=a-s/2`, and `delta=(w-w^2)(b-c)/2`, avoiding a general division
by three; the omitted inverse normalization can be folded into the eventual
CRT/carry scale.  An independent direct cyclic square in each quadratic field
validates all outputs.  The all-field kernel uses 40 registers, the M61-only
kernel 36, and the M31+q kernel 32, with no stack or spills.

```text
3M exact Good-Thomas radix-3/square/inverse PASS
all fields: 23.90--24.10 us
M61 only:   13.70--13.82 us
M31+q:      13.70--13.82 us
```

This is a favorable fused lower bound: it includes one coalesced read and write
of the three channels and the pointwise square.  A separate global radix-3 pass
would add traffic and is rejected; an implementation must absorb the edge into
the adjacent power-of-two tile/square kernels.  Field-separated scheduling is
preferred because keeping all three fields live saves less than 4 us while
raising the live state.

The q modulus also permits Harvey redundant Montgomery values in `[0,2q)`.
The same architecture-sized tile now contains a complete lazy q alternative,
with exact direct-convolution validation:

```text
ordinary 3M three-field tile: 112.06--113.66 us
Harvey-lazy q tile:           115.68--117.54 us
lazy/ordinary:                  1.027--1.034
```

Both variants use 40 registers without spills.  Although the lazy butterfly
removes canonical correction steps, it carries wider ranges through the
quadratic Karatsuba products and does not reduce M31/M61 work or memory traffic.
On this GPU that trade is negative.  **Do not repeat the lazy q tile for this
modulus.**

#### Exact 121-bit three-field CRT/carry gate

The registry contained exact carry implementations for the 93-bit
M31+two-Riesel path and the 154-bit four-field path, but not for this 121-bit
combination.  [`src/cuda/m31_m61_q3m_carry_bench.cu`](src/cuda/m31_m61_q3m_carry_bench.cu)
and `make m31-m61-q3m-carry-bench` add the missing production-sized paired
gate.  It reconstructs in the favorable order M31*q first: that product fits
below M61, so the final Garner digit uses the existing cheap M61 fold.  The
centered coefficient is carried in two 64-bit limbs, multiplied by three, and
propagated with the exact p136 mixture of 43- and 44-bit words.  Its roughly
80-bit carry is retained exactly rather than truncated to 64 bits.

The kernel also includes the field-boundary work omitted by a bare CRT timing:
q is read in Montgomery form, and each centered output digit is forward
weighted into M31, M61, and q Montgomery residues.  The paired current kernel
does the analogous M31/M61 reconstruction, 32/33-bit propagation, and forward
weighting.  Both sides therefore read and write 48 MiB of field state.  An
independent signed `__int128` reference validates every one of the 3,145,728
candidate digits and all 393,216 terminal chain carries over random values
across the complete balanced CRT range.

```text
3M M31+M61+q exact CRT/carry/forward-weight validation PASS
current 4M carry       median 28.160 us, mean 28.295 us, min 27.200 us
candidate 3M carry     median 32.320 us, mean 32.649 us, min 31.808 us
candidate/current ratio 1.148
```

Both kernels use 40 registers with no stack or spills.  The candidate pays only
about 4.2 us over the equal-payload current control: the 25% smaller population
largely compensates for the extra q normalization/output, 128-bit coefficient,
and wider carry.  This passes the arithmetic carry gate for the measured
modulus, but that modulus fails the weight-existence gate below.  It does not include the
small cross-chain carry correction or adjacent inverse/forward transform edges,
so it is not an end-to-end result; those costs must be included in the complete
standalone recurrence.  Unlike the rejected 2M four-field prototype, this gate
uses a real three-field CRT and exact wide carry rather than a placeholder.

#### Corrected whole-algorithm gate: the first q has no 3M weight

The prior gates checked primality, CRT range, transform roots, radix 3, and
carry, but omitted the final algebraic condition for a weighted Mersenne
convolution.  For `N=3*2^20` real words, each field must contain `theta` with
`theta^N=2`.  In a finite field this exists only if

```text
2^((q^2-1)/gcd(N,q^2-1)) = 1 (mod q).
```

For `q=738197503`, the test fails; equivalently, `2^((q-1)/3)` is
`665297573`, not one.  Thus 2 is not a cube in the relevant quadratic group,
and no arrangement of the otherwise valid radix-3 and power-of-two roots can
implement the required Crandall--Fagin weighting.  **Reject this q for the 3M
engine.**  This is an algebraic failure, so the favorable component timings
cannot rescue it.

A search over prime `k*2^s-1` fields with `s>=20` found a compatible existing
PRPLL modulus:

```text
q = 2141192191 = 1021*2^21 - 1
log2(M31*M61*q) = 122.995767
balanced CRT bits = 121.995767
theta = 756791506
theta^(3*2^20) mod q = 2
```

The conservative p150 coefficient estimate `2^118.17` leaves about 3.8 bits
below its balanced CRT limit.  Because `q == 1 (mod 3)` it also has a scalar
radix-3 root, avoiding an extension-field constant multiply in that edge.
This q is already one lane of the validated M31R2 implementation, so its
Montgomery arithmetic and order-`2^21` quadratic root are not new assumptions.
What has not been attempted is the one-q, 3M M31+M61 architecture.  The exact
tile, radix-3, and 121-bit carry gates must therefore be repeated with this
modulus before a production-shaped complete transform is justified.  A
standalone engine remains acceptable, but must validate the exact recurrence
and sustain at most 180 us for the full one-million-iteration run.

The three standalone gates are now parameterized by `Q3_MODULUS`, preserving
the rejected small-q binaries while building the admissible-q variants with:

```text
make m31-m61-q3m-viable-tile-bench
make m31-m61-q-viable-radix3-bench
make m31-m61-q3m-viable-carry-bench
```

All independent direct-convolution and full-range carry checks pass for
`q=2141192191`.  Three representative paired runs give:

```text
radix-16 current 4M M31+M61: 107.584--107.680 us
radix-16 viable 3M 3-field:  103.328--103.360 us  (ratio 0.9600--0.9604)

radix-3/square/inverse all:   23.968--24.288 us
radix-3 M61 path:             13.600--13.728 us
radix-3 M31+q path:           13.664--13.760 us

current 4M CRT/carry/weight:  27.936--28.128 us
viable 3M CRT/carry/weight:   33.184--33.632 us  (ratio 1.188--1.204)
```

The viable-q tile and carry kernels each use 40 registers without stack or
spills; the radix-3 resource counts are unchanged.  The larger q gives up the
Harvey range and reduces the fixed-tile lead from 4.9% to about 4.0%, while the
extra Mersenne fold in CRT raises the carry premium from about 4.2 to 5.3 us.
It nevertheless passes the corrected component gates: the exposed M61
population is still 25% shorter, the M31+q path is separately schedulable, and
the real wide carry is not by itself prohibitive.

An additional algebra audit rejected the initially proposed implementation of
that next gate before production code was changed.  Multiplying PRPLL's
power-of-two quadratic root by a scalar cube root does create an element of
order `3*2^k`, but the scalar cube root is fixed by quadratic conjugation and
is not norm-one.  Consequently `conjugate(root) != inverse(root)`, violating
the packing invariant used by the existing GF tail kernels.  A literal
`512:6:512` root-table extension would repeat the wrong-algebra mistake from
the first scalar-Goldilocks experiment.  **Do not implement a monolithic mixed
GF root or treat the existing radix-3 microbenchmark as proof of that layout.**

The valid mixed-radix realization is a Good--Thomas channel decomposition.
Write `N=3B`, `B=2^20`, map each natural coefficient index through the CRT
isomorphism `Z/N -> Z/3 x Z/B`, and run three independent packed, norm-one
power-of-two DGTs of `B` real words per field.  Apply the scalar radix 3 across
the three channel spectra, square pointwise, apply the inverse radix 3, and run
the three inverse DGTs.  Each radix-3 output channel still has base-field time
coefficients, so its power-of-two spectrum retains the Hermitian invariant
required by the mature tail kernel.  Natural-order weights and carry are
gathered/scattered through the same CRT index map.

This is not a renamed generic NTT: it retains PRPLL's packed M31/M61 DGT and
its field-local power-of-two roots, while exposing an explicit three-channel
batch dimension that a standalone CUDA engine may schedule directly.  The
next early-abandonment gate is therefore the measured cost of three
production-shaped `512:2:512` transform channels plus the already measured
radix-3 edge, not a naïve GF radix-6 kernel.  Only if that critical path leaves
credible room for the exact 121-bit joint carry should the q field and full
Good--Thomas permutation be integrated.

That corrected algebra and scheduling gate now passes.  The extended
[`src/riesel_algebra_test.cpp`](src/riesel_algebra_test.cpp) uses an independent
small `3*16` construction with the exact Good--Thomas CRT input/output map.  It
runs three norm-one quadratic power-of-two transforms, scalar radix 3,
pointwise square, inverse radix 3, and inverse power-of-two transforms.  The
result matches direct natural-order cyclic convolution at every coefficient,
and each post-radix channel independently retains Hermitian symmetry.  The
test is applied only when `3 | q-1`; this correctly admits the selected
`q=2141192191` and avoids treating the old lazy fields as 3M candidates merely
because they support the 4M power-of-two DGT.

The radix-3 CUDA benchmark now includes equal-payload square-only controls.
Across three fresh processes with the viable q:

```text
M61 radix3/square/inverse:       11.84--12.19 us
M61 square-only:                 11.68--11.81 us
incremental radix-3 edge:         0.13--0.38 us

M31+q radix3/square/inverse:     12.03--12.93 us
M31+q square-only:               11.68--12.42 us
incremental radix-3 edge:     about 0.2--0.6 us
```

The edge is bandwidth-hidden in this standalone gate; its arithmetic is not a
reason to reject the design.  It still must be embedded at the tail boundary,
because a separate pass would retain the measured 24-us all-field traffic.

Actual production `FFT3161 512:2:512` kernels were then profiled at the safe
exponent `30000001`.  With the normal two-queue schedule a complete independent
1M recurrence measured 67.4 us/iteration.  A single-queue profile separated
one channel's bottom half as follows:

```text
M61 middle-in + tail + middle-out:  9.8 + 14.3 + 9.7 = 33.8 us
M31 middle-in + tail + middle-out:  5.7 +  9.8 + 5.8 = 21.3 us
existing fused width/CRT/carry:                         20.1 us
single-queue complete iteration:                       81.3 us
```

Finally, three unrelated 1M `FFT3161` workers were run concurrently for 100k
iterations each.  All three completed with valid residues at 144.3, 144.6,
and 144.8 us per worker iteration.  For the proposed channel engine the
relevant batch latency is therefore about 144.8 us, not `3*67.4 = 202.2 us`.
This includes three redundant two-field carries; the real design substitutes
one previously validated 33.2--33.6-us 121-bit joint carry and adds one q
transform.  The remaining margin is narrow but credible enough that this gate
does **not** reject an end-to-end prototype.

An attempted `FFT31R2 512:2:512` profile is not used as q timing evidence.  It
repeated a deterministic ROE/check failure at iteration 2,000 even at safe
19- and 29-bit-per-word exponents; that older end-to-end mode was validated at
its 4M production geometry, not this 1M shape.  Nsight establishes aggregate
kernel cost but not a trustworthy recurrence.  The next implementation must
therefore make the three-channel layout and the single admissible q explicit,
rather than inferring correctness from the two-q mode.

#### Good--Thomas carry-layout gate: direct scatter fails, shared transpose passes

The experiment registry was checked before this implementation.  The existing
3M carry benchmark used natural-order scalar arrays and therefore omitted the
regular but nontrivial permutation between a natural coefficient `n` and its
three-channel Good--Thomas coordinates

```text
channel = n mod 3
base    = n mod 2^20.
```

No earlier experiment had measured this permutation for the admissible 3M
M31+M61+q path.  The prior Gold frequency-reversal bug is relevant as a warning
about factored indices, but it is a different transform and permutation.

[`src/cuda/m31_m61_q3m_carry_bench.cu`](src/cuda/m31_m61_q3m_carry_bench.cu)
now validates and times three storage choices with `q=2141192191`:

1. natural-order storage, the previous optimistic carry lower bound;
2. three globally planar transform channels, accessed directly through the
   CRT index map;
3. 512-word tile-interleaved channel planes, also accessed directly.

Both direct mapped layouts are exact, but they are performance no-gos.  A warp
visits all three planes, expanding memory transactions; putting the planes in
adjacent 512-word tiles does not change that cost.  Across fresh processes:

```text
natural-order exact carry:       32.99--33.98 us
direct channel-planar carry:     56.67--56.80 us
direct 512-tile planar carry:    56.70--56.86 us
direct/natural ratio:             1.67--1.72
```

This rejects direct Good--Thomas gather/scatter in the CRT thread.  **Do not
repeat it with a different plane distance or merely replace `% 3` by a lookup:**
the nearly identical global and tile-planar results identify coalescing, not
integer index arithmetic or long-distance locality, as the loss.

A fourth kernel implements the required memory architecture instead.  A
256-word base-index block coalesces all three fields from all three channel
planes into 12 KiB of shared memory.  Ninety-six workers then process three
sets of contiguous eight-word natural carry segments; across the three natural
`2^20`-word slabs those values consume every channel entry exactly once.  The
block finally writes each transform channel coalesced.  An independent signed
`__int128` oracle validates every digit and all 393,216 terminal segment
carries.  The kernel uses 40 registers, 12,288 bytes of shared memory, one
barrier resource, and no stack or spills.

```text
shared-transpose GT carry:       38.27--38.56 us
shared/natural ratio:             1.14--1.16
absolute permutation premium:     4.5--5.4 us
```

This passes the carry-layout gate.  It also fixes the standalone architecture:
keep field data in transform-friendly channel planes, and fuse a three-channel
shared transpose into the inverse-edge/CRT/carry/forward-edge boundary.  Do not
materialize a natural-order 48-MiB intermediate and do not issue direct
permuted global loads.  The benchmark still omits the small cross-segment carry
correction and the actual width transforms, as did its natural-order control;
both must be included in the complete recurrence.  Thus this is not an
end-to-end speedup, but the measured permutation cost is small enough that the
3M standalone route remains viable.

#### Full-kernel prototype selection and transform-core bring-up

The user explicitly permits a new CUDA program rather than a PRPLL-compatible
implementation.  The implementation audit nevertheless favors reusing PRPLL
as a kernel library: its field-separated power-of-two DGT, transposes, tail
square, queue scheduler, PRP state, and Gerbicz checks are already the hard
parts a clean-room driver would have to reproduce.  The new representation is
therefore opt-in and is not constrained to the old flat carry abstraction.
This choice is about shortening the route to an exact whole recurrence, not
preserving the existing API; a separate executable remains acceptable if the
eventual shared boundary cannot fit the current scheduler.

`GOOD_THOMAS3=1` now brings up the production-sized transform core for explicit
`FFT31R2 512:6:512`:

- middle slots are interpreted as `slot = 2*channel + binary_middle`, with
  three independent `512:2:512` packed DGTs;
- the middle-in kernel applies three radix-2 butterflies and the forward
  scalar radix 3, while middle-out applies the inverse scalar radix 3 and then
  the three radix-2 butterflies;
- the inverse radix 3 is deliberately unnormalised.  Its factor of three will
  be reconstructed exactly and divided once by the joint CRT/carry, avoiding a
  generic inverse-three multiply in every field;
- middle and tail trig tables use the power-of-two `middle=2` roots and are
  padded into the existing `middle=6` allocations, rather than constructing
  the invalid non-norm-one mixed root;
- Hermitian tail ordinals are mapped to three independent 1,024-line spectra,
  so mates never cross a channel boundary;
- only M31, M61, and the admissible q1 field are scheduled.  The old q0 buffer
  is still allocated by this first prototype but its transform kernels and
  traffic are skipped.

The runtime CUDA compiler accepts the complete square kernel set.  Generated
resource use is:

| Kernel | M31 | M61 | q1 |
|---|---:|---:|---:|
| middle in | 37 regs | 48 regs | 38 regs |
| middle out | 37 regs | 60 regs | 38 regs |
| tail square | 62 regs, no local | 96 regs, 24 B local | 62 regs, no local |

The core reached the 2,000-iteration checkpoint, proving that all required
kernels launch at the 3M geometry.  It is **not yet exact**: initial
premultiplication and carry still use the old flat/four-field FFT31R2 mapping,
q0 is intentionally absent, and the on-load residue is consequently
`10410840686ccc5c` instead of `3`.  A `GOOD_THOMAS3`-only on-load check bypass
is present solely to permit transform bring-up and must be removed before any
correctness or end-to-end claim.  The deterministic incorrect checkpoint was
`1e155a6a6d64f851`.  After compilation, the incompatible scaffold measured
about 254.8 us/iteration with its legacy long carry; that number is diagnostic,
not a candidate timing.

The next early-abandonment check measured the overlapped transform critical
path rather than extrapolating the individual kernel sums.  This check was not
present elsewhere in the experiment registry.  An Nsight Systems trace groups
each steady square from the first `fftMiddleIn` among the active field streams
through the last corresponding `fftMiddleOut`; checkpoint `tailMul` cycles and
compile/startup gaps are excluded.  The existing production trace is evaluated
by the identical script and is therefore the relevant control:

| Path | steady cycles | transform-core median | p10--p90 | iteration median |
|---|---:|---:|---:|---:|
| production 4M M31+M61 | 12,852 | 129.632 us | 118.016--130.688 us | 205.664 us |
| provisional 3M M31+M61+q GT | 5,987 | 144.000 us | 130.176--147.104 us | 262.944 us |

The candidate core's additional field and radix-3 work cost only about 14.4 us
over the mature two-field core when the three field queues overlap.  That is a
useful scheduling result, but it fails the complete-engine budget.  Adding the
already validated shared-transpose CRT/carry lower bound of 38.27--38.56 us to
the 144.00-us transform core gives **182.27--182.56 us before either inverse or
forward width edge is performed**.  The actual edge must also execute both
width transforms, the small cross-segment carry correction, and synchronization.
The legacy boundary consumes the remaining 118.5 us in the traced scaffold;
fusion can remove its intermediate traffic, but cannot make those omitted
operations free.

Decision: **the current production-kernel 3M Good--Thomas integration fails the
180-us lower-bound gate, so do not implement its exact legacy-style carry or
mistake a standalone host driver for the missing speedup.**  A standalone CUDA
program remains fully permitted, but it must change the transform/boundary
architecture enough to lower the 144-us core or overlap a demonstrably faster
complete edge; merely copying these kernels into a new executable repeats the
same failed bound.  The temporary `GOOD_THOMAS3` on-load residue bypass remains
prototype-only and cannot support a correctness claim.

Do not repeat a monolithic GF radix-6 root implementation or try to repair this
checkpoint by changing root orientation.  The remaining known boundary is the
natural-word/Good--Thomas shared transpose plus exact three-field CRT/carry,
whose standalone form has already passed above.

An exhaustive 32-bit Riesel search checked whether the failed budget could be
rescued by a different third field.  For every `q=k*2^s-1`, `s>=21`, below
`2^31` and above the capacity floor, it tested deterministic primality and the
exact 3M weight condition

```text
2^((q^2-1)/gcd(3*2^20,q^2-1)) == 1 (mod q).
```

There are 128 compatible `(k,s)` representations (including duplicates of the
same prime), so algebraic existence is not scarce.  For example,
`q=870318079=415*2^21-1` is independently prime, has a scalar cube root of 2,
passes the 3M weight gate, and gives 120.697 balanced CRT bits, about 2.53 bits
over the conservative p150 bound.  It is also below `2^30`.

This does **not** reopen implementation: the registry already contains the
architecture-sized sub-`2^30` Harvey/canonical q comparison.  Lazy arithmetic
made the three-field tile 2.7--3.4% slower, while the complete M31R2 use of the
same arithmetic class improved only 0.6%.  All compatible alternatives retain
the same 32-bit plane, generic Montgomery root multiplies, and carry width; the
near-`2^31` pseudo-Mersenne fold was also already about 2x slower than
Montgomery.  A new prime value without a new reduction/dataflow would therefore
repeat those measured experiments rather than remove the missing width-edge
budget.

#### Clean-room single-field M127 gate

The source and experiment registry were checked before implementation.  M89,
M31 times M89, M31 prime powers, Goldilocks, and multi-prime 32-bit fields had
all been tested, but no single `M127=2^127-1` field was present.  It is a
materially different design:

- use one `GF(M127^2)` field at `3*2^20` real words;
- retain norm-one power-of-two packing and a Good--Thomas scalar radix 3;
- obtain 126 balanced bits, comfortably covering the conservative p150 bound;
- store one 32-byte quadratic value for each of 1.5M pairs, exactly 48 MiB; and
- eliminate multi-field CRT at the carry boundary.

Its risk is wide arithmetic.  [`src/cuda/m127_limb_bench.cu`](src/cuda/m127_limb_bench.cu)
implements exact two-limb M127 addition, subtraction, four-product schoolbook
multiplication, Mersenne folding at bit 127, and quadratic Karatsuba.  Device
results after every timed chain agree with an independent bit-serial host
multiplier.  The M127 kernel uses 48 registers, no stack, and no spills; the
M61 control uses 32 registers.  At the actual architecture populations of
1,572,864 M127 versus 2,097,152 M61 quadratic values:

| chained quadratic multiplies | M61 | M127 | architecture-sized M127/M61 |
|---:|---:|---:|---:|
| 1 | 15.584 us | 28.224 us | 1.811x |
| 2 | 18.688 us | 48.608 us | 2.601x |
| 4 | 28.352 us | 87.072 us | 3.071x |
| 8 | 50.368 us | 161.664 us | 3.210x |
| 16 | 91.520 us | 313.792 us | 3.429x |
| 32 | 175.424 us | 617.920 us | 3.522x |

Three fresh 31-sample chain-32 processes reproduced the result at
`3.506--3.534x` architecture-sized M127/M61, so the rejection is not based on
a single clock or warm-up state.

The implementation uses twelve 64-bit wide products per quadratic multiply.
An ideal two-limb scalar Karatsuba rewrite could reduce that to nine, only a
25% reduction in the dominant products; applying that impossible best case to
the dense result still leaves about 2.64x M61 at architecture size.  A single
M127 transform therefore cannot approach the complete 180-us budget, much less
beat the field-overlapped M31/M61 production core.

Decision: **reject the single-field M127 architecture at the arithmetic gate.**
Do not implement roots, radix 3, carry, or a standalone driver unless a new
M127 multiplication primitive first improves this gate by well over 2x.

#### 3M FP32+M31+M61 Good--Thomas hybrid (active)

The experiment registry was checked before implementation.  The prior
FP32-estimated exact CRT architecture used the ordinary 4M power-of-two shape
and measured 258.7 us; the prior 3M Good--Thomas experiment used
M31+M61+generic-q and failed its lower-bound gate.  The combination of a 3M
FP32 estimator with the two production Mersenne fields was not attempted.

The active prototype uses FFT shape `4:512:6:512:202`.  FP32 retains the native
radix-6 middle, while M31 and M61 view the six slots as three independent
radix-2 1M channels followed by a scalar radix 3.  Three queues overlap the
FP32, M31, and M61 field work.  All transform kernels compile and launch; their
steady critical path, measured by the same timeline script as the controls, is:

| Path | transform-core median | p10--p90 | kernels/cycle |
|---|---:|---:|---:|
| production 4M M31+M61 | 129.632 us | 118.016--130.688 us | 6 |
| 3M M31+M61+q | 144.000 us | 130.176--147.104 us | 9 |
| active 3M FP32+M31+M61 | **127.648 us** | 119.680--135.967 us | 9 |

With the existing fused FP32-estimated CRT/carry used only as a launch and time
scaffold, three fresh uninstrumented processes warmed to 185.6--185.7
us/iteration at the 2,000-iteration checkpoint.  This is the first new complete-
budget estimate close enough to the requested 180 us to justify finishing the
boundary.  It is not yet a correctness result: the deterministic residue is
wrong because the NTT fields are left in channel-major Good--Thomas order where
the fused carry expects natural coefficient pairs.

All six supported transform variants (`000`, `101`, `202`, `010`, `111`, and
`212`) were then run as fresh 4,000-iteration scaffold processes.  Their warmed
2,000-iteration reports were 185.6--185.9 us, with no useful variant-dependent
headroom.  This check had not previously been performed for the odd 3M hybrid;
it rules out recovering the missing boundary allowance merely by reusing a
different 4M tuning variant.

An initially proposed six-register shuffle at `fftMiddleIn/Out` was rejected
after deriving the complete physical index map; it has **not** been implemented.
For base real indices `b0=2b`, `b1=b0+1`, channel `a` does satisfy
`a=(t+b) mod 3`, so at one fixed logical base coordinate the components obey

```text
P[t].x = V[(t+b0) mod 3].x
P[t].y = V[(t+b1) mod 3].y
```

However, PRPLL's physical natural pair index is
`j=x*(6H)+middle*H+y`, while a channel's power-of-two base index is laid out
with stride `2H`.  Substituting the exact CRT map also gives
`x_channel=(3*x+t) mod 512` (with component-dependent channel selection).
Consequently the mapping crosses width lanes as well as the six middle
registers.  A local shuffle alone would remain deterministically wrong.

The next exact boundary must therefore live at the width edge, where a
workgroup already owns the complete 512-point line in registers/shared memory,
or use the already measured shared three-channel transpose.  The preferred
gate is to absorb `x -> 3x+t` into that width edge and retain the existing exact
FP32/M31/M61 reconstruction; an extra pair of full-state permutation kernels is
unlikely to fit the remaining 5.6-us budget.  The inverse scalar radix 3 must
also be normalized by `1/3`, and the carry's NTT scale must use the 1M channel
length rather than treating middle 6 as middle 16.  Only after a matching
checkpoint will the prototype's timing be eligible for the 180-us gate and a
full 1M run.

The user explicitly permits this route, or any later route, to become a clean-
room CUDA program rather than a PRPLL integration.  That freedom changes host
and scheduling choices but not the acceptance test: exact p136 recurrence,
matching residue, at most 180 us/iteration sustained over 1,000,000 iterations.

#### 3M M31+M61+M19 architecture (component pass, complete-budget no-go)

The experiment registry was checked before implementation.  The only prior
small Mersenne-labelled arithmetic control used `2^23-1`, which is composite
(`47*178481`) and was explicitly rejected as a transform field.  M19 had not
been tested.  The prior 24-bit correction field was a generic Riesel field at
4M and did not combine the shorter 3M transform, M31+M61 capacity, or Mersenne
reduction.

For `N=3*2^20` p136 words, write

```text
E = 43*N + R,  R = 1,013,537.
```

The earlier uniform bound `6*N*2^86` assumes every balanced digit can have
magnitude `2^43`.  In fact the 44-bit words form the mechanical/Beatty set
`I[j]=floor((j+1)R/N)-floor(jR/N)`.  Because `gcd(R,N)=1`, multiplication by R
permutes indices; the maximum cyclic overlap of this set with its reflection
is exactly R.  If small-word magnitude is normalized to `1/2` and big-word
magnitude to `1`, the exact worst convolution sum is therefore

```text
2^86 * (N + 3R) / 4.
```

After the same factor-six square/MUL3 bound and the fixed point of the
base-`2^43` carry recurrence, the bound is `2^109.145617`.  The balanced
`M31*M61*M19` range exceeds it by **0.854380 bit** (factor 1.808).  This is a
p136-specific capacity proof; it does not cover p145--p150 at the same 3M
length.  Larger exponents can retain the production 4M path unless a separate
shape proof is established.

Algebraic gates pass:

- `M19=524287` is prime and `M19 == 1 (mod 3)`, so the radix-3 edge is scalar;
- `v2(M19^2-1)=20`, sufficient for a packed 1M-real-word channel whose
  quadratic transform has `2^19` values; and
- `2^((M19^2-1)/gcd(3*2^20,M19^2-1)) == 1 (mod M19)`, so the exact weighted
  transform root exists.

[`src/cuda/m31_m61_q3m_tile_bench.cu`](src/cuda/m31_m61_q3m_tile_bench.cu)
now has an M19-specialized path using ordinary Mersenne folding rather than
Montgomery form; `make m31-m61-m19-3m-tile-bench` builds it.  An independent
direct quadratic cyclic convolution validates every field.  It uses 38
registers with no stack or spills.  Three fresh paired runs measured:

```text
4M M31+M61 control:       113.088--113.984 us
3M generic-M19 control:   108.352--109.568 us
3M M19-special:           101.664--101.792 us
M19-special/current:      0.8930--0.8990
```

This passes the arithmetic tile gate and is a much stronger distinction than
merely selecting another 25--31-bit Montgomery prime.  The next required gates
are an exact M19 scalar radix-3/square/inverse edge and a full-range
M31*M61*M19 CRT/carry/forward-weight kernel, including Good--Thomas layout.

Both follow-up gates now pass.  The existing radix-3 harness was rebuilt with
M19 and independently validated the complete forward radix 3, quadratic
square, inverse radix 3 sequence in all fields:

```text
all fields radix3/square/inverse: 23.872 us
M61 path:                         12.096 us
M31+M19 path:                     12.352 us
M31+M19 square-only:              11.712 us
```

The M19 radix in this gate still uses the generic Montgomery helper, so the
0.64-us M31+M19 incremental edge is a conservative timing; the exact transform
will use Mersenne folding.

The carry benchmark was extended rather than creating a new untracked harness.
For M19 it now reduces M31 into M19 by folding, multiplies the first Garner
digit with M19 reduction, reconstructs the final digit through the existing
M61 fold, and writes a nonzero representative M19 power-of-two weight as a
rotation.  Random coefficients span the complete balanced 110-bit CRT range;
all 3,145,728 output digits, wide carries, and Good--Thomas maps match an
independent signed-`__int128` reference.  The specialized kernels use 40--42
registers, 12 KiB shared memory in the transpose version, and no spills:

```text
4M M31+M61 natural carry control: 28.192 us
3M M31+M61+M19 natural:           32.192 us
3M M31+M61+M19 GT shared:         38.400 us
```

This is essentially the same exact shared-layout cost as the large-q path, but
the transform tile is about 7--8 us faster.  M19 therefore passes the component
gates and should proceed to a production-shaped transform-core measurement.
It is not yet an end-to-end result: the exact width edge and full recurrence
remain mandatory.

That production-shaped gate is now implemented behind `M19_FIELD=1`.  The
source/registry audit above was performed before editing: M19 had not previously
been a transform field, while the direct-prime, generic-q Good--Thomas, M23
mislabel, 24-bit correction field, standalone scheduler, and large shared-tile
designs had all already been attempted.  The new opt-in path reuses only the q1
storage/scheduling lane and adds:

- ordinary M19 host roots and Crandall--Fagin weights, including cache-key
  separation from generic q1 tables;
- device `2^19-1` folding and power-of-two rotations instead of Montgomery
  multiplication; and
- the exact M19 radix-3 constants in the existing three-channel scheduler.

It compiles without changing the normal engine.  The current Good--Thomas
boundary is still the deliberately incompatible legacy carry scaffold, so its
residues are not correctness evidence.  Nsight Systems was therefore used only
for the same `first fftMiddleIn -> last fftMiddleOut` transform-core gate as the
earlier controls.  Fresh paired traces on the current source measured:

```text
3M M31+M61+generic-q core:  median 139.776 us, p10--p90 125.984--142.624
3M M31+M61+M19 core:        median 137.728 us, p10--p90 124.192--140.480
M19 improvement:                   2.048 us
```

The M19 stream itself is faster: its middle/tail/width/forward kernels save
about 6 us in their individual summed durations.  Most of that gain is hidden
because M19 overlaps the already-dominant M31 and M61 queues and competes for
the same integer/cache/power resources.  This is the same forecasting lesson as
the earlier q31 campaign: a faster isolated field does not translate
proportionally through a saturated multi-field schedule.

The complete-budget lower bound is now

```text
137.728 us  measured transform core
 38.400 us  exact M31*M61*M19 shared-layout CRT/carry
-----------
176.128 us  before inverse width and forward width
```

Only 3.872 us remains under the 180-us gate.  Both width transforms are
mandatory in PRPLL or a clean-room CUDA engine; changing the host scheduler
cannot remove their arithmetic.  For scale, one required M31 width kernel alone
is about 11 us in the paired trace, and the mature production fused
width/CRT/carry region is about 70.5 us for only two fields.  A new three-field
edge would have to fit both width transforms, the Good--Thomas permutation,
110-bit CRT, carry propagation, weighting, and stores into less time than the
standalone carry kernel alone.  The standalone-program permission therefore
does not rescue this representation.

Decision: **reject M31+M61+M19 as a route to 180 us in this transform
architecture.**  Retain the field specialization and exact component harnesses
as useful evidence, but do not implement the exact width boundary or run a
million incorrect scaffold iterations.  Reopen it only if a materially
different transform eliminates at least roughly 25--30 us at the complete-edge
level, not merely a few microseconds from the M19 field.

An isolated production regression check after adding the opt-in field still
matches the normal p=2,000 (`05d6515c416b83e2`) and p=10,000
(`52316d51aa52e6b7`) checkpoints; its warmed 2,000--10,000 interval was
197.0 us.  Thus the conditional M19 code does not change the default path.

## Experimental 2M M31+M61+q0+q1 architecture

The strongest local adaptation of the forum proposal is not a direct forum
claim: retain production M31 and M61, add the two already-integrated Riesel
fields, and halve the transform from 4M to 2M.  Its exact modulus product is

```text
P = (2^31-1)(2^61-1)(2090860543)(2141192191)
  = 22168704719850339839374414813441924894304501761
log2(P) = 153.957216843
```

At 2M words, `136279841 / 2^21 = 64.983...`, so centered digits use 64 or
65 bits and have magnitude below `2^64`.  A deliberately conservative cyclic
square bound is `2*NWORDS*2^128`; multiplication by three raises this to
`2^151.584963`.  The balanced CRT limit is `P/2 = 2^152.957217`, leaving
1.372 bits before the much smaller incoming carry.  Thus coefficient range is
not the immediate blocker.  The state is five 32-bit limbs by 2M, about 40 MiB
instead of 48 MiB, and first-order limb-stage work is approximately `0.80x`
the current transform.

### Exact 154-bit CRT/carry gate

`src/cuda/rns31_crt_bench.cu` now contains an exact full-range prototype.  It
uses a hierarchical Garner reconstruction:

1. Reconstruct canonical `x3` modulo `M31*q0*q1` with the existing q-pair plus
   M31 path.
2. Reduce that 93-bit value modulo M61 with Mersenne folds.
3. Compute the remaining M61 Garner digit with one general M61 multiply.
4. Form and balance the 154-bit result in a five-limb value.

Random values across the complete balanced range, including zero, both signs,
`floor(P/2)`, and boundary values, match the host reference.  A four-word
sequential carry prototype then multiplies by three, adds a signed 96-bit
incoming carry, extracts centered 64/65-bit digits, and propagates the resulting
approximately 89-bit carry.  It also matches the host reference.

At architectural sizes, 4M current coefficients versus 2M four-field
coefficients, median CUDA-event timings were:

| Gate | Median | Registers | Local stack |
|---|---:|---:|---:|
| Current 4M M31+M61 CRT | 0.024 ms | 14 | 0 |
| New 2M four-field CRT | 0.023 ms | 22 | 0 |
| New 2M four-field CRT + four-word carry | 0.020 ms | 46 | 0 |

The apparent carry timing is lower because it writes a narrower result and
reuses warm cache state; it is not an end-to-end speedup.  The useful result is
that exact 154-bit arithmetic does not spill and is not by itself a no-go.

### Four-plane transform-throughput gate

An opt-in `RIESEL_FOUR=1` prototype adds the production M61 plane to FFT type
54 while retaining separate M31, q0, and q1 planes.  Data and trig offsets are
independent, and no field tile is fused with another.  The old three-field
carry was temporarily used only as a timing placeholder; therefore residues
from these runs are intentionally invalid and the normal source again rejects
the mode at its on-load correctness check.

The first three-queue assignment put M61 and q0 together and measured about
236--242 us/iteration without profiling.  Sequential profiling showed why:
M61's transform region costs roughly 69 us, while M31, q0, and q1 each cost
about 45--46 us.  Reassigning the queues to

```text
queue 0: M31 + q0
queue 1: M61
queue 2: q1
```

reduced a 10k, profiling-disabled run to **225.0 us/iteration**.  A matched
profiling run was 240.4 us/iteration; profiling itself is material here.  A
single-queue control was 308.6 us/iteration, and `1024:4:256` was slower than
the `512:4:512` shape.

This is the first shorter-transform prototype with plausible headroom, but it
has not achieved a speedup.  The 225.0 us result still uses a cheaper,
incorrect three-field carry and is 9.5% slower than the 205.5 us production
path.  A proper four-field carry cannot improve that number by merely replacing
the placeholder.  It would need narrow inverse/forward edge fusion to recover
more than 19.5 us plus the extra exact-carry cost.  Consequently, a full
checkpoint-capable integration is deferred until an isolated narrow-edge
prototype demonstrates that recovery; a giant four-field fused kernel is
specifically contraindicated by both the local spill result and the literature.

## Sub-2^30 Harvey-field experiment

Harvey's modified Shoup/Montgomery butterflies require `p < beta/4` and keep
values in redundant `[0,2p)` or `[0,4p)` ranges.  With `beta=2^32`, the exact
threshold is `p < 2^30`.  They remove a final Montgomery correction and defer
two reductions in the inverse butterfly
([Harvey 2014](https://arxiv.org/pdf/1205.2926)).  The current q0/q1 near
`2^31` cannot use these ranges; the verified approximately 1.03-billion
candidates above can.

This now outranks more full-field fusion.  Its gate is one range-annotated
`GF(q^2)` radix tile, not an end-to-end prime substitution.  It must:

- prove every lazy input/output range, including complex multiplication;
- materially reduce the current q-tail counts of roughly 533 `IMAD`, 258
  `ISETP`, and 395 `VIMNMX` instructions;
- beat the canonical q tile in elapsed time; and
- retain the exact 90.8938-bit capacity proof for the 4M M31+qA+qB path.

If it passes, q fields remain separate grid blocks and only the narrow
inverse/forward edges are candidates for fusion.  On-the-fly twiddles and
unrelated-exponent batching remain secondary: published 30-bit versus 60-bit
GPU comparisons differ by only about 5% once the extra prime count is included,
and batching mainly helps transforms too small to saturate the GPU
([Kim et al. 2020](https://arxiv.org/abs/2012.01968)).  Dense Tensor Core NTT
reformulations remain rejected by their limb-expansion cost and size limit,
consistent with the local 42% tile deficit
([TensorFHE](https://arxiv.org/abs/2212.14191)).

### Implemented range and capacity gates

The experiment is now complete.  [`src/cuda/riesel_lazy_tile_bench.cu`](src/cuda/riesel_lazy_tile_bench.cu)
contains the isolated CUDA arithmetic gate, and opt-in `RIESEL_LAZY=1` selects
the same fields in complete FFT type 54 PRPLL runs.  The selected ascending
lane order is:

| Lane | Modulus | Form | exact order-`2^21` norm-one root | root-of-two `theta` |
|---:|---:|---:|---:|---:|
| 0 | 1,031,798,783 | `123*2^23-1` | `(451776209,678011953)` | 286,716,770 |
| 1 | 1,038,090,239 | `495*2^21-1` | `(219392826,888955106)` | 571,191,171 |

For each field, transform values stay in `[0,2q)`.  Addition and subtraction
reduce modulo `2q`; because `4q < 2^32`, they cannot overflow a `uint32`.
For Montgomery inputs `a,b < 2q`, `a*b < 4q^2 < 2^62`; the REDC sum fits in
64 bits and its uncorrected high word is below `2q`.  The final `result >= q`
correction is therefore omitted throughout the transform and performed only
when carry/CRT needs a canonical residue.  Complex Karatsuba multiplication
preserves the same invariant.

This is the `[0,2q)` portion of Harvey's redundant-range design.  It does not
implement the complete `[0,4q)` inverse-butterfly schedule: the generic
quadratic-field Karatsuba multiply needs both components normalized below
`2q`, and keeping `[0,4q)` complex values would either add normalizations or
replace three products with four.  The remaining opportunity is too small to
change the end-to-end decision below.

The exact combined modulus is

```text
P = (2^31-1) * 1031798783 * 1038090239
  = 2300170260959993715375472639
log2(P)   = 90.8937992168
log2(P/2) = 89.8937992168
```

For `N=2^22`, the conservative post-`*3` coefficient bound
`6*N*2^64` has logarithm 88.5849625007.  Adding the fixed point of the
base-`2^32` incoming-carry recurrence changes the margin negligibly, leaving
1.3088367158 bits.  `riesel-algebra-test` now verifies both primes, the exact
roots and weights, Hermitian convolution, `q<2^30`, and this integer capacity
inequality.

### Tile and SASS results

The CUDA gate processes 2,097,152 `GF(q^2)` values with a radix-8 transform
and pointwise square.  A one-round tile is bandwidth-bound and shows no gain:
the representative lazy/canonical median ratio is 1.007.  Repeating the
arithmetic four times in registers gives a median ratio of about 0.954 across
six process runs, or a 4.6% improvement.  Both variants use 40--44 registers
and no stack.

Final SASS for complete production-shaped q kernels shows that the arithmetic
change is real but leaves multiplication and traffic untouched:

| q kernel | Canonical instructions | Lazy instructions | Canonical/lazy `VIMNMX` | `IMAD` change |
|---|---:|---:|---:|---:|
| middle-in | 1,008 | 840 | 157 / 84 | 230 / 230 |
| tail multiply | 8,032 | 7,184 | 1,207 / 827 | 1,231 / 1,231 |
| middle-out | 1,008 | 848 | 157 / 84 | 230 / 230 |
| width | 776 | 712 | 106 / 77 | 105 / 105 |
| tail square | 2,856 | 2,552 | 395 / 267 | 533 / 533 |

The tail also drops from 264 to 246 `ISETP` instructions and from 61 to 60
registers, with no stack in either version.  The separate carry kernel remains
424 instructions and 38 registers; the prime substitution does not simplify
the cross-coefficient dependency or number of residue planes.

### End-to-end result

The integrated path uses the real roots/trig tables, weights, scale constants,
91-bit balanced Garner reconstruction, and canonicalizes lazy residues at the
carry boundary.  Packed, four-field, and giant-fused modes are deliberately
incompatible with `RIESEL_LAZY=1`; this test retains independent field queues
and the split carry.

Correctness results for exponent `136279841` were:

- 10k and 100k lazy runs matched the canonical M31R2 residue at every reported
  checkpoint.
- The 1M run produced `52b03a7cc55e677d`, matching production and the previous
  M31R2 run, with no carry or residue errors.
- The 1M lazy time was 285.4 us/iteration.

Two same-session, profiling-disabled 100k controls avoid comparing different
thermal/clock states:

| Path | Samples (us/iteration) | Mean/median of pair |
|---|---:|---:|
| Original near-`2^31` M31R2 | 291.3, 287.4 | 289.35 |
| Sub-`2^30` lazy M31R2 | 290.1, 285.2 | 287.65 |

The measured improvement is 1.70 us/iteration, or **0.59%**.  Instrumented
10k profiles show reductions in the transform regions—aggregate tail square
80.3 to 75.4 us, width 31.0 to 28.3 us, middle-in 39.1 to 37.3 us, and
middle-out 28.4 to 27.2 us—but these averages include the unchanged M31 field,
and profiling materially changes total time.

Decision: retain the mode and benchmark as evidence, but reject this prime
substitution as the main speedup route.  It validates the literature-derived
range idea and gives a small M31R2 improvement, yet remains about 80 us slower
than the 205.5 us production M31+M61 iteration.  Applying it to the 225.0 us
four-field scheduling probe can recover only a few microseconds and cannot
make that already-optimistic placeholder-carry result beat production.

## Field-specialized fused edge and 2M Goldilocks lower bound

[`src/cuda/warp_specialized_edge_bench.cu`](src/cuda/warp_specialized_edge_bench.cu)
tests whether the mature M31/M61 fusion can be reorganized without retaining
both field vectors in every thread.  Its production-shaped tile has 64 lanes,
eight quadratic elements per lane, two twiddled radix-8/shared-memory stages
plus a final radix-8, a scalar cross-field bridge, and a second transform edge.
It compares:

- a 64-thread block that keeps M31 and M61 vectors live per thread, matching
  the production fusion's register/dataflow shape; and
- a 128-thread block in which one 64-thread subgroup owns M31 and the other
  owns M61, with a shared-memory handoff at the scalar bridge.

At 2,097,152 field elements, five repeat processes measured 0.077 ms for the
serial form and 0.090 ms for the field-specialized form.  Specialized/serial
ratios were 1.163, 1.170, 1.166, 1.167, and 1.164: a stable **16.6% slowdown**.
At the 1,048,576 elements relevant to a 2M-word architecture, the slowdown was
about 9--10%.

The expected register separation did not occur.  Final SASS/resources were:

| Edge organization | Registers/thread | Shared/block | Static instructions | Barrier instructions |
|---|---:|---:|---:|---:|
| Serial per-thread fields | 90 | 9,216 B | 6,144 | 16 |
| Field-specialized block | 98 | 13,312 B | 6,376 | 11 |

CUDA assigns a kernel-wide register footprint based on the largest control-flow
path; logical field ownership does not create independently allocated register
files.  The specialized form therefore doubles threads per tile, increases
registers and shared memory, and adds field-selection/handoff work without an
occupancy benefit.  Inline PTX can alter individual instructions but cannot
change this resource-allocation model.  Warp specialization is rejected for
this fused edge.  A cluster/two-kernel form would restore separate compilation
but would communicate through global or distributed shared state and recreate
the split-path traffic that costs about 38 us end to end.

The same gate measures the remaining forum-derived
`M31*M61*Goldilocks` 2M proposal with the correct general quadratic arithmetic.
At 1,048,576 elements, five process medians were 0.029--0.030 ms for an M61
edge pair and 0.034 ms for Goldilocks.  Goldilocks/M61 ratios were 1.145--1.166,
with a median near 1.15.  This agrees with the earlier exact unit-norm radix
tile and shows that shortening the transform does not make the Goldilocks
field cheaper than M61.

The measured four-field schedule has a roughly 90 us critical queue
(`M31+q0`) and a 69 us M61 field.  Replacing both q fields by one Goldilocks
field estimated at `1.15*69 = 79.4 us` reduces the critical schedule by only
about 10.6 us.  Applied optimistically to the 225.0 us placeholder-carry result,
the lower bound is about **214 us/iteration before a correct 156-bit carry**.
It is already slower than production, so a full Goldilocks integration is
rejected without risking the main code path.

## Compact two-phase middle/width/carry fusion

Before implementing this route, the experiment registry was checked against the
earlier giant M31R2 fused kernel, field-specialized edge, carry transfer scan,
and cluster/DSM proposals.  The untested distinction was a two-kernel grid-wide
boundary at a *compact integer representation*, while retaining the incumbent
M31/M61 fields.  It was therefore sufficiently different to merit one complete
gate.

### Production-sized proxy

[`src/cuda/warp_specialized_edge_bench.cu`](src/cuda/warp_specialized_edge_bench.cu)
was extended to cover 2,097,152 quadratic values, the coefficient count of the
4M PRPLL transform.  It modeled

```text
inverse middle -> inverse width -> scalar bridge
grid-wide boundary
scalar bridge -> forward width -> forward middle
```

and compared the chain with the same stages launched separately.  Initial
results were unusually promising:

| Proxy bridge | Separate chain | Two-phase fused | Change |
|---|---:|---:|---:|
| Whole 96 KiB field tile, trivial bridge | 0.162--0.163 ms | 0.121--0.122 ms | about -25% |
| Two full 64-bit words/value | 0.159--0.160 ms | 0.140--0.141 ms | about -12% |
| 12 bytes/value | about 0.159 ms | about 0.138 ms | about -13.5% |
| Two u32 words plus two ballot bitplanes | 0.159--0.160 ms | 0.128--0.129 ms | about -19.5% |

At 32.49 bits/word, the first sloppy-carry word is an unsigned 33-bit value and
the second is a signed 33-bit value.  Their low words consume 8 bytes per
quadratic coefficient; the two 33rd bits consume two warp-ballot planes.  The
exact bridge is therefore 8.25 bytes/value, or 16.5 MiB, rather than a 32 MiB
pair of `int64` values.  The proxy passed its internal output comparison and its
two kernels used about 66/70 registers without local spills.

### Exact PRPLL integration and correctness

The opt-in implementation is in
[`src/cl/carrymiddle.cl`](src/cl/carrymiddle.cl).  Supporting changes teach the
CUDA kernel wrapper to request dynamic shared memory.  `MIDCARRY_FUSED=1` selects
the complete middle/width/carry path; `MIDCARRY_FUSED=2` retains the ordinary
middle kernels and substitutes only the split carry, which was used as a
diagnostic.  Both modes are restricted to CUDA, FFT3161, in-place
`512:8:512`, short carry, and no parity mode.  The default production path is
unchanged.

Two correctness defects were found and isolated:

1. The default width shared-memory padding overlapped the adjacent M31/M61
   tiles.  The new kernels compile with `LDSPAD_W=0` and explicitly request
   98,304 bytes of dynamic shared memory.
2. The forward half initially reused inverse weight counters adjusted for the
   `2*NWORDS` transform scale.  Production `carryFused` deliberately restores
   the saved pre-normalization counters before rebuilding the forward fields.
   Restoring those counters fixed both the split diagnostic and full fusion.

The full-64-bit diagnostic ruled out ballot packing as a source of error.  After
the counter fix, the compact representation reproduced every production
checkpoint through 10,000 iterations, including:

```text
iteration 1,000:  ec8a091a87311fac
iteration 2,000:  05d6515c416b83e2
iteration 10,000: 52316d51aa52e6b7
```

### Why the proxy gain reversed

Exact profiling gives the decisive comparison:

| Path/region | Time per call |
|---|---:|
| Production `carryFused` | about 68 us |
| Split compact carry, inverse half | 64.8 us |
| Split compact carry, forward half | 56.2 us |
| Full middle/carry fusion, inverse half | 114.2 us |
| Full middle/carry fusion, forward half | 95.6 us |

The split bridge therefore costs about 121 us before the ordinary middle stages,
versus 68 us for production's single stairway kernel.  Folding in the middle
stages raises the exact pair to about 210 us.  The generated full kernels use
98 and 90 registers/thread, 512 threads, 96 KiB dynamic shared memory, and no
local spills.  This admits only one large block per SM and has much less latency
hiding than production's 128-thread fused carry.  The exact inverse/forward
weighting, M31/M61 CRT, carry formation, and final carry propagation also create
long dependency chains that the proxy's simple scalar bridge did not model.

A profiling-disabled 10k run of the compact exact path measured 311.0
us/iteration and matched the expected residue.  A nearby production control
measured 230.3 us/iteration under the same short-run conditions.  Earlier
matched full-word tests were similarly unfavorable: about 327.3 us fused and
302.5 us split-only versus 228.5 us production.  Compact packing helps the
experimental path, but cannot recover the occupancy and extra-grid-boundary
cost.

Decision: reject this fusion architecture for performance.  Keep the exact
opt-in implementation and proxy as diagnostic evidence, but do not enable it by
default or infer a speedup from the proxy.  The new lesson is that a fusion proxy
must include the incumbent's exact CRT/carry arithmetic and resulting register
dependency graph before its memory-traffic result is predictive.

### Multi-worker confirmation, not a new experiment

After re-reading the registry, unchanged multi-worker execution was recognized
as already comprehensively tested above.  A short repeat was allowed to finish
but is not treated as a new route.  In the current thermal state, one worker
measured 199.4 us/iteration over 100k.  Two 50k workers measured 484.7/485.4 us,
or about 242.5 us per aggregate useful iteration.  `TAIL_KERNELS=3` improved
that only to 481.7/482.0 us, or 240.9 us aggregate.  This independently confirms
the existing conclusion that unchanged concurrent streams reduce throughput on
this power-limited GPU.

## From-scratch scope: mixed length and two-limb M61 gates

The performance gate is now explicitly implementation-independent: a standalone
program is acceptable if it performs the same exact PRP iteration and sustains
at most **180 us/iteration over 1,000,000 iterations** at exponent `136279841`,
with the expected final state/residue.  PRPLL compatibility is secondary to
that gate; this was explicitly reconfirmed after the folded-syndrome work, so a
clean-room CUDA scheduler/data layout is in scope if PRPLL's orchestration
becomes the limiting factor.  The same design must have a credible capacity path for the broader
140--150M exponent range.  The experiment registry and source tree are checked
before each new implementation so that a new harness is not merely a renamed
version of an earlier rejected route; any justified repeat must state the
architectural difference before code is written.

### 33/32-length FP32+M61 near-miss: algebraic no-go

The existing 4M `FP32+M61` path is fast enough at about 171 us/iteration but is
not exact at exponent `136279841`; its configured maximum is `132791664`.  A
transform only about 2.6% longer would nominally restore the missing bits while
leaving a plausible 180 us budget.  The attractive factorization is

```text
NWORDS = 33 * 2^17 = 4,325,376
ND     = 33 * 2^16 = 2,162,688
```

which is 33/32 times the current length.  This was checked against the logged
8M low-precision controls and the forum mixed-radix experiments; no 33/32
`FP32+M61` implementation had previously been attempted locally.

It fails before a kernel gate.  The quadratic DGT packs two real coefficients
only when the transform root is in the norm-one subgroup, so conjugation equals
inversion.  For M61 that subgroup has order

```text
M61 + 1 = 2^61
```

and therefore contains no order-3 or order-11 element.  M61 has odd roots in
its base-field subgroup, but multiplying one into the power-of-two root breaks
the conjugate-packing invariant.  Using the general quadratic field without
that packing doubles the M61 state/work; replacing M61 with a generic prime
whose `q+1` contains 33 loses the cheap Mersenne arithmetic.  Either change
erases the sub-180 projection.  PRPLL accordingly permits odd middle radices
only for FP64 and filters them out for every NTT type.  No mixed-radix kernel was
implemented after this mathematical gate failed.

#### Reopened clean-room gate: near-M61 norm-one radix-33 field

The standalone-engine permission was reaffirmed after the production-shaped
Goldilocks no-go.  It is a first-class acceptable outcome: PRPLL integration is
not required if a new exact program sustains the 180-us, one-million-iteration
gate and has a credible general path through the 140--150M range.  Reusing the
same production kernels under a new scheduler remains excluded by the measured
graph and multi-worker results.

The preceding radix-33 rejection was an unmeasured projection for a *generic*
prime, not an arithmetic experiment.  A source/registry search found no test of
a near-M61 field chosen jointly for a small pseudo-Mersenne complement and an
odd norm-one factor.  Deterministic 64-bit Miller--Rabin and `openssl prime`
both identify

```text
q = 2^61 - 58,327,041 = 2,305,843,009,155,366,911
q + 1 = 2 * (33 * 2^16) * 533,096,546,787.
```

Thus the norm-one subgroup of `GF(q^2)` contains the order-`33*2^16`
transform root, including the factor needed by the 4,325,376-word packed DGT,
while retaining a 61-bit modulus only 25.8 bits below `2^61`.  A wide product
can be reduced by folding at bit 61 twice: first multiply the high 61 bits by
58,327,041, then fold the at-most-26-bit high part once more.  This is
architecturally different from the rejected 31-bit Riesel shift/add test, which
used four serial shift/add folds, and from M31R2, which replaced one 64-bit
plane with two separately transformed 32-bit planes.  Here the intended full
architecture is the fast FP32 plus one packed quadratic 61-bit field at 33/32
the old length.

This entry reopens only an early-abandonment arithmetic/tile gate; it does not
claim a speedup.  The margin is severe: scaling the approximately 171-us 4M
FP32+M61 path by 33/32 leaves only about 4.6 us for all extra modular and
radix-33 cost.  The 33/32 length is enough for exponent 136279841 but not by
itself the entire 140--150M interval, so a successful implementation would
also need a related larger odd-factor shape.  Before any mixed-radix transform
is written, compare exact quadratic multiply chains against M61 and inspect
their SASS; reject the route if that component cannot plausibly preserve the
end-to-end budget.

That gate is now complete in
[`src/cuda/near_m61_radix33_bench.cu`](src/cuda/near_m61_radix33_bench.cu).
It validates every measured path against independent host `unsigned __int128`
arithmetic and compares two candidate reducers with the incumbent M61
quadratic multiply:

- the two-fold `2^61 == 58,327,041 (mod q)` reducer; and
- 64-bit Montgomery REDC, which omits the low word of `m*q` because its carry
  is exactly `productLo != 0`.

At 2,097,152 quadratic values, interleaved 21-sample medians were:

| Dependent complex-multiply chain | M61 | Two-fold q | Montgomery q |
|---:|---:|---:|---:|
| 4 | 0.0332 ms | 0.0482 ms (1.45x) | 0.0455 ms (1.37x) |
| 8 | 0.0542 ms | 0.0851 ms (1.57x) | 0.0815 ms (1.50x) |
| 16 | 0.0986 ms | 0.1526 ms (1.55x) | 0.1406 ms (1.43x) |
| 32 | 0.1804 ms | 0.2904 ms (1.61x) | 0.2693 ms (1.49x) |

All kernels have zero stack and spills; ptxas reports 32 registers for M61,
38 for the fold, and 36 for Montgomery.  In the chain kernel, final SASS has
39 IMAD-family instructions for M61 versus 102 for Montgomery.  The one-round
case is memory-bound and tied, but an NTT necessarily executes many arithmetic
rounds per loaded tile.  A 37--49% arithmetic penalty cannot fit in the roughly
4.6-us total allowance, even before the 3.125% extra data and radix-33
butterflies.  Decision: **reject the near-M61 radix-33 architecture at the
arithmetic gate**.  Do not implement its roots or mixed-radix NTT unless a new
one-wide-product modular reducer is first demonstrated on this GPU.

#### Batch-native M61 gate (distinct from multiple workers)

The registry contains independent-process/stream batching and a dense Tensor
Core batch tile, but not the proposed exponent-interleaved AoSoA implementation
from the 2026-08-09 follow-up list.  Two exponents of the same transform shape
can use one kernel instruction stream and one twiddle load for two independent
M61 values.  This is a real computational change, unlike `-workers 2`, and is
eligible because aggregate unrelated-exponent throughput is an accepted target.

The risk is that production M61 tiles already reuse small trig tables from
cache while their modular arithmetic saturates the integer pipelines.  Before
altering the PRPLL state layout, measure an exact paired kernel that deliberately
shares every twiddle and address calculation.  Its useful per-exponent time must
beat scalar execution by at least 12.4%; a smaller idealized gain cannot close
205.5 to 180 us and rejects the much more invasive whole-engine conversion.

Two exact gates are now implemented.  The simple
[`src/cuda/m61_batch_bench.cu`](src/cuda/m61_batch_bench.cu) applies one loaded
quadratic twiddle to two AoSoA exponent values.  Its useful throughput gain
falls rapidly as the arithmetic per load approaches an NTT tile:

| Dependent complex multiplies | Batch-2 throughput/scalar |
|---:|---:|
| 2 | 1.198x |
| 4 | 1.121x |
| 8 | 1.037x |
| 16 | 1.028x |
| 32 | 1.007x |

The stronger
[`src/cuda/m61_middle_batch_bench.cu`](src/cuda/m61_middle_batch_bench.cu)
models both production `MIDDLE=8` twiddle regions and the cheap M61 radix-8.
One thread owns the same coordinate from two exponents, so it shares the entire
`middleMul2` root-power chain as well as all direct `middleMul` root loads.  All
4,194,304 batch outputs exactly match the scalar path after every timed launch.
It deliberately omits the production transpose, making it an optimistic upper
bound.  ptxas reports no stack/spills, but registers rise from 66 to 98.

| Threads/block | Scalar | Batch time | Useful batch time | Throughput |
|---:|---:|---:|---:|---:|
| 64 | 0.0290 ms | 0.0524 ms | 0.0262 ms/exponent | **1.107x** |
| 128 | 0.0285 ms | 0.0525 ms | 0.0263 ms/exponent | 1.087x |
| 256 | 0.0294 ms | 0.0606 ms | 0.0303 ms/exponent | 0.969x |

The 64-thread result reduces this component's time by only 9.7%.  Applied to
the approximately 56-us production M61 middle-in/out pair, it saves about 5.4
us from the 205.5-us iteration.  Even unrealistically granting the same gain to
the complete approximately 110-us M61 transform chain projects about 195 us;
the roughly 74-us exponent-specific CRT/carry cannot share weights or state
between unrelated exponents.  A paired-lane version would retain scalar
registers but cannot reduce issued root arithmetic under SIMT: each warp covers
half as many coordinates and still issues the predicated instructions.

Decision: **reject batch-native M61 as a route to 180 us on this GPU**.  It is
a real component throughput gain at small blocks, unlike multiple workers, but
not an end-to-end architectural speedup of the required size.  Do not integrate
AoSoA PRPLL state unless a later design also shares or removes carry work.

#### Reopened direct-RNS gate: cache-fused scalar transform

The completion matrix's direct-three-prime no-go refers to the implemented
global-pass and generic shared radix paths.  The corresponding chronological
entry explicitly left a `512 x middle x 512`/cache-fused transform as future
work, and a source/registry search finds no such scalar implementation.  This
is therefore distinct from the measured 0.707--0.949-ms pipelines.  The earlier
exact radix-16 cyclic-square tile was 15.7--15.8% faster than the combined
M31/M61 tile, so the missing factorization is still an eligible clean-room gate.

The three existing scalar primes have a 92.639-bit product, slightly wider than
M31 times M61, and roots through at least order `2^22`.  They consequently retain
the exact coefficient-capacity route for the production 4M geometry and the
required 140--150M interval; no p136-only sparse correction is involved.  The
first prototype will keep one 4M scalar plane in Montgomery form and use three
large cache-fused kernels: high DIF, low DIF/square/DIT, and high DIT/scale.
This is the same whole-transform gate that made the scalar Goldilocks experiment
meaningful, but with 32-bit rather than generic 64-bit multiplication.

A one-plane cyclic square must approach 35--40 us.  Three planes at 40 us
already consume 120 us before a wider three-way CRT/carry, leaving at most 60 us
under the 180-us target.  If the exact one-plane implementation materially
misses that bound, neither plane fusion nor host scheduling can rescue it and a
three-plane integration should not be attempted.

[`src/cuda/rns31_shape_bench.cu`](src/cuda/rns31_shape_bench.cu) now closes
this exact gap.  It keeps the 4M state in Montgomery representation and factors
the transform into `1024 x 4096`: high DIF, a 12-stage shared-memory low
DIF/square/DIT kernel, and high DIT with inverse scaling.  The 16-KiB kernels use
38--40 registers and have no stack or spills.  A dense full-array round trip and
an independently checked sparse cyclic square both pass for every coefficient.

```text
4M one-plane q31 forward/square/inverse:
median 0.153088 ms, mean 0.153032 ms, minimum 0.150912 ms
```

This is about 32% faster than the prior 0.224-ms one-plane full-Montgomery
global-pass result, demonstrating that the cache-fused factorization was a real
missing experiment.  It is nevertheless almost four times the 40-us plane gate
and consumes 85% of the entire 180-us allowance by itself.  Three planes project
above 450 us before CRT/carry; fusing their indexing cannot supply the required
approximately 2.6x additional improvement.  Decision: **reject direct scalar
three-prime RNS even with the cache-fused transform**.  The completion matrix's
direct-RNS no-go now covers the previously unfinished production-shaped case.

### Exact radix-`2^31` M61 arithmetic

A distinct untested idea retained M61 itself, its roots, shifts, weights, CRT,
and 16-byte quadratic-field representation, but split each scalar residue as

```text
a = a0 + 2^31*a1,  a0 < 2^31, a1 < 2^30.
```

With `2^62 == 2 (mod 2^61-1)`, Karatsuba reduces a base product to three
32x32-to-64 products:

```text
t0    = a0*b0
t1    = a1*b1
cross = (a0+a1)(b0+b1) - t0 - t1
a*b   = t0 + 2*t1 + 2^31*cross  (mod M61).
```

[`src/cuda/m61_limb_bench.cu`](src/cuda/m61_limb_bench.cu) implements this
formula and compares it with the incumbent `a*b` plus `__umul64hi` Mersenne
fold.  Both paths use the same canonical `uint64`, the same three-product
quadratic multiply, identical input/output bytes, and 32 registers with no
spills.  All GPU outputs match each other and an independent `unsigned
__int128` host reference.

For 2,097,152 quadratic values:

| Chained quadratic multiplies | Current wide M61 | Radix-`2^31` limbs | Limb/wide |
|---:|---:|---:|---:|
| 1 | 0.01584 ms | 0.01594 ms | 1.006 |
| 2 | 0.01802 ms | 0.02384 ms | 1.323 |
| 4 | 0.02826 ms | 0.04038 ms | 1.429 |
| 8 | 0.04909 ms | 0.07322 ms | 1.492 |
| 16 | 0.09149 ms | 0.14032 ms | 1.534 |
| 32 | 0.17741 ms | 0.27565 ms | 1.554 |
| 64 | 0.34323 ms | 0.53786 ms | 1.567 |

Decision: reject before an NTT tile or PRPLL integration.  The nominal reduction
from four partial products to three is outweighed by limb extraction, the
Karatsuba recombination chain, two unequal-radix folds, and normalization.
Blackwell/ptxas already lowers the incumbent 64-bit expression efficiently.
This is the same forecasting lesson as pseudo-Mersenne shift/add q reduction:
count the complete reduction dependency graph, not only multiplier count.

### FP32+M61 stage-budget profile

The existing `FFT3261` correctness and shape experiments did not include a
kernel profile.  A profiling-only run was therefore made without changing or
retuning the algorithm.  At the target exponent it reproduced the deterministic
iteration-2000 failure, so a safe exponent (`120000007`) was used to let the
identical 4M `2:512:8:512:202` pipeline reach a steady profile:

| Main-queue FP32 region | Average kernel time | Profile share |
|---|---:|---:|
| fused width / reconstruction / carry / next width | 58.2 us | 36.89% |
| middle out | 51.5 us | 34.22% |
| tail square | 25.3 us | 16.80% |
| middle in | 15.6 us | 10.37% |
| periodic ROE fused carry | 55.1 us when invoked | 1.62% amortized |

The instrumented total was 182.8 us/iteration; the previous uninstrumented
target run warmed to about 171 us before its error.  The named middle/tail
regions are the FP32 cache group.  With `MULTI_Q=1` they execute concurrently
with the separate M61 cache group, whose corresponding bottom-half work is
roughly 102 us in the split M31/M61 profile.  Thus approximately 10 us of extra
FP32 bottom-half work can potentially hide behind M61, but extra work in the
58-us fused edge is exposed.  An exact-repair mechanism must fit that actual
schedule; the 9-us difference between 171 and the 180-us gate remains the
strict end-to-end allowance.

### Candidate gate: compressed compensated FP32 state

Registry/source search found no earlier compensated-float, float-expansion, or
quantized-residual transform experiment.  This differs from the rejected full
M31 sidecar and FP64-local-operation experiments: retain the missing mantissa
information as a narrow per-value residual throughout the FP32 transform,
rather than recomputing a complete residue field or making only the final
operation wider.  The objective is a custom FP32-plus-small-residual plane that
removes the approximately 87 wrong quotient directions per iteration while the
extra bottom-half arithmetic is hidden behind M61.

The first gate is a production-sized radix tile, including decomposition,
layout, both output components, and stores.  A Tensor Core version is relevant
because a radix-16 DFT is a dense 16-by-16 operation and published mixed-
precision work decomposes FP32 operands into multiple low-precision pieces to
recover FP32 accuracy.  This is not the already-rejected exact INT8 modular
Tensor tile: it has no four-limb modular product or reduction.  It is rejected
before integration unless it both materially improves accuracy over the SIMT
FP32 radix and has enough performance margin to pay for the remaining global
stages.  A tile win is only a gate, not a projected end-to-end speedup.

#### Tensor Core radix result

[`src/cuda/fp_compensated_tensor_bench.cu`](src/cuda/fp_compensated_tensor_bench.cu)
implements the gate over 2,097,152 complex values.  One tile is sixteen
independent radix-16 complex DFTs.  The timed compensated kernel includes TF32
high/low decomposition of both input and twiddle, high-high and both cross
products, layout, a `float2` high output, a scaled `half2` residual output, and
all loads/stores.  The ordinary comparator is a correct four-stage SIMT FP32
radix with the same input/output transform and bit-reversed order.

Against a host FP64 DFT over 32 checked tiles:

| Radix-16 path | RMS component error | Maximum complex error | Median at 2M values | Registers | Shared |
|---|---:|---:|---:|---:|---:|
| SIMT FP32 | `1.406e-7` | `7.812e-7` | 0.01619 ms | 16 | 0 |
| plain Tensor TF32 | `4.648e-4` | `1.901e-3` | 0.01440 ms | 64 | 4 KiB |
| compensated Tensor TF32 | `2.850e-7` | `2.004e-6` | 0.02179 ms | 118 | 8 KiB |

The compensated Tensor path is **1.35x slower** than SIMT and about **2x less
accurate** in RMS error.  Operand decomposition repairs TF32 truncation, but the
FP32 Tensor accumulator does not expose its internal summation error; retaining
only the cross-product residual therefore cannot exceed the factorized SIMT
radix's accuracy.  The dense DFT also expands a sparse butterfly into 24 MMA
operations and raises the live fragment footprint to 118 registers.  This
Tensor mapping is rejected before full-transform integration.  It is distinct
from, but reaches the same architectural conclusion as, the exact INT8 Tensor
NTT tile: peak matrix throughput does not compensate for the dense reformulation
and representation overhead in this radix.

The same harness also tested a sparse SIMT float-expansion butterfly to separate
the residual representation from the Tensor mapping.  Error-free FP32 products
and `TwoSum`-style additions propagate a `float` high part plus scaled `half`
residual.  It uses 24 registers with no spill and improves RMS error from
`1.406e-7` to `4.911e-8` (2.86x), but takes 0.02515 ms versus 0.01606 ms for the
ordinary radix, a **56.6% slowdown**.  Scaling that overhead across the 92-us
FP32 bottom half would add roughly 52 us, far beyond both the approximately
10-us M61 overlap slack and the 9-us end-to-end allowance.  Full-transform
compensated state is rejected.  Selective compensation would need evidence that
one stage alone causes the quotient errors; the existing FP64-local-operation
controls instead found that widening isolated multiplies/FMAs does not reduce
their population.

### New lead: sparse-error folded M31 syndrome

Registry and source search found no previous folded/aliased M31 transform,
sparse quotient-error syndrome, or subset-decoding experiment.  This is not a
shortened replacement transform and not the rejected full M31 sidecar.  It uses
the already-proved parity sidecar to reduce the remaining uncertainty to about
87 wrong `+/-1` choices among approximately 31,687 parity corrections per
iteration, then computes only enough exact M31 information to locate that sparse
set.

Let `b_j = w_j*a_j` be the M31-weighted input and `d=b*b` its length-`N`
cyclic convolution, so the coefficient consumed by carry is
`c_j = w_j^-1*d_j`.  For any `L` dividing `N`, folding before convolution gives

```text
B_r = sum_t b_(r+tL)
(B*B)_r = sum_t d_(r+tL) = sum_t w_(r+tL) c_(r+tL)  (mod M31).
```

Thus one length-`L` exact M31 square supplies one syndrome per alias bin.  The
FP32+M61 carry supplies its proposed coefficients and knows which positions had
a parity correction.  A wrong direction leaves a known signed error of two
M61 multiples, so each bin needs only choose a subset of its candidate positions
whose known weighted contributions equal the exact syndrome difference.

The production carry layout makes `L` especially favorable.  A thread already
owns `NW` coefficient pairs separated by `G_W*H` pairs; for the tuned 4M shape,
choosing the corresponding word stride groups the aliases that are simultaneously
live in its registers.  It can form both the folded next-input value and the
approximate output syndrome without atomics.  Candidate density is about 0.75%,
so an eight-alias bin normally has zero or one candidate; rare multi-candidate
bins can enumerate at most 256 subsets.  Only selected wrong pairs need their
carry reconstruction repeated.

This route has a credible sub-180 budget: a folded plane with one eighth of the
full M31 data has roughly one tenth of its stage work, while its transform can
run beside the approximately 102-us M61 bottom half.  The exposed costs are
fold construction, two M31 accumulators and candidate masks in fused carry, and
occasional reconstruction.  The gates, in order, are:

1. Prove the weighted folding identity and unique subset recovery against exact
   full coefficients, including Crandall--Fagin wrap and multiply-by-three.
2. Benchmark fold plus a production-shaped short M31 square; reject if the
   complete added path cannot fit the measured 9-us end-to-end allowance.
3. Integrate only after inspecting the fused carry's final registers, spills,
   and occupancy; proxy arithmetic that omits its live state is not predictive.
4. Validate every checkpoint through 10k, then 100k and the final 1M residue.

#### Algebra/decoder gate: pass

[`src/folded_syndrome_test.cpp`](src/folded_syndrome_test.cpp) and
`make CUDA=1 folded-syndrome-test` now validate the riskiest point before GPU
integration.  The test independently constructs M31 Crandall--Fagin weights,
checks that every direct product multiplier is exactly one or two, compares the
full weighted convolution with the direct coefficient formula, and proves that
fold-before-square equals folding the full weighted square.

For the real target geometry `N=2^22`, `L=2^19`, all eight alias weights are
distinct in every one of the 524,288 bins.  Also

```text
2*M61 == -1 (mod M31),
```

so a wrong correction contributes only a signed M31 power of two.  Any two
different subsets of the eight candidates have a nonzero signed sum of at most
eight distinct powers below M31; they therefore cannot collide modulo M31.
Subset decoding is deterministic, not probabilistic.  A target-density trial
with 31,687 candidate positions and 87 actual wrong directions decoded all 87
uniquely; the maximum was three candidates in a bin and 886 bins had more than
one candidate.

A standalone production-kernel scale probe forced an exact M31-only PRP at
512K words (`52:256:4:256:202`).  Its steady profile was:

| Short-plane region | Time |
|---|---:|
| forward weight + width | 5.8 us |
| middle in | 4.7 us |
| tail square | 6.9 us |
| middle out | 4.9 us |
| inverse width | 5.1 us |
| unrelated standalone carry A+B | 8.8 us |

The unfused whole small PRP measured 47.9 us/iteration, but the relevant
sidecar transform is the 27.4-us transform portion.  Its 16.5-us middle/tail
can run next to the approximately 102-us M61 bottom half, while its two width
edges and fold are candidates for the already-running fused carry.  This does
not prove the 9-us end-to-end gate, but it is sufficiently different from a
serial 27-us addition to justify a concurrency/resource probe before invasive
integration.

#### Concurrency/resource gate: pass at process level

The concurrency premise was tested before modifying the production carry.  A
safe `FFT3261` run at exponent `120000007` measured 174.1 us/iteration over
100,000 iterations.  While that run was active, a second process continuously
ran the complete 512K-word M31-only PRP (including its unnecessary standalone
carry); the main run measured 175.1 us and the side workload measured 26.6 us.
A subsequent uncontended control measured 175.8 us.  Within ordinary run-to-run
variation, the auxiliary M31 transform consumed otherwise-idle GPU resources
and did not reduce main-transform throughput.

This is deliberately stronger than the proposed sidecar workload and establishes
only a scheduling/resource fact, not end-to-end correctness.  The next guarded
implementation stages are: produce the folded M31 input in `carryFused` and
validate it against a standalone fold; attach the 512K cyclic M31 transform;
then consume the syndrome for sparse correction.  Before each stage, search
this registry and the source tree again so that no earlier prototype is
silently repeated.

#### Fused fold-construction gate: correctness pass

`FOLD_SYNDROME=2` now makes the real `FFT3261` `carryFused` kernel form the
eight-way folded, full-exponent M31-weighted next input while its eight aliases
are live in registers.  A deliberately separate validation kernel reads a
diagnostic dump of the final carried digits, recomputes every weight directly
from the logical word index (rather than reusing carry's incremental counter),
and folds them independently.  At exponent `120000007`, all 262,144 complex
bins / 524,288 word syndromes matched after the first fused iteration.

The 10k validation run measured 239.5 us/iteration, but mode 2 writes a 32-MiB
full-word diagnostic buffer on every iteration and is not a performance result.
`FOLD_SYNDROME=1` removes that dump and validator, retaining only two M31
accumulators, weight generation, and the 2-MiB folded output.  Its isolated
timing is the next gate, before attaching a shortened transform.

#### Fused fold-construction gate: performance pass; parity path exposed

Matched 100k runs on the same binary and safe exponent measured 208.3 us for
`PARITY_CORRECT=1` alone and 208.9 us with `FOLD_SYNDROME=1`.  The production
fold therefore costs only about **0.6 us/iteration**; its M31 rotations,
accumulation, and 2-MiB write do not reproduce the failed giant-fusion
register cliff.

The apparent gap from the earlier 174--176 us control is instead almost
entirely the existing parity-correction implementation.  It computes the
diagonal-square parity source mapping and performs a scattered parity load
inside each `weightAndCarryPairSloppy` call, extending the fused carry's
critical live range.  The folded sidecar requires this information but does
not require it to be produced there.  Before integrating the short NTT, move
expected-parity preparation to an independently schedulable kernel or write
the next iteration's expected parity directly in the final carry loop, then
measure whether the raw ~174-us core is recovered.

#### Prepared-parity experiment: partial speedup, still above gate

`PARITY_PREPARED=1` now permutes the exact square parity on a dedicated CUDA
queue while the FP32/M61 bottom half runs.  Both parity buffers use the
line-major physical carry layout, so carry's lookup is a direct coalesced load
instead of an exponent multiply, branch, address permutation, and scattered
load in each reconstruction.  The 10k residue at exponent `120000007` remained
`44c955af3fe955d9`, and the carry-ready kernel uses 12 registers.

This improved the 10k non-profiled parity path from about 208 us to 197.2 us.
The combined prepared-parity plus fused-fold run measured **196.1 us over
100k**, with the expected residue `e212dbb6397498a6`.  It is a real roughly
12-us recovery, but misses 180 us.  Profiling also shows why: while carry falls
to a direct lookup, the continuously overlapping permutation competes with the
main path (`tailSquare` rose from 25.8 to 34.0 us and `fftMiddleIn` from 15.6
to 17.6 us in profiling mode).  The next non-duplicate variant therefore
eliminates the extra kernel: final carry scatters the selected input-word
parity directly into the next iteration's carry-ready physical position.

The direct-scatter follow-up (`PARITY_PREPARED=2`) was exact at 10k but only
reached 195.5 us.  Moving the parity load inside the rare residual branch
(`PARITY_LAZY=1`) also failed to help: passing the sidecar pointer/index through
the already large inlined reconstruction raised normal/ROE carry spill from
0/8 to 8/24 bytes and measured 196.5 us.  This specific lazy-interface form is
rejected.  The useful remaining distinction is to retain a naturally packed,
coalesced parity output and derive/load it only for sparse candidates, without
the direct scatter's two inverse mappings per output pair.

#### Parity-free ternary folded decoder: algebraic gate pass

The more important conclusion is that parity need not be retained at all.
Instead of first using parity to choose an adjacent quotient and then using the
folded syndrome to flip the rare wrong directions, treat every residual-risk
position as having raw quotient delta `{-1,0,+1}`.  Its M31 syndrome
contribution is `delta*M61*w_j`; because `M61 mod M31 = 2^30-1`, these are
still signed power-of-two terms.

The expanded [`src/folded_syndrome_test.cpp`](src/folded_syndrome_test.cpp)
enumerates the real target's 248 distinct eight-alias weight patterns.  For
every nonempty candidate subset through all eight aliases, it enumerates all
ternary assignments.  All **63,240** candidate-subset assignment spaces were
injective modulo M31: zero ambiguous syndromes, including the full `3^8`
case.  Therefore one folded M31 value uniquely recovers whether each candidate
needs `-1`, zero, or `+1`; the 19--34 us parity machinery is mathematically
redundant.

The performance design is now based on the raw 174--176 us `FFT3261` core plus
the measured ~0.6-us fused fold.  Carry need emit only a sparse candidate mask
and enough signed-residual metadata for validation; the shortened M31 transform
and ternary decoder can overlap the main M61 path.  Before enabling correction,
an exact M31 oracle must confirm that every target-exponent quotient error is
in `{-1,+1}` and lies inside the chosen residual-candidate window.

#### Parity-free full-bin decoder: strengthened gate pass

The target-exponent exact `FFT323161` oracle invalidated the narrow ternary
assumption but enabled a stronger design.  Across 8,998 fused iterations at
`136279841`, approximately 35,741 raw quotient errors per iteration were
outside a `|residual|>=0.49` candidate window, so residual thresholding is not
usable.  About 0.086 coefficients/iteration had `|delta|>1` (490 total), but a
separate gate found **zero** cases with `|delta|>2`.  The required per-alias
alphabet is therefore `{-2,-1,0,+1,+2}`.

The host decoder now checks this full alphabet, not only observed sparse
subsets.  For each of the 248 real weight patterns it enumerates all `5^8 =
390,625` eight-alias error vectors.  Every mapping from the vector to the
folded M31 syndrome is injective: **zero ambiguous patterns**.  Thus the exact
folded value alone identifies all eight raw FP32/M61 quotient errors.  No
parity state, candidate mask, or residual metadata is required.

Efficient decoding need not search all `5^8` cases: most bins have zero error,
then single-error (32 possibilities) and pair-error (448 possibilities) tables
cover nearly all nonzero bins; a complete precomputed decoder or bounded
fallback handles the rest.  The GPU integration can now use the raw 174--176
us core plus the 0.6-us fold, eliminating every measured parity variant.

#### Existing-kernel folded transform: overlap implementation no-go

The registry and source were searched again before this implementation.  No
earlier experiment connected the fused eighth-length fold to a real shortened
M31 transform; the earlier full M31 sidecar and generic multi-stream tests have
different data volume and dependency structure.

`FOLD_TRANSFORM=1` now attaches a `512K`-word cyclic M31 square using the
production `52:256:4:256:202` kernels on an independent CUDA stream.  The first
width edge is a custom `foldP` because carry already emitted full-exponent
weighted values; carry also writes directly in its physical input layout, so no
transpose/copy kernel is present.  The generated kernels are well behaved:
`foldP/middle-in/tail/middle-out/width` use `39/32/48/34/40` registers and none
spills.  The side stream is event-chained from the producing carry and back to
the next carry, allowing all of the next FP32/M61 bottom half as overlap time.

This implementation nevertheless misses its performance gate.  In a matched
100k safe-exponent run in the current thermal state, the raw path measured
**193.1 us/iteration** and the folded-transform path measured **211.6 us**, both
ending at the expected residue `e212dbb6397498a6`.  Earlier matched 10k runs
were 188.8 us raw, 192.8 us with fold only, and 206.5 us with the transform.
CUDA's highest stream priority produced 204.9 us in a 10k check and did not
materially alter the result.

Nsight Systems establishes why.  In a representative steady cycle the main
carry occupied about 58--60 us.  The side `foldP` then began around +86 us and
the complete side chain ended around +189 us; the next carry began around
+190 us.  The shortened kernels themselves took roughly `9 + 12.5 + 13 + 3.3
+ 3.3 = 41 us` under contention, versus about 27.4 us in isolation.  They do
overlap most of the main bottom half, but compete with the simultaneously
running M61 kernels and expose roughly 9--19 us at the dependency boundary.
The earlier two-process concurrency result therefore did not predict a
dependency-closed single-iteration schedule: an unbounded independent process
can lag or be time-sliced, whereas every syndrome here must finish before its
consumer can overwrite the fold buffer.

This rejects the **unchanged production-kernel / generic stream** realization,
not the parity-free syndrome.  The next implementation must change at least one
architectural fact: fuse or shrink the first edge, use a carry-native short-NTT
factorization (for example a 64-point edge rather than repacking into 256-point
lines), or create a scheduler that explicitly places M31 work in M61-free
regions.  Repeating stream priority, a second generic queue, or the same
`256:4:256` chain is not justified.

#### Fine-grained side-stream scheduling: dependency removed, contention remains

The coarse event above protected the folded buffer until the complete side NTT
finished, although only `foldP` reads that buffer.  The implementation now
records a separate input-consumed event immediately after `foldP`; the following
main carry waits only for this event.  The remainder of the shortened M31 NTT can
pipeline across iteration boundaries on its private stream.  This removes the
artificial whole-transform lifetime dependency and is relevant to either an
integrated or clean-room scheduler.

It does **not** recover the target performance.  With high CUDA stream priority,
the side transform reliably finishes before the next carry, but it lengthens the
simultaneous M61 integer kernels: representative cycles remain about 188--194 us
in the profiler.  Giving the side stream low priority is better because M61 owns
the critical integer pipeline, but a 100k safe-exponent run still measured
**207.4 us/iteration**, compared with **193.1 us** for the matched raw path in
this thermal state.  The remaining approximately 14.3 us is resource contention,
not an event wait.  Three production-supported side shapes were also equivalent
at 10k (`202.1`, `202.2`, and `202.8 us`), and disabling `MULTI_Q` made absolute
performance worse (`197.3 us` raw, `206.6 us` with side transform).

This closes the generic-stream scheduling branch more strongly than the earlier
coarse experiment: neither priority, input-lifetime pipelining, transform shape,
nor the older single-queue mode makes the unchanged short M31 chain cheap enough.
Any future side syndrome must reduce its executed integer work, fuse a genuinely
cheap edge into carry, or use hardware resources that do not contend with M61.

#### Goldilocks N/16 folded syndrome: algebra pass, transform no-go

The experiment registry showed that the earlier Goldilocks rejection concerned
the full quadratic/Hermitian PRPLL transform: fast scalar roots violate that
packing and the correct unit-norm path was slow.  A side syndrome has no such
constraint, so a distinct untried architecture was tested: fold 16 aliases into
`N/16 = 262144` scalar coefficients and perform a direct base-field transform
modulo `2^64-2^32+1`.  Its state is still only 2 MiB, and a scalar root satisfying
`theta^(2^22) = 2 (mod q)` was constructed.

[`src/gold_fold_decoder_test.cpp`](src/gold_fold_decoder_test.cpp) performs an
exact meet-in-the-middle collision proof for the full error alphabet
`[-2,+2]^16`.  Equivalently it searches nonzero differences in `[-4,+4]^16`.
After correcting the test to cover actual representatives of all 16 fractional
weight classes, every class passed: the Goldilocks syndrome uniquely determines
all 16 correction values.  Thus this route passes its riskiest algebraic gate.

[`src/cuda/gold_fold_ntt_bench.cu`](src/cuda/gold_fold_ntt_bench.cu) then tests
the complete radix-16, direct scalar Goldilocks forward/inverse NTT at `2^18`
points.  Dense round-trip and sparse cyclic-square validation pass, but the
transform pair measures **0.043 ms** median/mean.  The matched existing generic
31-bit scalar full-root transform is **0.026 ms**, and the production shortened
M31 quadratic side transform is about **27.4 us** in isolation.  Goldilocks is
therefore approximately 65% slower than the 31-bit comparator before carry,
decoding, or resource contention.

This direct-scalar side use is now rejected on performance, despite its exact
decoder proof.  It should not be confused with, or repeated as, the earlier
unit-norm full-transform Goldilocks experiment.

#### Implementation-independent target (reaffirmed)

The user explicitly permits a program written from scratch.  The acceptance gate
is consequently the mathematical workload and measured throughput, not reuse of
PRPLL's classes or kernels: an exact, sustained exponent-`136279841` implementation
must reach at most **180 us/iteration over 1,000,000 iterations** and produce a
verifiable PRP/LL residue.  Existing PRPLL code remains valuable as a correctness
oracle and as the 205.5-us production reference, but architectural experiments
may replace its representation, transform schedule, carry scheme, and host
pipeline completely.  Standalone component timings are only early-abandonment
gates; they are not success until the complete exact recurrence meets this test.

#### Standalone Lucas--Lehmer recurrence: measured no-go

The registry was checked before this gate.  PRP/LL output was accepted above,
but no LL throughput measurement existed anywhere in the campaign.  LL is an
algorithm-level alternative, not a renamed host scheduler: it repeatedly
computes `s <- s^2-2 mod (2^p-1)` and does not require the PRP Gerbicz product or
proof pipeline.  PRPLL already has a dedicated `LL=1` fused carry, making it a
strictly better first gate than writing a duplicate standalone driver.

An isolated 10,000-iteration continuation initially printed 273.6 us/iteration,
and an Nsight run printed 283.1 us.  Those short figures include the synchronous
readback and 17-MiB LL checkpoint write before `IterationTimer` is reset, so
they are not the steady GPU cost.  The CUDA timeline shows consecutive
`carryFused -> middle in -> tail square -> middle out` cycles starting about
206 us apart.  Its kernel decomposition is:

```text
M31/M61 transform core: about 129.6 us (unchanged from production PRP)
dedicated LL carryFused: median 74.144 us
steady complete GPU cycle: about 206 us
```

The LL kernel uses 96 registers and has no separate width passes in steady
state; the width transforms and subtract-two normalization are already fused.
For comparison, the clean production PRP profile records about 70.5 us for its
fused edge.  LL therefore removes infrequent verification/proof work but does
not remove any transform and does not make carry cheaper.  CUDA graphs can save
only the already-measured roughly 1--2 us of launch overhead, far short of the
remaining gap.

Decision: **reject LL, including a clean-room LL driver using the same exact
M31/M61 representation, as a route to 180 us.**  A standalone implementation
would still need a materially faster multiplication representation; changing
only the recurrence, checkpointing, or scheduler cannot clear the gate.

#### Ballot-packed parity: real speedup, still not an exact target solution

The registry/source check found that existing “packed parity” stored two bits in
one **32-bit word**, not a bit-packed representation.  `PARITY_PACKED=1` now uses
CUDA warp ballots to store natural x/y parity as two bitplanes and prepared
expected parity as one bitplane.  Active parity state falls from 16 MiB to about
0.75 MiB (the two ping-pong natural buffers plus expected buffer allocate about
1.25 MiB total).  The preparation mapping is mathematically unchanged.

At exponent `120000007`, the first 10k run produced the known
`44c955af3fe955d9` residue.  A fresh 100k run produced the known
`e212dbb6397498a6` residue at **186.0 us/iteration**, versus **177.0 us** for a
matched raw FP32+M61 control.  Reordering the ballot output after the next-state
weight calculations shortened the fused carry live range: normal/ROE local
spill fell from `24/40` to `8/24` bytes per thread (raw is `0/8`).  The optimized
fresh 100k result is **185.4 us/iteration**.  `parityInit`/`parityPrepare` use
12/10 registers with no spill.  Giving the preparation stream low priority
regressed a 100k run from 186.0 to 187.2 us and was reverted.

At the actual exponent, the compact path reached **178.9 us** in its steady
2,000-iteration retry block, which demonstrates that its speed is relevant to
the 180-us gate.  It was nevertheless deterministically incorrect, producing
`ceb5b77e7fcd5f96` and an ROE failure instead of the exact production checkpoint
`05d6515c416b83e2`.  This agrees with the exact-oracle result: parity identifies
the quotient class, but about 87 coefficient pairs/iteration still choose the
wrong signed adjacent representative.  Packed parity is therefore a genuine
component speedup and useful building block, not an end-to-end solution.

#### Parity-assisted N/16 M31 syndrome: algebra pass, performance no-go

Before implementation, the registry confirmed that only the parity-free N/8
side transform had been connected to production.  The older parity-assisted
N/8 idea was conceptual; an N/16 implementation and full decoder proof had not
been attempted.

The extended [`src/folded_syndrome_test.cpp`](src/folded_syndrome_test.cpp)
proves the risky decoder point exactly.  For all **496** target weight patterns,
all `2^16 = 65,536` subsets of the sixteen aliases have distinct M31 sums.  The
result is invariant under the known sign of each candidate contribution, so a
single N/16 syndrome deterministically identifies every possible subset of
wrong parity corrections—there is no density or probabilistic assumption.

`FOLD_FACTOR=16` reuses carry's eight-alias partial fold, adds paired bins on the
side stream, transposes directly into a `256:2:256` layout, and executes a
256K-word M31 square.  Its compact/edge/middle/tail kernels use
`14/39/26/48/26/40` registers and no spill.  The safe 10k residue remains exact,
but packed parity plus this diagnostic side transform measures **202.8 us**,
versus **183.1 us** for a nearby packed-only continuation: approximately
19.7 us is exposed.  As with N/8, the loss is integer-pipeline contention, not
bad kernel resources or a coarse buffer-lifetime wait.

This rejects the unchanged production-kernel N/16 sidecar.  Even though its
decoder is exact, it consumes far more than the approximately 1-us target slack
seen in the real-exponent packed-parity run.  A future use would require a
fundamentally cheaper syndrome engine or an independently faster base transform;
simply shortening N/8 to N/16 is insufficient.

#### Correction-only 24-bit Riesel field: algebra pass, arithmetic no-go

The modulo-8 quotient observation opened a new clean-room formulation.  Given
an exact coefficient residue `n61` and `C = n61 + k*M61`, any independent small
modulus determines `k`; the conservative coefficient range makes the entire
quotient much smaller than a suitable correction modulus, not merely the
observed `[-2,+2]` error.  This could replace both parity and a sparse decoder.
The registry search confirmed that prior Riesel fields were all selected for
roughly 90-bit CRT capacity; no correction-only field had been tested.

The smallest useful candidate found is

```text
qC = 14,680,063 = 7*2^21 - 1  (< 2^24).
```

It is prime, `qC+1` contains the required order-`2^21` norm-one subgroup, and
the extended [`src/riesel_algebra_test.cpp`](src/riesel_algebra_test.cpp) finds
`root=(5118340,14210677)` and `theta=10048878`.  Root order, norm, the
`theta^(2^22)=2` weighting identity, and an independent Hermitian weighted
convolution all pass.

Two exact arithmetic engines were then added to
[`src/cuda/riesel_lazy_tile_bench.cu`](src/cuda/riesel_lazy_tile_bench.cu):

- Balanced FP32 residues use `fma(a,b,-high)` to retain the exact product error.
  Because centered inputs are below `qC/2`, every intermediate remainder stays
  within the exactly represented FP32 integer range.  All 2,097,152 quadratic
  values match Montgomery after one and four radix-8/square rounds.
- A prime-specific Montgomery radix `R=2^21` exploits `qC=7R-1`.  REDC becomes
  `high + 7*low` plus a bounded normalization, removing generic Montgomery's
  second full product.  It also matches the independent `R=2^32` result over
  all values and rounds.

Performance rejects both as speedup engines:

| qC quadratic tile, 4 repeated rounds | Median |
|---|---:|
| Generic Montgomery | 0.016 ms |
| Harvey lazy Montgomery | 0.016 ms |
| Balanced exact FP32/FMA | 0.020 ms |
| One-product `R=2^21` REDC | 0.016 ms |

The FP32 route is about 25% slower, while the special reducer is tied rather
than materially faster.  A one-round tile is memory-bound at about 0.010 ms in
all cases.  Thus reducing the modulus from ~30 bits to 24 bits does not reduce
the production-shaped GF(q^2) work; it still needs generic roots and a complete
extra transform.  Adding it to the 178.9-us approximate target path cannot fit
the roughly 1-us slack, and replacing M31 in the 205.5-us exact path would trade
the cheapest special field for a generic one.  The correction-field architecture
is rejected before end-to-end integration under the campaign's lower-bound rule.

##### Reopened hardware gate: three-limb Tensor qC sidecar

The rejection above covers CUDA-core Montgomery, Harvey, one-product REDC, and
balanced-FP32 arithmetic.  The registry has an exact Tensor Core experiment only
for a *scalar four-byte* 31-bit field; it has no quadratic three-byte correction
field or concurrent Tensor sidecar.  This is a material distinction rather than
a repeat.  For `qC < 2^24`, a base-field 16x16 matrix product needs nine INT8
MMAs instead of sixteen.  A quadratic-field matrix can use Karatsuba's three
base products, for 27 MMAs per radix tile; the corresponding old four-limb
complex construction would require 48.

The first gate is an exact standalone radix-16 comparison against the sparse
SIMT quadratic butterfly, including byte packing, modular reconstruction, and
stores.  Only if the Tensor tile wins substantially is a full side transform
or concurrency experiment justified.  This route targets the otherwise idle
dedicated Tensor pipeline, but it still risks dense-matrix excess work, register
pressure, and shared-memory synchronization seen in the earlier scalar test.

[`src/cuda/q24_tensor_bench.cu`](src/cuda/q24_tensor_bench.cu) now implements
that gate.  It constructs the order-16 root from the independently validated
order-`2^21` `GF(qC^2)` root, performs three Karatsuba matrix products with nine
INT8 MMAs each, and includes three-byte packing, shared partials, exact modular
reconstruction, and output stores.  The comparator is the actual sparse
four-stage radix-16 butterfly using Montgomery quadratic multiplication.  All
2,097,152 outputs match exactly after conversion from Montgomery form.

```text
three-limb Tensor GF(qC^2):  0.036256 ms
sparse SIMT GF(qC^2):        0.019872 ms
Tensor / sparse:             1.8245
```

ptxas reports 56 registers for Tensor and 16 for SIMT, with no stack or spills.
The 44% reduction in byte-limb MMAs is not enough: dense matrix work, packing,
shared synchronization, and modular recombination still dominate.  A complete
forward/inverse side transform would multiply this already slower stage and
cannot finish inside the FP32+M61 iteration merely by occupying dedicated MMA
units.  Decision: **reject the three-limb Tensor correction sidecar at the tile
gate**; do not implement its full NTT or concurrency path.

#### Sparse quotient repair fails the required 140--150M generality gate

The compact parity and folded-syndrome experiments above were motivated by the
very sparse quotient errors observed at exponent `136279841`. Before investing
in a more elaborate decoder, the exact three-field diagnostic path was run at
the upper end of the required exponent range with `QUOTIENT_STATS=1`, retaining
the production 4M transform shape. The requested exponents were adjusted to
valid nearby work values where noted:

| Requested exponent | Tested exponent | Maximum observed quotient error | Mean absolute scale |
|---:|---:|---:|---:|
| `140000003` | `139999991` | 7 | about 5.46 |
| `145000003` | `144999991` | 27--28 | about 21.9 |
| `150000001` | `150000001` | millions (4.77M--7.34M in short runs) | about 4.27M |

Thus the p136 `{-1,0,+1}`/parity phenomenon is a precision-cliff artifact, not a
general architectural property. Above p136 the errors rapidly become dense,
and by p150 a finite sparse alphabet or a few modular moments cannot repair
them. This closes sparse parity, subset syndrome, and multi-moment correction
as a credible route for the required 140--150M range, independently of their
p136 timing. The ballot-packed path remains a useful measured component, but
not a generally exact engine.

#### Exact production critical path and the limited TMA opportunity

Hardware performance counters are unavailable on this host
(`ERR_NVGPUCTRPERM`), but an exact M31/M61 Nsight Systems trace gives the kernel
timeline directly. Representative medians are:

| Production kernel | Median |
|---|---:|
| fused carry | 73.792 us |
| M61 tail/square | 56.128 us |
| M61 middle out | 28.832 us |
| M61 middle in | 26.976 us |
| M31 tail/square | 33.280 us |
| M31 middle out | 17.728 us |
| M31 middle in | 14.592 us |

The M31 work is substantially hidden by the M61 critical path. The approximately
110-us M61 transform chain plus approximately 74-us carry identifies where a
205.5-to-180-us design must save time; optimizing hidden M31 work alone cannot
do it. The observed payload/time ratios are already L2-class (for example,
roughly 64 MiB in 56 us for the tail and a read+write-scale 64 MiB in about
27 us for a middle kernel), while the special M61 path has high integer work.
Without counters there is no evidence for a 12% full-iteration gain from merely
replacing its existing cooperative loads with `cp.async`/TMA. TMA remains a
possible implementation tool for a new layout, but the unchanged-transform TMA
rewrite is low priority and must not be repeated as an assumed architectural
speedup.

##### Reopened fusion gate: M61-only middle/height resident tile

The registry has no field-only `middle-in -> tail -> middle-out` experiment.
The exact compact no-go below combines M31, M61, CRT/carry, and both directions
in a 512-thread, 96-KiB block; the DSM decision above rejects replacing short
carry polls, not transform residency.  A production-shaped M61 tile contains
`8*512` quadratic values, exactly 64 KiB, and can therefore retain both middle
and height transforms plus the pointwise square in one Blackwell block.

The first gate is a standalone exact 2D M61 transform comparing three kernels
(middle forward, height forward/square/inverse, middle inverse) with one
64-KiB resident kernel.  Both paths must execute identical quadratic arithmetic
and match every output.  This deliberately measures the maximum benefit of
removing the two intermediate global round trips before attempting the much
harder production transpose/twiddle mapping.  A gain substantially below the
roughly 25-us target closes this route; a large gain would justify deriving the
actual stripe mapping or a thread-block-cluster form.

Implemented [`src/cuda/m61_resident_tile_bench.cu`](src/cuda/m61_resident_tile_bench.cu)
and the `m61-resident-tile-bench` target.  It transforms 2,097,152 canonical
`GF(M61^2)` values in 512 independent `8x512` tiles.  Random forward/inverse
round trips, a sparse two-dimensional cyclic square, and the complete timed
outputs all agree exactly.  The resident kernel uses 82 registers, one barrier,
64 KiB dynamic shared memory, and has no stack or spills; the control uses
separate 96-register middle and 40-register height kernels.

Five 21-sample runs measured the following medians:

| Run | Three-kernel control | 64-KiB resident | Resident/control |
|---:|---:|---:|---:|
| 1 | 350.624 us | 360.160 us | 1.0272 |
| 2 | 345.664 us | 352.320 us | 1.0193 |
| 3 | 344.416 us | 352.576 us | 1.0237 |
| 4 | 344.480 us | 353.312 us | 1.0256 |
| 5 | 345.408 us | 353.728 us | 1.0241 |

Decision: **reject M61 field-only resident fusion.** It loses 6.7--9.5 us
instead of recovering the required roughly 25 us, even in the favorable proxy
that omits the production transpose/twiddle complication.  The experiment also
rules out register spilling as the explanation: occupancy/resource residency
and the low cost of the already cache-resident intermediate traffic are enough
to erase the saved global boundaries.  Do not derive the harder production
stripe mapping from this design unless a different representation reduces the
tile below 64 KiB or permits more than one resident block per SM.

##### New clean-room gate: 3M M61 times the M31 prime-power ring

The registry/source audit found no prime-power, Hensel, p-adic, dual-M31, or
`M31^2`-modulus experiment.  This is not the existing `GF(M31^2)` notation,
which means a two-component extension field over `M31`; the proposed coefficient
modulus is the integer prime power `M31*M31`, and its packed transform lives in
the quadratic ring `(Z/M31^2 Z)[i]`.

Use `3*2^20 = 3,145,728` real words and reconstruct coefficients modulo

```text
(2^61-1) * (2^31-1)^2,  log2(product) just below 123 bits.
```

The power-of-two DGT roots lift into the unramified quadratic ring, and the
odd radix-3 factor lives in the base-unit subgroup.  As in the independently
rejected M31 x M89 proposal, perform radix 3 in the base ring around packed
power-of-two subtransforms.  At exponent 150M, centered digits have magnitude
below `2^47`; the same conservative square-times-three estimate is about
`2^118.17`, leaving roughly 3.8 bits below the balanced 122-bit CRT limit.
This is a credible capacity path for the requested 140--150M range, subject to
a formal carry bound before integration.

The representation uses one 16-byte quadratic M61 value plus one 16-byte
quadratic prime-power value for each of 1.5M pairs, retaining the current
approximately 48-MiB state.  Its possible benefit is architectural rather than
a byte-count win: the exposed M61 population is 25% shorter, while the second
modulus avoids a generic Riesel field and has especially cheap base-`M31`
digits.  For `p=M31`, write `a=a0+p*a1`.  A scalar ring product is

```text
t  = a0*b0
c0 = t mod p
c1 = floor(t/p) + a0*b1 + a1*b0  (mod p).
```

Because `p=2^31-1`, both the remainder and `floor(t/p)` follow from one split
at bit 31 plus one bounded correction.  The first gate is an exact GPU
quadratic-multiply chain at the actual 3M population versus production M61 at
4M.  Reject before root lifting or mixed-radix code unless it is at most about
1.35x the M61 architecture-sized time; passing the arithmetic gate still does
not predict an end-to-end win.

[`src/cuda/m31_prime_power_bench.cu`](src/cuda/m31_prime_power_bench.cu) and
`make m31-prime-power-bench` implement that gate.  The device p-adic result is
checked against an independent host `unsigned __int128` product modulo
`M31^2` after every accumulated timing round; the M61 comparator has its own
independent host oracle.  The ring and M61 kernels use 30 and 32 registers,
respectively, with no stack or spills.  At their proposed architecture sizes:

| Quadratic multiply chain | 4M-word M61 population | 3M-word `M31^2` population | Ring/M61 |
|---:|---:|---:|---:|
| 1 | 15.872 us | 13.504 us | 0.851 |
| 2 | 18.112 us | 16.384 us | 0.905 |
| 4 | 28.128 us | 27.296 us | 0.970 |
| 8 | 48.672 us | 48.480 us | 0.996 |
| 16 | 91.456 us | 89.664 us | 0.980 |
| 32 | 175.392 us | 174.176 us | 0.993 |
| 64 | 343.424 us | 343.456 us | 1.000 |

The risky arithmetic gate therefore **passes**, including at arithmetic-dense
chain lengths.  This is not yet a speedup: the candidate must execute both a
shortened M61 transform and this new heavy transform, whereas production's
small M31 field is largely hidden.  The next early-abandonment gate is a real
lifted-root radix tile containing both proposed fields at 1.5M values, compared
with the current M31+M61 tile at 2M.  It must include the larger live state and
pointwise square before any 3M carry or mixed-radix integration is justified.

[`src/cuda/m31_prime_power_tile_bench.cu`](src/cuda/m31_prime_power_tile_bench.cu)
implements that next gate.  It Teichmuller-lifts the known order-`2^32` M31
quadratic generator into `(Z/M31^2 Z)[i]`, derives an order-16 root, and checks
`root^8=-1`, `root^16=1`.  Both the candidate and current control execute an
exact radix-16 forward transform, pointwise quadratic square, and inverse
transform.  The first candidate tile is also compared coefficient-for-
coefficient with an independent direct cyclic convolution.  All root and
square checks pass.

Five fresh 21-sample comparisons were tightly grouped:

| Exact architecture-sized radix-16 tile | Median range |
|---|---:|
| current 4M-word M31+M61 | 111.776--113.312 us |
| candidate 3M-word `M31^2`+M61 | 138.144--140.192 us |
| candidate/current | **1.2347--1.2406** |

Both kernels have no stack or spills and use only 36/38 registers, so this is
not another occupancy cliff.  Static SASS is 2,576 instructions for current
and 3,840 for the candidate; the p-adic path materially increases wide IMAD,
shift, compare, select, and shuffle work.  The earlier scalar-chain pass was
real but incomplete: one shortened prime-power field tied one full M61 field,
while the actual candidate must also execute its own shortened M61 field.

Decision: **reject the 3M M61 times M31-prime-power architecture at the exact
tile gate.** It is already about 24% slower before the radix-3 pass, 123-bit
CRT, carry, or weighting.  Those missing operations can only widen the gap,
so neither PRPLL integration nor a clean-room end-to-end driver is justified.
Keep both benchmarks as evidence; do not repeat the proposal as “two cheap
M31 planes,” because exact lifting requires the measured p-adic product carry.

##### New carry gate: quotient-domain normalization without `i96`

The registry/source audit found FP32-estimated CRT, block transfer functions,
compact carry bridges, and several wider CRT implementations, but no path that
eliminates the 96-bit coefficient representation from the existing M31+M61
carry itself.  This is exact algebra, not approximate quotient estimation.
After inverse weighting, production already obtains a balanced M61 quotient
`q` and an M31 residue `r` such that

```text
x = q*(2^31-1) + r.
```

After the PRP multiply-by-three and incoming carry, write

```text
A = 3*q
s = 3*r + incoming_carry - A
x = A*2^31 + s.
```

For a 32- or 33-bit output word, let `d=bits-31`, split
`A = hi*2^d + rem`, and set `t=rem*2^31+s`.  The centered output and next carry
are exactly

```text
word       = centered_low_bits(t, bits)
next_carry = hi + (t-word)/2^bits.
```

`rem` has at most two bits and `t` remains signed-64-bit throughout the current
coefficient/carry bounds.  This removes three-word construction, add-with-carry
chains, and 96-bit extraction while preserving the same CRT quotient.  The
first gate is a production-population CUDA comparison with eight coefficient
pairs per thread and chained carries.  Every digit and final carry must match
the current `i96` formulation; integration is justified only by a material
component gain because the complete iteration needs roughly 25 us.

[`src/cuda/m31_m61_direct_carry_bench.cu`](src/cuda/m31_m61_direct_carry_bench.cu)
and `make m31-m61-direct-carry-bench` implement the full 4M-word gate.  Each
thread reconstructs and normalizes eight coefficient pairs with chained signed
carries and the PRP multiply-by-three.  Random coefficients span the relevant
57-bit balanced CRT quotient range.  All 4,194,304 output words and all 262,144
final thread carries match the current three-word `i96` algorithm exactly.
Both kernels use 40 registers with no stack or spills.

Five fresh 31-sample runs measured:

| Carry arithmetic | Median range |
|---|---:|
| current `i96` | 16.384--16.928 us |
| direct quotient-domain | 19.488--19.776 us |
| direct/current | **1.163--1.198** |

The direct kernel has slightly fewer static SASS instructions (1,488 versus
1,536), but its signed quotient/remainder construction creates a longer serial
shift-and-add dependency path.  Blackwell executes the incumbent PTX
`add.cc/addc` chains efficiently enough that deleting the nominal 96-bit object
does not delete the expensive work.

Decision: **reject quotient-domain direct carry.** It loses 2.8--3.3 us in a
favorable standalone core instead of recovering any of the required 25 us.
Do not integrate it into `carryFused`; the current `i96` representation is not
the architectural bottleneck suggested by source-level operation width.

#### Exhaustive unsigned-32 Proth root search: no candidate

The earlier Proth search covered only primes below `2^31`. It has now been
extended over all unsigned 32-bit candidates `q = k*2^s + 1`, `s >= 22`,
`q > 2^31`, with the exact base-field condition
`2^((q-1)/2^22) = 1 (mod q)`. No prime passes. Consequently there is no
missing unsigned-32 direct scalar 4M field whose convolution-friendly root is
generated by 2; this closes the natural completion of the prior Proth search.

#### New untried clean-room candidate: mixed-length M31 x M89

The source and experiment registry contain no M89 implementation. A distinct
exact architecture is therefore eligible for an early-abandonment test:

- use the special fields `M31 = 2^31-1` and `M89 = 2^89-1`;
- use `3*2^20 = 3,145,728` real words (1,572,864 quadratic values), with the
  radix-3 dimension performed in the base field before independent power-of-two
  DGTs; and
- fuse only the mixed-field CRT/carry boundary, while keeping each transform's
  field-local state and schedule independent.

This is not the previously rejected monolithic odd-radix quadratic transform:
a base-field radix-3 decomposition preserves conjugate packing within each
power-of-two subtransform. At p150 the balanced digit magnitude is below
`2^47`; a conservative post-radix-3 convolution estimate is about `2^118.17`,
while the balanced half-range of `M31*M89` is `2^119`. The approximately
0.83-bit margin is narrow and requires a formal carry bound before integration.
The state is still about 48 MiB (`8+24` bytes per quadratic value), but the M89
critical transform has 25% fewer values than M61 at 4M.

The risk is M89 arithmetic, not PRPLL integration. Three radix-`2^30` limbs
need six 32x32 products for a scalar Karatsuba product and then a fold at bit 89.
Allowing for the shorter transform, the M89 quadratic multiply must be no more
than roughly 1.35x the production M61 quadratic multiply to remain plausible.
The next action is therefore a standalone, independently validated GF(M89^2)
tile benchmark. If it misses that gate, the whole M31xM89 design is rejected
without adding mixed-radix or carry code. If it passes, the next gate is the
formal p150 range proof and only then a clean-room end-to-end CUDA prototype.

The from-scratch permission is material here: success does not require this
engine to fit PRPLL's current transform classes. It still must produce an exact,
independently verifiable recurrence and sustain at most 180 us for 1M iterations;
a fast arithmetic tile or projected budget alone is not a success.

##### M89 arithmetic gate result: decisive no-go

[`src/cuda/m89_limb_bench.cu`](src/cuda/m89_limb_bench.cu) implements the gate
without any PRPLL dependencies. A canonical M89 value uses three radix-`2^30`
limbs; scalar multiplication uses six 32x32 Karatsuba products, exact carry
propagation, and a split/fold at bit 89. Quadratic multiplication uses the same
three-product complex formula as the M61 comparator. Sixteen random quadratic
lanes are independently checked after every timed/warm-up round with an 89-step
host double-and-add multiplier that shares no product or reduction logic with
the CUDA implementation.

The benchmark compares the actual proposed populations: 2,097,152 GF(M61^2)
values for 4M words against 1,572,864 GF(M89^2) values for 3M words. All tested
chains pass exact validation:

| Repeated quadratic multiplies | M61 median | M89 median | M89/M61 at architecture sizes | M89/M61 per value |
|---:|---:|---:|---:|---:|
| 1 | 0.015840 ms | 0.024224 ms | 1.529 | 2.039 |
| 2 | 0.018016 ms | 0.040384 ms | 2.242 | 2.989 |
| 4 | 0.028320 ms | 0.072608 ms | 2.564 | 3.418 |
| 8 | 0.050176 ms | 0.134528 ms | 2.681 | 3.575 |
| 16 | 0.091552 ms | 0.260160 ms | 2.842 | 3.789 |
| 32 | 0.175552 ms | 0.513600 ms | 2.926 | 3.901 |

The fresh 31-sample chain-8 means are 0.050046 and 0.135213 ms, respectively,
so the result is not a median outlier. Final `sm_120` code has no stack or
spills, but M89 needs 48 registers versus M61's 32. Static SASS is 448 versus
240 instructions; M89 replaces some wide multiplies with many more folds,
shifts, masks, and predicates (79 `SHF`, 55 `LOP3`, and 42 compare instructions
versus 9, 11, and 13 for M61). This is an intrinsic multi-limb cost rather than
a compiler spill accident.

The required arithmetic gate was at most about 1.35x per value, and even the
memory-dominated one-operation architecture-sized result is already 1.53x.
When a transform tile reuses values enough for arithmetic to matter, it is
2.7--2.9x slower despite the 25% shorter population. Radix-3 work and the wider
CRT/carry would only add cost. **M31 x M89 is therefore rejected before formal
carry proof or mixed-radix integration.** The standalone permission does not
rescue this design; it is the field arithmetic, not PRPLL scaffolding, that
fails.

#### New clean-room lead: full scalar Goldilocks plus packed M31

The registry check exposed a material distinction from both earlier Goldilocks
no-gos. The old fast scalar-root tile applied one scalar root independently to
the two components and then interpreted them as a quadratic-field convolution;
that does not pack consecutive coefficients of one real transform. The valid
unit-norm quadratic Goldilocks transform fixed the algebra but was 9--16% slower
than M61. Neither experiment tested a **full-length scalar** Goldilocks NTT.
The N/16 folded-syndrome benchmark is direct scalar, but only as an independent
small side transform with a generic global-stage implementation.

A 4M scalar Goldilocks plane still occupies 32 MiB, exactly the same as the 2M
`GF(M61^2)` plane. Store its even and odd coefficients in the two words of each
16-byte pair. Let their length-N/2 transforms be `E(k)` and `O(k)`, with a full
length root `omega`. The missing radix-2 coupling can be performed only at the
pointwise square:

```text
A(k)       = E(k) + omega^k O(k)
A(k + N/2) = E(k) - omega^k O(k)
E'(k)      = E(k)^2 + (omega^2)^k O(k)^2
O'(k)      = 2 E(k) O(k)
```

The two length-N/2 inverse transforms of `E'` and `O'` are exactly the even and
odd coefficients of the full scalar cyclic square. Thus all transform twiddles
are ordinary base-field Goldilocks multiplications; only the square needs four
scalar products per pair. Goldilocks also has the required full-length
Crandall--Fagin root-of-two weight (the existing folded decoder independently
checks `theta^N = 2`). Combining it with packed M31 provides about 95 bits and
keeps the total residue state at 48 MiB.

[`src/cuda/rns31_cuda_bench.cu`](src/cuda/rns31_cuda_bench.cu) now implements
this exact Cooley--Tukey coupling as `goldScalarPairRadix16Kernel`. The DIF output
index is bit-reversed when selecting `(omega^2)^k`. For every checked tile, the
result is compared with a direct 32-point scalar cyclic convolution after
interleaving the pair's even/odd inputs; all 64 checked coefficient pairs pass.
The M31 lane is simultaneously checked against its direct quadratic
convolution.

At 2,097,152 pairs and 21 samples on `sm_120`:

| Exact combined radix-16 forward/square/inverse tile | Median | Relative |
|---|---:|---:|
| M31 + M61 production arithmetic | 0.106 ms | 1.000 |
| M31 + scalar-root Goldilocks, old invalid quadratic square | 0.087 ms | 0.822 |
| **M31 + scalar Goldilocks, exact even/odd coupling** | **0.090 ms** | **0.845** |
| M31 + unit-norm quadratic Goldilocks | 0.116 ms | 1.095 |

The exact scalar-pair tile is therefore a real **15.5% component speedup**, not
the old wrong-algebra result. It is the first new representation in this phase
to clear a margin commensurate with the 180-us target.

The user has explicitly reaffirmed that this candidate may become a standalone
CUDA PRP engine rather than an `FFT3161` mode inside PRPLL.  The acceptance gate
is unchanged: the program must execute the complete exact recurrence, cover the
140--150M exponent range without relying on the p136-only sparse-error behavior,
and sustain at most 180 us/iteration for the full 1M-iteration p136 run.  This
permission changes where the scheduler, layout, and carry may live; it does not
turn the 0.090-ms tile result into an end-to-end result.  The immediate gate is
therefore still a complete production-shaped Gold scalar-pair transform, which
can be developed and validated independently before choosing PRPLL integration
or a clean-room driver.

It has not yet passed the whole-transform gate. The reusable but generic
global-stage direct Goldilocks harness takes 0.466 ms for a correct 4M
forward/square/inverse, far above the target; at 256K it was already much slower
than PRPLL's production-shaped shortened M31 transform. That baseline pays one
global radix pass per group of stages and a full-size root table, so it does not
measure the proposed 512x8x512/cache-fused layout, but it prevents treating the
tile ratio as an end-to-end forecast. The next gate is a production-shaped
multi-dimensional scalar-pair transform with compact roots and the exact
coupling. Only if that complete transform approaches the M61 critical-chain
budget should a Goldilocks/M31 CRT and carry be implemented.

##### Standalone three-pass scalar-pair transform: exact, but too slow

The source/registry audit found no prior full scalar-pair transform with fused
global stages, so the user's clean-room permission was exercised in
[`src/cuda/gold_pair_shape_bench.cu`](src/cuda/gold_pair_shape_bench.cu).  It
implements the complete 4M scalar cyclic square as two 2M even/odd transforms:

1. one kernel fuses the upper nine forward DIF stages;
2. one 64-KiB shared-memory kernel fuses the lower twelve forward stages, the
   exact `E^2 + root^k O^2, 2EO` coupling, and the lower twelve inverse stages;
3. one kernel fuses the upper nine inverse stages and normalization.

Dense random round-trip and the independently known sparse square of
`[1,2,3,4]` both pass.  The kernels use 38--40 registers with no local-memory
spill.  Fresh 21-sample timing is:

| Exact Gold scalar-pair region | Median |
|---|---:|
| upper forward | 0.065 ms |
| fused lower forward/coupling/inverse | 0.147 ms |
| upper inverse plus scale | 0.068 ms |
| **complete transform** | **0.274 ms** |

This is about 41% faster than the generic 0.466-ms direct-Gold harness, proving
that pass fusion matters, but it is still slower than the entire 0.2055-ms
production PRP iteration before adding M31, CRT, or carry.  The specific
4096-value/64-KiB shared tile is therefore rejected.  Its bottleneck is the
0.147-ms lower kernel: one resident block serializes 24 butterfly stages around
the coupling.  This result does **not** duplicate or invalidate the smaller
0.090-ms radix-16 component lead, nor does it yet measure PRPLL's mature
512x8x512 register/shared schedule.  Further work on scalar Gold is justified
only through that smaller-tile scheduler; another large shared-memory
clean-room transform should not be repeated.

Capacity is not the limiting issue.  The existing `1:512:8:512` `FFT3161`
configuration is rated to about 165.5M at 4M words, and M31 x Gold supplies
roughly three more CRT bits than M31 x M61.  The representation therefore has a
credible exact capacity path across the requested 140--150M range if its
transform and carry can meet the time gate.

##### M31 x Gold CRT/weight lower bound: passes, but is memory-hidden

Before a full transform port, the registry confirmed that no scalar-Goldilocks
CRT/carry bridge had been timed. [`src/cuda/m31_gold_crt_bench.cu`](src/cuda/m31_gold_crt_bench.cu)
adds an independent exact gate over all 4,194,304 scalar coefficients. The CRT
has unusually favorable constants:

```text
Gold = 2^64 - 2^32 + 1 = 3 (mod M31)
k    = (r31 - rGold) / 3 (mod M31)
x    = rGold + Gold*k
```

`Gold*k` is formed as shifts and adds (`k*2^64 - k*2^32 + k`), and the result
is balanced in a 96-bit container. All 4M random signed values drawn throughout
the balanced M31*Gold range match an independent host `unsigned __int128`
oracle. The CUDA CRT kernel uses 13 registers with no stack or spills.

The weighted bridge additionally performs one inverse and one forward generic
Goldilocks multiplication per coefficient. A generated-weight variant performs
two more runtime (non-constant-folded) multiplications to advance the weight
recurrences. Fresh 31-sample results are:

| 4M-scalar bridge | Median | Mean |
|---|---:|---:|
| M31 x Gold CRT only | 0.026176 ms | 0.026654 ms |
| Current M31+M61 weighted bridge | 0.129920 ms | 0.129498 ms |
| M31+Gold, fixed generic weights | 0.129536 ms | 0.129434 ms |
| M31+Gold, generated generic weights | 0.129760 ms | 0.129206 ms |

The fixed/generated ratios to current are 0.9970 and 0.9988. Both Gold weighted
kernels use 17 registers without spills versus 18 for the current comparator.
This gate therefore does **not** reject scalar Goldilocks: its extra modular
weight arithmetic is hidden by the synthetic pass's identical large input and
output traffic, and its CRT is as cheap as hoped.

It is deliberately not counted as a carry speedup. The production fused edge
keeps several coefficients live, propagates carry, and executes width transforms;
the earlier Riesel experience showed that arithmetic hidden in a flat memory
pass can become exposed there. A complete transform remains the next gate, and
an exact carry-shaped fused kernel is still required if it passes.

##### Production-shaped scalar-pair integration: in progress, not exact yet

The registry was checked again before this implementation.  This is the first
attempt to put the scalar Goldilocks even/odd representation through PRPLL's
actual `512x8x512` transform schedule; it is distinct from the rejected 64-KiB
standalone tile, the invalid old quadratic square, and the standalone CRT
lower bound above.  The path is opt-in with `GOLD_PAIR=1`; the normal M31/M61
engine remains unchanged.  A standalone driver remains equally acceptable if
this integration cannot meet the gate.

Implemented so far:

- scalar Goldilocks arithmetic and compact `{forward root, inverse root}` trig
  tables for width, middle, and height stages;
- exact even/odd tail square and multiply coupling;
- Crandall--Fagin forward/inverse weights and exact signed M31 x Gold CRT in the
  fused and split carry paths;
- separate multiplication for trig-record chaining, because transform data
  uses `{even, odd}` while roots use `{forward, inverse}`;
- scalar frequency reversal in the double-wide tail.  PRPLL uses a forward NTT
  for the inverse direction after reversing the spectrum; the quadratic M61
  pairing previously supplied this permutation implicitly.

The CUDA kernel set compiles, has no effect without the opt-in flag, and its
on-load product check still returns the exact initial residue `3`.  It is not
yet a valid PRP engine: before explicit frequency reversal the p=2000 residue
was the stable but incorrect `90cb5c98928b388c`; adding the reversal reaches a
different path but remains wrong at `64ffef9ab847341b`, versus the production
checkpoint `05d6515c416b83e2`.  The current ROE output is also intentionally
invalid because its diagnostic still interprets exact integer CRT values as
the old M61 roundoff estimate.

No speedup is claimed from this partial state.  The specialized square path
was about 252 us/iteration before the reversal, already over the 180-us gate;
the large incorrect ROE/check overhead makes the current short-run wall rate
unsuitable for comparison.  The next gate is an opt-in first-iteration residue
probe to locate whether the remaining mismatch is in the scalar transform/tail
convention or in weights/CRT/carry.  Do not repeat root-orientation changes or
the A-times-A `tailMul` oracle: both were already tested and produced stable,
incorrect checkpoints.

##### Production-shaped scalar-pair result: exact end to end, performance no-go

The remaining correctness failure was not CRT capacity or signed reconstruction.
PRPLL obtains the inverse NTT by feeding a frequency-reversed spectrum through
the forward-root kernels. For the factored scalar index
`k = line + H*column`, negation includes a borrow: line zero maps the column to
`-column mod L`, while every nonzero line maps it to `L-1-column`. The first
explicit reversal incorrectly used the line-zero rule everywhere and treated
the self-paired `H/2` line as line zero. That version matched iterations 1--13
and first diverged at iteration 14, exactly when the growing polynomial crossed
the first 512-word factor. Correcting the borrow in both `tailSquareGF61` and
`tailMulGF61` fixed the main recurrence and the independent Gerbicz product.
This is a reusable multidimensional-NTT lesson: a flat `k -> -k` permutation
cannot be implemented as independent negation of each factor coordinate.

The opt-in `GOLD_PAIR=1` prototype is now exact at the tested checkpoints:

- every diagnostic iteration from 1 through 20 matched production;
- p=200: `1173834bf6573c6e`;
- p=2,000: `05d6515c416b83e2`;
- p=10,000: `52316d51aa52e6b7`;
- the normal block-2,000/10,000 Gerbicz checks passed.

Temporary one-iteration block sizes and the check-bypass diagnostic used to
locate the factor-boundary bug were removed after validation. The ordinary
M31/M61 path remains unchanged, and the scalar-Gold engine remains opt-in.

Clean, sequential 10k measurements (after finding and terminating a stale
diagnostic process that had invalidated an earlier pair of timings) are:

| Exact 4M `1:512:8:512:202` path | p=10k interval |
|---|---:|
| production M31/M61 control | 202.7 us/iteration |
| M31/scalar-Gold, 96-register cap | 258.7 us/iteration |
| M31/scalar-Gold, 112-register cap | 251.5 us/iteration |
| M31/scalar-Gold, 128-register cap | 250.2 us/iteration |
| M31/scalar-Gold, compiler-default registers | **249.8 us/iteration** |
| M31/scalar-Gold, split/long carry | 277.7 us/iteration |
| M31/scalar-Gold, `1:1K:4:512` | 265.9 us/iteration |

All Gold rows reproduced the expected p=2,000 and p=10,000 residues. A clean
kernel profile explains the approximately 57-us default loss:

| Profiled region | M31/M61 | M31/Gold |
|---|---:|---:|
| fused width/CRT/carry | 70.5 us | 115.6 us |
| inverse middle | 49.8 us | 64.4 us |
| tail square | 39.3 us | 38.1 us |
| forward middle | 17.2 us | 16.6 us |

At the production 96-register cap, generated code reports 96 registers and no
local memory for M31/M61 carry, versus 96 registers plus 64 bytes/thread of
local memory for Gold. Removing that spill recovers only about 9 us. The
remaining cost is the generic 64-bit Goldilocks root and Crandall--Fagin weight
arithmetic in the fused width/carry and middle stages. The exact spectrum
reversal is not the bottleneck: the Gold tail is slightly faster than the
production tail in the clean profile.

Decision: **reject the 4M M31 plus scalar-Gold architecture for the 180-us
gate.** It is exact and supplies ample 140--150M capacity, but its best measured
rate is 23% slower than the already-over-gate production code. Separating the
carry and changing the factorization both make it worse. Further register
tuning, inline PTX, or Gold reduction micro-optimization cannot plausibly close
the roughly 70-us distance from 249.8 to 180, especially because eliminating
the entire measured Gold penalty would only return to the 202.7-us production
floor. Do not repeat this integration or the rejected 64-KiB standalone Gold
tile. The user's permission for a clean-room implementation remains applicable
to genuinely different representations or algorithms, not another scheduler
for this same 4M scalar-Gold design.

##### Standalone scheduler/graph gate and false 4.5M hybrid lead

The standalone permission also allows replacing PRPLL's host scheduler, but a
clean production A/B shows that launch orchestration is not the missing 22 us.
The exact p=10k M31/M61 control measured 202.7 us/iteration with CUDA graphs
disabled and 201.0 us/iteration with the existing per-iteration graph capture
enabled; both produced `52316d51aa52e6b7`. This 1.7-us (0.8%) gain agrees with
the earlier graph results for direct RNS and M31R2. Production already places
M31 and M61 bottom-half work on separate queues, and unchanged two-worker
execution is recorded above as reducing aggregate throughput. A clean-room
persistent/task-graph scheduler using the same kernels therefore has no
credible route to the 180-us gate and should not be implemented without a new
mechanism that removes arithmetic or memory work.

A seemingly untested 4.5M FP32+M61 point was also rejected at the geometry
audit, before implementation. M61 NTT configurations require a power-of-two
middle, so the hybrid path jumps from 4M directly to 8M; an odd-middle
`512x9x512` shape is invalid. The registry already contains the relevant longer
control: exact 8M FP32+M61 takes 0.416 ms. Do not propose a 4.5M hybrid or repeat
the 8M test. The fast 4M FP32+M61 result remains below the p136 capacity bound,
and its parity/syndrome repair variants are documented no-gos above.

##### Invalid M29 field gate and the next valid Mersenne choices

The standalone scope prompted a capacity screen of a 3M
`M31*M61*(2^29-1)` representation.  The registry/source check found no prior
M29 experiment; the old similarly labelled M23 control was already documented
as composite.  On paper the proposed product would have supplied about 120
balanced bits, enough for the conservative p150 coefficient bound, and its
48-MiB state would have matched production.  The existing exact 3M tile harness
was therefore parameterized for the new modulus before any PRPLL integration.

It fails the first required algebra gate:

```text
2^29 - 1 = 536,870,911 = 233 * 1103 * 2089.
```

The host search consequently cannot find an order-16 quadratic-field root;
there is no field in which to run the proposed NTT.  No GPU timing from that
binary is meaningful.  The temporary M29 build target was removed and the M19
harness now explicitly warns that 29 is not a Mersenne-prime exponent.

The next smaller genuine Mersenne primes are M17 and M13, not M29 or M23.
Supplying the p150 range with these adds enough quadratic planes to increase
the state and transformed limb count.  The only economical three-field fallback,
`M31*M61*M19` at `7*2^19` words, also repeats the already rejected mixed-radix-7
shape: the forum/local evidence shows that its 12.5% length reduction is lost
to radix-7 work even before adding M19.  Do not implement an M29 transform or
re-propose M23 as a prime field.

##### Standalone single-ring M31-cubed arithmetic gate

After closing M29, the experiment registry was checked for `M31^3`, cubic
M31 prime powers, and three-digit p-adic rings.  No such implementation or
measurement existed.  This is architecturally distinct from the rejected 3M
`M31^2+M61` design: use one 4M transform over
`(Z/(M31^3)Z)[i]`, replacing both production fields rather than placing the
prime-power ring beside a shortened M61 transform.

The proposal has a legitimate early gate:

- `M31^3` has just under 93 modulus bits, about one bit more than `M31*M61`, so
  it preserves every range supported by the production representation;
- one quadratic value is six 32-bit digits (24 bytes), exactly the same state
  width as one M31 quadratic value plus one M61 quadratic value;
- power-of-two roots and the root-of-two weight can be Hensel-lifted because
  the transform length is invertible modulo M31; and
- a single ring would eliminate the inter-field CRT at the carry boundary.

[`src/cuda/m31_cubic_power_bench.cu`](src/cuda/m31_cubic_power_bench.cu) and
`make m31-cubic-power-bench` implement the first exact gate at the real
2,097,152 packed-value population.  A scalar is stored in three canonical
base-M31 digits.  GPU multiplication uses `2^31=M31+1` to obtain both each
remainder and the next p-adic carry without integer division.  Quadratic
multiplication uses the usual three-product complex formula.  The independent
host oracle uses `unsigned __int128` schoolbook accumulators and ordinary
division/remainder, not the GPU fold.  Every measured chain matches it exactly.

Fresh 31-sample medians are:

| Dependent quadratic multiplies | Current M31/M61 | Single M31-cubed | Cubic/current |
|---:|---:|---:|---:|
| 1 | 73.728 us | 77.248 us | 1.048x |
| 2 | 74.688 us | 78.944 us | 1.057x |
| 4 | 76.704 us | 95.904 us | 1.250x |
| 8 | 78.400 us | 175.360 us | 2.237x |
| 16 | 140.448 us | 316.480 us | 2.253x |
| 32 | 261.184 us | 597.760 us | 2.289x |

Five separate chain-8 processes tightly reproduce `2.237--2.247x`.  Ptxas
reports 40 registers for M31-cubed versus 36 for the combined M31/M61 control,
with zero stack and zero spills in both.  Final `sm_120` SASS contains 304
instructions for the cubic kernel versus 160 for the comparator, including
100 versus 41 integer adds, 35 versus 10 compares, and 24 versus 9 funnel
shifts.  The loss is therefore the exact p-adic product/carry dependency graph,
not occupancy or PRPLL orchestration.

An NTT tile performs several dependent root products plus the pointwise square
per loaded value; the arithmetic-dense gate is the relevant limit, while the
near-tie at one operation is only equal memory traffic.  A 2.24x transform
arithmetic penalty cannot be repaid by deleting CRT from the approximately
205.5-us incumbent, and the candidate's arithmetic kernel alone is already
about 175 us at chain 8 before any transform additions, shuffles, weighting,
or carry.  Decision: **reject the standalone M31-cubed architecture before
root lifting or a transform tile.**  Keep the exact benchmark as evidence; do
not confuse this with the previously measured M31-squared-plus-M61 route.

## Production M31/M61 scheduling and remaining hardware audit

Before changing the production scheduler, the experiment registry and replay
code were checked.  `MULTI_Q=1` already puts the complete M31 bottom half on the
main queue and the complete M61 bottom half on an auxiliary queue.  The earlier
side-syndrome work tested priorities only for an added shortened transform, not
for prioritizing either field of the production transform itself.  The
production priority experiment was therefore new rather than a repeat.

The existing kernel profile also gives a hard budget for any scheduler that
retains these kernels.  Median constituent times are approximately:

| Production region | Time |
|---|---:|
| M61 middle in | 26.976 us |
| M61 tail | 56.128 us |
| M61 middle out | 28.832 us |
| complete M61 bottom-half sum | 111.936 us |
| M31 bottom-half sum | 65.600 us |
| fused width/CRT/carry edge | 73.792 us |
| concurrently scheduled transform core | about 129.6 us |

Even the unattainable ideal in which all M31 work is free and the M61 kernels
have no overlap penalty gives `111.936 + 73.792 = 185.728 us`, before residual
launch/dependency costs.  Thus stream phasing alone cannot satisfy the 180-us
gate.  It could at most supply a small gain to combine with a genuinely cheaper
edge.

An opt-in test recreated the M61 auxiliary stream at CUDA's highest priority.
Clean, alternating 30k runs at exponent `136279841` measured:

| Queue policy | 30k time | iteration-30,000 residue |
|---|---:|---:|
| equal priority, control 1 | 198.5 us/iteration | `9139db3046e846d4` |
| high-priority M61, trial 1 | 203.1 us/iteration | `9139db3046e846d4` |
| high-priority M61, trial 2 | 203.3 us/iteration | `9139db3046e846d4` |
| equal priority, control 2 | 199.3 us/iteration | `9139db3046e846d4` |
| low-priority M61, diagnostic | 199.7 us/iteration | `9139db3046e846d4` |

High priority is reproducibly about 4 us slower than its nearby controls; low
priority also supplies no gain.  CUDA priority operates at scheduling/block
boundaries and sacrifices complementary occupancy rather than eliminating the
integer-pipeline contention.  The experimental hook was removed after the
measurement.  Decision: **retain the incumbent equal-priority two-stream
schedule and do not pursue explicit phase barriers as an independent speedup.**

Two related ideas were closed by audit rather than duplicated implementation:

- Production's radix-8 local operations already use the cheap Mersenne
  rotations.  Recounting only generic-root products gives roughly 896 per
  512-point dimension for the incumbent, about 912 for conventional split
  radix, and about 864 for a radix-16-by-32 factorization.  The last is only a
  3.6% nominal reduction before extra exchanges/register pressure and cannot
  recover the measured roughly 25-us gap.  This agrees with the existing
  high-radix spill evidence; there is no new split-radix integration to run.
- The RTX PRO 6000 Blackwell exposes Tensor Core paths for FP16, BF16, TF32,
  INT8, and lower-precision formats, but not FP64 Tensor Core MMA.  The official
  [RTX PRO 6000 specifications](https://www.nvidia.com/en-us/products/workstations/professional-desktop-gpus/rtx-pro-6000/),
  [Blackwell professional architecture description](https://www.nvidia.com/content/dam/en-zz/Solutions/design-visualization/quadro-product-literature/NVIDIA-RTX-Blackwell-PRO-GPU-Architecture-v1.0.pdf),
  and [PTX ISA](https://docs.nvidia.com/cuda/archive/12.9.2/parallel-thread-execution/index.html)
  provide no FP64 Tensor MMA target for `sm_120`.  The dedicated-hardware route
  therefore remains the already measured INT8 decomposition, which was 42%
  slower than sparse SIMT; an FP64 Tensor rewrite is not available on this GPU.

These results do not change the overall outcome: no exact end-to-end speedup
has been achieved.  They narrow the remaining search to representations or
algorithms that remove work from the approximately 74-us fused edge and/or the
112-us M61 critical transform, rather than another host scheduler around the
same kernels.

## 3M FP32+M31+M61 factorization audit: the `1024x3x512` lead

Before extending the earlier `512x6x512` Good--Thomas prototype, the experiment
registry above was checked for alternative 3M factorizations and for a
middle-3 implementation.  Neither `1024x3x512` nor `512x3x1024` had previously
been attempted.  This matters more than changing a launch parameter: middle 3
turns each exact plane into three independent power-of-two channels with no
remaining binary middle stage.

The integrated experimental parser, trig generation, FP32 middle-3 transform,
M31/M61 scalar radix-3, and tail channel indexing were generalized for middle
3 and middle 6.  Two alternate middle-6 geometries were first rejected on the
existing deliberately incomplete scaffold:

| 3M shape | Wrong-scaffold whole iteration |
|---|---:|
| `4:1K:6:256:202` | about 224.3 us |
| `4:256:6:1K:202` | about 206.6 us |
| prior `4:512:6:512:202` | about 185.6 us |

Middle 3 gives a qualitatively different result.  Nsight Systems timelines
were grouped from the first `fftMiddleIn*` launch through the last
`fftMiddleOut*` launch in each cycle, excluding the known-invalid fused carry.
For 5,994 steady cycles with the stock width decomposition (`NW=4`):

| Shape | Mean | p10 | Median | p90 | Minimum | Maximum |
|---|---:|---:|---:|---:|---:|---:|
| `4:512:3:1K:202` | 120.134 us | 112.704 | 122.784 | 126.432 | 110.272 | 387.968 |
| `4:1K:3:512:202` | **88.039 us** | 83.616 | **88.480** | 92.288 | 82.336 | 439.871 |

The faster shape's median constituent kernels were:

| Plane/region | Middle in | Tail | Middle out | Sum |
|---|---:|---:|---:|---:|
| M61 | 13.600 us | 46.048 us | 11.488 us | 71.136 us |
| M31 | 11.008 us | 23.264 us | 11.680 us | 45.952 us |
| FP32 | 8.000 us | 11.552 us | 13.888 us | 33.440 us |

This is about 39--41 us below the earlier 3M/production transform critical
paths.  It is not a correctness result.  The old fused width/CRT/carry kernel
assumes a power-of-two flat transform order and normalization; with middle 3
the on-load check becomes zero.  The exact Good--Thomas mapping also crosses
the width coordinate.  For natural pair index represented by `(x,m,y)`, the
physical base coordinate is `(3*x+m) mod 1024`, while the two scalar residues
come from channels `(m+2*y) mod 3` and the following channel.  A dedicated
boundary transpose/CRT/carry is therefore mandatory.

An `NW=8` edge experiment exposed two implementation facts:

- `clDefines()` initialized five option members by reference while constructing
  `KernelCompiler`, but those members were declared after the compiler and all
  kernels.  Moving them before `KernelCompiler` removes this C++ lifetime/UB
  bug; it was real but not the CUDA failure.
- Compute Sanitizer localized the CUDA sticky illegal-address error to an
  out-of-bounds shared-memory access in the hybrid `fftP` at the padded
  `NW=8, SHUFL_BYTES_W=8` layout.  The dumped PTX compiles and loads in a fresh
  CUDA context.  `LDSPAD_W=0` avoids the invalid access and is the only valid
  timing used below.

With `NW=8, LDSPAD_W=0`, the fused carry median falls from about 111.616 to
105.216 us, but the transform-core median rises from 88.480 to 92.288 us.  The
known-wrong carry-to-carry scaffold has a 199.680-us median (5,992 steady
cycles).  Thus merely resizing the legacy fused edge does not turn the lead
into a result; it recovers only a few microseconds and remains incorrect.

The external Aevum/PrMers Type-4 PFA route was also audited before borrowing
its boundary.  A temporary, out-of-tree PFA3 extension launched end to end, but
Type-4 PFA3 disagreed with Type 1 even at exponent 120,000,007 where the paired
PFA3 range is safe.  The pre-existing Type-4 PFA9 path likewise disagreed with
Type 1 at exponent 175,000,039, so it is not a correctness oracle for this
work.  Its p136 PFA3 scaffold was also around 293.5 us.  No source from that
temporary branch was copied into this workspace.

Current decision: retain `4:1K:3:512:202` as the strongest unproven lead and
implement only a dedicated exact Good--Thomas edge for it.  Kill the route if
the inverse-width + mapping + exact CRT/carry + forward-width boundary cannot
fit within roughly 91.5 us (`180 - 88.5`), or if a full-run residue fails.
Until that gate passes, the campaign outcome remains **no exact end-to-end
speedup achieved**.

#### Middle-3 exact square recurrence: passes; split boundary remains slow

Before changing the boundary, the earlier Good--Thomas carry-layout benchmark
and every recorded 3M attempt above were re-read.  The first correctness
implementation deliberately used the split `fftW -> carry -> fftP` path, not a
new fused kernel: it is slower, but makes the component-wise permutation easy
to inspect and provides an oracle for a later fused edge.

The exact physical map is now implemented for the FP32+M31+M61 hybrid.  FP32
continues to use the ordinary natural-order 3M transform.  After inverse width,
a natural packed pair `(x,m,y)` gathers its exact components from

```text
x_channel = (3*x + m) mod 1024
channel_even = (m + 2*y) mod 3
channel_odd  = channel_even + 1 mod 3.
```

The forward edge applies the inverse map, gathering the two natural words that
form one packed channel value.  This direct global gather is a correctness
scaffold, not the intended final memory architecture; the earlier benchmark
already showed that an unfused direct channel gather is a coalescing no-go.

That bring-up found two independent bugs in the provisional middle-3 path:

1. `genMiddleTrigGF31/GF61(..., middle=1, ...)` returned one zero and the
   caller merely padded it to the physical middle-3 allocation.  Consequently
   `middleMul2` multiplied every exact value by a zero trig.  The one-channel
   generator now emits the real width and width-times-height root tables, and
   middle-3 reads them at the corresponding unpadded offsets.
2. The inverse radix 3 is normalized explicitly, leaving a `2^20`-word packed
   DGT per channel.  Its inverse output scale is still `2*2^20`, not `2^20`.
   The carry now removes `2^21`.  The old value made the deterministic
   recurrence `x -> 2*x^2` (`3 -> 18 -> 648`) rather than `x -> x^2`.

With those fixes, the complete square/carry/forward recurrence is exact at the
first normal checkpoint:

```text
p=30,000,001, iteration 2,000: 5b5052bd4e11a1c5
  (matches the trusted FP64 3M control)

p=136,279,841, iteration 2,000: 05d6515c416b83e2
  (matches the production 4M M31/M61 control)
```

The target run reports `Z ~= 13,358`, so FP32-assisted lifting has ample
distance from its decision boundary at this checkpoint.  The current run is
still labelled `EE` because the independent Gerbicz product uses `tailMul`,
whose Hermitian channel pairing has not yet been generalized; the main square
residue itself is the expected value.

The direct-gather, long/split implementation warms to about **265.3 us per
iteration**, so it is not a speedup.  Before profiling representative nonzero
data, the exact square recurrence temporarily made the 88.480-us core look like
an architectural lead.  The corrected trace below supersedes that provisional
classification; no fused edge or `tailMul` work was started from it.

That provisional performance conclusion was then invalidated by profiling the
corrected, nonzero path.  This check is important enough to supersede the last
paragraph: the 88.480-us core was measured while the middle-1 exact trig tables
were zero, so M31/M61 values collapsed to zero before the data-dependent
normalization/reduction paths.  Kernel topology and SASS were present, but the
executed modular work was not representative.

An Nsight Systems trace of 5,999 corrected square cycles, paired independently
across the FP32, M31, and M61 streams, gives:

```text
corrected transform-core mean:    133.721 us
corrected transform-core median:  136.160 us
p10 / p90:                        123.328 / 139.135 us
minimum / maximum (<1 ms):        103.808 / 465.727 us
```

The same trace's split-boundary kernel medians are approximately:

```text
inverse widths, overlapped:       max(17.1, 24.6, 29.9) ~= 29.9 us
direct mapped exact carry:        44.8 us
combined forward fftP:            60.6 us
carryB:                            2.3 us
```

Even an optimistic clean-room boundary that perfectly coalesces the channel
transpose, overlaps three separate forward width kernels, and reduces the
forward edge to roughly the inverse-width cost still has a lower bound near
`29.9 + 44.8 + 29.9 = 104.6 us`.  Added to the corrected 136.2-us core, that is
about 240.8 us.  More decisively, the core alone leaves only 43.8 us under the
180-us gate, less than the measured inverse width plus carry before any forward
width is done.

Decision: **reject the 3M FP32+M31+M61 middle-3 architecture for the 180-us
gate.**  Do not implement its fused edge or generalize `tailMul` merely to
obtain a complete slow engine.  The exact p136 square residue is a useful
validation of the Good--Thomas algebra and mapping, but there is no end-to-end
speedup.  The reusable performance lesson is that modular kernels with
predicated/data-dependent normalization must be benchmarked on representative
nonzero residues; a launch-complete zero-data scaffold can understate their
cost by tens of microseconds even when the generated instruction stream looks
plausible.

## Exact FP32-limb M31 resource-offload gate

Before implementation, the registry and source tree were searched for exact
FP32 M31, float-limb M31, double-single modular arithmetic, and related variants.
The earlier FP32 experiments were either approximate coefficient transforms,
compensated approximate FFT state, or an exact 24-bit generic Riesel field.
None retained the exact production M31 residue while moving its arithmetic to
the FP32 pipelines.

The motivation was the measured schedule rather than nominal FP32 throughput.
Production's approximately 65.6-us M31 bottom half increases the concurrently
scheduled M61 critical path from an ideal 111.9 us to about 129.6 us.  Conversely,
the approximate FP32 plane overlaps M61 much more effectively.  If an exact M31
engine could use FP32 resources and stay within the roughly 112-us M61 window,
the complete `M31*M61` capacity and its validity throughout 140--150M would be
preserved without another residue field.

[`src/cuda/m31_fp_limb_bench.cu`](src/cuda/m31_fp_limb_bench.cu) implements the
risky arithmetic gate at the production population of 2,097,152 quadratic
values.  Each M31 scalar is represented by exact FP32 integer limbs of
`11+11+9` bits.  A three-limb Karatsuba product needs six FP32 products, and its
largest operand product is

```text
(2*(2^11-1))^2 = 4094^2 < 2^24,
```

so every product and integer-valued intermediate is exactly representable.
The reducer normalizes in radix `2^11`, folds at bit 31 using
`2^31 == 1 (mod M31)`, and stays in FP32 limbs through repeated quadratic
products.  Only kernel input and output are converted to the incumbent `uint2`
format.  Every value agrees with an independent integer M31 path after each
tested chain length.

Three fresh processes tightly reproduce these medians:

| Dependent quadratic multiplies | Current integer M31 | Exact FP32 limbs | FP/integer |
|---:|---:|---:|---:|
| 1 | 8.010--8.011 us | 24.898--25.162 us | 3.108--3.141x |
| 2 | 10.270--10.371 us | 41.710--41.752 us | 4.024--4.065x |
| 4 | 16.190--16.274 us | 77.245--77.298 us | 4.747--4.774x |
| 8 | 28.445--28.488 us | 146.254--146.510 us | 5.135--5.151x |
| 16 | 52.333--52.360 us | 285.016--285.197 us | 5.444--5.450x |

Ptxas reports 40 registers for the eight-round FP32 kernel versus 20 for the
integer comparator, with zero stack, local memory, or spills in either.  Final
`sm_120` SASS has about 4,640 instructions for FP32 versus 928 for integer at
eight rounds.  The FP32 path includes 660 `FFMA`, 609 `FADD`, 492 `FMUL`, and
444 `FRND.FLOOR` instructions; exact radix normalization, not conversion or
memory traffic, dominates it.

Decision: **reject exact FP32-limb M31 offload before a transform or concurrency
integration.**  At transform-like arithmetic density it already takes about
146 us for a component whose entire M61 overlap window is roughly 112 us, and a
real NTT would add butterfly additions, shuffles, and three kernel-boundary
conversions.  FP32 pipeline complementarity cannot repay a fivefold expansion
of the M31 operation graph.  The reusable lesson is that exactness below the
24-bit mantissa is not sufficient: score normalization and limb maintenance
across the complete dependency chain.

## Exact FP32 q24 times M61 replacement gate

The failed 31-bit FP32-limb experiment prompted a deliberately narrower
registry check.  The earlier correction-field work had already validated

```text
qC = 7*2^21 - 1 = 14,680,063
```

with an exact balanced-FP32/FMA reducer, but only as an *additional* sidecar to
FP32+M61.  It had not tested replacing M31 outright or the intended FP32/M61
resource overlap.  This is distinct from repeating the correction sidecar:
the candidate has exactly two residue fields, the same 48-MiB packed state as
production, and no approximate coefficient plane.

### Capacity gate: narrow empirical pass through p150

The balanced qC*M61 range has 83.807354 bits.  This is much tighter than
production M31*M61, so the existing exact FP32+M31+M61 coefficient oracle was
extended with `COEFF_RANGE_STATS=24`.  It compares every reconstructed signed
128-bit coefficient directly against

```text
floor(qC*(2^61-1)/2) = 0xdffffefffffffff900000.
```

At exponent `150000001`, the exact 4M oracle produced the trusted checkpoints
through 100,000 iterations and observed no coefficient outside this range.  A
control against the smaller `(2^23-1)*M61` range triggered immediately and
heavily, while a `2^84` magnitude threshold did not trigger.  Thus qC lies in a
real, narrow capacity interval rather than passing because the diagnostic is
insensitive.

This is an empirical PRPLL-style capacity gate, not an adversarial all-input
proof.  A successful performance path would still require a 1M p150 range run
before acceptance.  It is adequate for early performance screening and covers
the requested upper endpoint rather than relying on p136's sparse FP errors.

### Production-population overlap gate: real but insufficient gain

[`src/cuda/q24_m61_overlap_bench.cu`](src/cuda/q24_m61_overlap_bench.cu)
compares the intended resource schedule at 2,097,152 quadratic values:

- incumbent exact integer M31 overlapped with exact integer M61; and
- exact centered-FP32 qC overlapped with the same M61 kernel.

The qC product uses rounded FP32 high product plus `fma(a,b,-high)` to recover
the exact low part, estimates the quotient, and returns an exactly represented
centered integer.  All qC and M31 outputs match independent host modular complex
products at every tested chain length.  Three processes measured:

| Quadratic multiply chain | M31 isolated | q24 FP isolated | M61 isolated | M31+M61 overlap | q24+M61 overlap | q24 overlap gain |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 10.15 us | 9.59--9.64 us | 17.16--17.17 us | 78.14--78.23 us | 77.64--77.94 us | 0.20--0.56 us |
| 4 | 16.98--17.07 us | 14.89--15.01 us | 29.66--29.78 us | 88.32--88.78 us | 87.08--87.49 us | 0.95--1.69 us |
| 8 | 31.21--31.24 us | 26.54--26.59 us | 54.28--54.38 us | 120.06--120.16 us | 115.86--116.09 us | **3.98--4.30 us** |
| 16 | 58.46--58.48 us | 48.96 us | 102.52--102.64 us | 193.99--194.09 us | 185.98--186.11 us | 7.89--8.11 us |

The absolute overlapped proxy times reflect severe simultaneous-kernel
contention and are not a production-iteration forecast.  The paired delta is
the relevant gate, and it is reproducible.  At the transform-like eight-product
density it recovers only about 4 us, versus 25.5 us required to move 205.5 to
180 us.

This proxy is deliberately optimistic for qC.  It compares generic complex
products, whereas production M31 uses especially cheap radix rotations and
power-of-two weights.  A complete qC transform must instead execute generic
roots; its earlier exact FP32 radix-8/square gate was 0.020 ms versus 0.015--
0.016 ms for integer qC.  Its edge would also replace the unusually cheap
M31/M61 CRT constants and M31 shift weights with generic qC multiplication.
None of those omitted costs can create the missing 21 us.

Decision: **reject q24*M61 before transform or carry integration.**  The
capacity observation is useful—the 4M p150 requirement is only just above an
84-bit product—but using an exact 24-bit field on FP32 pipelines does not supply
enough whole-core resource complementarity.  Do not repeat this as the old
sidecar proposal; both the two-field capacity and concurrent replacement gate
are now measured.

## Blackwell resource audit, power ceiling, and sparse Tensor Core gate

Three remaining hardware assumptions were checked before another implementation.

First, the RTX PRO 6000 Blackwell Max-Q was rerun at its available 300-W power
limit.  A fresh exact 100,000-iteration production run at exponent `136279841`
used `4M 1:512:8:512:202`, reproduced residue `52775eea4730be87`, and measured
201.7 us/iteration.  The device reports a 325-W maximum, but
`nvidia-smi -pl 325` fails with `Insufficient Permissions` in this container and
the enforced limit remains 300 W.  Raising the power limit is therefore an
unavailable control experiment, not a software speedup and not evidence for the
180-us gate.

Second, the independent Blackwell execution-unit analysis in
[Heinz et al.](https://arxiv.org/pdf/2507.10789) reports that the consumer GB203
SM's FP32 units are unified INT32/FP32 units: in a cycle they execute integer or
floating-point work rather than supplying two independent full-rate pipelines.
This explains why the exact q24 FP32/M61 overlap gate recovered only about 4 us.
It also corrects the initial resource model behind FP32 offload proposals.
FP32 can still change an instruction graph, but on this GPU it is not an
independent compute reservoir hidden beside the integer NTT.

Third, the experiment registry was searched for `mma.sp`, 2:4 sparsity, and a
sparse Tensor NTT before reopening Tensor Cores.  Only dense exact INT8 and TF32
experiments existed.  The current
[PTX ISA](https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#matrix-multiply-accumulate-operation-using-mma-sp-instruction-with-sparse-matrix-a)
does provide `mma.sp.m16n8k32` with unsigned or signed INT8 inputs and 2:4
structured sparsity, so the question was reduced algebraically before writing a
kernel.

The result is a no-go for the required fields:

- Dense unsigned-INT8 `mma.m16n8k16` performs 16 useful products per output.
- Sparse unsigned-INT8 `mma.sp.m16n8k32` doubles the logical K dimension to 32
  and stores 16 A values per row.  It also performs 16 useful products per
  output; the advertised sparse throughput compensates for doubled K.
- A direct radix-16 DFT has 16 nonzero coefficients per output.  Spreading its
  columns into a 2:4 16-by-32 matrix therefore leaves the MMA instruction count
  unchanged: two `n8` operations per limb pair, exactly as for dense
  `m16n16k16`.
- An 8-point transform over a packed quadratic field becomes a 16-scalar linear
  map.  Each real or imaginary output again depends on all 16 scalar inputs, so
  it has the same count.  A Karatsuba split into three 8-by-8 base-field products
  can pair two blocks in a sparse MMA, but its ideal reduction is only 25% and
  cannot erase the already measured 42% (q31) or 82% (q24 quadratic) dense
  Tensor deficits once packing and modular reduction are retained.
- A butterfly-stage matrix has only two useful coefficients per row, much less
  than 50% density.  The fixed sparse shape fills the remaining stored A slots
  with zeros and consequently performs more dummy work than SIMT butterflies.

Decision: **do not implement ordinary 2:4 `mma.sp` as another Tensor NTT tile.**
It does not meet the ledger's remaining condition of asymptotically fewer limb
products than the dense construction, and metadata plus expanded B packing can
only add work.  A future Tensor proposal must encode multiple independently
recoverable modular products in one accumulator or use a different sparse shape;
merely embedding the NTT matrix in 2:4 storage is now a closed duplicate.

## Negacyclic root-of-minus-two audit

An untried-looking route was to replace the Crandall--Fagin condition
`theta^N = 2` with a negacyclic transform and `theta^N = -2`, possibly admitting
cheaper Proth fields.  The registry had no implementation of that exact idea.
However, a negacyclic NTT also needs a twist `psi` with `psi^N = -1`.  Therefore

```text
(theta / psi)^N = (-2) / (-1) = 2.
```

Any field supporting the proposed construction already supports the original
Nth root of 2; the negacyclic change cannot enlarge the eligible-prime set.
An exhaustive search of sub-2^32 Proth primes with enough 2-adicity found no
counterexample, as the invariant predicts.  Decision: **reject the negacyclic
prime-search route algebraically and do not build a transform for it.**

## Tensor-Core Toeplitz M61 constant-multiplication gate

The dense and 2:4 Tensor NTT rejections do not cover a different mapping aimed
at the production bottleneck.  M61 middle transforms multiply many independent
values by repeated twiddle constants.  Sixteen 61-bit values can be written as
the rows of a base-256 limb matrix, while the eight limbs of one shared twiddle
form a 16-by-16 Toeplitz matrix.  One unsigned-INT8 MMA then computes all 15
convolution columns for 16 exact 61-by-61-bit products.  Three such products
implement the production Karatsuba quadratic-root multiplication.  This mapping
does not densify the NTT itself and had not appeared in the registry.

[`src/cuda/m61_tensor_constant_bench.cu`](src/cuda/m61_tensor_constant_bench.cu)
implements the production-population gate.  Its timed Tensor path includes:

- staging sixteen input residues as eight exact base-256 limbs;
- a genuine `u8*u8->s32` MMA;
- storing all accumulator columns;
- base-256 carry propagation into the exact 122-bit product;
- folding at bit 61 using `2^61 == 1`; and
- output stores.

The quadratic kernel repeats this for the same three operands/constants used by
the incumbent complex multiply.  Every one of 2,097,152 base and quadratic
outputs matches both the scalar GPU path and an independent host `unsigned
__int128` reference.  The initial one-warp CTA was corrected before judging the
idea: the final kernel places eight independent warp tiles in each CTA, reducing
the grid from 131,072 to 16,384 blocks.  Ptxas reports 26 registers for the base
Tensor kernel and 36 for quadratic, with no stack or spills.  Final SASS contains
two `IMMA.16816.U8.U8` instructions per logical 16-by-16 MMA, confirming native
Tensor execution.

Three fresh processes measured:

| Exact operation over 2,097,152 values | Scalar M61 | Tensor Toeplitz | Tensor/scalar |
|---|---:|---:|---:|
| one shared-constant base product | 9.47--9.57 us | 26.56--28.03 us | 2.78--2.96x |
| one shared quadratic-root product | 15.68--15.81 us | 70.18--71.26 us | 4.47--4.54x |

The MMA itself is not the limiting abstraction.  Turning its 15 exact
convolution columns back into a modular scalar requires accumulator traffic,
byte carries, 128-bit assembly, and Mersenne folds for every value.  The
incumbent instead obtains the low and high halves of a 64-bit product directly
and folds them with a short dependency graph.  Even perfect root reuse cannot
remove the roughly 55-us quadratic deficit at this population.

Decision: **reject Tensor-Core Toeplitz M61 multiplication before a middle-NTT
integration.**  It is a distinct and exact use of the dedicated units, but it
is much slower than the production primitive it would replace.  Do not repeat
it as a dense-DFT or sparse-MMA proposal; the conversion from MMA columns back to
M61 is the architectural boundary.

## Production M31/M61 transform-factorization confirmation

Before changing the transform/fusion boundary, the experiment registry was
checked for the equivalent 4M factorizations.  The `1024x4x512` and
`512x4x1024` shapes had already been timed for the slower M31/Riesel and
M31/Gold architectures, and old generic tuning data also contained
`1024x4x512` and `1024x2x1024`.  There was no matched, exact production
M31/M61 comparison, so this was treated as a narrow confirmation gate rather
than a new design.

Four fresh, sequential p=136279841 runs used the same binary, CUDA settings,
4M word count, NTT fields, and FFT variant.  Every shape reproduced the exact
iteration-2,000 residue `05d6515c416b83e2` and iteration-10,000 residue
`52316d51aa52e6b7`:

| Exact 4M M31/M61 shape | 10k interval |
|---|---:|
| production `1:512:8:512:202` | **197.8 us/iteration** |
| `1:512:4:1024:202` | 219.8 us/iteration |
| `1:1024:4:512:202` | 222.2 us/iteration |
| `1:1024:2:1024:202` | 250.9 us/iteration |

These layouts change no modulus, coefficient capacity, or transform length, so
the negative result applies directly across the requested 140--150M range.
Moving one binary stage out of the middle costs about 22--24 us; moving two
costs about 53 us.  Thus the incumbent symmetric `512x8x512` decomposition is
already the favorable fusion/cache boundary for this exact architecture.

Decision: **reject further equivalent 4M factorization sweeps.**  Do not profile
or tune these slower layouts, and do not present stage redistribution as the
missing architectural win.  A successful path must remove arithmetic or state,
not merely move the same stages between production kernels.

## Near-`2^53` Riesel replacement for M61

The registry was searched for 52/53-bit fields, medium-width Riesel primes, and
an M31 plus approximately 53-bit architecture before implementation.  None had
been tested.  This differs from the near-M61 radix-33 field: its purpose is not
an odd transform length, but to make the critical field's operands small enough
for three native 32-bit limb products while retaining the current 4M geometry.

A deterministic unsigned-64 prime search found

```text
q = 2^53 - (67*2^21 + 1)
  = 9,007,199,114,231,807.
```

Here `q+1` has 2-adicity 21, the quadratic field has the required order-`2^21`
norm-one subgroup, and the exact root-of-two exponent gate passes.  The
complement is

```text
67*2^21 + 1 = 2^27 + 2^22 + 2^21 + 1,
```

so it is unusually favorable for shift/add folding.  `M31*q` has almost exactly
84 bits, 0.193 bits more than the q24*M61 product that passed the empirical
p150 coefficient-range run.  It therefore has a credible, although narrow,
140--150M capacity path if performance succeeds.  It retains one 8-byte M31
quadratic value plus one 16-byte q53 value per pair, hence the production
48-MiB state size.

[`src/cuda/q53_m61_bench.cu`](src/cuda/q53_m61_bench.cu) implements two exact
q53 reducers and the M61 comparator over all 2,097,152 packed values:

- a native wide product followed by three folds at bit 53; and
- a radix-`2^27` Karatsuba product using three 32-bit products before the same
  folds.

Both q53 paths agree over every GPU output and with an independent host
`unsigned __int128` oracle.  Ptxas reports 32 registers for M61, 36 for the
wide q53 reducer, and 40 for the limb reducer, with no stack or spills.  Fresh
31-sample medians were:

| Dependent quadratic multiplies | M61 | q53 wide | wide/M61 | q53 limbs | limbs/M61 |
|---:|---:|---:|---:|---:|---:|
| 1 | 19.680 us | 18.016 us | 0.915x | 26.176 us | 1.330x |
| 2 | 28.096 us | 29.024 us | 1.033x | 39.200 us | 1.395x |
| 4 | 31.744 us | 49.312 us | 1.553x | 73.728 us | 2.323x |
| 8 | 52.480 us | 91.072 us | 1.735x | 134.464 us | 2.562x |
| 16 | 96.608 us | 162.848 us | 1.686x | 244.864 us | 2.535x |
| 32 | 179.360 us | 312.288 us | 1.741x | 470.016 us | 2.621x |

The one-round apparent win is a memory-bound consequence of identical 16-byte
input/output traffic.  Once values are reused as in an NTT tile, forming and
folding the 106-bit product dominates.  The low-Hamming complement still needs
multiple cross-word shifts/adds and repeated folds; explicit limb Karatsuba
adds even more recombination dependency than it removes multiply work.

Decision: **reject M31*q53 before roots, transform, or carry integration.**
The exact prime and range observation are reusable, but a modulus width just
below the FP64 mantissa is not itself a GPU arithmetic advantage.  Do not repeat
this with another near-`2^53` prime unless it has a fundamentally one-fold
reducer; changing only the complement cannot erase the measured 55--74%
arithmetic-dense deficit.

## Trivial Hensel-lift audit for the M31 generator

One possible way around the rejected serialized `M31^3` ring was checked
algebraically before another GPU implementation.  If the production Gaussian
generator `(7735,748621)` already had order `2^32` modulo `M31^3`, or if its
Teichmuller lift had zero/sparse higher p-adic digits, root multiplication could
have acted almost independently on three M31 planes.

Neither condition holds.  The original pair has order `2^32` only modulo M31;
its `2^32` power is nontrivial modulo both `M31^2` and `M31^3`.  Its unique
Teichmuller lift modulo `M31^3` has dense digits:

```text
real: (7735,   142746575, 1291122155)_M31
imag: (748621, 623754907, 1392908670)_M31.
```

Thus every general root product mixes the p-adic planes and retains the
dependency graph measured by the existing M31-cubed benchmark.  Decision:
**do not reopen M31-cubed as independent residue planes.**

## CUDA persisting-L2 trig-table gate

The registry contained L2-aware queue scheduling and rejected explicit L2
striping, but no CUDA access-policy-window measurement.  The wrapper contained
an unfinished, unused persisting-L2 helper, so this was completed as the opt-in
`L2_PERSIST=1` mode rather than inferred from cache-size specifications.

The mode queries the actual CUDA limits, reserves persisting L2, and marks the
largest real read-only allocation on every active stream.  On this run that is
the 6-MiB combined height/tail trig table; the device accepted a 6-MiB window
and 6-MiB set-aside with a 100% requested hit ratio.  Four sequential exact
p=136279841 runs all reproduced iteration 30,000 residue
`9139db3046e846d4`:

| Policy | 30k time |
|---|---:|
| ordinary control A | 202.0 us/iteration |
| persisting L2 A | 202.3 us/iteration |
| persisting L2 B | 202.9 us/iteration |
| ordinary control B | 202.6 us/iteration |

Control and trial means are both 202.3--202.6 us within noise, with no favorable
trial.  The table is already naturally cache-resident in the 128-MiB L2, and
the workload is integer/power limited.  Decision: **leave the mode opt-in as a
diagnostic, but reject persisting L2 as a speedup.**  Do not combine it with
graphs or attribute sub-microsecond variation to cache residency.

## Persistent radix-`2^31` M61 representation

Before repeating the prior INT32 M61 no-go, its exact scope was checked.  The
old benchmark retained canonical `u64` values, split them into limbs for every
product, and reassembled the result.  It did not test a memory-neutral transform
whose state and roots remain two `u32` digits throughout.  That distinction was
therefore eligible for one stronger gate.

[`src/cuda/m61_persistent_limb_bench.cu`](src/cuda/m61_persistent_limb_bench.cu)
stores each M61 scalar as canonical radix-`2^31` low/high digits.  Its quadratic
value remains 16 bytes, exactly matching production.  Addition, subtraction,
Karatsuba multiplication, the `2^62 == 2` fold, and every dependent product stay
in limbs; there is no per-product extraction or assembly.  All 2,097,152 GPU
outputs agree with the wide path, and 4,096 independently recomputed host
chains pass.  Ptxas reports 38 limb registers versus 32 wide, with no stack or
spills.

| Dependent quadratic multiplies | Wide M61 | Persistent limbs | Limb/wide |
|---:|---:|---:|---:|
| 1 | 17.696 us | 25.856 us | 1.461x |
| 2 | 19.072 us | 43.296 us | 2.270x |
| 4 | 27.968 us | 75.200 us | 2.689x |
| 8 | 48.736 us | 138.272 us | 2.837x |
| 16 | 91.168 us | 267.296 us | 2.932x |
| 32 | 177.216 us | 519.392 us | 2.931x |

Persistent storage removes the conversion concern but makes the result worse.
The incumbent can add/subtract and range-fold a scalar in a short 64-bit graph;
the limb path must propagate mixed-radix carries and normalize both digits after
every butterfly/product.  Three nominal 32-bit partial products do not offset
that work, and the RTX Blackwell compiler already lowers the wide M61 product
efficiently.

Decision: **reject persistent INT32 M61 before an NTT tile or integration.**
This closes the meaningful representation-level distinction left by the first
limb benchmark.  A future 32-bit proposal must avoid limb normalization across
several butterflies, not merely retain the same limbs in memory.

## Composite two-prime NTT ring audit

Before implementing another two-Riesel path, the registry was checked for the
stronger algebraic variant: choose compatible primes `q0` and `q1`, form
`Q=q0*q1`, and execute one NTT directly over `Z/QZ[i]`.  This is algebraically
valid by CRT and is different from the rejected packed-q kernel, which still
executes two independent modular transforms.  It would be useful only if `Q`
also admitted a reduction nearly as cheap as M61.

A search among compatible sub-`2^31` factors found no product close enough to a
power of two for a one-fold reducer.  Extending the factor search through
`2^40` found the best product near `2^61` as

```text
Q = 239075327 * 9644802047
  = 2305834203236794369
2^61 - Q = 8805976899583                 (44 bits)
```

The best sub-`2^31` products had roughly 50-bit complements near `2^61`, and
unsigned-32-bit pairs near `2^64` still had roughly 57-bit complements.  These
are not Mersenne-like reducers: they require either a generic 61-bit Montgomery
multiply or several wide folds.  The existing near-M61 benchmark already loses
by about 1.43--1.60x with a much more favorable 26-bit complement, so a 44-bit
complement cannot credibly replace production M61 arithmetic.

Decision: **reject the composite-ring implementation gate.**  Reopen it only if
a root-compatible product is found with a genuinely one-fold complement; CRT
isomorphism alone does not make the composite modulus cheap on this GPU.

## Alternate M61 pair-square identity

The experiment registry and ledger contained no prior run of
`ENABLE_BETTER_ONEPAIRSQ`.  The dormant path in
[`src/cl/tailsquare.cl`](src/cl/tailsquare.cl) replaces the independent
`a^2`, `b^2*t^2`, and `2ab` products by

```text
2ab = (a+b)^2 - a^2 - b^2,
```

saving one wide M61 multiply.  It was exposed as an opt-in `-use` key and first
validated through iteration 30,000; every trial reproduced residue
`9139db3046e846d4`.

Final `sm_120` SASS explains why the lower multiply count does not win.  At the
production 96-register cap, the control uses 614 `IMAD` instructions with no
spill.  The identity uses 580 `IMAD` instructions, but adds a single 8-byte
spill load/store pair and slightly more shift, logic, and predicate work.  More
importantly, computing `(a+b)^2-a^2-b^2` puts the `2ab` result after the two
input squares on the dependency graph; the production formula exposes its
three products independently.

Nsight Systems measured 2,000 M61 tail calls for each 1,000-iteration trace:

| M61 tail form | Registers/local | Average call | Median call |
|---|---:|---:|---:|
| Production identity, 96-register cap | 96 / 0 B | 52.664 us | 51.520 us |
| Alternate identity, 96-register cap | 96 / 8 B | 52.986 us | 51.840 us |
| Alternate identity, compiler-selected | 80 / 0 B | 52.745 us | 51.616 us |
| Production identity, compiler-selected | 80 / 0 B | 53.168 us | 52.128 us |

Thus removing the spill recovers most of the loss but still does not beat the
production identity.  Unprofiled exact 30k runs were `201.9` and `202.8`
us/iteration for two production controls, `202.7` for the capped alternate,
and `202.8` for the spill-free alternate.  The equal final pair is noise, not a
speedup; the isolated kernel medians consistently favor production.

Decision: **leave the alternate identity disabled.**  A lower-live-range source
rewrite is not justified because the compiler-selected spill-free form was
already measured and remained slower.  This is a dependency/ILP failure, not
merely the historical compiler-spill failure noted in the source comment.

## Whole-backend FP64 and mixed-radix-15 screens

The experiment registry contained FP32/FP64 hybrids and local FP64 arithmetic
tests, but no complete run of the existing pure-FP64 backend at the current
exponent.  An exact 8M `512:16:512:202` run was therefore made before considering
another FP64 architecture.  It reproduced the trusted residues through
iteration 10,000, but measured **937.5 us/iteration**.  This is 4.6x the
pre-existing M31/M61 result and more than five times the 180-us gate.  The pure
FP64 backend is a correctness reference, not a performance lead on this GPU.

A 3.75M (`15*2^18`) M61 plus q31 design was also screened algebraically before
implementation.  The initially attractive existing field
`q=2141192191=1021*2^21-1` has the radix factors, but the exact
Crandall--Fagin weight equation for that transform length is insoluble.  Other
Riesel primes do pass it; for example
`q=2116812799=8075*2^18-1`.  However, 3.75M removes only 6.25% of state and
first-order transform work.  Even a zero-overhead radix-15 and zero-cost generic
q field could not provide the required 12.4% complete-iteration improvement,
while all measured generic q31 fields are materially slower than M61.  This
route was rejected at the capacity/budget gate rather than repeating the
already measured mixed-radix and generic-field integrations.

## One-warp M61 512-point tail

The registry was checked before implementation.  The earlier warp-shuffle test
was a generic three-prime radix-256 transform; there was no test retaining the
production `GF((2^61-1)^2)` field and assigning one complete 512-point critical
tail line to one warp.  The new exact gate is
[`src/cuda/m61_warp_tail_bench.cu`](src/cuda/m61_warp_tail_bench.cu).

One physical lane owns 16 quadratic values, so a warp owns all 512 values.  The
first implementation uses nine radix-2 stages; the second reproduces the
production radix-8 factorization and replaces its two shared-memory transposes
in each direction with exact register/warp exchanges.  Both execute forward
NTT, pointwise square, and unnormalised inverse NTT over all 2,097,152 quadratic
values.  Random dense round trips and a sparse cyclic-square oracle pass.

Five independent 31-sample runs gave run-median medians:

| M61 512-point architecture | Median |
|---|---:|
| Register/warp radix 2 | 193.824 us |
| Register/warp radix 8, exact exchange | 110.528 us |
| Radix-8 arithmetic with all exchange deleted | 96.480 us |
| Production M61 tail (Nsight reference) | **56.128 us** |

The exact radix-8 kernel uses 148 registers, no stack, no spills, and no
barriers.  Its SASS has 1,422 `IMAD` and 1,536 `SHFL` instructions.  The
radix-2 form uses 127 registers and is spill-free, but has 3,655 `IMAD` and 960
`SHFL`.  Radix 8 therefore confirms why the incumbent's specialized transform
matters, but the one-warp ownership still loses badly.

The zero-exchange build is an intentionally impossible lower bound: it retains
the same radix-8 butterflies, twiddle products, loads, square, inverse, and
stores while making both required transposes free.  At 96.480 us it is already
72% slower than the complete production tail.  The gap reflects both the
standalone gate's canonical arithmetic versus production's carefully delayed
M61 ranges and the loss of thread-level parallelism from assigning twice as
many values to each lane.  Porting all production range-specialized arithmetic
cannot turn the ownership change into a 25.5-us whole-iteration saving: the
exact warp exchanges themselves add about 14 us over the impossible lower
bound, whereas the incumbent shared exchange is already included in 56.128 us.

Decision: **reject the one-warp 512-point tail.**  Do not implement a larger
bit-transpose network or integrate it into PRPLL.  Shared memory is not an
accidental bottleneck here; it efficiently supplies the cross-register
all-to-all exchange while the 64-thread-per-line layout preserves more
parallelism and the mature range-aware arithmetic.

## Folded Tensor sidecar concurrency audit

The exact folded-M31 correction design still had one apparently distinct
hardware option: tolerate a slower Tensor transform if it uses dedicated MMA
units and therefore avoids the 14--19-us integer-pipeline contention measured
for the SIMT side stream.  The registry contained exact Tensor radix tiles but
no dependency-closed concurrent sidecar.  Before implementing a full transform,
the existing exact quadratic Tensor kernel was rerun at the real folded-plane
population of 262,144 quadratic values.

The three-byte q24 version takes **9.408 us for one radix-16 group**, versus
5.504 us for sparse SIMT.  It uses 27 INT8 MMAs: three quadratic Karatsuba
products times nine byte-limb products.  M31 requires four bytes, or 48 MMAs
for the corresponding group.  A 262,144-point quadratic transform has eighteen
binary stages; forward plus inverse therefore needs four radix-16 groups and
one radix-4-equivalent group in each direction—ten groups before layout,
pointwise work, decoding, or synchronization.  Even scaling only the measured
MMA work places M31 near or beyond the complete approximately 102-us M61
overlap window.  Global exchanges and Tensor/M61 board-power contention can
only increase it.

Decision: **reject a full folded-M31 Tensor sidecar before integration.**  This
is not inferred only from the old 42% full-population tile deficit: the exact
side population and the number of mandatory quadratic limb MMAs were measured
and counted.  It cannot turn the p136-only folded syndrome into a general
140--150M solution either; the observed quotient alphabet grows beyond the
single-folded-residue information bound there.

## Single composite `M31*M61` transform ring

The registry contained a composite-q-pair search but no direct transform over
the product of the two production fields themselves.  This is algebraically
valid and retains the full general 4M capacity.  Let `B=2^31` and
`Q=(2^31-1)(2^61-1)`.  Then

```text
2Q = B^3 - B^2 - 2B + 2
B^3 = B^2 + 2B - 2  (mod Q)
B^4 = 3B^2 - 2       (mod Q)
```

A Q scalar has three radix-B digits and a quadratic value is 24 bytes, exactly
matching the separate M31 plus M61 state.  A six-product pairwise-Karatsuba
multiplier and the sparse reduction were first checked over 100,000 independent
host products.  The architecture-population CUDA gate is
[`src/cuda/m31_m61_composite_bench.cu`](src/cuda/m31_m61_composite_bench.cu).
Every one of 2,097,152 quadratic outputs is independently reduced modulo M31
and M61 and agrees with the existing separate-field multiplication path.

The first canonical implementation was 12--13x slower because it repeated a
three-limb normalization eight times per complex multiply.  It was not used as
the decision result.  The strengthened form instead computes four raw
three-limb products, combines their polynomial coefficients, and performs only
two final sparse reductions.  This is the direct-ring analogue of production's
delayed range handling.  Final 21-sample medians were:

| Dependent quadratic products | Separate M31+M61 | Composite Q | Q/separate |
|---:|---:|---:|---:|
| 1 | 28.928 us | 77.024 us | 2.663x |
| 2 | 29.280 us | 88.672 us | 3.028x |
| 4 | 44.544 us | 163.008 us | 3.659x |
| 8 | 81.312 us | 301.888 us | 3.713x |
| 16 | 147.072 us | 574.560 us | 3.907x |
| 32 | 281.856 us | 1,107.104 us | 3.928x |

Both kernels are spill-free.  The composite path uses 56 registers and 1,104
static SASS instructions versus 38 registers and 336 instructions for the
separate comparator.  Its complex product needs four six-partial-product
three-limb products, followed by signed carry/fold work.  In contrast, CRT form
lets every root product use the native special M31 and M61 reducers.  Combining
the fields deletes CRT only at one edge but makes every transform twiddle several
times slower.

Decision: **reject the single composite-ring NTT before roots or integration.**
The direct product ring is mathematically neat but RNS is already the efficient
factorization of its arithmetic on this GPU.  Do not confuse this with the
earlier search for a different q0*q1 product near `2^61`; this test uses exactly
the production modulus product and identical persistent bytes.

## 2025--2026 architecture literature refresh

A focused refresh was made after the composite-ring gate so that newer results
would not be mistaken for untried local ideas:

- [Zhang and Franchetti's 2025 MoMA/NTTX work](https://users.ece.cmu.edu/~franzf/papers/CGO_2025_Naifeng.pdf) generates strong multiword GPU
  NTTs, but its headline regime is a roughly 1,024-point transform that fits
  wholly in shared memory.  The local M89, M127, M31-cubed, M31-prime-power,
  persistent-limb, and composite-Q gates measure the corresponding multiword
  arithmetic at the actual 2--4M population; each loses before a large NTT.
- [Ozcan, Javeed, and Savas's 2025 four-step GPU NTT](https://research.sabanciuniv.edu/51889/1/High-Performance.pdf) improves large global-pass
  transforms through factorization, block shape, and memory locality.  PRPLL's
  production `512x8x512` schedule already embodies this class, and the exact
  alternate-factorization, resident-tile, direct-RNS cache-fusion, and
  Goldilocks three-pass tests cover its relevant choices locally.
- [New 2026 Ozaki FP4 work](https://arxiv.org/abs/2608.06812) can emulate FP64 efficiently on Tensor Cores for
  very large dense GEMMs.  Its own FFT analysis identifies reconstruction as
  the limiter for small inner factors.  PRPLL's radix factors are only 8--16;
  the exact local INT8, q24, compensated-TF32, Toeplitz-M61, and side-population
  Tensor gates include that packing/reconstruction cost and all lose to sparse
  SIMT.  Even an ideal 2x FP4 improvement leaves the measured 4.5x Toeplitz-M61
  deficit above the scalar path.
- Multi-GPU NTT work attacks communication across several devices.  This host
  exposes one RTX PRO 6000, and splitting one sequential exponent over PCIe
  would add communication to every recurrence.  It does not alter the accepted
  strategy of assigning unrelated exponents to separate GPUs when available.

The refresh supplies no architecture whose lower bound removes the required
roughly 22--26 us.  In particular, it does not reopen dense Tensor NTTs,
multiword composite fields, another four-step layout, or launch-only fusion.

## Production regression checkpoint after the architectural gates

After adding the benchmark-only warp-tail and composite-ring gates, the CUDA
PRPLL binary was rebuilt and the unchanged production configuration was run
for 10,000 exact iterations at exponent 136279841:

```text
FFT: 4M 1:512:8:512:202
iteration 2,000:  05d6515c416b83e2
iteration 10,000: 52316d51aa52e6b7
final reported interval: 197.4 us/iteration
```

Both residues match the trusted reference.  The 197.4-us final interval is
faster than the historical 205.5-us campaign baseline, but it comes from the
same production algorithm and configuration; none of the rejected experimental
kernels is linked into this path.  It must therefore be classified as
run-to-run/compiler-clock/environment variation, not as a speedup achieved by
an architectural change.  The accepted best remains the pre-existing
M31/M61 implementation, and the 180-us gate remains unmet.

## Reopened 3.5M M31+M61+M19 radix-7 gate

Before implementation, the registry was checked for `3.5M`, radix 7, and the
three-Mersenne combination.  The only prior entry was a rejection inferred
from forum radix-7 measurements; no local capacity proof, exact radix-7 edge,
or production-population scheduling measurement exists.  This candidate is
also distinct from the locally rejected 3M M31+M61+M19 engine: it uses seven
independent `2^19`-word Good--Thomas channels, rather than three `2^20`-word
channels, specifically to restore range through the requested p150 endpoint.

For `N=7*2^19=3,670,016` real words, write `p=b*N+R`.  The same exact
mechanical-word/Beatty overlap argument already validated for the 3M engine
bounds a weighted cyclic square, including the conservative factor six, by

```text
2^(2b) * (N + 3R)/4 * 6.
```

The balanced M31*M61*M19 range has 109.999997 bits.  For relevant exponents:

| Exponent | `b` | `R` | Bound bits | Balanced headroom |
|---:|---:|---:|---:|---:|
| 136,279,841 | 37 | 489,249 | 96.877671 | 13.122326 |
| 139,999,991 | 38 | 539,383 | 98.919298 | 11.080699 |
| 145,000,003 | 39 | 1,869,379 | 101.730369 | 8.269629 |
| 150,000,001 | 40 | 3,199,361 | 104.246421 | **5.753576** |

Each listed remainder is coprime to N, so the exact maximum-overlap argument
applies.  Capacity therefore passes generally through p150; unlike the 3M M19
route, this is not a p136-only range result.  All three base fields also contain
scalar seventh roots: seven-way Good--Thomas coupling preserves the packed
power-of-two DGT inside each channel.

The performance premise is narrower.  Counting a 64-bit field as two 32-bit
limbs, the binary-stage work is exactly tied at first order:

```text
production: (M31 + two-limb M61) * 2^21 packed values * 21 stages
candidate:  (M31 + two-limb M61 + M19) * (7*2^18) values * 18 stages
ratio:       1.000
```

The only plausible gain is scheduling: the M61 stream becomes about 25%
shorter, while M31 and cheap M19 can occupy another stream.  The extra state is
56 MiB versus 48 MiB, and radix 7 can erase that scheduling benefit.  The first
gate is therefore an exact architecture-population forward-radix-7, pointwise
square, inverse-radix-7 kernel, compared with an equal-payload square-only
control and timed both fused and as concurrent `M61 | M31+M19` streams.  Use a
Rader/Good--Thomas six-point convolution (ten scalar constant products per
component), not a deliberately slow 49-product direct DFT.  Do not build a
complete 3.5M engine unless this edge leaves a credible 18--22-us whole-cycle
opportunity.

That exact gate is now implemented in
[`src/cuda/m31_m61_m19_radix7_bench.cu`](src/cuda/m31_m61_m19_radix7_bench.cu).
It uses generator 3 for Rader's permutation.  Each length-six convolution is
diagonalized as coprime radix 2 by radix 3, so a scalar component needs two
radix-3 constant products in the forward transform, six frequency products,
and two radix-3 products in the inverse.  Both seven-point transforms are
left uniformly scaled; an independent direct cyclic-convolution oracle checks
the resulting factor `6^3*7=1512`.  Sixty-four random groups in all three
fields pass.  The M61 and M31+M19 kernels use 48 and 40 registers respectively,
with no stack, spills, shared memory, or barriers.

Five fresh 31-sample processes gave run-median medians:

| 3.5M edge gate | Median |
|---|---:|
| All fields, radix-7/square/inverse | 87.712 us |
| All fields, equal-payload square only | 72.896 us |
| Concurrent M61 and M31+M19 radix-7 paths | **86.592 us** |
| Concurrent equal-payload square-only paths | **76.352 us** |

Thus the intended split schedule pays about **10.24 us** for both radix-7
directions.  This is much less than a direct DFT forecast and passes the first
gate.  The forum radix-7 loss was directionally useful, but was not sufficient
evidence to skip this locally distinct three-Mersenne schedule.

The next measurement first raised the optional multi-worker policy limit from
four to eight and ran seven exact 512K M31/M61 PRP workers.  Their final 10K
intervals were 171.1--187.1 us/worker; the useful batch time is one such worker
latency, not seven times it.  However, every worker redundantly performs its own
CRT/carry, and an Nsight trace shows extensive overlap.  This establishes that
seven channel populations fit and execute concurrently, but it is not a native
transform-core time.

`GOOD_THOMAS7=1,M19_FIELD=1` therefore adds the actual native core gate at
`54:512:7:512`.  It runs seven independent 512K packed power-of-two channels
through the production M31, M61, and M19 middle/tail kernels on three queues.
The exact radix-7 coupling remains omitted intentionally and is charged from
the independently validated gate above.  The legacy boundary is deliberately
not a correctness implementation; as with the earlier 3M core gates, only the
representative nonzero transform timeline is admissible evidence.

The decisive timeline result over 5,900 steady cycles is:

```text
3.5M GT7 binary-channel core median: 152.368 us
p10 / p90:                           136.735 / 155.391 us
exact 4M production core, same span: 129.631 us
plus measured radix-7 increment:      10.240 us
dependency-complete candidate core:  162.608 us
```

The span is measured from the first of the three field `fftMiddleIn` launches
to the last `fftMiddleOut` completion.  Merely summing per-kernel medians would
incorrectly predict about 101 us: queue gaps, inter-field contention, cache
phasing, and the 300-W power limit add roughly 51 us to the actual dependency
span.  The production trace was recomputed with the identical interval method
and reproduces the established 129.6-us core, validating the analysis.

Even reusing the smaller 3M engine's **38.4-us** exact shared-layout
M31*M61*M19 CRT/carry as an impossibly favorable cost gives
`162.608+38.400=201.008 us` before either required width boundary.  The real
3.5M carry has 16.7% more coefficients.  It therefore cannot meet 180 us, and
an exact width/carry integration would only quantify an already decisive loss.

Decision: **reject the 3.5M radix-7 engine at the native transform-core gate.**
Keep the exact radix-7 benchmark, capacity proof, optional seven-worker
measurement, and `GOOD_THOMAS7` timing scaffold as evidence.  Do not implement
the end-to-end carry or mistake equal first-order limb-stage work for equal GPU
time: the extra field and three-way scheduling make the measured core 25.4%
slower than production after charging radix 7.

After these opt-in changes, a clean default 4M M31/M61 regression reproduced
`05d6515c416b83e2` at iteration 2,000 and `52316d51aa52e6b7` at iteration
10,000.  Its final interval was 197.5 us.  The normal production path is
therefore unchanged in both exact output and current-session performance.

## Reopened 4.5M FP32+M61 Good--Thomas radix-9 gate

The registry was checked for `4.5M`, radix 9, and FP32+M61 before changing
source.  There is no local implementation or timing of this combination.  The
earlier 33/32-length FP32+M61 entry rejected an odd packed-M61 root because the
norm-one subgroup has power-of-two order.  That reasoning predates the later
validated Good--Thomas construction: keep independent power-of-two Hermitian
channels and apply the odd base-field DFT across channel spectra.  The forum's
corrected radix-9 result is negative background evidence, but it does not time
this low-precision two-field schedule on the present GPU.

Choose

```text
N = 9*2^19 = 4,718,592 real words (4.5 Mi)
packed M61 population = 9*2^18 = 2,359,296
```

Nine divides `M61-1`, so the scalar radix-9 edge exists.  The exact weighted
convolution also passes:

```text
gcd(N, M61-1) = 18
2^((M61-1)/18) mod M61 = 1
N^-1 mod 61 = 30
theta = 2^30, theta^N mod M61 = 2
```

The last identity is particularly favorable: the M61 Crandall--Fagin weight is
still a cheap Mersenne rotation rather than a generic product.  A channel has
`2^18` packed values and eighteen binary stages.  Relative to production, the
M61 binary-stage count is

```text
(9*2^18*18) / (2^21*21) = 0.9642857.
```

FP32 replaces the hidden M31 plane and has previously overlapped M61 more
effectively, while radix 9 factors into two radix-3 transforms.  Persistent
state grows 12.5%, from 48 to 54 MiB.

Precision, not modular range, is the risky correctness point.  Scaling the
existing 4M FP32+M61 configured maximum `132,791,664` by 9/8 predicts only
`149,390,622`; exponent 150,000,001 is **0.408%** above it.  This is close
enough for an early transform-core gate, but no result is generally exact until
an exact M31 oracle measures quotient errors at p140, p145, and p150 and a
bounded correction scheme is proved.  The first implementation is therefore a
native binary-channel core timing scaffold.  It must beat the production
129.6-us core with enough room for its radix-9 increment and a complete exact
edge before any correction integration.

### Native radix-9 timing result: reject

The transform-core gate is complete.  `src/cl/fft9.cl` now has the missing
FP32 Nussbaumer nine-point DFT, copied operation-for-operation from the existing
FP64 routine with binary32 constants.  The modular channel edge in
`src/cl/fft-middle.cl` uses a `9=3*3` Cooley--Tukey factorization: six
radix-three butterflies and four scalar twiddles.  Its forward and inverse
primitive ninth roots were checked against direct modular DFTs on 2,000 random
vectors.  This is already the Winograd multiplicative count of ten scalar
products per coordinate for a length-nine DFT
([operation table](https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/iet-com.2017.0837));
the two coordinates of one `GF(M61^2)` value therefore execute twenty M61
products per direction.  A different radix-nine spelling cannot supply the
missing factor-of-two improvement.

The first trace deliberately omitted the M61 odd edge while retaining nine
independent power-of-two channels.  This is an optimistic lower bound, not a
valid recurrence:

```text
                                         4M FP32+M61 control   4.5M GT9 lower bound
transform-core median                         116.704 us             115.424 us
carryFused median                              59.488 us              65.056 us
carry-start period median                     178.176 us             182.272 us
```

Thus the shorter binary transforms save only 1.28 us on the overlapped core.
The 12.5% larger coefficient population adds 5.57 us to the already fused
carry and puts even the impossible zero-cost-radix lower bound above 180 us.

Adding the actual register-resident M61 radix-nine edge raises
`fftMiddleInGF61` from 26.624 to 30.464 us and `fftMiddleOutGF61` from 21.248
to 24.864 us.  Both kernels use 64 registers with zero local memory or spills.
The dependency-complete measurements become:

```text
transform-core median                 124.479 us
carryFused median                      66.528 us
carry-start period median             192.543 us
plain second 2k interval          192.9--194.2 us
```

The inverse leaves the factor nine in the modular residue for this gate; an
exact carry could absorb its inverse into a CRT constant.  The scaffold still
uses the legacy flat boundary and consequently produces a deterministic wrong
checkpoint (`78fb759db933e74d` at p150 after 2,000 iterations).  That expected
failure does not make the timing optimistic enough: the measured arithmetic
edge alone misses the final gate by 12.5 us before any FP32 precision repair,
Good--Thomas boundary permutation, or exact correction is charged.

A bounded tuning sweep could not recover the gap:

| Variant | second 2k interval |
|---|---:|
| 512x9x512, current schedule | **192.9--194.2 us** |
| CUDA graph | 196.6 us |
| tail organizations 0/1/2/3 | 204.5 / 208.5 / 196.5 / 197.7 us |
| L1 carve-outs 0/1/2/3 | 195.4 / 194.6 / 250.1 / **193.4 us** |
| compiler-default registers (`NOREG=1`) | 192.9 us |
| alternate carry-shuttle or FFT cache policies | 195.2--196.3 us |
| 256x9x1024 | 223.9 us |
| 1024x9x256 | 230.1 us |

Decision: **reject the 4.5M FP32+M61 Good--Thomas radix-nine architecture for
the 180-us gate.**  Do not repeat it as a standalone scheduler, a different
radix-nine flow graph, or another rectangular factorization.  The opt-in
`GOOD_THOMAS9` code is a timing scaffold and is not an exact PRPLL mode.  Its
useful lesson is that reducing binary-stage count is insufficient when a
larger fused carry and a multiplicatively minimal odd edge consume all of the
gain.

After the opt-in changes, a clean default 4M M31/M61 p136 regression reproduced
`05d6515c416b83e2` at iteration 2,000 and `52316d51aa52e6b7` at iteration
10,000, with a final 197.7-us interval.  `git diff --check` passes, and the
unbuffered `setvbuf(stdout, nullptr, _IONBF, 0)` call remains in `src/main.cpp`.

### Exact end-to-end radix-nine prototype: correct, but 296.2 us

This architecture was reopened at the user's request despite the preceding
timing-gate rejection.  Before changing it, the earlier section above was
checked: it had implemented and timed the odd edge, but had *not* implemented
the Good--Thomas boundary permutation, inverse normalization, correct
initialization/checkpoint loading, or the radix-nine `tailMulGF61` path.  The
new work completes those missing pieces rather than repeating the old timing
scaffold.  It lives on branch `prototype/good-thomas9-exact`.

For `N = 9 * 2^19` scalar words, the modular transform is stored as nine
independent `2^19` planes.  The exact natural-to-plane map follows the CRT
isomorphism `Z/N -> Z/9 x Z/2^19`, rather than treating the nine planes as nine
contiguous natural chunks.  At the `512 x 9 x 512` boundary, for storage
channel `c`, row `y`, and natural plane number `m`, the implemented forward
gather is

```text
m_even = (y + 4*c) mod 9
m_odd  = (m_even + 5) mod 9
x_even = (p - m_even) * 57 mod 512
x_odd  = (p - m_odd ) * 57 mod 512
```

where `57 = 9^-1 mod 512`.  The inverse carry-side gather uses

```text
x_storage  = 9*x + m mod 512
c_even     = 2*(y - m) mod 9
c_odd      = c_even + 1 mod 9.
```

The M61 Crandall--Fagin weight uses root exponent 30 because
`30*N = 1 (mod 61)`, while the unnormalized inverse radix-nine result is
multiplied by `9^-1 mod M61 = 2049638230412172401` at the exact reconstruction
boundary.  `tailMulGF61` now forms Hermitian pairs independently within each
of the nine channels, so Gerbicz multiplication and checks use the same exact
layout as the main recurrence.  Normal modular initialization and persisted
state validation are enabled.  The `512x9x512` BPW limit was also admitted
through 150M after the exact checkpoint tests below.

The current short fused carry cannot be made correct merely by changing an
index: one thread owns one M61 width line, while reconstructing a natural
complex value now gathers its two components from different Good--Thomas
channels.  The exact prototype therefore deliberately selects the existing
long/split carry.  A future fast implementation would require a genuinely new
multi-channel fused boundary (or an equivalent frequency shear in the middle
kernels), not reuse of the incorrect flat fused boundary.

Correctness was compared against the production M31/M61 implementation:

| Requested exponent | PRPLL exponent | iteration 2,000 | iteration 10,000 | result |
|---:|---:|---:|---:|---|
| 136279841 | 136279841 | `05d6515c416b83e2` | `52316d51aa52e6b7` | exact; Gerbicz OK |
| 140000011 | 139999991 | `0e6dbe7b2f54cbb4` | `bd59105fe62925d2` | exact; Gerbicz OK |
| 145000003 | 144999991 | `51c2eb2d052424bc` | `b24d260a279f9575` | exact; Gerbicz OK |
| 150000007 | 150000007 | `ab2b98556e771a7e` | `fdc46a772cb01070` | exact; Gerbicz OK |

The p136 checkpoint was also reloaded with ordinary on-load validation and
continued exactly through iterations 12,000 and 22,000.  Finally, a fresh
one-million-iteration p136 run produced the known exact residue
`52b03a7cc55e677d` and reported **296.2 us/iteration**.  Therefore the exact
implementation is 90.7 us (44.1%) slower than the historical 205.5-us
production result and misses the 180-us gate by 116.2 us (64.6%).  The earlier
192.9--194.2-us figure must not be quoted as end-to-end performance: it measured
an incorrect flat-boundary fused scaffold.

After committing, a detached clean checkout of the branch was built from
scratch with `make CUDA=1 -j4` and rerun at p136.  That committed executable
again produced `05d6515c416b83e2` at iteration 2,000, confirming that the
prototype does not depend on uncommitted experiment files or stale build
artifacts.

A four-point exact p136 tuning sweep, all matching the 10,000-iteration
checkpoint, gave 283.4, 288.1, 289.7, and 291.8 us.  Their arithmetic mean is
288.25 us and median is 288.9 us.  The best short-run setting was
`MULTI_Q=1,L1CUDA=3` at 283.4 us, but the much longer run above is the
authoritative sustained measurement.

Nsight profiling of that best exact organization attributes the main kernel
medians as follows:

| kernel | median |
|---|---:|
| split `carry` | 61.216 us |
| `tailSquareGF61` | 61.760 us |
| `fftP` | 48.896 us |
| `fftMiddleInGF61` | 34.912 us |
| `fftMiddleOutGF61` | 32.576 us |
| `fftWGF61` | 26.320 us |
| FP32 tail / width / middle-in / middle-out | 24.480 / 19.776 / 16.544 / 15.904 us |
| `carryB` | 3.904 us |

Decision: the completed prototype is a correctness reference and a useful
starting point for a redesigned fused Good--Thomas boundary, but it is **not a
speedup**.  Do not repeat the flat short-carry experiment, omit the cross-plane
permutation, or compare the 192.9-us scaffold against production.  Any follow-up
must first recover roughly 91 us by removing the split width/`fftP` boundary;
small kernel tuning cannot close that structural gap.

## Post-radix-nine architecture audit

The exact radix-nine result was used as a new lower-bound gate before attempting
another fused boundary.  Its earlier zero-cost-odd-edge timing was already
182.272 us, and the real radix-nine edge raised that to 192.543 us while still
using the wrong flat boundary.  Thus even a free Good--Thomas permutation would
not meet 180 us.  A multi-channel fused implementation could make the 296.2-us
reference much faster, but it cannot turn this transform into the requested
speedup.  No second fused radix-nine kernel was implemented.

### Dynamic-register warp specialization: stronger hardware no-go

The prior field-specialized edge used ordinary divergent field ownership, so
CUDA allocated every thread for the largest path.  Blackwell's architecture-
specific `setmaxnreg` instruction appeared to offer a materially different
test: decrease the M31 warpgroup's allocation and transfer registers to the
M61 warpgroup at a synchronized handoff.  The existing production-population
[`src/cuda/warp_specialized_edge_bench.cu`](src/cuda/warp_specialized_edge_bench.cu)
was extended rather than creating another proxy.

The dynamic kernel owns two tiles per 256-thread CTA and uses
`setmaxnreg.dec.sync.aligned` / `setmaxnreg.inc.sync.aligned` around the
field-specialized phases.  Its output matches the 64-thread per-thread-field
reference exactly.  A fresh 2,097,152-value run measured:

```text
serial per-thread fields:       0.079 ms
ordinary field-specialized:     0.090 ms  (1.146x serial)
dynamic-register specialized:   0.137 ms  (1.734x serial)
```

The final cubin reports 128 registers, a 176-byte stack, and 25,600 bytes of
shared memory for the dynamic kernel, versus 90 registers, no stack, and 9,216
bytes shared for the serial reference.  Dynamic redistribution does not shorten
the source-level live ranges or prevent ptxas from materializing the union of
the phases; it adds synchronization and crosses the spill boundary.  This
strengthens the previous warp-specialization rejection.  Do not repeat it with
different `setmaxnreg` limits unless the algorithm first removes live state.

### Exact-FP32 small-field 3M architecture: algebraic no-go

A new resource-complementary proposal was screened before GPU implementation:
shorten the critical M61 transform to `3*2^20` real words and replace the added
integer q field with enough sub-25-bit fields to execute their modular products
exactly through FP32/FMA.  Centered residues for `q < 2^25` fit exactly in one
FP32 value, and an error-free product can recover the rounded product residual;
the intent was to move the extra range work off the saturated integer pipelines.
This differs from the rejected 4M q24 replacement and from approximate FP32
coefficient repair.

The required field set does not exist.  An exhaustive deterministic search over
all primes `q = k*2^s +/- 1`, `19 <= s <= 24`, `q < 2^25`, applied all of:

1. a radix-three base-field edge (`3 | q-1`);
2. the power-of-two order for each packed 1M-word Good--Thomas channel; and
3. the exact Crandall--Fagin weight condition
   `2^((q^2-1)/gcd(3*2^20,q^2-1)) = 1 (mod q)`.

Only `q=2^19-1=524287` passes.  It supplies 19 bits, far short of the several
distinct fields needed beside M61 for the conservative p150 range.  In
particular, the previously useful q24 field `14680063=7*2^21-1` returns
`9024597`, not one, in the 3M weight test.

The scalar-NTT alternative was also exhausted over the same exact-FP32 range.
The only prime below `2^25` with `3*2^20 | q-1` is
`q=28311553=27*2^20+1`; its scalar weight test returns `2^9=512`, not one.
Consequently neither quadratic packing nor a full scalar small-field NTT can
supply multiple FP32-exact 3M residue planes.  Decision: reject this
representation algebraically; do not build a multi-field transform from q24
primes that individually fail the target weight equation.

### Independent native-PFA cross-check: also slower

After completing and measuring the local exact prototype, the current PrMers
tree was used as an independent architectural cross-check.  This was not used
as a substitute for the local correctness tests above: it is a separately
developed implementation with an Aevum native Good--Thomas path.  The tested
tree was `cherubrock-seb/PrMers` commit `d1c2e07`, and the forced plan was
`pfa9:1:512:9:512:202` at exponent 136279841 on the same RTX PRO 6000 Blackwell
Max-Q GPU.

After its OpenCL cache and checkpoint were warm, the representative interval
reported 2,883.91 iterations/s, or **346.75 us/iteration**.  This is 141.25 us
slower than the historical 205.5-us production path and misses the 180-us gate
by 166.75 us.  Its launch trace also exposes separate GF31/GF61 width,
carry/carryB, and forward-transform kernels.  In other words, this independent
native-PFA implementation reaches the same architectural boundary as the local
exact prototype: the odd-radix transform is feasible, but its cross-plane digit
map prevents reuse of the exceptionally cheap fused M31/M61 carry boundary.

For additional context, the same PrMers tree's independent Marin integer-IBDWT
backend selected an 8,388,608-word transform and reported 1,999.08
iterations/s, approximately **500.23 us/iteration**.  Neither external backend
is a speedup on this machine.  The local 296.2-us exact prototype remains the
faster of the completed radix-nine implementations, but all measured exact
paths are decisively behind production.

## Complex-point Toom Tensor correction tile

The correction-only q24 field was reconsidered only through an untried Tensor
mapping, after checking the earlier `q24_tensor_bench` result.  That benchmark
used schoolbook three-byte multiplication: nine unsigned-INT8 MMAs per dense
base-field matrix product and three Karatsuba products per quadratic value.
The new path instead treats each 24-bit operand as
`d0+d1*x+d2*x^2` and evaluates it at `0`, infinity, `1`, `-1`, and
`i`.  The product at `i` takes three real matrix products, for seven FP16 MMAs
per base-field product rather than nine INT8 MMAs.

This is exact floating-point integer arithmetic, not an approximate transform.
All three input digits are bytes; the largest real evaluation is 765, which is
exact in FP16.  Every 16-term matrix dot product is below `2^24`, so each exact
FP16 product and its FP32 accumulation is integral without rounding.  The
interpolation is particularly cheap:

```text
c2 = (C(1)+C(-1))/2 - c0 - c4
c1+c3 = (C(1)-C(-1))/2
c1-c3 = Im(C(i)).
```

[`src/cuda/q24_tensor_bench.cu`](src/cuda/q24_tensor_bench.cu) now validates
both Tensor forms against the sparse Montgomery radix-16 transform over all
2,097,152 `GF(qC^2)` values.  Every output agrees.  The Toom kernel uses 40
registers with no stack or spill, versus 56 for the prior INT8 kernel and 16
for sparse SIMT.  Final SASS contains native `HMMA.16816.F32` instructions, so
the result is not a compiler fallback.

```text
schoolbook INT8 Tensor, 27 logical MMAs:   36.000 us
complex-Toom FP16 Tensor, 21 logical MMAs: 87.776 us
sparse SIMT radix-16:                      19.872 us
Toom / INT8:                                2.438x
Toom / SIMT:                                4.417x
```

FP16 Tensor throughput and the evaluation/interpolation data path overwhelm
the nominal 22% reduction in dense matrix products.  This is a stronger failure
than the original q24 Tensor tile, not a surviving component lead.  Decision:
**reject complex-point Toom and do not build a full q24 correction sidecar from
it.**  A future Tensor proposal must avoid dense radix matrices or encode more
than one independently recoverable modular product per accumulator; merely
changing limb multiplication from schoolbook to Toom is closed.

## Exact FP64-FMA offload inside M61

The pure-FP64 backend no-go above does not answer a narrower hardware question:
whether one of the three independent scalar products in a quadratic M61
Karatsuba multiply can move to the otherwise separate FP64 pipeline while the
other two remain on integer units.  The experiment registry had no exact FP64
arithmetic inside the production M61 representation, so this was eligible for
an early primitive gate.

[`src/cuda/m61_limb_bench.cu`](src/cuda/m61_limb_bench.cu) now includes an
error-free FP64 product.  Each M61 scalar is split into 31- and 30-bit limbs.
For two exactly represented 32-bit integers, `high=a*b` followed by
`error=fma(a,b,-high)` gives an exact nonoverlapping product; converting and
adding the two integral doubles recovers every product bit.  Three such products
feed the same independently validated radix-`2^31` Karatsuba fold as the integer
limb control.  The mixed quadratic multiply uses incumbent integer M61 products
for `ac` and `bd` and the FP64 product only for `(a+b)(c+d)`.

All 2,097,152 values match the incumbent after every timed chain.  Both FP64
kernels use 40 registers and no stack or spills; the incumbent uses 32.  With
eight dependent quadratic products, 31-sample medians are:

```text
incumbent 64x64 M61:             50.144 us
integer radix-2^31 limbs:        73.504 us   (1.466x)
all FP64-FMA limb products:    1034.016 us  (20.621x)
two integer + one FP64 product: 346.976 us   (6.920x)
```

The FP64 multiply/FMA/conversion sequence has far lower throughput on this
`sm_120` GPU than the compiler's 32-bit-IMAD lowering of the wide integer
product.  Independent product scheduling cannot hide a roughly twentyfold
primitive deficit.  Decision: **reject exact FP64 offload inside M61 before an
NTT tile or production integration.**  This also rules out a hybrid FP64
pipeline as the missing complement to the power-limited integer transform.

## Fused two-dimensional lazy M61 limbs

The persistent-radix-`2^31` rejection above ended with one explicit reopening
condition: a future limb representation must avoid normalizing each scalar
product independently across the quadratic-field multiply.  That condition had
not been tested.  The new gate combines the radix-`2^31` and Gaussian algebra
before reduction instead of repeating the old limb multiplier.

Write a centered M61 scalar as `a0 + B*a1`, `B=2^31`.  Because
`B^2=2 (mod M61)`, an unreduced scalar product is the coefficient pair

```text
(a0*b0 + 2*a1*b1, a0*b1 + a1*b0).
```

Balanced digits bound every coefficient in signed 64 bits.  The first exact
form uses three 32-bit Karatsuba products for each of `ac`, `bd`, and
`(a+b)(c+d)`, combines the three raw pairs into the real and imaginary results,
and only then performs two final Mersenne normalizations.  A single carry fold
is sufficient: after folding the low coefficient, the high digit is in
`[0,B)` and the remaining low digit is within roughly `(-1.5B,2.5B)`, so its
signed assembly fits 64 bits and needs at most one M61 correction.

The strengthened representation stores both quadratic components persistently
as four balanced 32-bit digits, retaining the incumbent 16-byte state size.  A
second variant removes modular normalization of the complex Karatsuba sums:
`ac` and `bd` use three products each, while the unnormalised joined term uses
four schoolbook products.  This **ten-product** form eliminates two balanced
add/recenter paths and is the fastest result.  It is materially different from
both earlier limb gates, which reduced three scalar M61 products separately.

[`src/cuda/m61_limb_bench.cu`](src/cuda/m61_limb_bench.cu) checks all
2,097,152 quadratic values after the accumulated warm-up/timed chains against
the incumbent, and independently recomputes the first complete chain with host
`unsigned __int128`.  All variants pass.  The best
persistent-ten-product kernel uses 34 registers with no stack or spills, versus
32 for the incumbent.  Fresh 31-sample medians are:

| Dependent quadratic products | incumbent wide M61 | persistent fused ten-product | candidate/incumbent |
|---:|---:|---:|---:|
| 1 | 15.712 us | 16.448 us | 1.047x |
| 8 | 50.400 us | 56.800 us | 1.127x |
| 32 | 175.648 us | 208.256 us | 1.186x |

This is substantially better than the old persistent-limb result, which was
about 2.9x at long chains, and confirms that fusing both algebra dimensions was
the right strengthening.  It still loses before charging limb-wise butterfly
addition, subtraction, rotation, or transform integration.  Final SASS has only
ten signed `IMAD.WIDE` products in the loop but pays a longer serial network of
signed shifts, carry folds, and 64-bit additions; reducing the multiply count
does not reduce the critical dependency graph enough.

Decision: **reject the persistent fused-limb M61 representation at the exact
arithmetic gate.**  Do not repeat separately-normalized limb transforms.  A
future reopening would need to remove the remaining cross-limb normalization
from an entire radix group, not merely from one quadratic product; the current
candidate is already 12--19% behind before those group operations.

## Exact q24 FP32 short-quotient reduction

The q24 architecture was already rejected above after its compensated FP32
reducer recovered only about 4 us when overlapped with M61.  Before reopening
it, the registry and both q24 benchmarks were checked: every previous version
formed a quotient-error correction with one extra FMA, multiply, and add.  The
following shorter quotient selection was therefore new rather than a repeat.

For `q = 14680063` and centered integral operands `|a|,|b| <= 7340031`, let
`T=a*b`, `high=RN(T)`, `low=fma(a,b,-high)`, and let `inverseQ` be the binary32
value of `1/q`.  The old reducer used `low` to correct `high*inverseQ` before
rounding the quotient.  It is not necessary for this particular q.  Direct
error bounds give

```text
|high - T|                                  <= 2^21
|inverseQ - 1/q| = 1.5950962535145825e-15
rounding error in high*inverseQ              <= 1/8
|RN(high*inverseQ) - T/q|                    < 0.353795.
```

Thus `n=nearbyint(high*inverseQ)` guarantees `|T/q-n| < 0.853795`.
The exact remainder has magnitude below 12,533,760, and the largest
`high-n*q` FMA intermediate is below 14,630,912.  Both are comfortably below
`2^24`, so the error-free `low` term is needed only when forming the final
remainder:

```c++
high = a*b;
low = fma(a,b,-high);
n = nearbyint(high*inverseQ);
r = fma(-n,q,high) + low;
```

One q correction then returns the exact centered result.  This is a proof for
this q and centered range, not a license to use the shortcut for arbitrary
24-bit primes.

[`src/cuda/riesel_lazy_tile_bench.cu`](src/cuda/riesel_lazy_tile_bench.cu)
now retains corrected and short-quotient variants side by side.  The short
variant passed all 132,120,567 scalar pairs formed by enumerating every
centered `a` against both endpoints, endpoint neighbors, `-1`, `0`, `1`, and
two independent affine-hash `b` values.  It also matched Montgomery for every
output of the 2,097,152-value radix-8/square tests through four rounds.

Fresh production-population tile medians were:

| Four radix-8/square rounds | median |
|---|---:|
| corrected exact FP32 | 20 us |
| short-quotient exact FP32 | **18 us** |
| canonical integer qC Montgomery | 16 us |
| Harvey integer qC Montgomery | 15 us |

The shortcut is a real 10% FP32 primitive improvement but still does not beat
the integer qC tile.  Both FP32 tiles use 40 registers without stack or spills.

The stronger test replaced the compensated reducer in
[`src/cuda/q24_m61_overlap_bench.cu`](src/cuda/q24_m61_overlap_bench.cu) and
reran three fresh processes.  At the transform-like chain length of eight,
the isolated q24 kernel improved from `26.590--26.595 us` to
`22.197--22.226 us`.  Its q24+M61 population time was
`110.016--110.328 us`, versus `119.787--120.048 us` for M31+M61: a repeatable
`9.459--9.962 us` advantage.  The old corrected q24 reducer had recovered only
`4.000--4.320 us` in the same gate.  At chain 16, the new isolated q24 time was
`39.888--39.950 us` and its overlap advantage was `18.197--18.251 us`.

For chain eight, the short kernel uses 19 registers versus 22 for the
corrected kernel, with no spills in either.  Final SASS shrinks from 944 to 800
instructions, including 48 fewer `FFMA` and 24 fewer `FADD` instructions.

This improves an exact reusable primitive, but it does **not** reverse the
architecture decision.  The relevant eight-product overlap gain is about
9.7 us, only 38% of the 25.5 us needed to move 205.5 to 180 us.  Moreover this
proxy omits q24's generic transform roots, generic Crandall--Fagin weights, and
less favorable CRT, while production M31 obtains unusually cheap rotations,
weights, and carry constants.  All omitted whole-transform effects are costs,
not sources for the missing 16 us.  Decision: **retain the short exact reducer
as a component result, but keep q24 times M61 rejected for the 180-us gate.**
Do not repeat the compensated q24 quotient or infer a full-engine speedup from
the improved isolated FP32 timing.

## Two-block DSM M61 resident-transform gate

The earlier M61 resident experiment above was checked before reopening this
route.  Its one 256-thread block retained an entire `8x512` quadratic tile in
64 KiB of shared memory, but serialized all eight height lines and measured
2--3% slower than three cache-resident kernels.  That result explicitly left
one architectural distinction untested: split the tile across a thread-block
cluster so that each block uses less shared memory and multiple line groups
execute concurrently, while distributed shared memory carries only the
factorization boundary.

[`src/cuda/m61_resident_tile_bench.cu`](src/cuda/m61_resident_tile_bench.cu)
now includes this exact DSM design.  A two-block cluster assigns four complete
512-point lines and 32 KiB of shared memory to each 256-thread block.  The first
forward and last inverse radix-8 middle stages exchange the opposite four lines
through DSM; every 512-point forward/square/inverse transform remains local.
This removes the same two global intermediate boundaries as the rejected
one-block kernel without making one block serialize all eight lines.

Random round trips, the independently known sparse square of `[1,2]`, and every
value after the accumulated timed squares match the three-kernel exact M61
control.  The clustered square kernel uses 80 registers, 32 KiB shared per
block, one CTA barrier resource, and no stack or spills.  The one-block resident
kernel uses 82 registers and 64 KiB; the separate middle/height controls use
`96/40/94` registers.

Five fresh 21-sample processes gave:

| Exact 2,097,152-value M61 core | Median range | Relative to separate |
|---|---:|---:|
| three-kernel control | 344.288--345.632 us | 1.000 |
| one-block 64-KiB resident | 351.968--353.824 us | 1.021--1.025 |
| **two-block 32-KiB DSM cluster** | **306.240--307.328 us** | **0.886--0.890** |

The cluster reproducibly saves `37.728--39.264 us`, or about **11%**, whereas
the old one-block fusion remains slower.  This is the first evidence here that
Blackwell DSM can profitably remove PRPLL-like transform boundaries: residency
alone was insufficient; residency plus restored line-level parallelism is the
mechanism.

Two strengthenings delimit the useful layout.  A four-block cluster gives each
block two lines and 16 KiB but adds a second cross-block middle stage; it is
exact, spill-free at 46 registers, and measures 381.888 us versus a nearby
348.608-us control (**1.096x**, a decisive loss).  Splitting the two-block
cluster across the 512-point height rather than across the eight middle lines
keeps radix 8 local and sends only the first/last height stage through DSM.  It
is also exact and spill-free at 76 registers, but is effectively tied with the
winning layout: 311.584 versus 310.688 us in the same run.  Thus two blocks are
the useful cluster size; merely increasing cluster parallelism is not a win.

This is a **passing component gate, not yet an end-to-end speedup**.  A direct
ratio forecast would reduce the measured 111.9-us M61 bottom-half sum by about
12 us.  That is large enough to justify a production-shaped gate because the
ideal `M61 core + fused edge` floor would move from 185.7 us toward 174 us, but
there is little allowance for lost M31 overlap.  The next implementation must
therefore include production M61 range-redundant arithmetic, middle/height
twiddles, transposes, the Hermitian pointwise square, and simultaneous M31
work.  Do not count the 11% proxy result as a PRPLL speedup until that exact
dependency-closed path is measured.

### Hermitian-pair strengthening: reject production DSM residency

The required follow-up exposed a dependency omitted by the passing gate above.
Its independent cyclic square keeps one `8x512` tile live (64 KiB), whereas
PRPLL's `onePairSq` couples each spectrum line to its reversed conjugate line.
A resident implementation must therefore hold **two** complete tiles, 128 KiB,
not one.  This was checked before attempting runtime cluster-launch integration.

The benchmark now ports the canonical exact type-zero `onePairSq` algebra:
`X2(a,conj(b))`, two field squares, multiplication by `t_squared`, `2ab`, the
second `X2`/conjugation, and the production component swap.  Tile B is reversed
in shared memory before pairing and reversed back afterwards, matching the
actual `revCrossLine` dependency.  A spill-free 512-thread control processes
two conjugate 512-point lines concurrently in 16 KiB, with the middle forward
and inverse remaining separate.  The resident candidate uses a four-block
cluster: ranks 0/1 retain tile A and ranks 2/3 tile B, each block owning four
height lines in 32 KiB.  Radix eight crosses the two blocks of each tile and
the pair square crosses corresponding blocks through DSM.

Every value from the separate and resident paths matches after both an isolated
application and the accumulated timed applications.  The double-wide control
uses 48 registers; the clustered kernel uses 80 registers and one barrier.
Neither has a stack frame or spills.  Across five fresh 21-sample processes:

| Exact 2,097,152-value Hermitian-pair core | Median range | Relative |
|---|---:|---:|
| separate middle / double-wide tail / middle | **238.016--239.680 us** | 1.000 |
| four-block, 128-KiB paired DSM residency | 366.144--368.416 us | **1.536--1.547** |

The dependency-correct resident layout loses `128.000--130.272 us`, roughly
54%.  The control exposes all 2,048 line pairs as independent blocks; the
cluster instead serializes four 512-point lines per block and adds cluster-wide
synchronization.  Eliminating two global tile boundaries does not compensate.
A two-block version with one 64-KiB tile per block would serialize all eight
lines per block and inherits the already measured loss of the one-block
single-tile layout, so it is not a stronger alternative.

Decision: **reject full M61 transform residency for the production Hermitian
tail.**  The earlier 11% result remains a valid single-tile DSM primitive, but
its direct 12-us PRPLL forecast is invalid because it omitted the conjugate
tile.  Do not integrate cluster launch support or production lazy arithmetic
for this layout, and do not quote the single-tile timing as an engine speedup.
Future DSM work would need a substantially smaller live representation or a
way to preserve the existing line-pair block parallelism; merely distributing
the two complete tiles cannot pass the gate.

## Two-61-bit-field 3M architecture

After the Hermitian DSM rejection, the registry was searched for a 3M design
using two quadratic 61-bit fields rather than the previously tested
`M31+M61+q31` three-field system.  This is a real representation-level
distinction: two 16-byte quadratic values over 1.5M packed positions retain the
same 48-MiB state, supply roughly 122 CRT bits, and avoid transforming two
separate 32-bit fields.  No earlier local benchmark included its shorter
population or its two-stream overlap against production.

A deterministic 64-bit prime search jointly required:

1. `q+1` divisible by `2^19` for the three packed 1M-word channels;
2. `3 | q-1` for the scalar Good--Thomas radix three;
3. `2^((q^2-1)/gcd(3*2^20,q^2-1)) = 1 (mod q)` for the exact
   Crandall--Fagin weight; and
4. the smallest pseudo-Mersenne complement to make reduction as favorable as
   possible.

The first useful candidate is

```text
q = 2^61 - 72,351,745 = 2,305,843,009,141,342,207.
```

Deterministic Miller--Rabin and `openssl prime` both classify it as prime;
`q+1` is divisible by `2^20`, `q-1` is divisible by three, and the weight test
returns one.  `M61*q/2` has just under 121 usable balanced bits, versus the
existing conservative p150 3M coefficient bound of about 118.17 bits, so this
is a general 140--150M-capable range rather than a p136 sparse-error trick.

The early arithmetic test first compared exact canonical pseudo-Mersenne fold
and radix-`2^64` Montgomery reducers over 1,572,864 quadratic values.  At eight
dependent products the fold is 1.565x M61 and Montgomery is 1.494x; Montgomery
is therefore the candidate representation.  It uses 36 registers versus 32
for M61, with no stack or spills.

[`src/cuda/m61_near61_overlap_bench.cu`](src/cuda/m61_near61_overlap_bench.cu)
then performs the decisive matched gate.  It compares concurrent M31/M61
kernels at the real 4M population (2,097,152 quadratic values per field) with
concurrent M61/q kernels at the shorter 3M population (1,572,864 per field).
The q stream remains in Montgomery form.  Every tested output agrees with an
independent host `unsigned __int128` calculation.  Five fresh 31-sample
processes gave:

| Dependent quadratic products | 4M M31+M61 | 3M M61+q | Candidate/control |
|---:|---:|---:|---:|
| 2 | 77.770--78.536 us | 77.061--77.662 us | 0.981--0.999 |
| 4 | 84.379--85.752 us | 86.458--87.979 us | 1.014--1.030 |
| 8 | 116.464--117.587 us | 129.045--130.178 us | **1.100--1.111** |
| 16 | 195.208--196.240 us | 219.514--220.549 us | **1.121--1.127** |

The one-load/two-product case nearly ties because the candidate moves 25% less
data.  As soon as arithmetic reuse approaches an NTT tile, two wide modular
streams contend for the same integer pipelines and the generic field loses
10--13% despite its shorter population.  A full transform would add the radix-3
edge, generic q roots/weights, and a wider 122-bit CRT/carry; none can reverse
an already negative transform-density gate or produce the required 12.4%
whole-iteration gain.

Decision: **reject the two-61-bit 3M architecture before a radix tile or PRPLL
integration.**  Its algebra and capacity are sound, but replacing the cheap
M31 stream by a second wide field removes the incumbent's complementary
arithmetic scheduling.  Do not repeat this design with a larger complement;
the selected q is already the most reduction-favorable compatible candidate
found, and generic Montgomery was faster than its pseudo-Mersenne fold.

## Width ownership reduction at fixed `512x8x512`

After the two-61-bit rejection, the registry was checked for an experiment
that changed the number of width values owned by each thread without changing
the exact production transform shape.  The existing shape sweeps changed both
the width factorization and ownership together; the one-warp tail work instead
increased ownership.  There was no prior test of fewer values per thread in
the production `512x8x512` M31/M61 pipeline.

The first scaffold overrode `NW=8` with `NW=4`, giving 128 threads and four
values per thread.  It exposed two assumptions that made its early timings
invalid:

1. when width and height are both 512, the production radix-eight width table
   aliases the larger height/tail trig table; a different width radix needs a
   standalone table and different per-field offsets; and
2. the radix-four padded shared-memory formulas were sized for the original
   256-point case.  Compute Sanitizer found out-of-bounds shared reads in
   `fftP`; disabling padding removed the memory error.

Even after both layout issues were removed, radix four is not an exact
factorization of a 512-point transform.  The generic loop performs five
radix-four levels, algebraically describing 1024 points.  Its apparently fast
output was therefore an incorrect timing scaffold and must not be quoted as a
candidate result.  A correct mixed-radix implementation would require a real
binary coupling stage, not just a smaller `NW` define.

The decisive implementation instead adds a true radix-two width transform:
`NW=2`, `G_W=256`, two values per thread, nine exact binary levels, and the
compact unpadded shuffle layout.  It retains the production 4M coefficient
population, M31/M61 fields, middle/height transforms, Hermitian square, CRT,
and carry.  It matched the production residues at every checked point:

```text
iteration  2,000: 05d6515c416b83e2
iteration 10,000: 52316d51aa52e6b7
iteration 20,000: 6a5c8b8989125413
iteration 40,000: b597cca6031938ff
iteration 50,000: cce5a14b7c17aecb
```

A matched pair of fresh 50,000-iteration runs on the RTX PRO 6000 Blackwell
Max-Q measured:

| Exact p136 path | Final timing | Relative |
|---|---:|---:|
| production radix eight, eight values/thread | **202.3 us** | 1.000 |
| radix two, two values/thread | 264.8 us | **1.309** |

Event profiling explains the regression.  The kernels outside the fused edge
were essentially unchanged, but `kCarryFused` increased from 71.0 to 136.1
us/call.  Reducing live vectors does lower per-thread state, but it doubles the
thread population and replaces three radix-eight levels with nine
shuffle/barrier/twiddle levels inside the already critical fused carry kernel.
Changing `WMUL` from two to one did not recover the loss (259.7 versus 260.6
us in matched 10,000-iteration smoke runs).

Decision: **reject reduced width ownership for this production shape.**  The
exact radix-two prototype is useful evidence that register ownership is not
the limiting resource: the incumbent's coarse radix and low synchronization
count are substantially more valuable.  Do not implement the missing
radix-four-by-radix-two mixed path unless a new design removes rather than
adds fused-edge synchronization; its best plausible behavior lies between
the rejected binary path and the already faster radix-eight incumbent.
