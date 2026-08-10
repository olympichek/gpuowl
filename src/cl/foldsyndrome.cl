// Diagnostic support for the sparse-error folded M31 sidecar.
// This kernel is intentionally independent of carryFused's incremental
// weight counter: it recomputes each logical word's weight from its index.

#include "base.cl"
#include "math.cl"
#include "weight.cl"

#if FOLD_SYNDROME_VALIDATE

Z31 foldWeightWord(Word word, u32 wordIndex) {
  const u32 log2RootTwo = (u32) (((1ULL << 30) / NWORDS) % 31);
  const u32 bigwordShift = (NWORDS - EXP % NWORDS) * log2RootTwo % 31;
  const u32 shiftStep = (bigwordShift + 30) % 31;
  u64 combo = comboFracBits(wordIndex) +
              make_u64(wordIndex * shiftStep, 0xFFFFFFFF);
  return shl(make_Z31(word), hi32(combo) % 31);
}

KERNEL(256) foldValidate(CP(Word2) words, CP(GF31) folded,
                         P(u32) mismatches) {
  u32 const basePair = get_global_id(0);
  if (basePair >= NWORDS / 16) return;

  GF31 expected = U2((Z31)0, (Z31)0);
  for (u32 alias = 0; alias < 8; ++alias) {
    u32 const pair = basePair + alias * (NWORDS / 16);
    Word2 const value = words[pair];
    expected = add(expected,
                   U2(foldWeightWord(value.x, 2 * pair),
                      foldWeightWord(value.y, 2 * pair + 1)));
  }

  GF31 const actual = folded[basePair];
  if (get_Z31(expected.x) != get_Z31(actual.x) ||
      get_Z31(expected.y) != get_Z31(actual.y)) {
    atomic_add(mismatches, 1u);
  }
}

#endif
