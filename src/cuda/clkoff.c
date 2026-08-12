// clkoff.c — query/set the NVML SM/graphics clock V-f offset (device 0, P0).
//
// Usage: ./clkoff              query only: print supported types + allowed ranges
//        sudo ./clkoff <MHz>   set the offset (0 restores stock)
// Build: gcc -I/usr/local/cuda/include -o clkoff src/cuda/clkoff.c -lnvidia-ml
//
// Campaign context: on the RTX PRO 6000 Blackwell SERVER Edition the vBIOS
// reports GRAPHICS P0 range [0, 0] (SM type unsupported) — overclocking is
// fused off on that SKU.  Workstation-edition GB202 boards allow offsets;
// +60 MHz was projected worth ~-3.6 us/iter on the production workload.

#include <stdio.h>
#include <stdlib.h>
#include <nvml.h>

#define CK(x) do { nvmlReturn_t r = (x); if (r != NVML_SUCCESS) { \
  fprintf(stderr, "FAIL %s: %s\n", #x, nvmlErrorString(r)); nvmlShutdown(); exit(1); } } while (0)

static const char* clkName(nvmlClockType_t t) {
  return t == NVML_CLOCK_SM ? "SM" : t == NVML_CLOCK_GRAPHICS ? "GRAPHICS" : "?";
}

int main(int argc, char** argv) {
  int queryOnly = (argc < 2);
  int want = queryOnly ? 0 : atoi(argv[1]);

  CK(nvmlInit_v2());
  nvmlDevice_t dev;
  CK(nvmlDeviceGetHandleByIndex_v2(0, &dev));

  // Try SM first, fall back to GRAPHICS (boards expose one or the other).
  nvmlClockType_t types[2] = {NVML_CLOCK_SM, NVML_CLOCK_GRAPHICS};
  int rc = 1;
  for (int i = 0; i < 2; ++i) {
    nvmlClockOffset_t info = {0};
    info.version = nvmlClockOffset_v1;
    info.type = types[i];
    info.pstate = NVML_PSTATE_0;
    nvmlReturn_t r = nvmlDeviceGetClockOffsets(dev, &info);
    if (r != NVML_SUCCESS) {
      fprintf(stderr, "get %s offsets: %s\n", clkName(types[i]), nvmlErrorString(r));
      continue;
    }
    printf("%s P0: current %+d MHz, allowed [%d, %d]%s\n",
           clkName(types[i]), info.clockOffsetMHz, info.minClockOffsetMHz,
           info.maxClockOffsetMHz,
           (info.minClockOffsetMHz == 0 && info.maxClockOffsetMHz == 0)
             ? " (offset LOCKED by vBIOS)" : "");
    if (queryOnly) { rc = 0; continue; }
    if (want < info.minClockOffsetMHz || want > info.maxClockOffsetMHz) {
      fprintf(stderr, "requested %+d outside allowed range\n", want);
      nvmlShutdown();
      return 1;
    }
    info.clockOffsetMHz = want;
    r = nvmlDeviceSetClockOffsets(dev, &info);
    if (r != NVML_SUCCESS) {
      fprintf(stderr, "set %s offset: %s\n", clkName(types[i]), nvmlErrorString(r));
      continue;
    }
    nvmlClockOffset_t back = {0};
    back.version = nvmlClockOffset_v1;
    back.type = types[i];
    back.pstate = NVML_PSTATE_0;
    CK(nvmlDeviceGetClockOffsets(dev, &back));
    printf("%s P0 offset now %+d MHz\n", clkName(types[i]), back.clockOffsetMHz);
    nvmlShutdown();
    return 0;
  }
  nvmlShutdown();
  return rc;
}
