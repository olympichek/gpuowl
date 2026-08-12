#!/usr/bin/env python3
"""Carve cubin ELFs out of driver JIT-cache blobs and nvdisasm them.

Usage: e2carve.py <e2-sass-dir>
For each jit-<tag>/ subdir: finds \\x7fELF magics in every cache file, carves
[magic..next magic) slices, nvdisasm's each, names output by the kernel entry
functions found in the disassembly -> sass-<tag>/<kernels>.sass.txt
"""
import os, re, subprocess, sys, glob

NVDISASM = "/usr/local/cuda/bin/nvdisasm"

def carve(blob: bytes):
    pos, out = [], []
    i = blob.find(b"\x7fELF")
    while i != -1:
        pos.append(i)
        i = blob.find(b"\x7fELF", i + 4)
    for n, p in enumerate(pos):
        end = pos[n + 1] if n + 1 < len(pos) else len(blob)
        out.append(blob[p:end])
    return out

def main(root):
    for jc in sorted(glob.glob(os.path.join(root, "jit-*"))):
        tag = os.path.basename(jc)[4:]
        outdir = os.path.join(root, f"sass-{tag}")
        os.makedirs(outdir, exist_ok=True)
        n_ok = n_fail = 0
        for path in glob.glob(os.path.join(jc, "**", "*"), recursive=True):
            if not os.path.isfile(path):
                continue
            blob = open(path, "rb").read()
            for k, cub in enumerate(carve(blob)):
                tmp = os.path.join(outdir, f"_tmp{k}.cubin")
                with open(tmp, "wb") as f:
                    f.write(cub)
                r = subprocess.run([NVDISASM, "-c", tmp], capture_output=True, text=True)
                os.remove(tmp)
                if r.returncode != 0:
                    n_fail += 1
                    continue
                sass = r.stdout
                names = sorted(set(re.findall(r"\.text\.(\w+)", sass)))
                if not names:
                    n_fail += 1
                    continue
                name = names[0] if len(names) == 1 else f"{names[0]}_and_{len(names)-1}more"
                out = os.path.join(outdir, f"{name}.sass.txt")
                with open(out, "w") as f:
                    f.write(sass)
                n_ok += 1
        print(f"{tag}: {n_ok} cubins disassembled, {n_fail} skipped")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "e2-sass")
