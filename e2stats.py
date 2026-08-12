#!/usr/bin/env python3
"""E2 SASS-stats: per-kernel opcode-class histograms + .reuse density across
NVRTC versions, diffed along the E1 clock-cost table's axes.

Usage: e2stats.py <e2-sass-dir> [kernelFilter]
Reads sass-<tag>/ subdirs; prints per-kernel per-version rows and deltas.
"""
import os, re, sys, glob, collections

# Class bins, ordered roughly by E1 relevance (int pipes first).
CLASSES = [
    ("IMAD.WIDE", r"^U?IMAD\.WIDE"),
    ("IMAD.MOV",  r"^U?IMAD\.MOV"),
    ("IMAD",      r"^U?IMAD(?!\.WIDE|\.MOV)"),
    ("IADD64",    r"^U?IADD\.64|^U?IADD3\.64"),
    ("IADD3",     r"^U?IADD3(?!\.64)"),
    ("LOP3",      r"^U?LOP3"),
    ("SHF",       r"^U?SHF"),
    ("PRMT",      r"^U?PRMT"),
    ("SEL",       r"^U?SEL"),
    ("MOV",       r"^U?MOV"),
    ("FFMA",      r"^U?FFMA|^HFMA2"),
    ("LDS",       r"^LDS"),
    ("STS",       r"^STS"),
    ("LDG",       r"^LDG|^LD\.E"),
    ("STG",       r"^STG|^ST\.E"),
    ("LDCx",      r"^U?LDCU?\b|^U?LDC\b|^U?LDC\."),
    ("ISETP",     r"^U?ISETP"),
    ("BRA",       r"^BRA|^JMP|^RET|^EXIT|^BSSY|^BSYNC"),
    ("NOP",       r"^NOP"),
    ("BAR",       r"^BAR|^DEPBAR"),
]
CRE = [(n, re.compile(p)) for n, p in CLASSES]

def analyze(path):
    ops = collections.Counter()
    reuse = 0
    uniform = 0
    total = 0
    for line in open(path):
        m = re.search(r"/\*[0-9a-f]{4,}\*/\s+(@!?U?P\w+\s+)?([A-Z][A-Z0-9.]*)", line)
        if not m:
            continue
        op = m.group(2)
        total += 1
        if op.startswith("U") and not op.startswith("UN"):  # uniform-datapath op (crude)
            uniform += 1
        reuse += line.count(".reuse")
        for name, cre in CRE:
            if cre.match(op):
                ops[name] += 1
                break
        else:
            ops["other:" + op.split(".")[0]] += 1
    return total, ops, reuse, uniform

def main(root, filt=None):
    dirs = sorted(glob.glob(os.path.join(root, "sass-*")))
    tags = [os.path.basename(d)[5:] for d in dirs]
    kernels = collections.defaultdict(dict)  # kernel -> tag -> stats
    for d, tag in zip(dirs, tags):
        for f in glob.glob(os.path.join(d, "*.sass.txt")):
            k = os.path.basename(f)[:-9]
            if filt and filt not in k:
                continue
            kernels[k][tag] = analyze(f)
    classcols = [n for n, _ in CLASSES if n not in ("NOP", "BAR", "BRA", "ISETP", "LDCx")]
    for k in sorted(kernels):
        print(f"\n== {k}")
        rows = kernels[k]
        base_tag = tags[0] if tags[0] in rows else sorted(rows)[0]
        for tag in tags:
            if tag not in rows:
                continue
            total, ops, reuse, uniform = rows[tag]
            cols = " ".join(f"{c}:{ops.get(c,0)}" for c in classcols if ops.get(c, 0))
            others = sum(v for kk, v in ops.items() if kk.startswith("other:"))
            print(f"  {tag:5s} n={total:5d} reuse={reuse:4d} uni={uniform:4d} other={others:3d}  {cols}")
        if len(rows) > 1 and base_tag in rows:
            bt, bops, bre, bun = rows[base_tag]
            for tag in tags:
                if tag == base_tag or tag not in rows:
                    continue
                t, o, r, u = rows[tag]
                dc = {c: o.get(c, 0) - bops.get(c, 0) for c in set(list(o) + list(bops))}
                dd = " ".join(f"{c}:{v:+d}" for c, v in sorted(dc.items(), key=lambda x: -abs(x[1])) if v)
                print(f"  d({tag}-{base_tag}): n={t-bt:+d} reuse={r-bre:+d} uni={u-bun:+d}  {dd[:200]}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "e2-sass",
         sys.argv[2] if len(sys.argv) > 2 else None)
