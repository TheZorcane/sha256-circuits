"""Converters used for the ABC control comparators (abc_control.py); not part of any proof.

xag_to_aag   : XAG -> ASCII AIGER (XOR -> 3 ANDs); unused by abc_control.py
xag_to_blif  : XAG -> BLIF (.names, own writer; use this one, not verifier/xag.py write_blif)
blif_to_xag  : ABC-mapped BLIF (.gate AND2/XOR2/XNOR2/INV/BUF/ZERO/ONE) or Yosys/ABC .names BLIF -> XAG
"""
import os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "verifier"))
from xag import XAG, AND, XOR, CircuitError  # noqa: E402


def xag_to_aag(x: XAG) -> str:
    n = x.n_inputs
    amap = [0] * x.n_vars()           # xag var -> aiger positive literal
    for i in range(n):
        amap[i + 1] = 2 * (i + 1)
    ands, nxt = [], n + 1

    def new(a, b):
        nonlocal nxt
        l = 2 * nxt
        nxt += 1
        ands.append((l, a, b))
        return l

    for k, (op, a, b) in enumerate(x.nodes):
        la = amap[a >> 1] ^ (a & 1)
        lb = amap[b >> 1] ^ (b & 1)
        if op == AND:
            amap[n + 1 + k] = new(la, lb)
        else:
            n1 = new(la, lb ^ 1)
            n2 = new(la ^ 1, lb)
            amap[n + 1 + k] = new(n1 ^ 1, n2 ^ 1) ^ 1
    outs = [amap[o >> 1] ^ (o & 1) for o in x.outputs]
    M = nxt - 1
    L = [f"aag {M} {n} 0 {len(outs)} {len(ands)}"]
    L += [str(2 * (i + 1)) for i in range(n)]
    L += [str(o) for o in outs]
    L += [f"{l} {a} {b}" for l, a, b in ands]
    L += [f"i{i} i{i}" for i in range(n)] + [f"o{i} o{i}" for i in range(len(outs))]
    return "\n".join(L) + "\n"


def xag_to_blif(x: XAG, model="xag") -> str:
    def nm(v):
        return "c0" if v == 0 else (f"i{v-1}" if v <= x.n_inputs else f"n{v}")
    out = [f".model {model}", ".inputs " + " ".join(f"i{i}" for i in range(x.n_inputs)),
           ".outputs " + " ".join(f"o{i}" for i in range(len(x.outputs))), ".names c0"]
    for k, (op, a, b) in enumerate(x.nodes):
        v = x.n_inputs + 1 + k
        pa, pb = ("0" if a & 1 else "1"), ("0" if b & 1 else "1")
        out.append(f".names {nm(a>>1)} {nm(b>>1)} {nm(v)}")
        if op == AND:
            out.append(f"{pa}{pb} 1")
        else:
            qa, qb = ("1" if pa == "0" else "0"), ("1" if pb == "0" else "0")
            out += [f"{pa}{qb} 1", f"{qa}{pb} 1"]
    for i, o in enumerate(x.outputs):
        out.append(f".names {nm(o>>1)} o{i}")
        out.append("0 1" if o & 1 else "1 1")
    out.append(".end")
    return "\n".join(out) + "\n"


def _idx(name, pref):
    m = re.fullmatch(pref + r"(\d+)", name)
    return int(m.group(1)) if m else None


def blif_to_xag(text: str) -> XAG:
    """Combinational BLIF with gates of <=2 inputs. Inputs/outputs ordered by trailing index
    of names i<k>/o<k> (our naming); raises if names do not follow it."""
    text = re.sub(r"\\\n", " ", text)
    lines = [l.split("#")[0].strip() for l in text.split("\n")]
    lines = [l for l in lines if l]
    inputs, outputs, defs = [], [], {}
    i = 0
    while i < len(lines):
        t = lines[i].split()
        if t[0] == ".inputs":
            inputs += t[1:]
        elif t[0] == ".outputs":
            outputs += t[1:]
        elif t[0] == ".gate":
            g = t[1]
            pins = dict(p.split("=") for p in t[2:])
            o = pins.get("O") or pins.get("Y") or pins.get("out")
            ins = [pins[k] for k in sorted(pins) if k not in ("O", "Y", "out")]
            if g in ("AND2", "XOR2", "XNOR2"):
                defs[o] = (g, ins[0], ins[1])
            elif g in ("INV", "BUF"):
                defs[o] = (g, ins[0])
            elif g in ("ZERO", "ONE"):
                defs[o] = (g,)
            else:
                raise CircuitError(f"unsupported gate {g}")
        elif t[0] == ".names":
            sig = t[1:]
            cubes = []
            while i + 1 < len(lines) and not lines[i + 1].startswith("."):
                i += 1
                cubes.append(lines[i].split())
            defs[sig[-1]] = ("TT", sig[:-1], cubes)
        elif t[0] in (".model", ".end"):
            pass
        else:
            raise CircuitError(f"unsupported BLIF directive {t[0]}")
        i += 1
    n = len(inputs)
    pos = {}
    for s in inputs:
        k = _idx(s, "(?:i|pi)")
        if k is None or k >= n:
            raise CircuitError(f"unexpected input name {s}")
        pos[s] = k
    if sorted(pos.values()) != list(range(n)):
        raise CircuitError("input names not a permutation of i0..")
    x = XAG(n)
    lit = {s: x.input(k) for s, k in pos.items()}

    def tt_gate(ins, cubes):
        nin = len(ins)
        if nin > 2:
            raise CircuitError(">2-input .names")
        rows = []
        for v in range(1 << nin):
            bits = [(v >> (nin - 1 - j)) & 1 for j in range(nin)]
            val = 0
            for c in cubes:
                pat, o = (c[0], c[1]) if nin else ("", c[0])
                if all(p == "-" or int(p) == b for p, b in zip(pat, bits)):
                    val = 1 if o == "1" else 0
            rows.append(val)
        return rows

    def resolve(sig):
        if sig in lit:
            return lit[sig]
        stack = [sig]
        expanded = set()
        while stack:
            s = stack[-1]
            if s in lit:
                stack.pop()
                continue
            if s not in defs:
                raise CircuitError(f"undriven signal {s}")
            d = defs[s]
            if d[0] == "TT":
                deps = d[1]
            elif d[0] in ("AND2", "XOR2", "XNOR2", "INV", "BUF"):
                deps = d[1:]
            else:
                deps = []
            miss = [q for q in deps if q not in lit]
            if miss:
                # a dependency that is expanded but unresolved lies on the current DFS path: a cycle.
                # (A dependency merely pushed by a sibling is not; reconvergent fan-in is allowed.)
                for q in miss:
                    if q in expanded:
                        raise CircuitError("combinational cycle")
                expanded.add(s)
                stack += miss
                continue
            if d[0] == "AND2":
                r = x.AND(lit[d[1]], lit[d[2]])
            elif d[0] == "XOR2":
                r = x.XOR(lit[d[1]], lit[d[2]])
            elif d[0] == "XNOR2":
                r = x.XOR(lit[d[1]], lit[d[2]]) ^ 1
            elif d[0] == "INV":
                r = lit[d[1]] ^ 1
            elif d[0] == "BUF":
                r = lit[d[1]]
            elif d[0] == "ZERO":
                r = 0
            elif d[0] == "ONE":
                r = 1
            else:
                ins, cubes = d[1], d[2]
                rows = tt_gate(ins, cubes)
                nin = len(ins)
                if nin == 0:
                    r = rows[0]
                elif nin == 1:
                    a = lit[ins[0]]
                    r = {(0, 0): 0, (1, 1): 1, (0, 1): a, (1, 0): a ^ 1}[tuple(rows)]
                else:
                    a, b = lit[ins[0]], lit[ins[1]]
                    pc = sum(rows)
                    if pc == 0:
                        r = 0
                    elif pc == 4:
                        r = 1
                    elif rows == [0, 1, 1, 0]:
                        r = x.XOR(a, b)
                    elif rows == [1, 0, 0, 1]:
                        r = x.XOR(a, b) ^ 1
                    elif pc == 1 or pc == 3:
                        inv = pc == 3
                        k = rows.index(0 if inv else 1)
                        ba, bb = (k >> 1) & 1, k & 1
                        r = x.AND(a ^ (1 - ba), b ^ (1 - bb)) ^ (1 if inv else 0)
                    else:
                        r = {(0, 0, 1, 1): a, (1, 1, 0, 0): a ^ 1, (0, 1, 0, 1): b, (1, 0, 1, 0): b ^ 1}[tuple(rows)]
            lit[s] = r
            stack.pop()
        return lit[sig]

    outs = {}
    for s in outputs:
        k = _idx(s, "(?:o|po)")
        if k is None:
            raise CircuitError(f"unexpected output name {s}")
        outs[k] = resolve(s)
    if sorted(outs) != list(range(len(outputs))):
        raise CircuitError("output names not a permutation of o0..")
    x.outputs = [outs[k] for k in range(len(outputs))]
    x.validate()
    return x
