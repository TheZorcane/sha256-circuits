"""Circuit IR, parsers, validator and writers.

Two layers:

RawCircuit  -- exactly what a file says: named wires, gates with their original
               type (AND, XOR, INV, XNOR, NOT, MAND, EQ, EQW, ...).  Used for raw
               (pre-canonical) gate counts and structural validation.
XAG         -- canonical IR: 2-input AND / XOR nodes over literals with a
               complement bit.  NOT is free (a complemented edge), constants are
               free.

Literal encoding (AIGER-like): lit = 2*var + neg.  var 0 is constant FALSE
(lit 0 = FALSE, lit 1 = TRUE); vars 1..n_inputs are the primary inputs in
order; var n_inputs+1+k is node k.  Nodes are stored in topological order:
every operand of node k refers to a var < n_inputs+1+k.
"""

import re
from dataclasses import dataclass, field


class CircuitError(ValueError):
    pass


# --------------------------------------------------------------------------- #
# XAG
# --------------------------------------------------------------------------- #
AND, XOR = "AND", "XOR"


@dataclass
class XAG:
    n_inputs: int
    nodes: list = field(default_factory=list)    # (op, lit_a, lit_b)
    outputs: list = field(default_factory=list)  # lits

    # ---- builder helpers (no simplification; see cost.canonicalize) --------
    def input(self, i):
        if not 0 <= i < self.n_inputs:
            raise CircuitError(f"input {i} out of range")
        return 2 * (i + 1)

    def _add(self, op, a, b):
        top = 2 * (self.n_inputs + 1 + len(self.nodes))
        if not (0 <= a < top and 0 <= b < top):
            raise CircuitError("operand refers to undefined node")
        self.nodes.append((op, a, b))
        return top

    def AND(self, a, b):
        return self._add(AND, a, b)

    def XOR(self, a, b):
        return self._add(XOR, a, b)

    @staticmethod
    def NOT(a):
        return a ^ 1

    def OR(self, a, b):
        return self.AND(a ^ 1, b ^ 1) ^ 1

    def n_vars(self):
        return self.n_inputs + 1 + len(self.nodes)

    def counts(self):
        a = sum(1 for n in self.nodes if n[0] == AND)
        return {"AND": a, "XOR": len(self.nodes) - a, "total": len(self.nodes)}

    def validate(self):
        base = self.n_inputs + 1
        for k, (op, a, b) in enumerate(self.nodes):
            if op not in (AND, XOR):
                raise CircuitError(f"node {k}: unknown op {op!r}")
            lim = 2 * (base + k)
            for x in (a, b):
                if not isinstance(x, int) or not 0 <= x < lim:
                    raise CircuitError(f"node {k}: operand {x} not defined before use")
        lim = 2 * self.n_vars()
        for o in self.outputs:
            if not isinstance(o, int) or not 0 <= o < lim:
                raise CircuitError(f"output literal {o} undefined")
        return True


# ---- native text format (.xag) ---------------------------------------------
# xag 1
# inputs <n>
# nodes <m>
# A <a> <b> | X <a> <b>      (one per node, literals)
# outputs <k>
# <lit> ...                  (one per line)

def write_xag(x: XAG) -> str:
    out = ["xag 1", f"inputs {x.n_inputs}", f"nodes {len(x.nodes)}"]
    out += [f"{'A' if op == AND else 'X'} {a} {b}" for op, a, b in x.nodes]
    out.append(f"outputs {len(x.outputs)}")
    out += [str(o) for o in x.outputs]
    return "\n".join(out) + "\n"


def read_xag(text: str) -> XAG:
    toks = text.split()
    try:
        if toks[:2] != ["xag", "1"]:
            raise CircuitError("not an xag v1 file")
        i = 2
        assert toks[i] == "inputs"; n_in = int(toks[i + 1]); i += 2
        assert toks[i] == "nodes"; m = int(toks[i + 1]); i += 2
        x = XAG(n_in)
        for _ in range(m):
            t, a, b = toks[i], int(toks[i + 1]), int(toks[i + 2]); i += 3
            if t == "A":
                x.AND(a, b)
            elif t == "X":
                x.XOR(a, b)
            else:
                raise CircuitError(f"unknown node tag {t!r}")
        assert toks[i] == "outputs"; k = int(toks[i + 1]); i += 2
        x.outputs = [int(t) for t in toks[i:i + k]]
        if len(x.outputs) != k or len(toks) != i + k:
            raise CircuitError("bad output section")
    except (IndexError, AssertionError, ValueError) as e:
        raise CircuitError(f"malformed xag file: {e}") from None
    x.validate()
    return x


# --------------------------------------------------------------------------- #
# RawCircuit
# --------------------------------------------------------------------------- #
# arity table: type -> (n_in, n_out) ; MAND is variable (checked separately)
ARITY = {
    "AND": (2, 1), "XOR": (2, 1), "XNOR": (2, 1), "OR": (2, 1),
    "NAND": (2, 1), "NOR": (2, 1),
    "INV": (1, 1), "NOT": (1, 1), "EQW": (1, 1),
    "EQ": (0, 1),   # constant: value stored in Gate.const
}


@dataclass
class Gate:
    type: str
    ins: list
    outs: list
    const: int = None   # only for EQ


@dataclass
class RawCircuit:
    inputs: list          # wire names, in file order
    outputs: list         # wire names, in file order
    gates: list           # Gate
    fmt: str = ""
    input_groups: list = None    # sizes of input groups as declared
    output_groups: list = None
    meta: dict = field(default_factory=dict)

    def raw_counts(self):
        """Counts by original gate type.  MAND counted both as #MAND lines and
        as the number of 2-input ANDs they contain (MAND_ands)."""
        c = {}
        for g in self.gates:
            c[g.type] = c.get(g.type, 0) + 1
            if g.type == "MAND":
                c["MAND_ands"] = c.get("MAND_ands", 0) + len(g.outs)
        return dict(sorted(c.items()))

    # ---- validation --------------------------------------------------------
    def validate(self):
        """Structural checks; returns a topological order of gate indices.
        Rejects: unknown types, wrong fan-in/out, multiply-driven wires,
        unbound (undriven) wires, undefined outputs, cycles."""
        driver = {}
        for w in self.inputs:
            if w in driver:
                raise CircuitError(f"input wire {w} declared twice")
            driver[w] = -1
        for k, g in enumerate(self.gates):
            if g.type == "MAND":
                n = len(g.outs)
                if n < 1 or len(g.ins) != 2 * n:
                    raise CircuitError(f"gate {k}: MAND needs 2n inputs / n outputs")
            elif g.type in ARITY:
                if (len(g.ins), len(g.outs)) != ARITY[g.type]:
                    raise CircuitError(f"gate {k}: {g.type} has wrong fan-in/fan-out "
                                       f"{len(g.ins)}/{len(g.outs)}")
                if g.type == "EQ" and g.const not in (0, 1):
                    raise CircuitError(f"gate {k}: EQ constant must be 0/1")
            else:
                raise CircuitError(f"gate {k}: unknown gate type {g.type!r}")
            for w in g.outs:
                if w in driver:
                    raise CircuitError(f"wire {w} driven more than once")
                driver[w] = k
        for k, g in enumerate(self.gates):
            for w in g.ins:
                if w not in driver:
                    raise CircuitError(f"gate {k}: input wire {w} is unbound")
        for w in self.outputs:
            if w not in driver:
                raise CircuitError(f"output wire {w} is undefined")
        # Kahn topological sort over gates
        indeg = [0] * len(self.gates)
        users = {}
        for k, g in enumerate(self.gates):
            for w in g.ins:
                d = driver[w]
                if d >= 0:
                    indeg[k] += 1
                    users.setdefault(d, []).append(k)
        ready = [k for k in range(len(self.gates)) if indeg[k] == 0]
        ready.reverse()
        order = []
        while ready:
            k = ready.pop()
            order.append(k)
            for u in reversed(users.get(k, [])):
                indeg[u] -= 1
                if indeg[u] == 0:
                    ready.append(u)
        if len(order) != len(self.gates):
            raise CircuitError("circuit contains a cycle")
        return order

    # ---- translation to XAG (no simplification) ----------------------------
    def to_xag(self) -> XAG:
        order = self.validate()
        x = XAG(len(self.inputs))
        lit = {w: x.input(i) for i, w in enumerate(self.inputs)}
        for k in order:
            g = self.gates[k]
            t, i = g.type, [lit[w] for w in g.ins]
            if t == "AND":
                r = [x.AND(i[0], i[1])]
            elif t == "XOR":
                r = [x.XOR(i[0], i[1])]
            elif t == "XNOR":
                r = [x.XOR(i[0], i[1]) ^ 1]
            elif t in ("INV", "NOT"):
                r = [i[0] ^ 1]
            elif t == "EQW":
                r = [i[0]]
            elif t == "EQ":
                r = [g.const]
            elif t == "OR":
                r = [x.OR(i[0], i[1])]
            elif t == "NAND":
                r = [x.AND(i[0], i[1]) ^ 1]
            elif t == "NOR":
                r = [x.OR(i[0], i[1]) ^ 1]
            elif t == "MAND":
                n = len(g.outs)
                r = [x.AND(i[j], i[n + j]) for j in range(n)]
            else:  # unreachable: validate() rejects unknown types
                raise CircuitError(t)
            for w, v in zip(g.outs, r):
                lit[w] = v
        x.outputs = [lit[w] for w in self.outputs]
        x.validate()
        return x


# --------------------------------------------------------------------------- #
# Parsers
# --------------------------------------------------------------------------- #
def _int_lines(text):
    return [ln.split() for ln in text.splitlines()]


def _parse_bristol_gates(lines, start, ngates):
    gates = []
    for ln in lines[start:]:
        if not ln:
            continue
        t = ln[-1]
        try:
            nin, nout = int(ln[0]), int(ln[1])
            nums = ln[2:-1]
            if t == "EQ":
                # "1 1 c w EQ": constant c assigned to wire w
                if (nin, nout, len(nums)) != (1, 1, 2):
                    raise CircuitError(f"bad EQ line: {' '.join(ln)}")
                gates.append(Gate("EQ", [], [int(nums[1])], const=int(nums[0])))
                continue
            if len(nums) != nin + nout:
                raise CircuitError(f"gate line arity mismatch: {' '.join(ln)}")
            vals = [int(v) for v in nums]
        except CircuitError:
            raise
        except ValueError:
            raise CircuitError(f"non-integer field in gate line: {' '.join(ln)}") from None
        gates.append(Gate(t, vals[:nin], vals[nin:]))
    if len(gates) != ngates:
        raise CircuitError(f"header says {ngates} gates, found {len(gates)}")
    return gates


def parse_bristol(text: str, fmt: str = "auto") -> RawCircuit:
    """Bristol circuits.

    fmt="old": original Bristol / SCAPI format (also used by Goldfeder):
        ngates nwires
        n_in1 n_in2 n_out
        <gates>   "nin nout in.. out.. TYPE"
      inputs = wires 0 .. n_in1+n_in2-1, outputs = last n_out wires.
    fmt="fashion": Bristol Fashion:
        ngates nwires
        niv i_0 .. i_{niv-1}
        nov o_0 .. o_{nov-1}
        <gates>  (adds INV, EQ, EQW, MAND)
      inputs = first sum(i) wires, outputs = last sum(o) wires.
    fmt="auto": Bristol Fashion iff line 3 consists only of integers and has
      exactly int(first token)+1 tokens; otherwise old format.
    """
    lines = _int_lines(text)
    nonblank = [i for i, ln in enumerate(lines) if ln]
    if len(nonblank) < 2:
        raise CircuitError("file too short")
    l0, l1 = lines[nonblank[0]], lines[nonblank[1]]
    try:
        ngates, nwires = int(l0[0]), int(l0[1])
        if len(l0) != 2:
            raise ValueError
    except (ValueError, IndexError):
        raise CircuitError("bad header line 1") from None
    if fmt == "auto":
        l2 = lines[nonblank[2]] if len(nonblank) > 2 else []
        is_fashion = (l2 and all(re.fullmatch(r"\d+", t) for t in l2)
                      and len(l2) == int(l2[0]) + 1
                      and len(l1) == int(l1[0]) + 1)
        fmt = "fashion" if is_fashion else "old"
    try:
        if fmt == "old":
            if len(l1) != 3:
                raise CircuitError("old Bristol: line 2 must be 'n_in1 n_in2 n_out'")
            n1, n2, no = map(int, l1)
            in_groups, out_groups = [n1, n2], [no]
            gate_start = nonblank[1] + 1
        elif fmt == "fashion":
            l2 = lines[nonblank[2]]
            niv, nov = int(l1[0]), int(l2[0])
            in_groups = [int(v) for v in l1[1:]]
            out_groups = [int(v) for v in l2[1:]]
            if len(in_groups) != niv or len(out_groups) != nov:
                raise CircuitError("Bristol Fashion: bad input/output header")
            gate_start = nonblank[2] + 1
        else:
            raise CircuitError(f"unknown Bristol format {fmt!r}")
    except ValueError:
        raise CircuitError("non-integer header field") from None
    gates = _parse_bristol_gates(lines, gate_start, ngates)
    nin, nout = sum(in_groups), sum(out_groups)
    if nin + nout > nwires:
        raise CircuitError("more I/O wires than wires")
    for g in gates:
        for w in g.ins + g.outs:
            if not 0 <= w < nwires:
                raise CircuitError(f"wire {w} outside 0..{nwires - 1}")
    rc = RawCircuit(inputs=list(range(nin)), outputs=list(range(nwires - nout, nwires)),
                    gates=gates, fmt=f"bristol_{fmt}", input_groups=in_groups,
                    output_groups=out_groups, meta={"ngates": ngates, "nwires": nwires})
    return rc




def _expand_names(tokens):
    out = []
    for tok in tokens:
        m = re.fullmatch(r"([A-Za-z_]+)(\d+):([A-Za-z_]+)?(\d+)", tok)
        if m:
            p, a, p2, b = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
            if p2 not in (None, p) or b < a:
                raise CircuitError(f"bad range {tok!r}")
            out += [f"{p}{k}" for k in range(a, b + 1)]
        elif re.fullmatch(r"[A-Za-z_]+\d+", tok):
            out.append(tok)
        else:
            raise CircuitError(f"bad wire name {tok!r}")
    return out


SLP_ARITY = {"AND": 2, "XOR": 2, "XNOR": 2, "OR": 2, "NAND": 2, "NOR": 2, "NOT": 1}


def parse_slp(text: str) -> RawCircuit:
    """NIST circuit-complexity straight-line program (.slp).

    '#' comments; header lines 'Inputs: ...', 'Outputs: ...', 'Internal: ...'
    with ranges like x1:x512; 'GateSyntax: GateName Output Inputs'; gates
    between 'begin SLP' and 'end SLP', one per line: 'OP out in1 [in2]'.
    """
    inputs = outputs = internal = None
    gates, in_slp, seen_end, meta = [], False, False, {}
    for raw in text.splitlines():
        ln = raw.strip()
        if ln.startswith("#"):
            if ln.startswith("# Tally:") or ln.startswith("# Depth") or ln.startswith("# Description"):
                meta.setdefault("header_comments", []).append(ln[2:])
            continue
        if not ln:
            continue
        if in_slp:
            if ln == "end SLP":
                in_slp, seen_end = False, True
                continue
            tok = ln.split()
            op = tok[0]
            if op not in SLP_ARITY:
                raise CircuitError(f"unknown SLP gate {op!r}")
            if len(tok) != 2 + SLP_ARITY[op]:
                raise CircuitError(f"wrong fan-in: {ln!r}")
            gates.append(Gate("NOT" if op == "NOT" else op, tok[2:], [tok[1]]))
            continue
        if ln.startswith("Inputs:"):
            inputs = _expand_names(ln.split()[1:])
        elif ln.startswith("Outputs:"):
            outputs = _expand_names(ln.split()[1:])
        elif ln.startswith("Internal:"):
            internal = _expand_names(ln.split()[1:])
        elif ln.startswith("GateSyntax:"):
            if ln.split()[1:] != ["GateName", "Output", "Inputs"]:
                raise CircuitError(f"unsupported GateSyntax {ln!r}")
        elif ln == "begin SLP":
            in_slp = True
        elif ln.startswith("begin circuit"):
            meta["name"] = ln[len("begin circuit"):].strip()
        elif ln == "end circuit":
            pass
        else:
            raise CircuitError(f"unexpected SLP line {ln!r}")
    if inputs is None or outputs is None or not seen_end:
        raise CircuitError("SLP missing Inputs/Outputs or 'end SLP'")
    declared = set(inputs) | set(outputs) | set(internal or [])
    for g in gates:
        for w in g.ins + g.outs:
            if w not in declared:
                raise CircuitError(f"undeclared wire {w}")
    return RawCircuit(inputs=inputs, outputs=outputs, gates=gates, fmt="nist_slp",
                      input_groups=None, output_groups=[len(outputs)], meta=meta)


def load(path: str, fmt: str = "auto") -> RawCircuit:
    with open(path) as f:
        text = f.read()
    if fmt == "slp" or (fmt == "auto" and path.endswith(".slp")):
        return parse_slp(text)
    if fmt == "xag" or (fmt == "auto" and path.endswith(".xag")):
        raise CircuitError("use read_xag for .xag files")
    return parse_bristol(text, "auto" if fmt == "auto" else fmt)


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #
def write_bristol_fashion(x: XAG, in_groups=None, out_groups=None) -> str:
    """Bristol Fashion text.  Complemented edges become INV gates (one per
    complemented var), constants become EQ gates.  Inputs are wires
    0..n_in-1, outputs are copied (EQW / INV) into the last wires."""
    in_groups = in_groups or [x.n_inputs]
    out_groups = out_groups or [len(x.outputs)]
    assert sum(in_groups) == x.n_inputs and sum(out_groups) == len(x.outputs)
    lines = []
    wire_of_var = {v: v - 1 for v in range(1, x.n_inputs + 1)}
    nxt = [x.n_inputs]
    inv_wire, const_wire = {}, {}

    def new():
        w = nxt[0]; nxt[0] += 1
        return w

    def w_of(l):
        v, neg = l >> 1, l & 1
        if v == 0:
            if neg not in const_wire:
                const_wire[neg] = new()
                lines.append(f"1 1 {neg} {const_wire[neg]} EQ")
            return const_wire[neg]
        if not neg:
            return wire_of_var[v]
        if v not in inv_wire:
            inv_wire[v] = new()
            lines.append(f"1 1 {wire_of_var[v]} {inv_wire[v]} INV")
        return inv_wire[v]

    for k, (op, a, b) in enumerate(x.nodes):
        wa, wb = w_of(a), w_of(b)
        w = new()
        wire_of_var[x.n_inputs + 1 + k] = w
        lines.append(f"2 1 {wa} {wb} {w} {op}")
    srcs = [w_of(o) for o in x.outputs]
    for s, o in zip(srcs, x.outputs):
        lines.append(f"1 1 {s} {new()} EQW")
    head = [f"{len(lines)} {nxt[0]}",
            " ".join(map(str, [len(in_groups)] + in_groups)),
            " ".join(map(str, [len(out_groups)] + out_groups)), ""]
    return "\n".join(head + lines) + "\n"


def write_blif(x: XAG, model="xag") -> str:
    def nm(v):
        if v == 0:
            return "c0"
        return f"i{v - 1}" if v <= x.n_inputs else f"n{v}"

    out = [f".model {model}",
           ".inputs " + " ".join(f"i{i}" for i in range(x.n_inputs)),
           ".outputs " + " ".join(f"o{i}" for i in range(len(x.outputs))),
           ".names c0"]
    for k, (op, a, b) in enumerate(x.nodes):
        v = x.n_inputs + 1 + k
        pa, pb = "0" if a & 1 else "1", "0" if b & 1 else "1"
        out.append(f".names {nm(a >> 1)} {nm(b >> 1)} {nm(v)}")
        if op == AND:
            out.append(f"{pa}{pb} 1")
        else:
            qa, qb = "1" if pa == "0" else "0", "1" if pb == "0" else "1"
            out += [f"{pa}{qb} 1", f"{qa}{pb} 1"]
    for i, o in enumerate(x.outputs):
        out.append(f".names {nm(o >> 1)} o{i}")
        out.append("0 1" if o & 1 else "1 1")
    out.append(".end")
    return "\n".join(out) + "\n"
