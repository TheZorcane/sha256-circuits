def parse(path):
    t = open(path).read().split('\n')
    it = iter([l.split() for l in t if l.strip()])
    assert next(it) == ['xag', '1']
    n = int(next(it)[1]); m = int(next(it)[1]); nodes = []
    for j in range(m):
        k, a, b = next(it); a = int(a); b = int(b); assert k in ('A', 'X')
        assert a >> 1 < n + 1 + j and b >> 1 < n + 1 + j, "not topological"
        nodes.append((k, a, b))
    k = int(next(it)[1]); outs = [int(next(it)[0]) for _ in range(k)]
    for o in outs:
        assert o >> 1 < n + 1 + m
    return n, nodes, outs


def write(path, n, nodes, outs):
    with open(path, 'w') as f:
        f.write("xag 1\ninputs %d\nnodes %d\n" % (n, len(nodes)))
        for k, a, b in nodes:
            f.write("%s %d %d\n" % (k, a, b))
        f.write("outputs %d\n" % len(outs))
        for o in outs:
            f.write("%d\n" % o)
