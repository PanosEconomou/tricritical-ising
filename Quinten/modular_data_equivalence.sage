"""
Group modular data (S, T) by permutation equivalence.

Two pairs (S1, T1) and (S2, T2) are equivalent if there is a permutation
matrix P with

        P * S1 * P^T == S2     and     P * T1 * P^T == T2.

Equivalently, there is a bijection sigma of {0..n-1} with
        S2[sigma(i), sigma(j)] == S1[i, j]   and   T2[sigma(i), sigma(j)] == T1[i, j]
for all i, j.  (T is normally diagonal, but this code does not assume it.)

Method
------
1. Each pair is encoded as a complete directed graph on n vertices whose edge
   (i, j) is coloured by the pair (S[i,j], T[i,j]).  Entries are mapped to
   small integers through one shared dictionary, so all later comparisons are
   integer comparisons.
2. Colour refinement (1-dimensional Weisfeiler-Leman) gives every pair a
   permutation-invariant fingerprint.  Pairs are bucketed by fingerprint; pairs
   in different buckets are never equivalent.
3. Inside a bucket, each pair is tested against one representative of every
   class found so far.  The test is a backtracking search: individualise a
   vertex of the smallest non-trivial colour class, try each possible image,
   refine again, and recurse.  A candidate permutation is always verified
   against the full matrices, so a "match" is exact.  Equivalence is an
   equivalence relation, so comparing with one representative per class is
   enough.

Entries are compared exactly (they must be hashable and compare with ==),
which is the right thing for exact rings such as QQ, number fields,
CyclotomicField(N) or QQbar.  Numerical (RR/CC) entries are not suitable;
convert them to an exact ring first.

Example (Sage)
--------------
    K.<z> = CyclotomicField(16)
    s2 = z^2 + z^-2                                  # sqrt(2)
    S = matrix(K, [[1, s2, 1], [s2, 0, -s2], [1, -s2, 1]]) / 2
    T = diagonal_matrix(K, [1, z, -1])
    P = Permutation([2, 1, 3]).to_matrix()           # swaps labels 0 and 1
    pairs = [(S, T), (P*S*P.transpose(), P*T*P.transpose()), (S, diagonal_matrix(K, [1, -1, z]))]
    modular_data_equivalence_classes(pairs)          # -> [{0, 1}, {2}]
"""

from collections import Counter, defaultdict


# --------------------------------------------------------------------------
# internal helpers (pure python on integer data)
# --------------------------------------------------------------------------

def _encode(S, T, palette, ring=None):
    """Return the n x n table E with E[i][j] = id of the pair (S[i,j], T[i,j])."""
    if ring is not None:
        S = S.change_ring(ring)
        T = T.change_ring(ring)
    n = S.nrows()
    if S.ncols() != n or T.nrows() != n or T.ncols() != n:
        raise ValueError("S and T must be square matrices of the same size")
    return [[palette.setdefault((S[i, j], T[i, j]), len(palette))
             for j in range(n)] for i in range(n)]


def _refine(E, colors, table):
    """Colour refinement to a stable colouring.

    `table` is a dictionary shared by everything that is ever compared, so
    colour ids are comparable between different matrices.
    """
    n = len(E)
    ncls = len(set(colors))
    while True:
        new = []
        for i in range(n):
            Ei = E[i]
            nbrs = tuple(sorted((Ei[j], E[j][i], colors[j])
                                for j in range(n) if j != i))
            new.append(table.setdefault((colors[i], Ei[i], nbrs), len(table)))
        new_ncls = len(set(new))
        colors = new
        if new_ncls == ncls:
            return colors
        ncls = new_ncls


def _search(EA, EB, cA, cB, table, depth=0):
    """Find sigma with EB[sigma(i)][sigma(j)] == EA[i][j], or return None.

    cA, cB are colourings of A and B that any isomorphism must respect.
    """
    n = len(EA)
    cA = _refine(EA, cA, table)
    cB = _refine(EB, cB, table)
    if Counter(cA) != Counter(cB):
        return None

    classes = defaultdict(list)
    for i, c in enumerate(cA):
        classes[c].append(i)

    if len(classes) == n:                       # discrete colouring: forced map
        where = {c: j for j, c in enumerate(cB)}
        sigma = [where[cA[i]] for i in range(n)]
        for i in range(n):
            Ei, Fi = EA[i], EB[sigma[i]]
            for j in range(n):
                if Ei[j] != Fi[sigma[j]]:
                    return None
        return sigma

    # individualise a vertex of a smallest non-trivial class of A
    c = min((c for c, v in classes.items() if len(v) > 1),
            key=lambda c: (len(classes[c]), c))
    v = classes[c][0]
    mark = table.setdefault(('individualised', depth), len(table))
    for w in range(n):
        if cB[w] != c:
            continue
        a = list(cA)
        b = list(cB)
        a[v] = mark
        b[w] = mark
        sigma = _search(EA, EB, a, b, table, depth + 1)
        if sigma is not None:
            return sigma
    return None


class _Encoded(object):
    __slots__ = ("E", "colors", "fingerprint")


def _prepare(S, T, palette, table, ring=None):
    obj = _Encoded()
    obj.E = _encode(S, T, palette, ring)
    n = len(obj.E)
    obj.colors = _refine(obj.E, [0] * n, table)
    obj.fingerprint = (n, tuple(sorted(Counter(obj.colors).items())))
    return obj


# --------------------------------------------------------------------------
# public interface
# --------------------------------------------------------------------------

def modular_data_equivalence_classes(pairs, ring=None):
    """Group modular data by permutation equivalence.

    INPUT
        pairs : list of (S, T), Sage matrices over an exact ring
        ring  : optional common ring to which all matrices are first converted
                (useful if the inputs live in different but compatible rings)

    OUTPUT
        A list of sets of indices into `pairs`.  Two indices lie in the same
        set iff there is a permutation matrix P with P S P^T = S' and
        P T P^T = T'.  Sets are ordered by their smallest element.
    """
    palette, table = {}, {}
    data = [_prepare(S, T, palette, table, ring) for (S, T) in pairs]

    buckets = defaultdict(list)
    for idx, d in enumerate(data):
        buckets[d.fingerprint].append(idx)

    classes = []
    for members in buckets.values():
        reps = []                               # list of (representative index, set)
        for idx in members:
            for rep, cls in reps:
                if _search(data[rep].E, data[idx].E,
                           data[rep].colors, data[idx].colors, table) is not None:
                    cls.add(idx)
                    break
            else:
                reps.append((idx, {idx}))
        classes.extend(cls for _, cls in reps)

    classes.sort(key=min)
    return classes


def find_modular_data_permutation(S1, T1, S2, T2, ring=None):
    """Return a permutation matrix P with P*S1*P^T == S2 and P*T1*P^T == T2,
    or None if the two pairs are not equivalent."""
    from sage.all import matrix, ZZ
    palette, table = {}, {}
    A = _prepare(S1, T1, palette, table, ring)
    B = _prepare(S2, T2, palette, table, ring)
    if A.fingerprint != B.fingerprint:
        return None
    sigma = _search(A.E, B.E, A.colors, B.colors, table)
    if sigma is None:
        return None
    n = len(sigma)
    P = matrix(ZZ, n, n)
    for i in range(n):
        P[sigma[i], i] = 1
    return P
