r"""
Branching functions (coset characters) for WZW cosets  L_0 / (L_1 x ... x L_n).

Idea
----
Write every character as  q^{val} * sum_{k>=0} (Laurent polynomial)_k q^k  and, for a
fixed numerator weight, use *relative* integer orders k.  With D_a = prod_i chi_{lambda_i^a}
and  b_a = q^{val_0 - val(D_a)} * sum_{j>=0} beta_a[j] q^j,  the branching rule

        chi_{lambda_0}(q, P y) = sum_a b_a(q) D_a(q, y)

becomes, at relative order k,   N[k] = sum_a sum_{j+m=k} beta_a[j] D_a[m].

Instead of comparing whole multivariate Laurent polynomials, compare only the coefficients of
the N monomials y^{nu^b}, where nu^b = classical (Dynkin-label) part of the weights of combo b.
Because the factors live in disjoint variables, the coefficient of y^{nu^b} in D_a is a product
of single-factor weight multiplicities,

        F_ab(q) = prod_i  X^{(i)}_{lambda_i^a, nu_i^b}(q),
        X^{(i)}_{lambda, nu}(q) = [y_i^nu] chi_lambda(q, y_i),

and the equations read   n_b(q) = sum_a beta_a(q) F_ab(q).   F_ab(0) = prod_i mult of the weight
nu_i^b in the finite-dimensional irrep lambda_i^a, which is unitriangular for the dominance
order, so the system is solved order by order with exact integer arithmetic.

Everything below the Sage layer (``_coset_core``) is plain python/numpy.
"""
import math
from collections import defaultdict
from fractions import Fraction

import numpy as np

# Use fast int64 arithmetic only if a (conservative) bound shows it cannot overflow.
_INT64_LIMIT = 2 ** 62


# =============================================================================
# Sage-agnostic helpers
# =============================================================================
def _frac(x):
    """Sage Rational / Integer / int / Fraction  ->  fractions.Fraction."""
    if isinstance(x, Fraction):
        return x
    num, den = x.numerator, x.denominator
    if callable(num):                      # Sage objects expose methods, python ints attributes
        num, den = num(), den()
    return Fraction(int(num), int(den))


class _Series:
    """
    Relative-integer view of a Puiseux series   q^val * sum_k c_k q^k.

    val   : Fraction, lowest exponent present
    prec  : Fraction or None (None = exact).  Exponents < prec are correct.
    terms : dict k -> (monomials (m x rank int64 array), coefficients (m,) object array)
            (raw Sage Laurent polynomials are accepted and parsed lazily)
    """
    __slots__ = ("val", "prec", "terms", "rank")

    def __init__(self, val, prec, terms, rank):
        self.val, self.prec, self.terms, self.rank = val, prec, terms, rank

    def n_valid(self):
        """Number of correct relative orders k = 0, 1, ..., n_valid-1  (None = exact)."""
        if self.prec is None:
            return None
        return max(0, math.ceil(self.prec - self.val))

    def get(self, k):
        t = self.terms.get(k)
        if t is None:
            return None
        if not isinstance(t, tuple):       # raw Sage polynomial -> parse once
            t = _parse_poly(t, self.rank)
            self.terms[k] = t
        return t


def _parse_poly(poly, rank):
    """Sage Laurent polynomial -> (int64 array of exponent vectors, object array of int coeffs)."""
    try:
        d = poly.dict()
    except AttributeError:                 # a plain integer coefficient
        d = {(0,) * rank: poly}
    if not d:
        return np.zeros((0, rank), dtype=np.int64), np.zeros(0, dtype=object)
    rows = []
    for key in d.keys():
        if hasattr(key, "__len__"):        # ETuple / tuple (multivariate)
            rows.append(tuple(int(t) for t in key))
        else:                              # plain int (univariate Laurent ring)
            rows.append((int(key),))
    mons = np.array(rows, dtype=np.int64)
    if mons.shape[1] != rank:
        raise ValueError(
            "A character has %d Cartan variables but %d were expected "
            "(check the projection matrices / weight keys)." % (mons.shape[1], rank))
    coefs = np.empty(len(rows), dtype=object)
    for i, v in enumerate(d.values()):
        coefs[i] = int(v)
    return mons, coefs


def _components(pattern):
    """Connected components of a boolean adjacency matrix (the system decouples along them)."""
    N = pattern.shape[0]
    try:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components
    except ImportError:                    # pragma: no cover
        return np.zeros(N, dtype=int)
    return connected_components(csr_matrix(pattern), directed=False)[1]


def _layers(F0):
    """
    Topological layers of the (unit-)triangular matrix F0:  b depends on a (a != b)
    iff F0[a, b] != 0.  Layer 0 = maximal elements, etc.
    """
    g = F0.shape[0]
    B = (F0 != 0)
    np.fill_diagonal(B, False)
    indeg = B.sum(axis=0).astype(np.int64)
    alive = np.ones(g, dtype=bool)
    layers, done = [], 0
    while done < g:
        layer = np.flatnonzero(alive & (indeg == 0))
        if layer.size == 0:
            raise RuntimeError(
                "The leading-order matrix is not triangular: the weight combinations probably "
                "contain duplicates or the weight keys / projections are inconsistent.")
        layers.append(layer)
        alive[layer] = False
        indeg -= B[layer].sum(axis=0)
        done += layer.size
    return layers


def _trunc_mul(A, B, n):
    """Product of two q-series stored as arrays of shape (n, ...), elementwise in the rest."""
    out = np.zeros_like(A)
    for m in range(n):
        acc = A[0] * B[m]
        for j in range(1, m + 1):
            acc = acc + A[j] * B[m - j]
        out[m] = acc
    return out


def _solve_block(Xs, la_b, sb_b, nb_blk, n, dtype):
    """Solve  n_b(q) = sum_a beta_a(q) F_ab(q)  on one connected block.  Returns beta (n, g)."""
    F = None
    for Xi, l, s in zip(Xs, la_b, sb_b):
        G = Xi[:, l[:, None], s[None, :]]                    # (n, g, g)
        F = G if F is None else _trunc_mul(F, G, n)
    F0 = F[0]
    g = F0.shape[0]

    layers = _layers(F0)
    blocks, done = [], np.empty(0, dtype=np.intp)
    for t, L in enumerate(layers):
        blocks.append(None if t == 0 else (done, F0[np.ix_(done, L)]))
        done = np.concatenate([done, L])

    def solve_row(r):                      # beta @ F0 = r,   F0 unit-triangular
        beta = np.zeros(g, dtype=dtype)
        for L, blk in zip(layers, blocks):
            if blk is None:
                beta[L] = r[L]
            else:
                d, B = blk
                beta[L] = r[L] - np.dot(beta[d], B)
        return beta

    beta = np.zeros((n, g), dtype=dtype)
    for k in range(n):
        r = nb_blk[k].copy()
        for j in range(k):
            bj = beta[j]
            if bj.any():
                r = r - np.dot(bj, F[k - j])
        beta[k] = solve_row(r)
    return beta


def _orders(num, dens, combos, max_ord):
    """Number n of reliable relative orders and the leading exponents shift[a] of every b_a."""
    f = len(dens)
    shift = [num.val - sum((dens[i][tuple(c[i])].val for i in range(f)), Fraction(0))
             for c in combos]
    avail, seen = [], set()
    nn = num.n_valid()
    if nn is not None:
        avail.append(nn)
    for c in combos:
        for i in range(f):
            key = (i, tuple(c[i]))
            if key not in seen:
                seen.add(key)
                m = dens[i][key[1]].n_valid()
                if m is not None:
                    avail.append(m)
    n_av = min(avail) if avail else None
    if max_ord is not None:
        n = max(0, max(math.ceil(max_ord - s) for s in shift))
        if n_av is not None and n > n_av:
            raise ValueError(
                "max_ord=%s needs more q-orders than the input characters provide; "
                "the largest possible max_ord is %s." % (max_ord, min(shift) + n_av))
    else:
        if n_av is None:
            raise ValueError("All input series are exact; please specify max_ord.")
        n = n_av
    return n, shift


def _candidate_combos(num, den_keys, P_arrs, n):
    """
    All combinations of denominator weights that can possibly occur in the decomposition of
    ``num`` up to relative order n.  Necessary condition (all multiplicities are >= 0): if
    beta_a[j] != 0 then the classical weight nu^a of combo a occurs in the projected numerator
    at grade j.  So we collect the dominant projected monomials that are integrable weights.
    """
    f = len(den_keys)
    by_cl = [{tuple(int(t) for t in w[1:]): tuple(w) for w in keys} for keys in den_keys]
    ranks = [len(next(iter(d))) for d in by_cl]
    cuts = np.concatenate([[0], np.cumsum(ranks)]).astype(int)
    P_totT = np.vstack(P_arrs).T
    found = set()
    for k in range(n):
        t = num.get(k)
        if t is None or len(t[1]) == 0:
            continue
        proj = t[0] @ P_totT
        proj = proj[(proj >= 0).all(axis=1)]
        if not len(proj):
            continue
        for row in np.unique(proj, axis=0).tolist():
            parts = [tuple(row[cuts[i]:cuts[i + 1]]) for i in range(f)]
            if all(parts[i] in by_cl[i] for i in range(f)):
                found.add(tuple(by_cl[i][parts[i]] for i in range(f)))
    return found


def _conv_ints(a, b, n):
    out = [0] * n
    for i, x in enumerate(a):
        if x:
            for j in range(n - i):
                out[i + j] += x * b[j]
    return out


# =============================================================================
# The Sage-free core
# =============================================================================
def _coset_core(num, dens, combos, P_arrs, max_ord=None, check=True, den_keys=None):
    """
    num      : _Series, numerator character (Cartan variables x_1..x_r0)
    dens     : list (one entry per denominator factor) of dicts  affine-weight -> _Series
    combos   : list of tuples (wt_1, ..., wt_n) of denominator affine weights
               (affine Dynkin labels; the classical Dynkin labels are wt[1:])
    P_arrs   : list of integer arrays, P_i of shape (rank_i, r0)
    den_keys : if given (list of all integrable weights of every factor), ``combos`` is
               completed with every channel that can occur in the numerator's decomposition.
    Returns (beta, shift, n, combos): b_a = q^{shift[a]} * sum_{j<n} beta[j, a] q^j (+ O(...)),
    where ``combos`` is the (possibly completed) list the columns refer to.
    """
    f = len(dens)
    if f == 0:
        raise ValueError("At least one denominator factor is required.")
    if len(P_arrs) != f:
        raise ValueError("Need one projection matrix per denominator factor.")
    combos = [tuple(tuple(w) for w in c) for c in combos]

    if den_keys is not None:                           # add every channel that can occur
        have = set(combos)
        n, _ = _orders(num, dens, combos, max_ord)
        while True:
            new = sorted(c for c in _candidate_combos(num, den_keys, P_arrs, n) if c not in have)
            if not new:
                break
            combos += new
            have.update(new)
            n, _ = _orders(num, dens, combos, max_ord)
    N = len(combos)

    # ---- bookkeeping of weights --------------------------------------------------
    nu = [[tuple(int(t) for t in combos[a][i][1:]) for a in range(N)] for i in range(f)]
    lam_list, nu_list, la, sb = [], [], [], []
    for i in range(f):
        lam_idx, nu_idx = {}, {}
        la_i = np.empty(N, dtype=np.intp)
        sb_i = np.empty(N, dtype=np.intp)
        for a in range(N):
            la_i[a] = lam_idx.setdefault(tuple(combos[a][i]), len(lam_idx))
            sb_i[a] = nu_idx.setdefault(nu[i][a], len(nu_idx))
        lam_list.append(list(lam_idx))
        nu_list.append(list(nu_idx))
        la.append(la_i)
        sb.append(sb_i)

    pivots = [tuple(t for i in range(f) for t in nu[i][a]) for a in range(N)]
    if len(set(pivots)) != N:
        raise ValueError("weight_combos contains duplicate entries for one numerator weight.")
    pivot_index = {key: b for b, key in enumerate(pivots)}

    ranks = [len(nu_list[i][0]) for i in range(f)]
    for i in range(f):
        if P_arrs[i].shape[0] != ranks[i]:
            raise ValueError("Projection matrix %d has %d rows but factor %d has rank %d."
                             % (i, P_arrs[i].shape[0], i, ranks[i]))
    P_tot = np.vstack(P_arrs)
    P_totT = P_tot.T
    if num.rank != P_tot.shape[1]:
        raise ValueError("Projection matrices have %d columns but the numerator has rank %d."
                         % (P_tot.shape[1], num.rank))

    # ---- orders -----------------------------------------------------------------
    n, shift = _orders(num, dens, combos, max_ord)
    if n == 0:
        return np.zeros((0, N), dtype=object), shift, 0, combos

    # ---- numerator: multiplicity of every pivot weight at every order -------------
    nb = np.zeros((n, N), dtype=object)
    dimN = [0] * n
    for k in range(n):
        t = num.get(k)
        if t is None or len(t[1]) == 0:
            continue
        mons, coefs = t
        dimN[k] = int(coefs.sum())
        proj = mons @ P_totT
        mask = (proj >= 0).all(axis=1)                 # pivots are dominant weights
        if mask.any():
            for row, c in zip(proj[mask].tolist(), coefs[mask]):
                b = pivot_index.get(tuple(row))
                if b is not None:
                    nb[k, b] += c

    # ---- denominators: weight multiplicities X^{(i)}_{lambda, nu}[k] --------------
    X, dimD = [], []
    for i in range(f):
        nu_index = {w: j for j, w in enumerate(nu_list[i])}
        Xi = np.zeros((n, len(lam_list[i]), len(nu_list[i])), dtype=object)
        Di = [[0] * n for _ in lam_list[i]]
        for l, lam in enumerate(lam_list[i]):
            s = dens[i][lam]
            for k in range(n):
                t = s.get(k)
                if t is None or len(t[1]) == 0:
                    continue
                mons, coefs = t
                if mons.shape[1] != ranks[i]:
                    raise ValueError("Factor %d: character has %d variables, expected %d."
                                     % (i, mons.shape[1], ranks[i]))
                Di[l][k] = int(coefs.sum())
                mask = (mons >= 0).all(axis=1)
                if mask.any():
                    for row, c in zip(mons[mask].tolist(), coefs[mask]):
                        j = nu_index.get(tuple(row))
                        if j is not None:
                            Xi[k, l, j] += c
        X.append(Xi)
        dimD.append(Di)

    for a in range(N):                                   # sanity check of the weight convention
        for i in range(f):
            if X[i][0, la[i][a], sb[i][a]] != 1:
                raise ValueError(
                    "Factor %d, weight %s: the classical part %s does not occur with "
                    "multiplicity 1 at lowest order.  Weight keys must be affine Dynkin "
                    "labels (a_0, a_1, ..., a_r) with classical part (a_1, ..., a_r)."
                    % (i, combos[a][i], nu[i][a]))

    # ---- choose integer type -----------------------------------------------------
    Gmax = [max(1, abs(int(Xi.max())), abs(int(Xi.min()))) for Xi in X]
    nmax = max(1, abs(int(nb.max())), abs(int(nb.min())))
    Fbound = max(1, n) ** (f - 1)
    for gm in Gmax:
        Fbound *= gm
    dtype = np.int64 if Fbound * nmax * N * n < _INT64_LIMIT else object
    X = [Xi.astype(dtype) for Xi in X]
    nb = nb.astype(dtype)

    # ---- the system decouples along connected components -------------------------
    pattern = np.ones((N, N), dtype=bool)
    for i in range(f):
        nz = (X[i] != 0).any(axis=0)
        pattern &= nz[la[i][:, None], sb[i][None, :]]
    labels = _components(pattern)

    beta = np.zeros((n, N), dtype=dtype)
    order = np.argsort(labels, kind="stable")
    bounds = np.flatnonzero(np.diff(labels[order])) + 1
    for idx in np.split(order, bounds):
        beta[:, idx] = _solve_block([Xi for Xi in X],
                                    [l[idx] for l in la], [s[idx] for s in sb],
                                    nb[:, idx], n, dtype)

    # ---- consistency checks --------------------------------------------------------
    if check:
        if (beta < 0).any():
            raise RuntimeError(
                "Negative branching coefficients found: the weight combinations are incomplete "
                "or the projection matrices / characters are inconsistent.")
        tot = [0] * n
        for a in range(N):
            dD = [1] + [0] * (n - 1)
            for i in range(f):
                dD = _conv_ints(dD, dimD[i][la[i][a]], n)
            ba = [int(x) for x in beta[:, a]]
            for k, c in enumerate(_conv_ints(ba, dD, n)):
                tot[k] += c
        for k in range(n):
            if tot[k] != dimN[k]:
                raise RuntimeError(
                    "Dimension check failed at relative order %d (%d != %d): channels are "
                    "missing (a list of coset primaries modulo field identification is NOT a "
                    "complete list of channels of a numerator -- use complete=True) or P is "
                    "wrong." % (k, tot[k], dimN[k]))
    return beta, shift, n, combos


# =============================================================================
# Sage layer
# =============================================================================
def _series_from_puiseux(s, rank, Infinity):
    try:
        exps = [_frac(e) for e in s.exponents()]
        coefs = list(s.coefficients())
    except AttributeError:
        Ls, e = s.laurent_part(), int(s.ramification_index())
        exps = [Fraction(int(a), e) for a in Ls.exponents()]
        coefs = list(Ls.coefficients())
    if not exps:
        raise ValueError("A zero character was passed as input.")
    val = min(exps)
    terms = {}
    for e, c in zip(exps, coefs):
        d = e - val
        if d.denominator != 1:
            raise ValueError("Exponents of a character differ by a non-integer (%s); "
                             "these do not look like WZW characters." % d)
        terms[int(d)] = c
    prec = None
    for name in ("precision_absolute", "prec"):
        fn = getattr(s, name, None)
        if callable(fn):
            p = fn()
            if p != Infinity:
                prec = _frac(p)
            break
    return _Series(val, prec, terms, rank)


def _make_puiseux(R_q, coeffs, shift, n, max_ord):
    from math import lcm
    from sage.all import QQ, ZZ
    prec = (shift + n) if max_ord is None else _frac(max_ord)
    e = lcm(shift.denominator, prec.denominator)       # so that e*exponent and e*prec are integers
    base = int(shift * e)
    data = {base + j * e: ZZ(int(c)) for j, c in enumerate(coeffs) if c != 0}
    s = R_q(data, e=e) if data else R_q(0)
    return s.add_bigoh(QQ(prec.numerator) / QQ(prec.denominator))


class _LazyDens(dict):
    """weight -> _Series, parsed from the Sage series on first access."""
    def __init__(self, loader):
        super().__init__()
        self._loader = loader

    def __missing__(self, key):
        val = self[key] = self._loader(key)
        return val


def compute_coset_characters(num_chars, den_chars_list, weight_combos, P_list,
                             max_ord=None, check=True, complete=True):
    r"""
    Branching functions (coset characters) of a WZW coset L_0 / (L_1 x ... x L_n).

    Parameters
    ----------
    num_chars : dict
        numerator affine highest weight -> Puiseux series in q (Laurent-polynomial coefficients
        in the Cartan variables, exponents = Dynkin labels), e.g. from ``affine_char_dict``.
    den_chars_list : list of dict
        one such dictionary per denominator factor (all integrable weights of the factor).
    weight_combos : list of tuples
        (numerator weight, weight_1, ..., weight_n): the combinations whose branching functions
        are returned.  This may be just a set of representatives modulo field identification.
    P_list : list of matrices
        P_i maps Dynkin labels of L_0 to Dynkin labels of L_i (shape rank(L_i) x rank(L_0)).
    max_ord : rational, optional
        If given, every returned series is correct for all exponents < max_ord (i.e. is
        ``... + O(q^max_ord)``).  If None, as many terms as the input precision allows are returned.
    check : bool
        Verify non-negativity and the dimension identity chi_0(q, y=1) = sum_a b_a D_a(q, 1).
    complete : bool
        The decomposition of a numerator character involves *all* channels that occur, not only
        the requested ones.  If True (default) the remaining channels are found automatically
        (from the projected support of the numerator) and solved for together with the
        requested ones; only the requested combinations are returned.

    Returns
    -------
    dict : combo -> Puiseux series in q with integer coefficients, with its precision set.
        The keys appear in exactly the order of ``weight_combos``, so ``list(result.values())``
        is aligned with ``weight_combos`` (and with S and T matrices built from the same list).
    """
    from sage.all import Matrix, ZZ, PuiseuxSeriesRing, Infinity

    f = len(den_chars_list)
    if len(P_list) != f:
        raise ValueError("P_list and den_chars_list must have the same length.")
    P_arrs = [np.array([[int(x) for x in row] for row in Matrix(ZZ, P).rows()], dtype=np.int64)
              for P in P_list]
    if len({P.shape[1] for P in P_arrs}) != 1:
        raise ValueError("All projection matrices must have rank(L_0) columns.")

    by_num = defaultdict(list)
    ordered = []
    for combo in weight_combos:
        combo = tuple(tuple(w) for w in combo)
        if len(combo) != f + 1:
            raise ValueError("Each weight combination needs 1 + %d weights." % f)
        by_num[combo[0]].append(combo)
        ordered.append(combo)

    cache = {}

    def make_loader(i):
        def load(lam):
            key = (i, tuple(lam))
            if key not in cache:
                cache[key] = _series_from_puiseux(den_chars_list[i][lam], len(lam) - 1, Infinity)
            return cache[key]
        return load

    den_keys = [[tuple(w) for w in d.keys()] for d in den_chars_list] if complete else None
    R_q = PuiseuxSeriesRing(ZZ, 'q')
    result = {}
    for wt0, combos in by_num.items():
        num = _series_from_puiseux(num_chars[wt0], P_arrs[0].shape[1], Infinity)
        dens = [_LazyDens(make_loader(i)) for i in range(f)]
        beta, shift, n, used = _coset_core(num, dens, [c[1:] for c in combos], P_arrs,
                                           max_ord=max_ord, check=check, den_keys=den_keys)
        index = {c: a for a, c in enumerate(used)}
        for combo in combos:
            a = index[combo[1:]]
            result[combo] = _make_puiseux(R_q, beta[:, a], shift[a], n, max_ord)
    return {combo: result[combo] for combo in ordered}          # keys in the order of weight_combos
