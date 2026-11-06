r"""
Fast non-specialised characters of integrable highest-weight modules of affine Lie algebras.

Same mathematics as the original ``affine_char_dict`` (Kac-Peterson decomposition

    ch_Lambda = sum_{mu in dominant maximal weights} q^{m(mu)} S_mu(q) * (sum over the classical
                orbit of mu, modulo k*M, of the lattice theta functions Theta_v),

with S_mu the string function), but organised for speed:

* the lattice points of every theta function are enumerated with a *vectorised* Fincke-Pohst
  search (numpy) and all exact arithmetic (weights, exponents) is done in integers; the only
  floating point step is a safety-margin superset that is then filtered exactly;
* theta functions are cached per *coset* of the lattice k*M*Z^r (not per representative vector),
  and only enumerated up to the order that can actually contribute to the truncated result;
* the classical Weyl orbit modulo k*M is generated with integer arithmetic;
* the product  q^m * S_mu(q) * Theta_orbit(q)  is done directly on Laurent polynomials, once per
  (theta exponent, string index), instead of through repeated Puiseux-series additions;
* Sage objects (Laurent polynomials, the final Puiseux series) are only built at the very end.

Everything up to ``_assemble`` is plain python / numpy and can be tested without Sage.
"""
import math
import time
from collections import deque
from fractions import Fraction
from math import gcd as _igcd

import numpy as np

_BLOCK = 1 << 21          # maximal number of partial lattice points kept in memory at once
_CACHE_LIMIT = 30_000_000  # maximal number of cached theta-function points


def _lcm(*xs):
    out = 1
    for x in xs:
        out = out * x // _igcd(out, x)
    return out


def _frac(x):
    """Sage Rational / Integer / int / Fraction  ->  fractions.Fraction."""
    if isinstance(x, Fraction):
        return x
    num, den = x.numerator, x.denominator
    if callable(num):
        num, den = num(), den()
    return Fraction(int(num), int(den))


# =============================================================================
# Lattice data and the vectorised theta-function enumeration
# =============================================================================
class _Lat:
    """
    All lattice data of ``affine_char_dict`` as exact integers / floats.

    cartan, M, inv_C, inv_M, quadF : square lists of Fractions (the objects of the original code)
    qF_gcd, ram_ord                : Fractions
    """

    def __init__(self, cartan, M, inv_C, inv_M, quadF, qF_gcd, ram_ord, level):
        r = len(cartan)
        level = int(level)                  # may arrive as a Sage Integer, which does not mix with Fraction/float
        self.r, self.level = r, level
        self.cartan = self._ints(cartan, "Cartan matrix")
        self.M = self._ints(M, "M = k*C")
        ram_ord = Fraction(ram_ord)
        if ram_ord.denominator != 1:
            raise ValueError("ram_ord is not an integer: %s" % ram_ord)
        self.ram_ord = int(ram_ord)

        # key(vec) = qF_gcd * vec^T quadF vec  is an integer;  we use  s*qF_gcd*quadF  integral
        Qf = [[Fraction(qF_gcd) * Fraction(x) for x in row] for row in quadF]
        self.s = _lcm(*[x.denominator for row in Qf for x in row])
        self.Qs = np.array([[int(x * self.s) for x in row] for row in Qf], dtype=np.int64)

        # centre of the search ellipsoid  u0 = -(1/k) inv_C v   and Cholesky factor of A = M^T quadF M
        self.u0_mat = np.array([[float(-Fraction(x) / level) for x in row] for row in inv_C])
        A = self.M.T.astype(float) @ np.array([[float(x) for x in row] for row in quadF]) @ self.M.astype(float)
        self.V = np.linalg.cholesky(A).T                    # upper triangular, A = V^T V

        # canonical coset representative of  v  modulo  M Z^r   (integer version of the original)
        Dm = _lcm(*[Fraction(x).denominator for row in inv_M for x in row])
        self.Dm = Dm
        self.Madj = [[int(Fraction(x) * Dm) for x in row] for row in inv_M]
        self.cart_cols = [[self.cartan[j, i].item() for j in range(r)] for i in range(r)]

        self.cache = {}
        self.cached_points = 0
        self.t_enum = 0.0

    @staticmethod
    def _ints(Mx, what):
        out = []
        for row in Mx:
            out.append([])
            for x in row:
                x = Fraction(x)
                if x.denominator != 1:
                    raise ValueError("%s has a non-integer entry %s" % (what, x))
                out[-1].append(int(x))
        return np.array(out, dtype=np.int64)

    # ---------------------------------------------------------------- orbit
    def canon(self, v):
        return tuple(sum(a * b for a, b in zip(row, v)) % self.Dm for row in self.Madj)

    def orbit(self, v):
        """Classical Weyl orbit of the weight v (Dynkin labels) modulo M Z^r, one vector per coset."""
        r = self.r
        start = tuple(int(x) for x in v)
        orb, seen, queue = [start], {self.canon(start)}, deque([start])
        while queue:
            cur = queue.popleft()
            for i in range(r):
                ci = cur[i]
                if ci == 0:
                    continue
                col = self.cart_cols[i]
                nxt = tuple(cur[j] - ci * col[j] for j in range(r))
                c = self.canon(nxt)
                if c not in seen:
                    seen.add(c)
                    orb.append(nxt)
                    queue.append(nxt)
        return orb

    # ---------------------------------------------------------------- enumeration
    def _expand(self, u0, R2):
        """Generator of integer arrays U (m x r) containing all u with |V(u-u0)|^2 <= R2 (superset)."""
        r, V = self.r, self.V

        def rec(U, rem, k):
            if k < 0:
                yield U
                return
            m = U.shape[0]
            S = (U[:, k + 1:] - u0[k + 1:]) @ V[k, k + 1:] if k + 1 < r else np.zeros(m)
            b = np.sqrt(np.maximum(rem, 0.0)) + 1e-7
            Vkk = V[k, k]
            lo = np.ceil((-b - S) / Vkk + u0[k]).astype(np.int64)
            hi = np.floor((b - S) / Vkk + u0[k]).astype(np.int64)
            cnt = np.maximum(hi - lo + 1, 0)
            if int(cnt.sum()) == 0:
                return
            csum = np.cumsum(cnt)
            a = 0
            while a < m:                                    # split parents into memory-bounded blocks
                base = int(csum[a - 1]) if a > 0 else 0
                e = max(int(np.searchsorted(csum, base + _BLOCK, side="right")), a + 1)
                e = min(e, m)
                c = cnt[a:e]
                tot = int(c.sum())
                if tot:
                    idx = np.repeat(np.arange(a, e), c)
                    off = np.arange(tot) - np.repeat(np.cumsum(c) - c, c)
                    uk = lo[idx] + off
                    Un = U[idx]
                    Un[:, k] = uk
                    term = Vkk * (uk - u0[k]) + S[idx]
                    yield from rec(Un, rem[idx] - term * term, k - 1)
                a = e

        yield from rec(np.zeros((1, r), dtype=np.int64), np.array([R2]), r - 1)

    def _enumerate(self, v, keymax):
        """All (key, vec) with vec in v + M Z^r and key = qF_gcd vec^T quadF vec <= keymax; sorted by key."""
        r = self.r
        v = np.asarray(v, dtype=np.int64)
        R2 = 2.0 * self.level * keymax / self.ram_ord
        R2 = R2 * (1 + 1e-9) + 1e-9
        u0 = self.u0_mat @ v.astype(float)
        keys, vecs = [], []
        for U in self._expand(u0, R2):
            vec = U @ self.M.T + v
            key = ((vec @ self.Qs) * vec).sum(axis=1)
            if self.s > 1:
                if (key % self.s).any():
                    raise ValueError("qF_gcd*vec^T quadF vec is not an integer for some lattice point")
                key = key // self.s
            keep = key <= keymax                            # exact filter
            if keep.any():
                keys.append(key[keep])
                vecs.append(vec[keep].astype(np.int32))
        if not keys:
            return np.zeros(0, dtype=np.int64), np.zeros((0, r), dtype=np.int32)
        key = np.concatenate(keys)
        vec = np.concatenate(vecs)
        order = np.argsort(key, kind="stable")
        return key[order], vec[order]

    def theta(self, v, keymax):
        """Points of Theta_v with key <= keymax (cached per coset of M Z^r)."""
        ck = self.canon(v)
        ent = self.cache.get(ck)
        if ent is None or ent[0] < keymax:
            t0 = time.perf_counter()
            keys, vecs = self._enumerate(v, keymax)
            self.t_enum += time.perf_counter() - t0
            if ent is not None:
                self.cached_points -= len(ent[1])
            if self.cached_points + len(keys) > _CACHE_LIMIT:
                self.cache.clear()
                self.cached_points = 0
            ent = (keymax, keys, vecs)
            self.cache[ck] = ent
            self.cached_points += len(keys)
        _, keys, vecs = ent
        cut = int(np.searchsorted(keys, keymax, side="right"))
        return keys[:cut], vecs[:cut]


# =============================================================================
# Assembly of one character from its maximal weights
# =============================================================================
def _assemble(lat, maxws, p, make_poly):
    """
    maxws     : list of (v_dom, m, s):  dominant classical weight of a maximal weight (Dynkin labels),
                its modular characteristic m (Fraction) and the string function s (list of ints)
    p         : Fraction, the precision of the result (all terms have exponent < p)
    make_poly : callable (m x r int array of exponent vectors) -> Laurent polynomial with all
                coefficients 1; the returned objects must support  int*poly  and  poly+poly

    Returns (terms, e_tot): terms maps t -> Laurent polynomial, the exponent being t/e_tot.
    """
    p = Fraction(p)
    e_tot = _lcm(lat.ram_ord, p.denominator, *[Fraction(m).denominator for _, m, _ in maxws])
    step = e_tot // lat.ram_ord
    terms = {}
    for v_dom, m, s in maxws:
        m = Fraction(m)
        keymax = math.ceil((p - m) * lat.ram_ord) - 1        # theta exponent < p - m  (n >= 0)
        if keymax < 0:
            continue
        parts = [lat.theta(v, keymax) for v in lat.orbit(v_dom)]
        keys = np.concatenate([k for k, _ in parts])
        if not len(keys):
            continue
        vecs = np.concatenate([x for _, x in parts])
        order = np.argsort(keys, kind="stable")
        keys, vecs = keys[order], vecs[order]
        bounds = np.flatnonzero(np.diff(keys)) + 1
        starts = np.concatenate([[0], bounds])
        ends = np.concatenate([bounds, [len(keys)]])
        base = int(m * e_tot)
        for a, b in zip(starts.tolist(), ends.tolist()):
            t0 = base + int(keys[a]) * step
            if t0 >= p * e_tot:
                break
            poly = make_poly(vecs[a:b])
            for n, sn in enumerate(s):
                t = t0 + n * e_tot
                if t >= p * e_tot:
                    break
                if sn == 0:
                    continue
                term = poly if sn == 1 else sn * poly
                terms[t] = term if t not in terms else terms[t] + term
    return terms, e_tot


# =============================================================================
# Helpers to build the lattice data from a Cartan matrix without Sage (used by the tests)
# =============================================================================
def _finv(A):
    n = len(A)
    M = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[c])]
    return [row[n:] for row in M]


def _fmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def _ft(A):
    return [list(r) for r in zip(*A)]


def _qgcd(xs):
    xs = [Fraction(x) for x in xs if x != 0]
    num = 0
    for x in xs:
        num = _igcd(num, x.numerator)
    return Fraction(num, _lcm(*[x.denominator for x in xs]))


def _lattice_from_cartan(cartan, d, level):
    """The quantities of the original code (Fractions), computed from the Cartan matrix and symmetrizer d."""
    A = [[Fraction(x) for x in row] for row in cartan]
    r = len(A)
    dm = max(d)
    d = [Fraction(x) / dm for x in d]
    inv = _finv(A)
    G = [[d[i] * A[i][j] for j in range(r)] for i in range(r)]
    quadF = _fmul(_fmul(_ft(inv), G), inv)
    C = [[A[i][j] / d[j] for j in range(r)] for i in range(r)]
    M = [[level * x for x in row] for row in C]
    qF_gcd = 1 / _qgcd([x for row in quadF for x in row])
    return _Lat(A, M, _finv(C), _finv(M), quadF, qF_gcd, 2 * level * qF_gcd, level)


def _comarks(cartan_aff, zero_pos=0):
    """
    Comarks a_i^vee of an affine Cartan matrix  A_ij = <alpha_i^vee, alpha_j>  (Kac convention, the
    one used throughout the original code): the left null vector  x A = 0, scaled so that the
    comark of node 0 is 1.  Exact (Fractions); returns a list of positive python ints.
    """
    n = len(cartan_aff)
    B = [[Fraction(cartan_aff[j][i]) for j in range(n)] for i in range(n)]       # B = A^T, solve B x = 0
    pivots, row = [], 0
    for col in range(n):
        piv = next((i for i in range(row, n) if B[i][col] != 0), None)
        if piv is None:
            continue
        B[row], B[piv] = B[piv], B[row]
        pv = B[row][col]
        B[row] = [x / pv for x in B[row]]
        for i in range(n):
            if i != row and B[i][col] != 0:
                f = B[i][col]
                B[i] = [x - f * y for x, y in zip(B[i], B[row])]
        pivots.append(col)
        row += 1
    free = [c for c in range(n) if c not in pivots]
    if len(free) != 1:
        raise ValueError("the affine Cartan matrix does not have a one-dimensional null space")
    f = free[0]
    x = [Fraction(0)] * n
    x[f] = Fraction(1)
    for r, c in enumerate(pivots):
        x[c] = -B[r][f]
    x = [v / x[zero_pos] for v in x]
    if any(v.denominator != 1 or v <= 0 for v in x):
        raise ValueError("comarks are not positive integers: %s" % x)
    return [int(v) for v in x]


# =============================================================================
# Sage layer
# =============================================================================
def affine_char_dict(cartan_type, level, max_ord=15, verbose=False):
    r'''
    Input:
      - cartan_type: Cartan type string (e.g., "A1", "A4") or CartanType object
      - level: positive integer (level k)
      - max_ord: maximum expansion order for q
      - verbose: print a timing breakdown
    Output:
      - Dictionary mapping highest weight vectors to their non-specialised characters
    '''
    from sage.all import (CartanType, RootSystem, QQ, ZZ, vector, diagonal_matrix, gcd,
                          LaurentPolynomialRing, PuiseuxSeriesRing)
    from sage.combinat.root_system.integrable_representations import IntegrableRepresentation

    t_start = time.perf_counter()
    level = int(level)
    ct = CartanType(cartan_type)
    if ct.is_affine():
        ct = ct.classical()

    cartan_mat = ct.cartan_matrix()
    inv_cartan_mat = cartan_mat.inverse()

    d = vector(QQ, ct.symmetrizer())
    d /= max(d)

    G = diagonal_matrix(QQ, d) * cartan_mat
    quadF = inv_cartan_mat.transpose() * G * inv_cartan_mat

    ct_aff = ct.affine()
    index = list(ct_aff.index_set())
    WL = RootSystem(ct_aff).weight_lattice(extended=True)
    Lambda = WL.fundamental_weights()
    rnk = len(index) - 1

    # Simple coroots in fundamental weight coordinates
    C = cartan_mat * diagonal_matrix(QQ, [1 / d[i] for i in range(rnk)])
    M = level * C
    inv_M = M.inverse()
    inv_C = C.inverse()

    qF_gcd = 1 / gcd(quadF.list())
    ram_ord = 2 * level * qF_gcd

    mat = lambda X: [[_frac(x) for x in row] for row in X.rows()]
    lat = _Lat(mat(cartan_mat), mat(M), mat(inv_C), mat(inv_M), mat(quadF),
               _frac(qF_gcd), _frac(ram_ord), int(level))

    var_names = ['x%d' % (i + 1) for i in range(rnk)]
    Rx = LaurentPolynomialRing(ZZ, var_names)
    P = PuiseuxSeriesRing(Rx, 'q')
    q = P.gen()

    def make_poly(rows):
        if rnk == 1:
            return Rx(dict.fromkeys(rows[:, 0].tolist(), 1))
        return Rx(dict.fromkeys(map(tuple, rows.tolist()), 1))

    # Integrable weights setup.  The comarks (level of Lambda = sum_i comark_i * a_i) are computed
    # exactly from the affine Cartan matrix instead of from the dual root lattice's null root,
    # which does not give the integer comarks for some non-simply-laced types (e.g. G2).
    cart_aff = [[_frac(x) for x in row] for row in ct_aff.cartan_matrix().rows()]
    comarks = dict(zip(index, _comarks(cart_aff, index.index(0))))

    int_weights_coeffs = []

    def find_weights(idx, current_sum, current_dict):
        if idx == rnk + 1:
            if current_sum == level:
                int_weights_coeffs.append(current_dict.copy())
            return
        i = index[idx]
        max_val = (level - current_sum) // comarks[i]
        for val in range(max_val + 1):
            current_dict[i] = val
            find_weights(idx + 1, current_sum + val * comarks[i], current_dict)
            del current_dict[i]

    find_weights(0, 0, {})
    int_weights = [sum(Lambda[i] * int_weights_coeffs[j][i] for i in index)
                   for j in range(len(int_weights_coeffs))]

    t_setup = time.perf_counter() - t_start
    t_sage = t_asm = t_final = 0.0

    # Character Assembly
    char_dict = {}
    for iwght in int_weights:
        t0 = time.perf_counter()
        int_rep = IntegrableRepresentation(iwght)
        md_weights = int_rep.dominant_maximal_weights()
        maxws = []
        for wght in md_weights:
            v_dom = tuple(int(x) for x in wght.to_classical().to_vector())
            mod_char = int_rep.modular_characteristic(wght)
            str_func = [int(c) for c in int_rep.string(wght, depth=max_ord)]
            maxws.append((v_dom, mod_char, str_func))

        # precision of the result, obtained from Sage's own arithmetic on empty series
        # (identical to what the term-by-term construction of the original code yields)
        probe = P(0).add_bigoh(max_ord)
        for _, mod_char, str_func in maxws:
            probe += (q ** mod_char) * P(str_func) * P(0).add_bigoh(max_ord)
        prec = getattr(probe, "precision_absolute", None) or probe.prec
        p = _frac(prec())
        t_sage += time.perf_counter() - t0

        t0 = time.perf_counter()
        terms, e_tot = _assemble(lat, [(v, _frac(m), s) for v, m, s in maxws], p, make_poly)
        t_asm += time.perf_counter() - t0

        t0 = time.perf_counter()
        char = P({int(t): poly for t, poly in terms.items()}, e=ZZ(e_tot))
        char = char.add_bigoh(QQ(p.numerator) / QQ(p.denominator))
        char_dict[tuple(iwght.to_vector()[:-1])] = char
        t_final += time.perf_counter() - t0

    if verbose:
        print("comarks:", comarks)
        print("affine_char_dict(%s, %s, max_ord=%s): total %.2fs = setup %.2fs + Sage "
              "IntegrableRepresentation %.2fs + lattice enumeration %.2fs + assembly %.2fs + "
              "Sage series construction %.2fs"
              % (cartan_type, level, max_ord, time.perf_counter() - t_start, t_setup, t_sage,
                 lat.t_enum, t_asm - lat.t_enum, t_final))
    return char_dict
