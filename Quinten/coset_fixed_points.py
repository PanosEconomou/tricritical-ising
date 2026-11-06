# =============================================================================
#  Coset primaries + modular S and T matrices, INCLUDING fixed-point resolution
#
#  Replaces  find_primaries  and  get_modular_matrices.
#
#  Main entry point:
#       primaries, S, T = get_coset_modular_data(coset)
#
#  Reference: Fuchs, Schellekens, Schweigert, "A matrix S for all simple current
#  extensions", hep-th/9601078  (formula (5.1), treated as a fact), together with
#  the orbit-Lie-algebra table of section 6.3 and the Kac-Peterson formula of 6.4.
#
#  Needs from the notebook:  get_prop_of_lie_algebra, find_weights   (unchanged)
#  Denominators with U(1) factors are NOT supported (a ValueError is raised).
# =============================================================================
from fractions import Fraction
from math import gcd, isqrt
from sage.rings.universal_cyclotomic_field import UniversalCyclotomicField
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
import itertools as itr
import numpy as np

#from coset_helper_functions import get_prop_of_lie_algebra,find_weights

UCF = UniversalCyclotomicField()

CHECK_CONSISTENCY = True   # cheap internal sanity checks (set to False to skip)

# NB: this file is written so that it behaves identically whether it is run as plain python
# (load('coset_fixed_points.py')) or pasted into a preparsed Sage notebook cell: rationals are python
# Fractions built only from python ints, and there are no float / complex literals.
_TOL = float('1e-9')
_TOL8 = float('1e-8')
_TOL6 = float('1e-6')
_HALF = float('0.5')
_TWO_PI = float(np.pi) + float(np.pi)


# -----------------------------------------------------------------------------
#  Small helpers
# -----------------------------------------------------------------------------
def _c(x):
    """Numerical value (python complex) of an exact or numerical number."""
    try:
        return complex(x.n())
    except AttributeError:
        return complex(x)


def _is_zero(x):
    return abs(_c(x)) < _TOL


def _frac(x):
    """Sage rational / integer / python Fraction  ->  python Fraction."""
    return Fraction(str(x))


def _Q(a, b=1):
    """Fraction from (possibly Sage) integers."""
    return Fraction(int(a), int(b))


def _is_int(q):
    return q.denominator == 1


def _det(rows):
    """Determinant of a square matrix (list of rows) over UCF."""
    n = len(rows)
    if n == 1:
        return rows[0][0]
    if n == 2:
        return rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]
    try:
        return matrix(UCF, rows).determinant()
    except Exception:          # fall back to Laplace expansion
        tot = UCF(0)
        for j in range(n):
            minor = [r[:j] + r[j + 1:] for r in rows[1:]]
            tot = tot + (rows[0][j] if j % 2 == 0 else -rows[0][j]) * _det(minor)
        return tot


def _mod1(q):
    """Fractional part of a Fraction, in [0,1)."""
    return q - (q.numerator // q.denominator)


_E_cache = {}


def _E(q):
    """exp(2 pi i q) in the universal cyclotomic field, q a Fraction."""
    q = _mod1(Fraction(q))          # q is already a Fraction
    if q not in _E_cache:
        _E_cache[q] = UCF.gen(int(q.denominator)) ** int(q.numerator)
    return _E_cache[q]


_sqrt_cache = {}


def _ucf_sqrt(N):
    """Exact sqrt(N) in UCF for a positive integer N, via the quadratic Gauss sum
    G(4N) = sum_k exp(2 pi i k^2/(4N)) = 2(1+i) sqrt(N)."""
    N = int(N)
    if N not in _sqrt_cache:
        g = UCF.gen(4 * N)
        G = sum((g ** ((k * k) % (4 * N)) for k in range(4 * N)), UCF(0))
        _sqrt_cache[N] = G / (2 * (1 + UCF.gen(4)))
    return _sqrt_cache[N]


def _perm_orbits(perm):
    """Orbits of a permutation (tuple i -> perm[i]), each sorted, ordered by smallest element."""
    seen, orbits = set(), []
    for i in range(len(perm)):
        if i in seen:
            continue
        orb, j = [i], perm[i]
        while j != i:
            orb.append(j)
            j = perm[j]
        seen.update(orb)
        orbits.append(tuple(sorted(orb)))
    return orbits


def _perm_order(perm):
    order = 1
    for orb in _perm_orbits(perm):
        order = order * len(orb) // gcd(order, len(orb))
    return order


# -----------------------------------------------------------------------------
#  Twisted orbit algebra  Btilde_n^(2) = A_{2n}^(2): Kac-Peterson S-matrix
#
#  Horizontal subalgebra C_n, level p, h^vee = 2n+1, lattice M = root lattice of B_n.
#  In an orthonormal basis (long root of C_n has norm 4) the Weyl sum over the
#  hyperoctahedral group is a determinant:
#       S = i^{n(n-2)} N^{-n/2} det_{jk}[ E_N^{x_j y_k} - E_N^{-x_j y_k} ],   N = p+2n+1
#  with x_j = sum_{i>=j} (c_i+1) for horizontal weight c=(c_1..c_n), E_N=exp(2 pi i/N).
#  For ODD p this equals the S-matrix of C_n^(1) at level (p-1)/2 (Ginory, Prop. 4.1).
#  For EVEN p it is a genuinely different matrix (NOT C_n at level p/2 with signs).
# -----------------------------------------------------------------------------
def _twisted_S_entry(n, p, ca, cb):
    N = p + 2 * n + 1
    xa = [sum(ca[i] + 1 for i in range(j, n)) for j in range(n)]
    xb = [sum(cb[i] + 1 for i in range(j, n)) for j in range(n)]
    g = UCF.gen(N)
    rows = [[g ** ((xa[j] * xb[k]) % N) - g ** ((-xa[j] * xb[k]) % N) for k in range(n)]
            for j in range(n)]
    det = _det(rows)
    phase = UCF.gen(4) ** ((n * (n - 2)) % 4)
    return phase * det / (_ucf_sqrt(N) ** n)


# -----------------------------------------------------------------------------
#  Orbit Lie algebras (table on p.39 of the paper)
#
#  Given a simple factor (letter, rank n) and a non-trivial diagram automorphism
#  `perm` (permutation of the affine nodes 0..n) this returns
#     kind   : 'untwisted' | 'twisted' | 'trivial'(= A_0^(1))
#     name   : Sage name of the untwisted orbit algebra (e.g. 'C3')
#     n_tw   : n of Btilde_n^(2)  (twisted case)
#     reps   : representative node of each orbit, in the order of the nodes of the
#              orbit algebra  (labels of the fixed point = labels at these nodes)
#     ratio  : level of the factor / level p of the orbit algebra
#     dpp    : shift delta per unit of p   (delta = dpp * p)
#  Node orderings assume Sage's conventions for the affine diagrams.
# -----------------------------------------------------------------------------
def _orbit_algebra(letter, n, perm, adjacency=None):
    orbits = _perm_orbits(perm)
    r = int(_perm_order(perm))

    def fail():
        raise NotImplementedError("Unrecognised diagram automorphism %s of %s%d (orbits %s)"
                                  % (perm, letter, n, orbits))

    def c_name(m):                     # C_m^(1); C_1 = A_1
        return 'A1' if m == 1 else 'C%d' % m

    if letter == 'A':
        ell = len(orbits)
        reps = [o[0] for o in orbits]
        if reps != list(range(ell)) or any(len(o) != r for o in orbits):
            fail()
        return dict(kind='untwisted' if ell > 1 else 'trivial', name='A%d' % (ell - 1),
                    reps=reps, ratio=r, dpp=_Q((n + 1) ** 2 - ell ** 2, 24 * ell))

    if letter == 'B':                  # B_n^(1) -> Btilde_{n-1}^(2),  p = level
        if orbits != [(0, 1)] + [(i,) for i in range(2, n + 1)]:
            fail()
        return dict(kind='twisted', n_tw=n - 1, reps=[0] + list(range(2, n + 1)),
                    ratio=int(1), dpp=_Q(2 * n + 1, 24))

    if letter == 'C':
        if orbits != [tuple(sorted({i, n - i})) for i in range(n // 2 + 1)]:
            fail()
        if n % 2 == 0:                 # C_{2m}^(1) -> Btilde_m^(2),  p = level
            return dict(kind='twisted', n_tw=n // 2, reps=list(range(n // 2 + 1)),
                        ratio=int(1), dpp=_Q(5 * (n // 2), 24))
        return dict(kind='untwisted', name=c_name(n // 2),   # C_{2m+1}^(1) -> C_m^(1), level 2p
                    reps=list(range(n // 2 + 1)), ratio=int(2), dpp=_Q(n, 8))

    if letter == 'D':
        vec = sorted([(0, 1), (n - 1, n)] + [(i,) for i in range(2, n - 1)])
        if sorted(orbits) == vec:      # vector current -> C_{n-2}^(1), level 2p
            return dict(kind='untwisted', name=c_name(n - 2),
                        reps=[0] + list(range(2, n - 1)) + [n - 1], ratio=int(2), dpp=_Q(1, 4))
        if n % 2 == 0 and r == 2:      # spinor currents of D_{2m} -> B_m^(1), level 2p
            m = n // 2
            if [o[0] for o in orbits] != list(range(m + 1)) or \
               any(len(o) != (1 if i == m else 2) for i, o in enumerate(orbits)):
                fail()
            return dict(kind='untwisted', name='B%d' % m, reps=list(range(m + 1)),
                        ratio=int(2), dpp=_Q(m, 8))
        if n % 2 == 1 and r == 4:      # spinor currents of D_{2m+1} -> C_{m-1}^(1), level 4p
            m = (n - 1) // 2
            if orbits != [(0, 1, n - 1, n)] + [tuple(sorted({i, n - i})) for i in range(2, m + 1)]:
                fail()
            return dict(kind='untwisted', name=c_name(m - 1), reps=[0] + list(range(2, m + 1)),
                        ratio=int(4), dpp=_Q(n + 2, 8))
        fail()

    if letter == 'E' and n in (6, 7):  # E6 -> G2 (level 3p), E7 -> F4 (level 2p)
        if adjacency is None:
            fail()
        orb_of = {i: k for k, o in enumerate(orbits) for i in o}
        order = [orb_of[0]]
        while len(order) < len(orbits):
            nxt = {orb_of[j] for i in orbits[order[-1]] for j in adjacency[i]} - set(order)
            if len(nxt) != 1:
                fail()
            order.append(nxt.pop())
        return dict(kind='untwisted', name='G2' if n == 6 else 'F4',
                    reps=[orbits[k][0] for k in order], ratio=r,
                    dpp=_Q(2, 3) if n == 6 else _Q(3, 8))
    fail()


# -----------------------------------------------------------------------------
#  One simple factor of the coset (numerator, or one denominator factor)
# -----------------------------------------------------------------------------
class _CosetFactor(object):

    def __init__(self, name, level, conj):
        self.name = name
        self.letter = name[0]
        self.rank = int(name[1:])
        self.level = int(level)
        self.conj = bool(conj)                       # denominators enter complex conjugated
        self.sign = int(-1) if conj else int(1)
        dim, dcox, cmks, theta, F, aut, v, mod = get_prop_of_lie_algebra(self.letter, self.rank, detail=True)
        n = self.rank
        self.dim, self.dcox = int(dim), int(dcox)
        self.cmks = tuple(int(x) for x in cmks)
        self.F = [[_frac(F[i, j]) for j in range(n)] for i in range(n)]
        self.aut = [tuple(int(x) for x in p) for p in aut]
        self.ident = tuple(range(n + 1))
        self.vac = (self.level,) + (0,) * n
        self.primaries = [tuple(int(x) for x in w) for w in find_weights(self.cmks, self.level)]
        self._rhoFrho = sum(sum(row) for row in self.F)
        self._fusion_data = None
        self._S1 = {}
        self._SJ = {}
        self._fixed = {}
        self._orbit = {}
        self._hcur = {}
        self._adj = None

    # ---- elementary data ---------------------------------------------------
    def act(self, perm, lab):
        return tuple(lab[perm[i]] for i in range(len(lab)))

    def h(self, lab):
        """conformal weight of the primary with affine Dynkin labels lab"""
        x = [_Q(l + 1) for l in lab[1:]]
        n = self.rank
        q = sum(x[i] * self.F[i][j] * x[j] for i in range(n) for j in range(n))
        return (q - self._rhoFrho) / _Q(2 * (self.level + self.dcox))

    def central_charge(self):
        return _Q(self.level * self.dim, self.level + self.dcox)

    def h_current(self, perm):
        """conformal weight of the simple current  perm.(vacuum)"""
        if perm not in self._hcur:
            self._hcur[perm] = self.h(self.act(perm, self.vac))
        return self._hcur[perm]

    def fixed(self, perm):
        if perm not in self._fixed:
            self._fixed[perm] = [p for p in self.primaries if self.act(perm, p) == p]
        return self._fixed[perm]

    # ---- ordinary (unitary) S-matrix of the factor --------------------------
    def _fusion(self):
        if self._fusion_data is None:
            fr = FusionRing(self.name, self.level)
            qd = fr.global_q_dimension().sqrt()
            if _c(qd).real < 0:
                qd = -qd
            self._fusion_data = (fr, qd)
        return self._fusion_data

    def S1(self, a, b):
        key = (a, b)
        if key not in self._S1:
            fr, qd = self._fusion()
            self._S1[key] = UCF(fr.s_ij(fr(a[1:]), fr(b[1:])) / qd)
        return self._S1[key]

    # ---- orbit algebra data and S^J -----------------------------------------
    def _adjacency(self):
        if self._adj is None:
            A = CartanType([self.letter, self.rank, 1]).cartan_matrix()
            n = self.rank
            self._adj = [[j for j in range(n + 1) if j != i and A[i, j] != 0] for i in range(n + 1)]
        return self._adj

    def orbit_data(self, perm):
        orbs = tuple(_perm_orbits(perm))
        if orbs not in self._orbit:
            od = _orbit_algebra(self.letter, self.rank, perm,
                                self._adjacency() if self.letter == 'E' else None)
            od['p'] = int(self.level // od['ratio']) if self.level % od['ratio'] == 0 else None
            od['delta'] = od['dpp'] * od['p'] if od['p'] is not None else None
            od['orbits'] = orbs
            od['fusion'] = None
            self._orbit[orbs] = od
        return self._orbit[orbs]

    def orbit_weight(self, od, lab):
        """Weight of the orbit algebra corresponding to the fixed point lab."""
        b = [lab[r] for r in od['reps']]
        if od['kind'] == 'twisted':                  # horizontal C_n weight (c_1..c_n) = (b_n..b_1)
            return tuple(reversed(b[:-1]))
        return tuple(b[1:])

    def _check_level(self, od, lab):
        b = [lab[r] for r in od['reps']]
        if od['kind'] == 'twisted':
            p = 2 * sum(b[:-1]) + b[-1]
        elif od['kind'] == 'trivial':
            p = b[0]
        else:
            lt, rk = od['name'][0], int(od['name'][1:])
            cm = [int(x) for x in get_prop_of_lie_algebra(lt, rk, detail=True)[2]]
            p = sum(c * x for c, x in zip(cm, b))
        if p != od['p']:
            raise RuntimeError("Level mismatch for fixed point %s of %s: %s != %s" % (lab, self.name, p, od['p']))

    def SJ(self, perm, a, b):
        """S^J_{a,b} of this factor (not conjugated), a and b fixed points of perm."""
        if perm == self.ident:
            return self.S1(a, b)
        od = self.orbit_data(perm)
        key = (od['orbits'], a, b)
        if key not in self._SJ:
            if self.act(perm, a) != a or self.act(perm, b) != b:
                return UCF(0)
            if CHECK_CONSISTENCY:
                self._check_level(od, a)
                self._check_level(od, b)
            wa, wb = self.orbit_weight(od, a), self.orbit_weight(od, b)
            if od['kind'] == 'trivial':
                val = UCF(1)
            elif od['kind'] == 'twisted':
                val = _twisted_S_entry(od['n_tw'], od['p'], wa, wb)
            else:
                if od['fusion'] is None:
                    fr = FusionRing(od['name'], od['p'])
                    qd = fr.global_q_dimension().sqrt()
                    if _c(qd).real < 0:
                        qd = -qd
                    od['fusion'] = (fr, qd)
                fr, qd = od['fusion']
                val = UCF(fr.s_ij(fr(wa), fr(wb)) / qd)
            self._SJ[key] = _E(-_Q(3) * od['delta']) * val          # exp(-6 pi i delta) * S-check
        return self._SJ[key]

    def check_SJ(self, perm):
        """(S^J T^J)^3 = (S^J)^2 and unitarity on the fixed points (condition {2},{3})."""
        fx = self.fixed(perm)
        if len(fx) == 0 or len(fx) > 80:
            return
        S = np.array([[_c(self.SJ(perm, a, b)) for b in fx] for a in fx])
        c24 = self.central_charge() / _Q(24)
        T = np.diag([np.exp(complex(0, _TWO_PI * float(self.h(a) - c24))) for a in fx])
        ST = S @ T
        if np.abs(S @ S.conj().T - np.eye(len(fx))).max() > _TOL8 or \
           np.abs(ST @ ST @ ST - S @ S).max() > _TOL8:
            raise RuntimeError("S^J for factor %s (level %d), automorphism %s fails unitarity or "
                               "(S^J T^J)^3 = (S^J)^2.  (If this factor is of type A, check that "
                               "FusionRing.s_ij uses the same convention as T = exp(2 pi i (h - c/24)).)"
                               % (self.name, self.level, perm))

    # ---- phase F(a,K,J) of this factor, as an exponent mod 1 ----------------
    def F_exp(self, a, K, J):
        """F(a,K,J) = exp(2 pi i * result);  a is fixed by K and J.
        Uses  F(a,K,J) = exp(-2 pi i Q_K(b))  for any b with S^J_{a,b} != 0."""
        if K == self.ident or J == self.ident:
            return _Q(0)
        qs = set()
        for b in self.fixed(J):
            if not _is_zero(self.SJ(J, a, b)):
                qs.add(_mod1(self.h(b) + self.h_current(K) - self.h(self.act(K, b))))
                if not CHECK_CONSISTENCY:
                    break
        if len(qs) != 1:
            raise RuntimeError("F(a,K,J) is not well defined for factor %s, a=%s: monodromy charges %s"
                               % (self.name, a, sorted(qs)))
        return _mod1(-qs.pop())


# -----------------------------------------------------------------------------
#  Finite abelian group helpers (group elements = tuples of permutations)
# -----------------------------------------------------------------------------
def _gmul(g, h):
    return tuple(tuple(q[p[i]] for i in range(len(p))) for p, q in zip(g, h))


def _closure(gens, ident):
    group, frontier = {ident}, [ident]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = _gmul(x, g)
                if y not in group:
                    group.add(y)
                    new.append(y)
        frontier = new
    return group


def _characters(elems, ident):
    """All characters of a finite abelian group: list of dicts  element -> Fraction in [0,1).
    The first one is the trivial character."""
    elems = list(elems)
    if len(elems) == 1:
        return [{ident: _Q(0)}]
    gens, sub = [], {ident}
    for g in elems:
        if g not in sub:
            gens.append(g)
            sub = _closure(gens, ident)
    expo = 1
    for g in gens:
        o, x = 1, g
        while x != ident:
            x = _gmul(x, g)
            o += 1
        expo = expo * o // gcd(expo, o)
    chars = []
    for xs in itr.product(range(expo), repeat=len(gens)):
        val, queue, ok = {ident: _Q(0)}, [ident], True
        while queue and ok:
            e = queue.pop()
            for g, x in zip(gens, xs):
                new, v = _gmul(e, g), _mod1(val[e] + _Q(x, expo))
                if new in val:
                    if val[new] != v:
                        ok = False
                        break
                else:
                    val[new] = v
                    queue.append(new)
        if ok:
            chars.append(val)
    if len(chars) != len(elems):
        raise RuntimeError("character enumeration failed (group not abelian?)")
    return chars


def _isqrt_exact(n):
    n = int(n)
    r = isqrt(n)
    if r * r != n:
        raise RuntimeError("|S_a|/|U_a| = %d is not a perfect square" % n)
    return r


# -----------------------------------------------------------------------------
#  Core: resolve fixed points, formula (5.1) of the paper.
#    factors : list of _CosetFactor (numerator first)
#    group   : list of group elements (identity included) -- the extension group G
#    fields  : all fields of the unextended theory that are neutral under G
# -----------------------------------------------------------------------------
def _resolve_fixed_points(factors, group, fields, check=True):
    nf = len(factors)
    ident = tuple(f.ident for f in factors)
    group = list(group)
    fields = list(fields)
    field_set = set(fields)

    def gact(g, a):
        return tuple(factors[i].act(g[i], a[i]) for i in range(nf))

    # ---- sanity checks on the extension group -------------------------------
    if check:
        for g in group:
            spin = sum(factors[i].sign * factors[i].h_current(g[i]) for i in range(nf))
            if not _is_int(spin):
                raise ValueError("Identification current %s has non-integer spin %s" % (g, spin))
        for g in group:
            for a in fields:
                Q = sum(factors[i].sign * (factors[i].h(a[i]) + factors[i].h_current(g[i])
                                           - factors[i].h(factors[i].act(g[i], a[i]))) for i in range(nf))
                if not _is_int(Q):
                    raise ValueError("Field %s has non-zero monodromy charge %s under %s: "
                                     "selection rules and identification group are inconsistent." % (a, Q, g))

    # ---- orbits and stabilisers (representative = first member in `fields`) --
    seen, reps, Sa = set(), [], {}
    for a in fields:
        if a in seen:
            continue
        orbit = {gact(g, a) for g in group}
        if not orbit <= field_set:
            raise ValueError("Orbit of %s leaves the list of selected fields" % (a,))
        seen |= orbit
        reps.append(a)
        Sa[a] = [g for g in group if gact(g, a) == a]

    # ---- F(a,K,J), untwisted stabiliser U_a, characters ----------------------
    def F_exp(a, K, J):
        tot = _Q(0)
        for i in range(nf):
            tot += factors[i].sign * factors[i].F_exp(a[i], K[i], J[i])
        return _mod1(tot)

    Ua, mult, chars = {}, {}, {}
    for a in reps:
        if len(Sa[a]) == 1:
            Ua[a] = Sa[a]
        else:
            Ua[a] = [J for J in Sa[a] if all(F_exp(a, K, J).numerator == 0 for K in Sa[a])]
        mult[a] = _isqrt_exact(len(Sa[a]) // len(Ua[a])) if len(Sa[a]) % len(Ua[a]) == 0 else None
        if mult[a] is None:
            raise RuntimeError("U_a does not divide S_a")
        chars[a] = _characters(Ua[a], ident)
        if check:
            for J in Sa[a]:                                     # F(a,J,K) = F(a,K,J)^*   (A.12)
                for K in Sa[a]:
                    if _mod1(F_exp(a, J, K) + F_exp(a, K, J)).numerator != 0:
                        raise RuntimeError("F(a,J,K) F(a,K,J) != 1 for a=%s" % (a,))

    # ---- check the factor matrices S^J for every non-trivial J that occurs ---
    if check:
        done = set()
        for a in reps:
            for J in Sa[a]:
                for i in range(nf):
                    if J[i] != factors[i].ident and (i, J[i]) not in done:
                        factors[i].check_SJ(J[i])
                        done.add((i, J[i]))

    # ---- primaries (fixed points are repeated |U_a| times) and T -------------
    c_tot = sum(f.sign * f.central_charge() for f in factors)
    primaries, offsets = [], {}
    for a in reps:
        offsets[a] = len(primaries)
        primaries += [a] * len(chars[a])
    N = len(primaries)

    T = matrix(UCF, N, N)
    for a in reps:
        expo = sum(factors[i].sign * factors[i].h(a[i]) for i in range(nf)) - c_tot / _Q(24)
        for k in range(len(chars[a])):
            T[offsets[a] + k, offsets[a] + k] = _E(expo)

    # ---- S, formula (5.1) -----------------------------------------------------
    def SJ_total(J, a, b):
        val = UCF(1)
        for i in range(nf):
            v = factors[i].SJ(J[i], a[i], b[i])
            if _is_zero(v):
                return UCF(0)
            val = val * (v.conjugate() if factors[i].conj else v)
        return val

    S = matrix(UCF, N, N)
    order_G = len(group)
    Uset = {a: set(Ua[a]) for a in reps}
    Sset = {a: set(Sa[a]) for a in reps}
    for a in reps:
        for b in reps:
            common = [J for J in Ua[a] if J in Uset[b]]
            if check:                                           # S^J_{ab} = 0 unless J in U_a and U_b   (5.6)
                for J in Sa[a]:
                    if J in Sset[b] and not (J in Uset[b] and J in Uset[a]):
                        if not _is_zero(SJ_total(J, a, b)):
                            raise RuntimeError("S^J_{a,b} != 0 although J not in U_a and U_b (a=%s, b=%s)" % (a, b))
            vals = [(J, SJ_total(J, a, b)) for J in common]
            pref = _Q(order_G, len(Ua[a]) * mult[a] * len(Ua[b]) * mult[b])
            prefU = UCF(int(pref.numerator)) / int(pref.denominator)
            for i, ca in enumerate(chars[a]):
                for j, cb in enumerate(chars[b]):
                    s = UCF(0)
                    for J, v in vals:
                        s = s + _E(ca[J] - cb[J]) * v
                    S[offsets[a] + i, offsets[b] + j] = prefU * s
    return primaries, S, T


# -----------------------------------------------------------------------------
#  Selection rules and joint centre: this is the code of find_primaries, unchanged,
#  except that the field identification step is dropped (it is now done above).
# -----------------------------------------------------------------------------
def _selection_rules(coset):
    ''' Returns (coset_prim, joint_centre):
          coset_prim   : all fields (lam, mu_1, ..., mu_n) allowed by the selection rules
                         (no field identification applied)
          joint_centre : non-trivial elements of the product of the outer automorphism groups that
                         act trivially on the branching (the field identification group, minus 1) '''

    # Get the basic info for the numerator Lie algebra
    nu_rnk = Integer( coset[1][0][1:] )
    nu_letter = coset[1][0][0]
    nu_level = coset[1][1]
    nu_dim, nu_dcox, nu_cmks, nu_theta, nu_F, nu_aut, nu_v, nu_mod = get_prop_of_lie_algebra(nu_letter, nu_rnk, detail=True)

    # Get the basic info for the denominator Lie Algebra
    num_de = len(coset[2])
    de_rnk = [Integer( coset[2][i][0][1:] ) for i in range(num_de)]
    de_level = [coset[2][i][1] for i in range(num_de)]
    de_dim=[]
    de_dcox=[]
    de_cmks=[]
    de_theta=[]
    de_F=[]
    de_aut=[]
    de_v=[]
    de_mod = []
    for i in range(num_de):
        letter = coset[2][i][0][0]
        de_dim0, de_dcox0, de_cmks0, de_theta0, de_F0, de_aut0, de_v0, de_mod0 = get_prop_of_lie_algebra(letter, de_rnk[i], detail=True)
        de_dim.append(de_dim0)
        de_dcox.append(de_dcox0)
        de_cmks.append(de_cmks0)
        de_theta.append(de_theta0)
        de_F.append(de_F0)
        de_aut.append(de_aut0)
        de_v.append(de_v0)
        de_mod.append(de_mod0)

    # Projection matrices as integer matrices
    proj_mat = [matrix(ZZ, coset[3][i][1]) for i in range(num_de)]

    # Get the primaries of the numerator simple Lie algebra
    nu_prim = find_weights(nu_cmks,nu_level)

    # Get the primaries of the denominator simple Lie algebras
    de_prim = [find_weights(de_cmks[i],de_level[i]) for i in range(num_de)]

    # Deal with the branching of the outer automorphism groups
    joint_centre = []
    trivial_aut =  tuple( [tuple(range(nu_rnk+1))] ) + tuple( tuple(range(de_rnk[i]+1)) for i in range(num_de) )
    for poss_aut in itr.product(nu_aut,*de_aut):
        if poss_aut == trivial_aut:
            continue
        for k in range(nu_rnk):
            nu_lam = vector(QQ, nu_rnk, {k:1})
            de_lam = [proj_mat[i]*nu_lam for i in range(num_de)]
            lhs = vector(ZZ, [1 if j==0 else 0 for j in poss_aut[0][1:]])*nu_F*nu_lam
            rhs = sum( vector(ZZ, [1 if j==0 else 0 for j in poss_aut[i+1][1:]])*de_F[i]*de_lam[i] for i in range(num_de) )

            if lhs-rhs not in ZZ:
                break
        else:
            joint_centre.append( poss_aut )

    # Find all coset primaries according to selection rules.
    # Selection rule:
    #   There exists phi in Q_nu such that
    #       P_i(phi) = P_i(lambda) - mu_i
    #   simultaneously for every denominator factor i.
    # The condition phi in Q_nu is imposed by
    #       nu_v * phi = 0 mod nu_mod.
    coset_prim = []
    total_proj = proj_mat[0]
    for i in range(1, num_de):
        total_proj = total_proj.stack(proj_mat[i])
    if nu_mod:
        full_mat = block_matrix(ZZ,
            [
                [ nu_v, -diagonal_matrix(ZZ, nu_mod) ],
                [ total_proj, zero_matrix(ZZ, total_proj.nrows(), nu_v.nrows()) ]
            ]
        )
    else:
        full_mat = total_proj

    for lam in nu_prim:
        nu_lam = vector(ZZ, lam[1:])
        total_de_lam = total_proj * nu_lam
        for mu in itr.product(*de_prim):
            total_mu = vector(ZZ, [ entry for mu_i in mu for entry in mu_i[1:] ])
            difference = total_de_lam - total_mu
            if nu_mod:
                rhs = vector(ZZ, [0] * nu_v.nrows() + list(difference))
            else:
                rhs = difference
            try:
                full_mat.solve_right(rhs, extend=False)
            except ValueError as error:
                if str(error) == "matrix equation has no solutions":
                    continue
                raise error

            coset_prim.append( tuple([lam]) + mu )

    coset_prim = [tuple(tuple(int(x) for x in lab) for lab in field) for field in coset_prim]
    joint_centre = [tuple(tuple(int(x) for x in p) for p in g) for g in joint_centre]
    return coset_prim, joint_centre


def _coset_factors(coset):
    if 'U(1)' in str(coset[0]) or any(str(alg[0]).startswith('U') for alg in coset[2]):
        raise ValueError("Cosets with U(1) factors in the denominator are not supported.")
    return [_CosetFactor(coset[1][0], coset[1][1], False)] + \
           [_CosetFactor(alg[0], alg[1], True) for alg in coset[2]]


# -----------------------------------------------------------------------------
#  MAIN FUNCTION
# -----------------------------------------------------------------------------
def get_coset_modular_data(coset, check=True):
    '''
    Input : coset - list in the format of an entry of all_target_c_cosets
            check - run internal consistency checks (monodromy charges, F well defined,
                    (S^J T^J)^3 = (S^J)^2 for every factor, ...)
    Output: primaries - list of the primaries given by their affine Dynkin labels
                        (numerator labels, denominator labels...). A fixed point of the field
                        identification group appears |U_a| times, once for each resolved field.
            S, T      - modular S and T matrices (over UCF) in the ordering of `primaries`.
                        T = exp(2 pi i (h - c/24)).
    '''
    factors = _coset_factors(coset)
    fields, joint_centre = _selection_rules(coset)
    ident = tuple(f.ident for f in factors)
    group = sorted(_closure(joint_centre, ident))
    return _resolve_fixed_points(factors, group, fields, check=check)


# -----------------------------------------------------------------------------
#  Extras
# -----------------------------------------------------------------------------
def get_conformal_weights(coset, primaries):
    '''Conformal weights h of the primaries returned by get_coset_modular_data.'''
    factors = _coset_factors(coset)
    out = []
    for p in primaries:
        h = sum(f.sign * f.h(p[i]) for i, f in enumerate(factors))
        out.append(QQ(int(h.numerator)) / QQ(int(h.denominator)))
    return out


def _numpy(M):
    n = M.nrows()
    return np.array([[_c(M[i, j]) for j in range(n)] for i in range(n)])


def check_modular_relations(S, T, verbose=True):
    '''Numerical check of: S unitary & symmetric, S^4 = 1, (ST)^3 = S^2.'''
    Sn, Tn = _numpy(S), _numpy(T)
    I = np.eye(Sn.shape[0])
    S2, ST = Sn @ Sn, Sn @ Tn
    res = dict(unitary=np.abs(Sn @ Sn.conj().T - I).max() < _TOL8,
               symmetric=np.abs(Sn - Sn.T).max() < _TOL8,
               S4=np.abs(S2 @ S2 - I).max() < _TOL8,
               ST3=np.abs(ST @ ST @ ST - S2).max() < _TOL8)
    if verbose:
        print('  '.join('%s: %s' % (k, v) for k, v in res.items()))
    return all(res.values())


def check_fusion_rules(S, verbose=True):
    '''Verlinde formula: are all fusion coefficients non-negative integers?
    (The identity is the unique primary whose row of S is real and positive.)'''
    Sn = _numpy(S)
    n = Sn.shape[0]
    cands = [i for i in range(n) if np.all(np.abs(Sn[i].imag) < _TOL) and np.all(Sn[i].real > _TOL)]
    if len(cands) != 1:
        if verbose:
            print('could not identify the vacuum uniquely: candidates', cands)
        return False
    i0 = cands[0]
    ok = True
    for i in range(n):                       # N_{ij}^k = sum_m S_im S_jm conj(S_km) / S_0m
        Ni = (Sn * (Sn[i] / Sn[i0])) @ Sn.conj().T
        Nr = Ni.real
        if np.abs(Ni.imag).max() > _TOL6 * (1 + np.abs(Nr).max()) or \
           np.any(np.abs(Nr - np.round(Nr)) > _TOL * (1 + np.abs(Nr))) or Nr.min() < -_HALF:
            ok = False
            break
    if verbose:
        print('fusion coefficients are non-negative integers: %s  (vacuum = index %d)' % (ok, i0))
    return ok
