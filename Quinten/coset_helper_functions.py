from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.extended_affine_weyl_group import ExtendedAffineWeylGroup
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix


def congruence_data(ct):
    """
    Compute the congruence invariants of P/Q directly from the
    Smith normal form of the Cartan matrix.

    In the fundamental-weight basis,

        P = Z^r
        Q = A Z^r.

    If U*A*V = S is the Smith normal form, then the quotient
    coordinates are

        U * lambda  (mod diagonal entries of S).

    Returns:

        (v, moduli)

    where v is a list of the independent 
    congruence vectors.
    """
    A = matrix(ZZ, ct.cartan_matrix())

    S, U, V = A.smith_form()

    rows = []
    moduli = []

    for i in range(A.nrows()):
        s = abs(ZZ(S[i,i]))

        if s > 1:
            rows.append(vector(ZZ, U.row(i)))
            moduli.append(s)

    return matrix(ZZ,rows), moduli

    
def additional_data(cartan_type):
    """
    Return (comarks, theta, F, v, modulus).

    comarks:
        the comarks of the Lie Algebra
    
    theta:
        highest root in fundamental-weight coordinates.

    quadF:
        Gram matrix of the fundamental weights, with the longest
        roots normalized to have squared norm 2.
    
    outer_aut:
        List of tuples of outer automorphisms

    v:
        congruence vectors for P/Q.

    modulus:
        order of simple abelian groups in P/Q
    """
    ct = CartanType(cartan_type)

    if not ct.is_finite():
        raise ValueError("Cartan type must be finite")

    # Cartan Matrix
    A = ct.cartan_matrix()

    # After normalization, d_i = (alpha_i, alpha_i)/2.
    d = vector(QQ, ct.symmetrizer())
    d /= max(d)

    # Gram matrix of simple roots.
    G = diagonal_matrix(QQ, d) * A

    if G != G.transpose():
        raise ValueError("DA is not symmetric")

    # Gram matrix of fundamental weights.
    quadF = A.transpose().inverse() * G * A.inverse()

    # Root system.
    R = ct.root_system().root_lattice()
    # Highest root (multiply by A to convert to fundamental weight coordinates).
    theta = A * vector(ZZ, R.highest_root())

    # Calculate the comarks via them being in the left kernel of the affine Cartan matrix
    Ahat = ct.affine().cartan_matrix()
    comarks = Ahat.left_kernel().basis()[0]
    comarks /= comarks[0] # Normalise the zeroeth comark to one

    # Now get the outer automorphism group
    Og = ExtendedAffineWeylGroup(ct.affine()).fundamental_group()
    outer_aut = []
    for s in Og.special_nodes():
        outer_aut.append( tuple(Og.action(s)(v) for v in Og.index_set()) )

    # Congruence vector is handled separately below.
    v, modulus = congruence_data(ct)

    return comarks, theta, quadF, outer_aut, v, modulus



def get_prop_of_lie_algebra(letter, rnk, detail = False):
    ''' 
    Function to compute the dimension and dual coxeter number of a Lie Algebra
    Input: letter (string) - cartan label of the lie algbra
            rnk (int) - rank of the lie algebra
            detail (bool) - whether or not to compute theta, F, v, modulus
    Output: dim (int) - dimension of the lie Algebra
            dcox (int) - dual Coxeter number
            *** Only if details=True *** 
            theta (vector) - highest root vector
            comarks (vector) - the affine comarks
            quadF (matrix) - quadratic form
            outer_aut (list) - list of outer automorphisms
            v (matrix) - congruence matrix
            modulus (list) - list of orders of the abelian groups in P/Q
    '''

    # Default values indicate errors
    dim = -1 # Dimension of the Lie Algebra
    dcox = -1 # Rank of the Lie Algebra

    # Go through all the Cartan types and compute the dimension and rnk
    if letter == 'A':
        dim = rnk**2 + 2*rnk
        dcox = rnk + 1
        
    elif letter == 'B':
        dim = 2*rnk**2 + rnk
        dcox = 2*rnk - 1

    elif letter == 'C':
        dim = 2*rnk**2 + rnk
        dcox = rnk + 1

    elif letter == 'D':
        dim = 2*rnk^2 - rnk
        dcox = 2*rnk - 2

    elif letter == 'E':
        if rnk == 8:
            dim = 248
            dcox = 30

        elif rnk == 7:
            dim = 133
            dcox = 18

        elif rnk == 6:
            dim = 78
            dcox = 12

    elif letter == 'F' and rnk == 4:
        dim = 52
        dcox = 9

    elif letter == 'G' and rnk == 2:
        dim = 14
        dcox = 4
    
    # If detail == true we compute more stuff
    if detail:
        comarks, theta, quadF, outer_aut, v, modulus = additional_data([letter,rnk])
        # if letter == 'D' and (rnk/2).denominator()==1:
        #     v = matrix(ZZ, [[0]*(rnk-2)+[1,1]])
        #     modulus = [2]
        return dim, dcox, comarks, theta, quadF, outer_aut, v, modulus

    return dim,dcox



def find_weights(cmks, lvl):
    """
    Finds all non-negative integer tuples `weight` of length len(cmks)
    such that sum(weight[i] * cmks[i]) == lvl.
    """
    n = len(cmks)

    # Non-negative bounds: 0 + 1*w_i >= 0
    ieqs = [[0] + [1 if j == i else 0 for j in range(n)] for i in range(n)]

    # Linear constraint: -lvl + sum(w_i * cmks[i]) == 0
    eqns = [[-lvl] + list(cmks)]

    # Construct the polyhedron and extract integer points
    P = Polyhedron(ieqs=ieqs, eqns=eqns)
    return [tuple(pt) for pt in P.integral_points()]