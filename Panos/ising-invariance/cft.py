#===========================================================#
#                        ┬┌─┐┬┌┐┌┌─┐                        #
#                        │└─┐│││││ ┬                        #
#                        ┴└─┘┴┘└┘└─┘                        #
#                                                           #
# A small library to calculate things about minimal models  #
# using allmighty sage.                                     #
#                                                           #
# The interpreter is located in:                            #
# /private/var/tmp/sage-10.7-current/local/bin/python3      #
#===========================================================#


#===========================================================#
# IMPORTS                                                   #
#===========================================================#

from sage.all import Field, ComplexField, SR
from sage.all import matrix
from sage.all import pi, gcd, sqrt, sin
from sage.all import latex

from IPython.display import display, Math
from inspect import signature

#===========================================================#
# S-Matrix Calculations                                     #
#===========================================================#

def get_S_matrix(*args,model='minimal',**kwargs):
    """Returns the smatrix of a Rational CFT as well as a list of indices 
    for the irreducible rerpesentations of the algebra

    Args:
        model (string, optional): The CFT you want. Defaults to 'minimal'.
            Other options include 'su(2)', 'parafermion', 'su(2)_diagonal_coset', 'fold_minimal'
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        k (int, optional): level of the WZW algebra. Defaults to 2.
        l (int, optional): level of the second WZW algebra in the diagonal su(2) coset. Defaults to 2.
        K (field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.
    """
    
    allowed_kwargs = {'p', 'q', 'k', 'l', 'K'}
    unknown = set(kwargs) - allowed_kwargs
    if unknown:
        raise TypeError(f'Unexpected keword argument(s) {unknown}')

    sig = signature(S_MATRIX_MODELS[model])
    arguments = sig.bind(*args,**kwargs)
    arguments.apply_defaults()

    return S_MATRIX_MODELS[model](*arguments.args,**arguments.kwargs)


def get_Kac_labels(p:int = 4, q:int = 3):
    """Return a list of pairs (s,r) of Kac labels for the minimal model p,q

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
    """
    if gcd(p,q) != 1: raise ValueError("p and q must be coprime")

    labels = set()

    for r in range(1, q):
        for s in range(1,p):
            if (r < q - r) or (r == q-r and s <= p - s):
                labels.add((r,s))
    
    return labels


def get_minimal_model_S_matrix(p:int = 4,q:int = 3, K = None):
    """Calculate the S-matrix of a minimal model

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        K (Field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.
    """
    if K == None: K = SR

    labels  = get_Kac_labels(p,q)
    n       = len(labels)
    S       = matrix(base_ring = K, nrows = n, ncols = n)

    for i,(ri,si) in enumerate(labels):
        for j,(rj,sj) in enumerate(labels):
            S[i,j] = (sqrt(8/p*q)*(-1)**(ri*sj + rj*si + 1)*sin(pi*ri*rj*p/q)*sin(pi*si*sj*q/p))
            if K == SR:
                S[i,j] = S[i,j].canonicalize_radical()
    return S, labels

def get_folded_minimal_model_S_matrix(p:int = 4,q:int = 3, K = None):
    """Calculate the S-matrix of a folded minimal model

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        K (Field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.
    """
    s = get_minimal_model_S_matrix(p,q,K)
    return s.kronecker_product(s) # type: ignore

#===========================================================#
# Sage Helpers                                              #
#===========================================================#

def show(obj=None, *args, **kwargs):
    """Replacement for Sage's show() that works on stupid VSCode
    """
    if obj is None:
        return

    try:
        display(Math(latex(obj))) 
    except Exception:
        print(obj)


#===========================================================#
# Globals                                                   #
#===========================================================#

S_MATRIX_MODELS = {
        'minimal': get_minimal_model_S_matrix,
        'folded_minimal': get_folded_minimal_model_S_matrix
}

#===========================================================#
#oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo#
#===========================================================#
