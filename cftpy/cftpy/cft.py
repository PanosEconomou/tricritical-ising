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
try:
    from sage.all import SR, QQ, QQbar
    from sage.all import matrix
    from sage.all import pi, gcd, sqrt, sin, exp, I 
    from sage.all import latex
except ImportError:
    raise RuntimeError("cftpy requires SageMath, but Sage cannot be imported.\nInstall SageMath and run with its Python interpreter.")

from .cythonft import cstring_function_su2 # type: ignore

from IPython.display import display, Math
from inspect import signature
from re import sub


#===========================================================#
# Module Parameters and Cache                               #
#===========================================================#

MAX_ORDER           = 1000          # Order to calculate q-expansions in
_STRING_FUNCTIONS   = {}            # Caches string functions once they're calculated
_CHARACTERS         = {}            # Caches charaacters once they're calculated


#===========================================================#
# Weights, Central Chrages, and other Mode Info             #
#===========================================================#

def info(*args, model:str ='minimal', **kwargs) -> dict:
    """Returns information like the labelling system for a specific model

    Args:
        model (string, optional): The CFT you want. Defaults to 'minimal'.
            Other options include:'{"', '".join(MODEL_NAMES)}'
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        k (int, optional): level of the WZW algebra. Defaults to 2.
        l (int, optional): level of the second WZW algebra in the diagonal su(2) coset. Defaults to 2.

    Returns:
        info (dict): Information with fields:
            info['model'] (str): model
            info['labels'] (set): Labels for the irreps of the model with a given structure. 
                For example minimal models returns a set of elements (r,s) corresponding to the Kaz indices
    """

    sig       = signature(MODEL_LABELS[model])
    arguments = sig.bind(*args,**kwargs)
    arguments.apply_defaults()
    
    labels  = MODEL_LABELS[model](*arguments.args,**arguments.kwargs)
    info    =  {
        'model' : model, 
        'labels': labels,
        **kwargs
        }
    return info


def h(*args, model:str = 'minimal', **kwargs):
    """Returns information like the labelling system for a specific model

    Args:
        model (string, optional): The CFT you want. Defaults to 'minimal'.
            Other options include:'{"', '".join(MODEL_NAMES)}'
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        k (int, optional): Level of the WZW algebra. Defaults to 2.
        l (int, optional): Level of the second WZW algebra in the diagonal su(2) coset. Defaults to 2.
        label (tuple): Label of irreducible representation. If unset it returns a list of coformal weights for this theory.

    Returns:
        rational: Conformal weight of module labeled by `label` or
        dict: A dictionary of conformal weights indexed by the theory's labels
    """
    sig       = signature(MODEL_WEIGHTS[model])
    arguments = sig.bind(*args,**kwargs)
    arguments.apply_defaults()

    return MODEL_WEIGHTS[model](*arguments.args,**arguments.kwargs)


def h_minimal(label:tuple = (), p:int = 4, q:int = 3):
    """Conformal weights for irreducible representations of minimal models

    Args:
        label (tuple, optional): (r,s) label of irreducible module. If empty a list is returned.
        p (int, optional): Kac-index. Defaults to 4.
        q (int, optional): Kac-index. Defaults to 3.

    Returns:
        rational: Conformal weight of module (r,s) or
        dict: A dictionary of conformal weights indexed by Kac-labels
    """
    if not label:
        labels = Kac_labels(p,q)
        return {l: h_minimal(l, p, q) for l in labels}
    
    return QQ(((p*label[0] - q*label[1])**2 - (p-q)**2)/(4*p*q))

def h_folded_minimal(label:tuple=(), p:int=4, q:int=3):
    """Conformal weights for irreducible representations of minimal models

    Args:
        label (tuple, optional): ((r,s),(r',s')) label of irreducible module. If empty a list is returned.
        p (int, optional): Kac-index. Defaults to 4.
        q (int, optional): Kac-index. Defaults to 3.

    Returns:
        rational: Conformal weight of module ((r,s),(r',s')) or
        dict: A dictionary of conformal weights indexed by Kac-labels
    """
    if not label:
        labels = folded_Kac_labels(p,q)
        return {l: h_folded_minimal(l, p, q) for l in labels}

    return h_minimal(label[0],p,q) + h_minimal(label[1],p,q)  # type: ignore

def c(*args, model:str = 'minimal',**kwargs) -> QQ:
    """Central charge for a spefic cft

    Args:
        model (string, optional): The CFT you want. Defaults to 'minimal'.
            Other options include:'{"', '".join(MODEL_NAMES)}'
        p (int, optional): Kac-index. Defaults to 4.
        q (int, optional): Kac-index. Defaults to 3.
        k (int, optional): Level of the WZW algebra. Defaults to 2.
        l (int, optional): Level of the second WZW algebra in the diagonal su(2) coset. Defaults to 2.

    Returns:
        rational: central charge 
    """
    sig       = signature(MODEL_CENTRAL_CHARGES[model])
    arguments = sig.bind(*args,**kwargs)
    arguments.apply_defaults()

    return MODEL_CENTRAL_CHARGES[model](*arguments.args,**arguments.kwargs)

def c_minimal(p:int=4, q:int=3) -> QQ:
    """Central charge for a minimal Model

    Args:
        p (int, optional): Kac-index. Defaults to 4.
        q (int, optional): Kac-index. Defaults to 3.

    Returns:
        rational: central charge 
    """
    return QQ(1 - 6*(p-q)**2/(p*q))

def c_folded_minimal(p:int=4, q:int=3) -> QQ:
    """Central charge for a folded minimal Model

    Args:
        p (int, optional): Kac-index. Defaults to 4.
        q (int, optional): Kac-index. Defaults to 3.

    Returns:
        rational: central charge 
    """
    return 2*c_minimal(p,q)

#===========================================================#
# S-Matrix Calculations                                     #
#===========================================================#

def S_matrix(*args, model:str = 'minimal', **kwargs) -> (matrix, set):
    """Returns the smatrix of a Rational CFT as well as a list of indices 
    for the irreducible rerpesentations of the algebra

    Args:
        model (string, optional): The CFT you want. Defaults to 'minimal'.
            Other options include:'{"', '".join(MODEL_NAMES)}'
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        k (int, optional): level of the WZW algebra. Defaults to 2.
        l (int, optional): level of the second WZW algebra in the diagonal su(2) coset. Defaults to 2.
        K (field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.
    
    Returns:
        S (matrix): The s matrix in the field K
        info (dict): Information about the S matrix
            info['labels'] (set): A set of kac-labels for each row/column of S
            info['model'] (str): The model you chose
    """
    
    allowed_kwargs = {'p', 'q', 'k', 'l', 'K'}
    unknown = set(kwargs) - allowed_kwargs
    if unknown:
        raise TypeError(f'Unexpected keword argument(s) {unknown}')

    sig = signature(MODEL_S_MATRIX[model])
    arguments = sig.bind(*args,**kwargs)
    arguments.apply_defaults()

    return MODEL_S_MATRIX[model](*arguments.args,**arguments.kwargs)


def Kac_labels(p:int = 4, q:int = 3) -> list:
    """Return a list of pairs (s,r) of Kac labels for the minimal model p,q

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
    
    Returns:
        list: the set of labels in (r,s) notation
    """
    if gcd(p,q) != 1: raise ValueError("p and q must be coprime")

    labels = []

    for r in range(1, q):
        for s in range(1,p):
            if (r < q - r) or (r == q-r and s <= p - s):
                if (r,s) not in labels:
                    labels.append((r,s))
    
    return labels

def Kac_labels_inverse(label = None, p:int = 4, q:int = 3):
    """Return a dictionary {i:(s,r)} that matches the Kac label (s,r) 
    to this module's internal index i, or it returns the corresponding index given a label 

    Args:
        label (tuple, optional): Kac-label to convert. If unset a dictionary is returned.
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
    
    Returns:
        dict: the map of labels and indices {i:(r,s)}
    """
    if gcd(p,q) != 1: raise ValueError("p and q must be coprime")

    if not label:
        return {i:label for i,label in enumerate(Kac_labels(p, q))}

    return int((label[0]-1)*q + label[1] - 1)

def folded_Kac_labels(p:int = 4,q:int = 3, K = None) -> list:
    """Return a list of pairs ((s,r),(s',r')) of Kac labels for the folded minimal model p,q

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
    
    Returns:
        set: the set of labels in ((r,s),(r',s')) notation
    """
    labels = Kac_labels(p,q)
    return [(l,m) for l in labels for m in labels]

def folded_Kac_labels_inverse(label = None, p:int = 4,q:int = 3, K = None):
    """eturn a dictionary {i:((s,r),(s',r'))} that matches the Kac label ((s,r),(s',r'))
    to this module's internal index i, or it returns the corresponding index given a label 

    Args:
        label (tuple, optional): Kac-label to convert. If unset a dictionary is returned.
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
    
    Returns:
        dict: the map of labels and indices {i:((r,s),(r',s'))}
    """
    if not label:
        return {i:label for i,label in enumerate(folded_Kac_labels(p, q))}

    return int(Kac_labels_inverse(label[0], p, q)*(p-1)*(q-1)//2 + Kac_labels_inverse(label[1], p, q))

def minimal_model_S_matrix(p:int = 4, q:int = 3, K = None) -> (matrix, set):
    """Calculate the S-matrix of a minimal model

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        K (Field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.
    
    Returns:
        S (matrix): The s matrix in the field K
        labels (set): the set of labels in (r,s) notation
    """
    if K == None: K = SR

    labels  = Kac_labels(p,q)
    n       = len(labels)
    S       = matrix(base_ring = K, nrows = n, ncols = n)

    for i,(ri,si) in enumerate(labels):
        for j,(rj,sj) in enumerate(labels):
            S[i,j] = (sqrt(8)/sqrt(p*q)*(-1)**(ri*sj + rj*si + 1)*sin(pi*ri*rj*p/q)*sin(pi*si*sj*q/p)) # type: ignore
            if K == SR:
                S[i,j] = S[i,j].canonicalize_radical()
    return S, labels

def folded_minimal_model_S_matrix(p:int = 4,q:int = 3, K = None) -> (matrix, set):
    """Calculate the S-matrix of a folded minimal model

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        K (Field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.

    Returns:
        S (matrix): The s matrix in the field K
        labels (set): the set of labels in ((r,s),(r',s')) notation
    """
    if K == None: K = SR
    
    s, labels   = minimal_model_S_matrix(p,q,K)
    idx         = {l:i for i,l in enumerate(labels)}
    labels      = [(l,m) for l in labels for m in labels]
    n           = len(labels)
    S           = matrix(base_ring = K, nrows = n, ncols = n)

    for i,l in enumerate(labels):
        for j,m in enumerate(labels):
            S[i,j] = s[idx[l[0]], idx[m[0]]]*s[idx[l[1]], idx[m[1]]]
    
    return S, labels


#===========================================================#
# T-Matrix Calculations                                     #
#===========================================================#

def T_matrix(*args, K = None, **kwargs) -> (matrix, set):
    """Calculate the T-matrix of a minimal model

    Args:
        model (string, optional): The CFT you want. Defaults to 'minimal'.
            Other options include:'{"', '".join(MODEL_NAMES)}'
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        K (Field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.
    
    Returns:
        T (matrix): The T matrix in the field K
        labels (set): the set of labels in the model's notation
    """
    if K == None: K = SR

    labels  = info(*args,**kwargs)['labels']
    cc      = c(*args,**kwargs)
    hh      = h(*args,**kwargs)
    n       = len(labels)
    T       = matrix(base_ring = K, nrows = n, ncols = n)

    for i,label in enumerate(labels):
        T[i,i] = K(exp(2*I*pi*(hh[label] - cc/24)))
    return T, labels


#===========================================================#
# Topological Defect Line Tools                             #
#===========================================================#

def verlinde_line(*args, model:str = 'minimal', **kwargs) -> matrix:
    """Calculates the matrix of how a Verlinde line in the theory acts on the primaries

    Args:
        model (string, optional): The CFT you want. Defaults to 'minimal'.
            Other options include:'{"', '".join(MODEL_NAMES)}'
        label (tuple, optional): Label of conformal defect. Defaults to identity.
        S (matrix, optional): S-matrix for model. If unset, the S-matrix is calculated.
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        k (int, optional): level of the WZW algebra. Defaults to 2.
        l (int, optional): level of the second WZW algebra in the diagonal su(2) coset. Defaults to 2.
        K (field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.
    
    Returns:
        matrix: The matrix L whose diagonal contains how the Verlinde line acts on the corresponding primary.
    """
    
    allowed_kwargs = {'label', 'S', 'p', 'q', 'k', 'l', 'K'}
    unknown = set(kwargs) - allowed_kwargs
    if unknown:
        raise TypeError(f'Unexpected keword argument(s) {unknown}')

    sig = signature(MODEL_VERLINDE[model])
    arguments = sig.bind(*args,**kwargs)
    arguments.apply_defaults()

    return MODEL_VERLINDE[model](*arguments.args,**arguments.kwargs)


def minimal_model_verlinde_line_matrix(p:int = 4, q:int = 3, label:tuple = (1,1), S:matrix = None, K = None) -> matrix:
    """Calculates the matrix of how a Verlinde line from a minimal model acts on the primaries of that theory.

    Args:
        p (int, optional): Kac p-index. Defaults to 4.
        q (int, optional): Kac q-index. Defaults to 3.
        label (tuple, optional): Kac label of conformal defect. Defaults to (1,1).
        S (matrix, optional): S-matrix for model. If unset, the S-matrix is calculated.
        K (Ring, optional): The Ring to calculate the matrix as. Defaults to Symbolic Ring (SR).

    Returns:
        matrix: The matrix L whose diagonal contains how the Verlinde line acts on the corresponding primary.
    """
    if not S: S, _ = minimal_model_S_matrix(p, q, K)
    if not K: K = SR
    i = Kac_labels_inverse(label, p, q)

    L = matrix(base_ring = K, nrows = S.nrows(), ncols = S.ncols())
    for j,Sij in enumerate(S[i]):
        L[j,j] = Sij/S[j,0]

    return L 

def folded_minimal_model_verlinde_line_matrix(p:int = 4, q:int = 3, label:tuple = ((1,1),(1,1)), S:matrix = None, K = None) -> matrix:
    """Calculates the matrix of how a Verlinde line from a folded minimal model acts on the primaries of that theory.

    Args:
        p (int, optional): Kac p-index. Defaults to 4.
        q (int, optional): Kac q-index. Defaults to 3.
        label (tuple, optional): Kac label of conformal defect. Defaults to (1,1).
        S (matrix, optional): S-matrix for model. If unset, the S-matrix is calculated.
        K (Ring, optional): The Ring to calculate the matrix as. Defaults to Symbolic Ring (SR).

    Returns:
        matrix: The matrix L whose diagonal contains how the Verlinde line acts on the corresponding primary.
    """
    if not S: S, _ = folded_minimal_model_S_matrix(p, q, K)
    if not K: K = SR
    i = folded_Kac_labels_inverse(label, p, q)

    L = matrix(base_ring = K, nrows = S.nrows(), ncols = S.ncols())
    for j,Sij in enumerate(S[i]):
        L[j,j] = Sij/S[j,0]

    return L 


#===========================================================#
# Characters (These will be cythonized soon)                #
#===========================================================#

def string_function_su2(l:int = -1, m:int = 0, k:int = 2, order:int = MAX_ORDER):
    """Generates the string function of su(2) at level k for highest weight l, and string weight m

    Args:
        l (int, optional): Highest weight. If left blank all the string functions for level k are calculated
        m (int, optional): String weight.
        k (int, optional): Level of su(2)_k. Defaults to 2.
        order (int, optional): maximum order to calculate the q expansion. Defaults to MAX_ORDER=1000, or whatever the user overrides thsi to 

    Returns:
        PowerSeriesRing: Containing the q expansion for the corresponding string function or
        dict: Containing all string functison if l is unspecified
    """
    
    if l == -1: 
        strings = {}
        # TODO: Optimize the range of the iterations here to account for multiplicites
        for l in range(k+1):
            for m in range(l & 2 + 1,l,2):
                strings[(l,m)] = string_function_su2(l,m,k,order)
        return strings

    if (l - m) & 1 != 0:
        raise ValueError("l - m must be an even integer")

    if 'su(2)' in _STRING_FUNCTIONS:
        if (l,m) in _STRING_FUNCTIONS['su(2)']:
            try:
                if order > _STRING_FUNCTIONS['su(2)'][(l,m)]['order']:
                    _STRING_FUNCTIONS['su(2)'][(l,m)] = {}
                else:
                    return _STRING_FUNCTIONS['su(2)'][(l,m)]['series'].add_bigoh(order)
            except KeyError:
                del _STRING_FUNCTIONS['su(2)'][(l,m)]
    else:
        _STRING_FUNCTIONS['su(2)'] = {}

    _STRING_FUNCTIONS['su(2)'][(l,m)] = {
            'order'  : order, 
            'series' : cstring_function_su2(l,m,k,order)
        }

    return _STRING_FUNCTIONS['su(2)'][(l,m)]['series']




#===========================================================#
# Sage Helpers                                              #
#===========================================================#

def showvs(obj=None, *args, **kwargs):
    """Replacement for Sage's show() that works on stupid VSCode
    """
    if obj is None:
        return
    try:
        display(Math(latex(obj))) 
    except Exception:
        print(obj)

def fill_docstring(func):
    """Fills docstrings for functions with dynamic options-- INTERNAL
    """

    def repl(match):
        expr = match.group(1).strip()
        return str(eval(expr))

    # match {anything-but-braces}
    func.__doc__ = sub(r"\{([^{}]+)\}", repl, func.__doc__)

def delete_cache():
    """Removes all the cached series from memory
    """
    _STRING_FUNCTIONS   = {}
    _CHARACTERS         = {}


#===========================================================#
# Globals                                                   #
#===========================================================#

SPECIAL_NAMES       = {
    'minimal': {
        'ising'                     : lambda f: lambda *args, **kwargs: f(p=4,q=3,*args,**kwargs),
        'tricritical_ising'         : lambda f: lambda *args, **kwargs: f(p=5,q=4,*args,**kwargs),
        'yang_lee'                  : lambda f: lambda *args, **kwargs: f(p=5,q=2,*args,**kwargs),
    },
    'folded_minimal': {
        'folded_ising'              : lambda f: lambda *args, **kwargs: f(p=4,q=3,*args,**kwargs),
        'folded_tricritical_ising'  : lambda f: lambda *args, **kwargs: f(p=5,q=4,*args,**kwargs),
        'folded_yang_lee'           : lambda f: lambda *args, **kwargs: f(p=5,q=2,*args,**kwargs),
    },
}

MODEL_LABELS        = {
    'minimal'                   : Kac_labels,
    'folded_minimal'            : folded_Kac_labels,
}

MODEL_WEIGHTS       = {
    'minimal'                   : h_minimal,
    'folded_minimal'            : h_folded_minimal,
}

MODEL_CENTRAL_CHARGES= {
    'minimal'                   : c_minimal,
    'folded_minimal'            : c_folded_minimal,
}

MODEL_S_MATRIX      = {
    'minimal'                   : minimal_model_S_matrix,
    'folded_minimal'            : folded_minimal_model_S_matrix,
}

MODEL_VERLINDE      = {
    'minimal'                   : minimal_model_verlinde_line_matrix,
    'folded_minimal'            : folded_minimal_model_verlinde_line_matrix,
}

DICTS               = [MODEL_LABELS, MODEL_WEIGHTS, MODEL_CENTRAL_CHARGES, MODEL_S_MATRIX, MODEL_VERLINDE]
DYNAMIC_DOCSTRINGS  = [info, h, c, S_matrix, T_matrix, verlinde_line]
MODEL_NAMES         = list(MODEL_LABELS.keys())

#===========================================================#
#oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo#
#===========================================================#
