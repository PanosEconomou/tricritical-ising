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

from sage.all import SR, QQ
from sage.all import matrix
from sage.all import pi, gcd, sqrt, sin
from sage.all import latex

from IPython.display import display, Math
from inspect import signature
from re import sub

#===========================================================#
# Weights, Central Chrages, and other Mode Info             #
#===========================================================#

def get_info(model='minimal', *args, **kwargs) -> dict:
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


def h(model='minimal', *args, **kwargs):
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


def h_minimal(label:tuple=(), p:int=4, q:int=3):
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
        labels = get_Kac_labels(p,q)
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
        labels = get_folded_Kac_labels(p,q)
        return {l: h_folded_minimal(l, p, q) for l in labels}

    return h_minimal(label[0],p,q) + h_minimal(label[1],p,q)  # type: ignore

def c(model='minimal',*args,**kwargs) -> QQ:
    """Central charge for a spefic cft

    Args:
        model: The CFT you want. Defaults to 'minimal'.
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

def get_S_matrix(*args,model='minimal',**kwargs) -> (matrix, set):
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


def get_Kac_labels(p:int = 4, q:int = 3) -> set:
    """Return a list of pairs (s,r) of Kac labels for the minimal model p,q

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
    
    Returns:
        set: the set of labels in (r,s) notation
    """
    if gcd(p,q) != 1: raise ValueError("p and q must be coprime")

    labels = set()

    for r in range(1, q):
        for s in range(1,p):
            if (r < q - r) or (r == q-r and s <= p - s):
                labels.add((r,s))
    
    return labels

def get_folded_Kac_labels(p:int = 4,q:int = 3, K = None) -> set:
    """Return a list of pairs ((s,r),(s',r')) of Kac labels for the folded minimal model p,q

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
    
    Returns:
        set: the set of labels in ((r,s),(r',s')) notation
    """
    labels = get_Kac_labels(p,q)
    return set([(l,m) for l in labels for m in labels])


def get_minimal_model_S_matrix(p:int = 4,q:int = 3, K = None) -> (matrix, set):
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

    labels  = get_Kac_labels(p,q)
    n       = len(labels)
    S       = matrix(base_ring = K, nrows = n, ncols = n)

    for i,(ri,si) in enumerate(labels):
        for j,(rj,sj) in enumerate(labels):
            S[i,j] = (sqrt(8)/sqrt(p*q)*(-1)**(ri*sj + rj*si + 1)*sin(pi*ri*rj*p/q)*sin(pi*si*sj*q/p)) # type: ignore
            if K == SR:
                S[i,j] = S[i,j].canonicalize_radical()
    return S, labels

def get_folded_minimal_model_S_matrix(p:int = 4,q:int = 3, K = None) -> (matrix, set):
    """Calculate the S-matrix of a folded minimal model

    Args:
        p (int, optional): Kac-index p. Defaults to 4.
        q (int, optional): Kac-index q. Defaults to 3.
        K (Field, optional): Sage field to calculate the matrix in. Defaults to SymbolicRing.

    Returns:
        S (matrix): The s matrix in the field K
        labels (set): the set of labels in ((r,s),(r',s')) notation
    """
    s, labels = get_minimal_model_S_matrix(p,q,K)
    labels    = set([(l,m) for l in labels for m in labels])

    return s.tensor_product(s), labels


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

#===========================================================#
# Globals                                                   #
#===========================================================#

SPECIAL_NAMES       = {
    'minimal': {
        'ising'                     : lambda f: lambda *args, **kwargs: f(p=4,q=3,*args,**kwargs),
        'tricritical_ising'         : lambda f: lambda *args, **kwargs: f(p=5,q=4,*args,**kwargs)
    },
    'folded_minimal': {
        'folded_ising'              : lambda f: lambda *args,**kwargs: f(p=4,q=3,*args,**kwargs),
        'folded_tricritical_ising'  : lambda f: lambda *args,**kwargs: f(p=5,q=4,*args,**kwargs)
    },
}

MODEL_LABELS        = {
    'minimal'                   : get_Kac_labels,
    'folded_minimal'            : get_folded_Kac_labels,
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
    'minimal'                   : get_minimal_model_S_matrix,
    'folded_minimal'            : get_folded_minimal_model_S_matrix,
}

DICTS               = [MODEL_LABELS, MODEL_WEIGHTS, MODEL_CENTRAL_CHARGES, MODEL_S_MATRIX]

# Populate the dictionaries with examples listed in SPECIAL_NAMES
for name in SPECIAL_NAMES:
    for labels in DICTS:
        for example in SPECIAL_NAMES[name]:
            labels[example] = SPECIAL_NAMES[name][example](labels[name])

MODEL_NAMES         = list(MODEL_LABELS.keys())
DYNAMIC_DOCSTRINGS  = [get_info, h, c, get_S_matrix]

# Assing the appropriate dynamic docstring
for func in DYNAMIC_DOCSTRINGS: fill_docstring(func)

#===========================================================#
#oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo.oOo#
#===========================================================#
