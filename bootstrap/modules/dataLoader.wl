(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
(* :!CodeAnalysis::Disable::UppercaseVariable:: *)

BeginPackage["dataLoader`"]

ModularData::usage                   = "Define a new conformal field theory";
ModularDataFold::usage               = "From a cft obtain the folded one";
ModularDataFoldZ2Fixed::usage        = "From a cft obtain the folded in a smaller VOA";

ModularDataMinimal::usage            = "Easy access to minimal models";
ModularDataFoldMinimal::usage        = "Easily fold minimal models";
ModularDataFoldZ2FixedMinimal::usage = "Fold minimal models in a smaller voa";

Begin["`Private`"]

(* ADE modular invariants *)
kacLabels[p_, pp_] := SortBy[
    Join[
        Flatten[Table[{r, s}, {r, 1, Floor[(pp - 1)/2]}, {s, 1, p - 1}], 1],
        If[EvenQ[pp], Table[{pp/2, s}, {s, 1, Floor[(p - 1)/2]}], {}]
    ], {((p #[[1]] - pp #[[2]])^2 - (p - pp)^2)/(4 p pp) &, Identity}
];

(* {b, ...} -> Sum_b (Sum_{a in b} e_a) (Sum_{a in b} e_a)^T *)
blocks[n_, bs_] := Total[Outer[Times, #, #] & /@
    Map[Table[Boole[MemberQ[#, a]], {a, n - 1}] &, bs]];

ade[n_, "A"] := IdentityMatrix[n - 1];
ade[n_, "D"] /; EvenQ[n] && n >= 6 := Table[
    If[Mod[n, 4] == 2,
        Boole[OddQ[a]] (Boole[a == b] + Boole[a + b == n]),   (* extension   *)
        If[OddQ[a], Boole[a == b], Boole[a + b == n]]],       (* automorphism *)
    {a, n - 1}, {b, n - 1}];
ade[12, "E"] := blocks[12, {{1, 7}, {4, 8}, {5, 11}}];
ade[18, "E"] := blocks[18, {{1, 17}, {5, 13}, {7, 11}, {3, 9, 15}}] - blocks[18, {{3, 15}}];
ade[30, "E"] := blocks[30, {{1, 11, 19, 29}, {7, 13, 17, 23}}];

MinimalModularInvariant[p_, pp_, invariant_] := Module[
    {
        labels = kacLabels[p, pp], mr, ms, lift, z
    },
    (* the D/E factor sits on the label with even Coxeter number *)
    {mr, ms} = If[EvenQ[pp],
        {ade[pp, invariant], IdentityMatrix[p - 1]},
        {IdentityMatrix[pp - 1], ade[p, invariant]}
        ];

    (* sum over both lifts of each label, halving the (r,s) ~ (pp-r,p-s) double count *)
    lift[{r_, s_}] := {{r, s}, {pp - r, p - s}};
    z[{r_, s_}, {a_, b_}] := mr[[r, a]] ms[[s, b]];
    Table[Sum[z[x, y], {x, lift[f]}, {y, lift[g]}]/2, {f, labels}, {g, labels}]
]

(* Extracts an association of the modules in a given theory *)
ModuleData[labels_, weights_, z_] := Module[
    {
        modules = Position[z, Except[0], {2}, Heads->False],
        moduleLables
    },
    moduleLables = labels[[#]] & /@ modules;
    <|
        "labels"        -> moduleLables,
        "index"         -> AssociationThread[moduleLables -> Range@Length[moduleLables]],
        "multiplicities"  -> Extract[z, modules],
        "weights"       -> (weights[[#]] & /@ modules)
    |>
];


ModularData[labels_, weights_, s_, t_, z_, theory_, algebra_, c_] := Module[
    {
        order = Ordering[Transpose[{weights, labels}]],
        l, w
    },
    l = labels[[order]];
    w = weights[[order]];
    <|
        "irreps" -> <|
            "labels"  -> l,
            "index"   -> AssociationThread[l -> Range@Length[l]],
            "weights" -> w
        |>,
        "S"       -> s[[order, order]],
        "T"       -> t[[order, order]],
        "Z"       -> z[[order, order]],
        "modules" -> ModuleData[l, w, z[[order, order]]],
        "theory"  -> theory,
        "algebra" -> algebra,
        "c"       -> c
    |>
];

(* Get the modular data of a minimal model in a way that everyone can understand  *)
ModularDataMinimal[p_:5, pp_:4, invariant_: "A"] := Module[
    {
        c       = 1 - 6 (p - pp)^2 / (p pp), 
        labels  = kacLabels[p, pp], 
        hf, sf, tf,
        h , s , t , z, theory, algebra
    },
        (* Conformal Weights *)
        hf[r_, s_] := ((p r - pp s)^2 - (p - pp)^2) / (4 p pp);
        h           = hf @@ # & /@ labels;

        (* Modular Data *)
        sf[a_, b_]   := 2 Sqrt[2 / (p pp)] (-1)^(1 + a[[2]] b[[1]] + a[[1]] b[[2]]) Sin[a[[1]] b[[1]] Pi p / pp] Sin[a[[2]] b[[2]] Pi pp / p];
        tf[a_, e_:1] := Exp[2 Pi I e (hf[a[[1]], a[[2]]] - c/24)];
        s = Outer[sf, labels, labels, 1];
        t = DiagonalMatrix[tf /@ labels];
        z = MinimalModularInvariant[p, pp, invariant];

        (* Some Metadata*)
        theory  = "M("   <> ToString[p] <>  ", "  <> ToString[pp] <> ")";
        algebra = "Vir(" <> ToString@InputForm[c] <> ")";

        ModularData[labels, h, s, t, z, theory, algebra, c]
];


(* Take any theory's modular data and fold it *)
ModularDataFold[md_] := ModularData[
    Tuples[md["irreps"]["labels"], 2],                  (* Labels        *)
    Total /@ Tuples[md["irreps"]["weights"], 2],        (* Weights       *)
    KroneckerProduct[md["S"], md["S"]],                 (* S-matrix      *)
    KroneckerProduct[md["T"], md["T"]],                 (* T-matrix      *)
    KroneckerProduct[md["Z"], md["Z"]],                 (* Z-matrix      *)
    md["theory"]  <> " x " <> md["theory"],             (* Theory        *)
    md["algebra"] <>  "x"  <> md["algebra"],            (* Algebra       *)
    md["c"] * 2                                         (* Center Charge *)
];

(* Data of folding a minimal model *)
ModularDataFoldMinimal[p_:5, pp_:4, invariant_: "A"] := ModularDataFold[
    ModularDataMinimal[p, pp, invariant]
];

ModularDataFoldZ2Fixed[md_ ] := Module[
    {
        labels = md["irreps"]["labels"],
        h      = md["irreps"]["weights"],
        s      = md["S"],
        c      = md["c"],
        n, vacuum, indices, hf, weights, twist, sf, z
    },
    n       = Length[labels];
    vacuum  = First @ FirstPosition[h, 0];

    (* Calculate the new indices *)
    indices = Join [
        Join @@ Table[{{i, i, 1}, {i, i, -1}}, {i, n}],
        {#[[1]], #[[2]], 0}& /@ Subsets[Range[n], {2}],
        Join @@ Table[{{i, i, 2}, {i, i, -1}}, {i, n}]
    ];

    (* New conformal weights  *)
    hf[{a_, b_, k_}] := If[Abs[k] == 2,
        h[[a]]/2 + c/16 - (k - 2)/8    + Boole[{a, b, k} === {vacuum, vacuum, -2}],
        h[[a]] + h [[b]] + k (k - 1)/2 + Boole[{a, b, k} === {vacuum, vacuum, -1}]
    ];
    weights = hf /@ indices;

    (* Now for the new S-matrix  *)
    twist = s . MatrixPower[md["T"], 2] . s;
    sf[{a_, b_, j_}, {x_, y_, k_}] := Which[
        j == 0      && k == 0,      s[[a, x]] s[[b, y]] + s[[a, y]] s[[b, x]],
        j == 0      && Abs[k] == 1, s[[a, x]] s[[b, x]],
        Abs[j] == 1 && k == 0,      s[[a, x]] s[[a, y]],
        Abs[j] == 1 && Abs[k] == 1, 1/2 s[[a, x]]^2,
        Abs[j] == 1 && Abs[k] == 2, 1/2 j s[[a, x]],
        Abs[j] == 2 && Abs[k] == 1, 1/2 k s[[a, x]],
        Abs[j] == 2 && Abs[k] == 2, 1/8 j k twist[[a, x]] Exp[Pi I (h[[a]] + h[[x]] - c/12)],
        True, 0
    ];

    (* The modular invariant *)
    z = Outer[(Boole[#1 === #2] + Boole[#1[[1]] === #2[[1]] && #1[[2]] === #2[[2]] && #1[[3]] === -#2[[3]]]) * Boole[Abs[#1[[3]]] != 2] Boole[Abs[#2[[3]]] != 2] &,states, states, 1];

    ModularData[
        {labels[[#[[1]]]], labels[[#[[2]]]], #[[3]]} & /@ indices,
        weights,
        Outer[sf, indices, indices, 1],
        DiagonalMatrix[Exp[2 Pi I (# - 2 c/24)] & /@ weights],
        z,
        md["theory"]  <> " x " <> md["theory"],
        md["algebra"] <> "x"   <> md["algebra"] <> "^Z2",
        2 * c
    ]
];

(* Do the same thing for the a minimal modela a bit faster *)
ModularDataFoldZ2FixedMinimal[p_: 5, pp_: 4] := ModularDataFoldZ2Fixed[
    ModularDataMinimal[p, pp, "A"]
];


(* Get the modular data of a theory and output a basis for its ishisbashi states *)
IshibashiBasis[modularData_] := Module[
    {
        m = modularData["modules"],
        spinless, modules, labels
    },
    spinless = Flatten@Position[m["weights"], {x_, x_}, {1}];
    modules  = Join @@ Table[
        {i, k},
        {i, spinless},
        {k, m["multiplicities"][[i]]}
    ];
    labels = {m["labels"][[First@#]], Last@#} & /@ modules;
    <|
        "labels"  -> labels,
        "index"   -> AssociationThread[labels -> Range@Length[labels]],
        "weights" -> (m["weights"][[First@#, 1]] & /@ modules)
    |>
];


End[]
EndPackage[]

(* :!CodeAnalysis::EndBlock:: *)
