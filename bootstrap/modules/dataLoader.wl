(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
(* :!CodeAnalysis::Disable::UppercaseVariable:: *)

BeginPackage["dataLoader`"]

(* A basis is an association that depends on a sector association defined like so
 *
 * sectors = <| "labels"  -> {...}, "index" -> AssociationThread[labels -> Range[n]],
 *              "weights" -> {...}, "S" -> ..., "T" -> ..., "order" -> <sort key> |>
 * 
 * basis   = <| "id"    -> "...", "sectors" -> <ref>, "invariant" -> z,
 *              "slots" -> {{a,b}, ...}, "index" -> AssociationThread[slots -> Range[n]],
 *              "degeneracy" -> {...}, "degPos" -> ..., "nondegPos" -> ...|>
 *)

FoldedMinimalModularData::usage = "Obtain parameters for folded minimal models";

Begin["`Private`"]

MinimalModularData[p_:5, pp_:4] := Module[{
        sectors = <| |>, c, h, s, t
    },
        c = 1 - 6 (p - pp)^2 / (p pp);

        (* Labels for each of the virasoro irreps *)
        sectors["labels"] = Join[Flatten[Table[{r,s}, {r, 1, Floor[(pp - 1)/2]}, {s, 1, p - 1}], 1], If[EvenQ[pp], Table[{pp/2,s}, {s, 1, Floor[(p - 1)/2]}], {}]];
        sectors["index"] = AssociationThread[sectors["labels"] -> Range@Length[sectors["labels"]]];

        (* Conformal Weights *)
        h[r_, s_] := h[r,s] = ((p r - pp s)^2 - (p - pp)^2) / (4 p pp);
        sectors["weights "] = hv @@ # & /@ sectors["labels"];

        (* Modular Data *)
        s[a_, b_] := s[a,b] = 2 Sqrt[2 / (p pp)] (-1)^(1 + a[[2]] b[[1]] + a[[1]] b[[2]]) Sin[a[[1]] b[[1]] Pi p / pp] Sin[a[[2]] b[[2]] Pi pp / p];
        t[a_, e_:1] := tv[a,e] = Exp[2 Pi I e (h[a[[1]], a[[2]]] - cv/24)];

        <| 
            "id" -> StringJoin["M(", ToString[p], ", ", ToString[pp], ")"],
            "sectors" -> sectors
        |>
]

FoldedMinimalFixedPointModularData[p_:5, pp_:4] := Module[{
        cv, hv, sv, tv, H, i, ho, t, sMaker, s, nH, HIndex, svMat, tv2, twist
    },
        cv = 1 - 6 (p - pp)^2 / (p pp);
        hv[r_, s_] := hv[r,s] = ((p r - pp s));
        sv[a_, b_] := sv[a,b] = 2 Sqrt[2 / (p pp)] (-1)^(1 + a[[2]] b[[1]] + a[[1]] b[[2]]) Sin[a[[1]] b[[1]] Pi p / pp] Sin[a[[2]] b[[2]] Pi pp / p];
        tv[a_, e_:1] := tv[a,e] = Exp[2 Pi I e (hv[a[[1]], a[[2]]] - cv/24)];
    ]

End[]
EndPackage[]

(* :!CodeAnalysis::EndBlock:: *)
