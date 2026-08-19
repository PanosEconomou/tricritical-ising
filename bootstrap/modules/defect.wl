(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
(* :!CodeAnalysis::Disable::UppercaseVariable:: *)

BeginPackage["defect`"]

(* Defect Mathematica Helpers  *)
ShowDefect::usage = "Display Defect as a collection of matrices over modules";

(* Defect Creation *)
DefectVariables::usage = "From the modules of a theory get the appropriate blank defect";
DefectVariablesWithTraces::usage = "Use the traces to construct a blank defect";
IdentityDefect::usage = "I mean, just get the identity";

(* Defect Algebra *)
DefectFusion::usage = "DefectFusion[a, b] fuses two defects blockwise. Also available as a \[CircleTimes] b.";
tr::usage = "tr[defect] returns a list of the traces on each subspace";
inv::usage = "inv[defect] inverts each subspace map if possible";
dual::usage = "dual[defect] gets the representation of the dual defect";

Begin["`Private`"]

(* Defect Mathematica Helpers  *)
ShowDefect[defect_] := MatrixForm[MatrixForm /@ defect];


(* Defect Creation *)
DefectVariables[modules_, var_Symbol] := MapIndexed[
    Function[{multiplicity, module},
        Table[var[First @ module, i, j], {i, multiplicity}, {j, multiplicity}]
    ],
    modules["multiplicities"]
];

DefectVariablesWithTraces[modules_, var_Symbol, traces_] := MapThread[
    Function[{trace, block},
        Module[
            {
                n = Length[block],
                matrix = block
            },
            matrix[[n, n]] = trace - Total[Diagonal[matrix][[;; n - 1]]];
            matrix
        ]
    ],
    {traces, DefectVariables[modules, var]}
];

IdentityDefect[modules_, var_Symbol] := IdentityMatrix /@ modules["multiplicities"];

(* Defect Algebra *)
DefectFusion[a_, b_] := MapThread[Dot, {a, b}];
Unprotect[CircleTimes];
SetAttributes[CircleTimes, {Flat, OneIdentity}];
a_List \[CircleTimes] b_List /; Length[a] === Length[b] := DefectFusion[a, b];
Protect[CircleTimes];

tr[defect_]   := Tr /@ defect;
inv[defect_]  := Inverse /@ defect;
dual[defect_] := ConjugateTranspose /@ defect;

End[]
EndPackage[]


(* :!CodeAnalysis::EndBlock:: *)
