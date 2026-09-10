(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

BeginPackage["debugTools`"]

Time::usage     = "Time[expr] evaluates expr and echoes how long it took.";
ToRoots::usage  = "Convert elements of a matrix to closest root approximant.";
ToFile::usage   = "Export lists of algebraic numbers to a file";

Begin["`Private`"]

(* A Timing Function *)
SetAttributes[Time, HoldFirst]
Time[x_] := Module[
    {
        res = AbsoluteTiming[x]
    },
    Echo[res[[1]], "Time: "];
    res[[2]]
    ];


(* A small unitility to efficiently convert to algebraic numbers *)
ToRoots[x_List, roundError_:10^-10]:=Module[
    {
        unique = DeleteDuplicates[Flatten@x,Equal@@Round[{#1,#2},roundError]&], 
        roots, rules
    },
    roots = (ToRadicals@RootApproximant[#] &) /@ unique;
    rules = AssociationThread[Round[unique, roundError], roots];

    x /. a_Real :> rules[Round[a, roundError]]
];

(* A way to send symbolic matrices to a reasonable C format *)
ToFile[m_List, filename_, precision_: 30] := Module[
    {
        dims = Dimensions[m], 
        vals = Union@Flatten[m],
        idx, PolyCoef, NumValue, ParseVal, output
    },
    idx = AssociationThread[vals -> Range @ Length[vals]];
    
    PolyCoef[x_] := CoefficientList[MinimalPolynomial[x, y], y]; 
    NumValue[x_] := Round[N[x, precision + 9] * (10^precision)];
    ParseVal[x_] := Module[
        {
            coefficients = PolyCoef[x]
        },
        {Length[coefficients] - 1} ~Join~ coefficients ~Join~ {NumValue[x]}
    ];

    output = Join[
        {Length[dims]},
        {dims},
        {{Length[vals], precision}},
        ParseVal /@ vals,
        idx /@ Flatten[m]
    ];

    Export[filename, output, "Table"];
    output
];

End[]
EndPackage[]

(* :!CodeAnalysis::EndBlock:: *)
