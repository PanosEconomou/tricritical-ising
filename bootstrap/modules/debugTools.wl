(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

BeginPackage["debugTools`"]

Time::usage     = "Time[expr] evaluates expr and echoes how long it took.";
ToRoots::usage  = "Convert elements of a matrix to closest root approximant.";

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
ToRoots[x_List,roundError_:10^-10]:=Module[
    {
        unique = DeleteDuplicates[Flatten@x,Equal@@Round[{#1,#2},roundError]&], 
        roots, rules
    },
    roots = (ToRadicals@RootApproximant[#] &) /@ unique;
    rules = AssociationThread[Round[unique, roundError], roots];

    x /. a_Real :> rules[Round[a, roundError]]
]

End[]
EndPackage[]

(* :!CodeAnalysis::EndBlock:: *)
