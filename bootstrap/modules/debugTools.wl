(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

BeginPackage["debugTools`"]

Time::usage = "Time[expr] evaluates expr and echoes how long it took.";

Begin["`Private`"]

(* A Timing Function *)
SetAttributes[Time, HoldFirst]
Time[x_] := Module[{res=AbsoluteTiming[x]},
    Echo[res[[1]], "Time: "];
    res[[2]]
    ];

End[]
EndPackage[]

(* :!CodeAnalysis::EndBlock:: *)
