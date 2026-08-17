(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)
(* :!CodeAnalysis::Disable::UppercaseVariable:: *)

BeginPackage["bootstraptdl`"]

ModularBootstrapTraces::usage = "Given some modular data obtain the traces of the tdls";

Begin["`Private`"]

(* Look I know this looks ugly. I copied it from my mathematica notebook.
 * It's supposed to rescale the inequalities so that there are no denominators.
 * Honestly, if you're reading this, I probably have messed up so please 
 * message me at po524 at nyu dot edu. *)
processInequalities[ineqs_List, vars_List] := Module[{c = 0},
    Flatten[
        Function[lhs,
            With[{coeffs = Append[Table[Coefficient[lhs, v], {v, vars}], lhs /. Thread[vars -> 0]]},
                If[AnyTrue[coeffs, ! IntegerQ[#] &],
                    With[{lcm = LCM @@ Denominator /@ coeffs, w = Symbol["w" <> ToString[++c]]},
                        {w >= 0, Expand[lcm*lhs] - lcm*w == 0}], lhs >= 0]]] /@ (ineqs[[All, 1]]),1]
];

(* Get some modular data of a rational cft and obtain a set of traces for the simple tdls *)
ModularBootstrapTraces[modularData_] := Module[
    {
        s       = modularData["S"],
        modules = modularData["modules"],
        vacuum, st, pos, X, m, coeff, A, eqn, multipliers, basis, traces
    },
    Needs["math4ti2`"];

    vacuum = First @ FirstPosition[modules["weights"], {0, 0}];
    st  = ConjugateTranspose[s];
    pos = modules["positions"];

    X = Array[\[FormalX], Length[pos]];
    m = st . ReplacePart[modularData["Z"], Thread[pos -> X]] . s // Flatten;

    (* Coefficient of each variable in the S-transformed mass matrix *)
    coeff = CoefficientArrays[m, X][[2]] // Normal // N // Chop;

    (* Pick a maximal independent set of equations and invert them *)
    A = Chop @ Inverse @ N @ coeff[[Flatten[
        FirstPosition[#, 1] & /@ DeleteCases[
            Chop @ RowReduce @ N @ Transpose @ coeff, 
            {0 ..}
        ]
    ]]];

    (* Re-express the constraints in that basis *)
    m   = Chop[N[coeff . A]] . X;
    eqn = DeleteDuplicates @ Rationalize @ Chop @ Expand[# >= 0 & /@ m];

    (* Rescale variables so the system has integer coefficients *)
    multipliers = Min @ Cases[eqn, (a_. * #) >= 0 :> Abs[a]]^-1 & /@ X;
    eqn = (#[[0]][Subtract @@ # // ExpandAll, 0]) & /@ (eqn /. Thread[X -> multipliers*X]);
    eqn = processInequalities[eqn, X];

    (* Hilbert basis of the cone: the candidate simple objects *)
    basis  = math4ti2`zsolve[eqn];
    basis  = basis[[2]][[All, Flatten[Position[basis[[-1]], _[_Integer]]]]];
    traces = Chop[A . (multipliers * #)] & /@ basis;

    <|
        "theory"  -> modularData["theory"],
        "algebra" -> modularData["algebra"],
        "modules" -> modules,
        "defects" -> SortBy[traces, #[[1]] &]
    |>
];

End[]
EndPackage[]


(* :!CodeAnalysis::EndBlock:: *)
