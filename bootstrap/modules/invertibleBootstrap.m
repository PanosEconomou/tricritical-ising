(* ::Package:: *)
(* :!CodeAnalysis::BeginBlock:: *)
(* :!CodeAnalysis::Disable::SuspiciousSessionSymbol:: *)

(* TODO: REDO THIS ENTIRE THING!!!! *)

BeginPackage["invertibleBootstrap`"]

InvertibleGroupCandidates::usage =
  "InvertibleGroupCandidates[bootstrap] returns the finite groups compatible \
with the g = 1 families of a modular bootstrap, as a list of associations.";
InvertibleBlocks::usage =
  "InvertibleBlocks[bootstrap, candidate] returns one defect per group element, \
each a list of matrices indexed by modules.";
DeterminedBlocks::usage =
  "DeterminedBlocks[bootstrap] gives, for each invertible family, the blocks that \
are already fixed by the traces alone, i.e. every multiplicity-1 module, plus any \
multiplicity-m module whose trace is +-m.";

Begin["`Private`"]

$gapScript = FileNameJoin[{DirectoryName[$InputFileName], "candidates.g"}];

$gapLastScript = None;
$gapLastResult = None;

(* Runs a GAP body and returns the parsed contents of the file it wrote.
   The body refers to its output path as $OUT$.

   Load-bearing details:
     - no -A: it disables package autoloading, and SmallGroup is in a package
     - "" as standard input: closes the pipe, so a break loop cannot block
     - --quitonbreak: errors exit instead of prompting, but still exit 0,
       so the check below reads stderr rather than the exit code
   Output goes to a file, not stdout: GAP wraps stdout at screen width and the
   backslash continuations are a syntax error on the Mathematica side. *)
runGap[bodyTemplate_String] := Module[{stem, script, outfile, result, text},
    stem = FileNameJoin[{$TemporaryDirectory,
             "invertible-" <> ToString[$SessionID] <> "-" <> ToString[Unique[]]}];
    script = stem <> ".g"; outfile = stem <> ".m";
    Quiet@DeleteFile[outfile];
    Export[script, "Read(\"" <> $gapScript <> "\");\n" <>
        StringReplace[bodyTemplate, "$OUT$" -> outfile] <> "\nQUIT;\n", "Text"];
    result = RunProcess[{"gap", "-q", "-b", "--quitonbreak", script}, All, ""];
    $gapLastScript = script; $gapLastResult = result;
    If[StringContainsQ[result["StandardError"], "Error"] || ! FileExistsQ[outfile],
        Message[runGap::err, result["StandardError"], script]; Return[$Failed]];
    text = Import[outfile, "Text"];
    DeleteFile /@ {script, outfile};
    ToExpression[text]
];
runGap::err = "GAP reported an error:\n`1`\nScript kept at `2`";
InvertibleGroupCandidates::gap    = "GAP failed: `1`";
InvertibleGroupCandidates::nogap  = "GAP was not found on PATH.";
InvertibleGroupCandidates::cyclo  = "Some characters are not rational; the GAP \
bridge only transports rationals, so those entries come back as $Failed.";

(* Get the invertible defects *)
invertibleClasses[bootstrap_] := Module[
    {vacuum = First @ FirstPosition[bootstrap["modules"]["weights"], {0, 0}]},
    Select[bootstrap["defects"], Chop[N[#[[vacuum]] - 1]] == 0 &]
];

(* --- the character data ---------------------------------------------------
   The invertible families are the CONJUGACY CLASSES of the invertible symmetry
   group G, not its elements: conjugate lines have equal traces. A module of
   multiplicity m carries an m-dimensional representation of G, and the family
   traces restricted to that module are its character. *)

characterData[bootstrap_] := Module[
    {mult = bootstrap["modules"]["multiplicities"], traces, identity, order},
    traces  = invertibleClasses[bootstrap];
    (* the identity line acts as the identity matrix on every module, so its
       trace vector is exactly the list of multiplicities *)
    identity = FirstPosition[traces, t_ /; Chop[N[t - mult]] === 0*mult, Missing[], 1];
    If[MissingQ[identity] || identity === Missing[],
        Message[InvertibleGroupCandidates::noid]; Return[$Failed]];
    order = Prepend[Delete[Range@Length[traces], identity], First[identity]];
    <|
        "classes"        -> order,                       (* into bootstrap["defects"] *)
        "nclasses"       -> Length[traces],
        "multiplicities" -> mult,
        "maxdeg"         -> Max[mult],
        (* one row per DISTINCT module character; duplicates carry no extra
           information and multiply the GAP search for nothing *)
        "characters"     -> DeleteDuplicates[Transpose[traces[[order]]]]
    |>
];
InvertibleGroupCandidates::noid = "No g = 1 family acts as the identity on every \
module; the traces are inconsistent with an invertible line.";

gapNumber[x_] := Which[
    IntegerQ[x],  ToString[x],
    Head[x] === Rational, ToString[Numerator[x]] <> "/" <> ToString[Denominator[x]],
    True, $Failed];

(* --- the search ----------------------------------------------------------- *)

Options[InvertibleGroupCandidates] = {"MaxOrder" -> Automatic};

InvertibleGroupCandidates[bootstrap_, OptionsPattern[]] := Module[
    {data, maxOrder, entries, body, out},
    If[FindFile["gap"] === $Failed && RunProcess[{"which", "gap"}]["ExitCode"] =!= 0,
        Message[InvertibleGroupCandidates::nogap]; Return[$Failed]];
    data = characterData[bootstrap];
    If[data === $Failed, Return[$Failed]];

    entries = Map[gapNumber, data["characters"], {2}];
    If[! FreeQ[entries, $Failed],
        Message[InvertibleGroupCandidates::cyclo]; Return[$Failed]];

    (* |G| = sum of squares of irrep degrees; there are `nclasses` of them and
       none exceeds the largest multiplicity, so this bound is safe. *)
    maxOrder = OptionValue["MaxOrder"] /. Automatic -> data["nclasses"] * data["maxdeg"]^2;

    body = StringJoin[
        "chars := [", StringRiffle[("[" <> StringRiffle[#, ","] <> "]") & /@ entries, ","], "];;\n",
        "res := CandidateGroups(", ToString[data["nclasses"]], ",",
            ToString[data["maxdeg"]], ", chars, ", ToString[maxOrder], ");;\n",
        "WriteCandidates(\"$OUT$\", res);\n"];

    out = runGap[body];
    If[out === $Failed, Return[$Failed]];

    (* Several class orderings usually match: that is the character table's own
       symmetry, not extra physics. Report distinct GROUPS, keeping the
       orderings so a later fusion check can pick among them. *)
    Association[
        "classes"    -> data["classes"],
        "candidates" -> (Module[{g = First[#]},
             <|"name" -> g["name"], "order" -> g["order"], "id" -> g["id"],
               "degrees" -> g["degrees"],
               "classOrderings" -> (#["classOrder"] & /@ #),
               "decomposition"  -> g["decomposition"]|>] & /@
            GatherBy[out, {#["order"], #["id"]} &])
    ]
];

(* --- what the traces already fix ------------------------------------------ *)

DeterminedBlocks[bootstrap_] := Module[
    {
        mult = bootstrap["modules"]["multiplicities"], 
        classes = invertibleClasses[bootstrap]
    },
    Function[cl, MapThread[
        Which[
            #2 === 1,          {{#1}},                       (* 1x1: block IS the trace *)
            #1 ==  #2,         IdentityMatrix[#2],           (* trace = +m on a unitary: +I *)
            #1 == -#2,        -IdentityMatrix[#2],           (* trace = -m: -I *)
            True,              Missing["Undetermined", #2]   (* needs the representation *)
        ] &, {cl, mult}]] /@ classes
];

(* --- assemble the lines --------------------------------------------------- *)

InvertibleBlocks[bootstrap_, candidate_] := Module[
    {report, mult = bootstrap["modules"]["multiplicities"], decomp, ordering, byModule},
    report = runGap["WriteGroupReport(\"$OUT$\", " <> ToString[candidate["order"]] <>
                    "," <> ToString[candidate["id"]] <> ");"];
    If[report === $Failed, Return[$Failed]];

    ordering = First[candidate["classOrderings"]];
    decomp   = candidate["decomposition"];

    (* Each module's block is the direct sum of the irreps in its decomposition,
       evaluated at the group element. This fixes the blocks only up to a change
       of basis inside each module -- the same frame choice the copy index in
       IshibashiBasis leaves open. *)
    byModule = Function[elt,
        MapThread[Function[{d, m},
            Module[{blocks = Join @@ MapThread[
                ConstantArray[report["matrices"][[#2, elt]], #1] &,
                {d, Range@Length[d]}]},
                If[blocks === {}, {{1}}, SparseArray@Join @@
                    MapIndexed[ArrayFlatten[{{#1}}] &, {BlockDiagonalCompose[blocks]}]]
            ]], {decomp, mult}]];

    <|"group" -> candidate["name"],
      "elementClass" -> report["elementClass"],
      "lines" -> (byModule /@ Range@Length[report["elementClass"]])|>
];

BlockDiagonalCompose[blocks_] := Module[{n = Total[Length /@ blocks], out, k = 0},
    out = ConstantArray[0, {n, n}];
    Do[With[{d = Length[b]},
        out[[k + 1 ;; k + d, k + 1 ;; k + d]] = b; k += d], {b, blocks}];
    out];

End[]
EndPackage[]

(* :!CodeAnalysis::EndBlock:: *)
