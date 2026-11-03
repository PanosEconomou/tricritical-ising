(* ::Package:: *)
(* SVirCharacters.wl
   Characters of the unitary N=1 super-Virasoro minimal models
     c = 3/2 (1 - 8/(m(m+2))),  m = 3, 4, 5, ...   (k = m - 2 in the GKO coset
     SU(2)_2 x SU(2)_k / SU(2)_{k+2}).

   Formula (Goddard-Kent-Olive; Meurman-Rocha-Caridi):
     chi_{r,s} = q^{h-c/24} P(q) Sum_k [ q^{a_k(r,s)-h0} - q^{a_k(r,-s)-h0} ],
     a_k(r,s) = ((2m(m+2)k + (m+2)r - m s)^2 - 4) / (8 m (m+2)),
     P_NS = Prod (1+q^(n-1/2))/(1-q^n),   P_R = Prod (1+q^n)/(1-q^n).
   NS: r-s even. R: r-s odd. Identification (r,s) ~ (m-r, m+2-s).

   Internally everything is an integer coefficient list in x = q^(1/2),
   built from PartitionsP / PartitionsQ and ListConvolve, so it is fast.

   Grading conventions (option "Grading"):
     "Full"     : the whole irrep including (-1)^F.  For R with h != c/24 the
                  ground space is 2-dim (G_0 pairs the two parities), so
                  Full = 2 x Ungraded.  This matches Z_R = 2|chi_{7/16}|^2 + ...
     "Ungraded" : the irrep of SVir alone, no (-1)^F (G_0 acts as a scalar on
                  the R ground state).  Same as Full in NS.
     "Even" / "Odd" : (-1)^F = +1 / -1 parts, ground state taken even.
                  NS: integer / half-integer levels.
                  R, h != c/24: both equal Ungraded.
                  R, h == c/24 (fixed point): (Ungraded +- q^0 term)/2.
     "Both"     : {Even, Odd}.

   Run SVirSelfTest[] once: it checks m = 3 against the Virasoro M(5,4)
   (tricritical Ising) characters computed independently.
*)

BeginPackage["SVirCharacters`"];

SVirCentralCharge::usage =
  "SVirCentralCharge[m] gives c = 3/2 (1 - 8/(m(m+2))).";
SVirWeight::usage =
  "SVirWeight[m, {r,s}] gives the conformal weight of the irrep (r,s).";
SVirKacTable::usage =
  "SVirKacTable[m] lists the irreps (one per identification class) as \
associations with keys \"Label\", \"Sector\", \"h\", sorted by sector then h.";
SVirCharacter::usage =
  "SVirCharacter[m, {r,s}, nmax, opts] gives the character of irrep (r,s) \
through level nmax above its ground state. Options: \"Output\" -> \"Series\" \
(SeriesData, default), \"Expression\" (plain sum), \"List\" (coefficients of \
q^(h-c/24+j/2), j=0..2 nmax); \"Grading\" -> \"Full\" | \"Ungraded\" | \"Even\" | \
\"Odd\" | \"Both\"; \"Variable\" -> q. SVirCharacter[\"k\" -> k, ...] uses m = k+2.";
SVirCharacters::usage =
  "SVirCharacters[m, nmax, opts] gives an Association {r,s} -> character for \
every irrep. Extra option \"Sector\" -> All | \"NS\" | \"R\". \
SVirCharacters[\"k\" -> k, nmax, opts] uses m = k+2.";
SVirProduct::usage =
  "SVirProduct[m, {{r1,s1}, {r2,s2}, ...}, nmax, opts] gives the product of \
the characters of the listed irreps of the same model, exact through level \
nmax above the product's ground state. SVirProduct[{{m1,{r1,s1}}, {m2,{r2,s2}}, \
...}, nmax, opts] allows factors from different models. Options as for \
SVirCharacter; \"Grading\" may be a single value or a list, one per factor \
(\"Both\" is not allowed). SVirProduct[\"k\" -> k, ...] uses m = k+2.";
SVirProductTable::usage =
  "SVirProductTable[m, nmax, opts] gives an Association {{r1,s1},{r2,s2}} -> \
product for every unordered pair of irreps of model m. Extra option \
\"Sector\" -> \"Matched\" (NS x NS and R x R only, default) | All.";
SVirDecompose::usage =
  "SVirDecompose[aList, bList] writes each character in aList as an integer \
combination of the characters in bList, a = K.b. Entries may be \"List\" \
associations (from SVirCharacter/SVirProduct with \"Output\" -> \"List\") or \
SeriesData. Returns an Association with \"K\" (rows: a, columns: b) and \
\"Rows\" (per-row diagnostics: \"Status\", \"Rank\", \"KnownThrough\", \
\"Untested\" = indices of b whose leading power lies beyond what a determines).";
SVirSelfTest::usage =
  "SVirSelfTest[nmax:20] checks m = 3 against independent Virasoro M(5,4) \
characters and the m = 10 Ramond fixed point. All entries should be True.";

Begin["`Private`"];

SVirCharacter::badlabel =
  "`1` is not a Kac label for m = `2`: need 1 <= r <= m-1, 1 <= s <= m+1.";
SVirCharacter::badgrading =
  "Unknown grading `1`. Use \"Full\", \"Ungraded\", \"Even\", \"Odd\" or \"Both\".";
SVirCharacter::badoutput =
  "Unknown output `1`. Use \"Series\", \"Expression\" or \"List\".";

(* ---------- basic data ---------- *)

SVirCentralCharge[m_] := 3/2 (1 - 8/(m (m + 2)));

SVirWeight[m_, {r_, s_}] :=
  (((m + 2) r - m s)^2 - 4)/(8 m (m + 2)) + If[OddQ[r - s], 1/16, 0];

SVirKacTable["k" -> k_Integer] := SVirKacTable[k + 2];
SVirKacTable[m_Integer /; m >= 3] := SVirKacTable[m] = Module[{labels},
   labels = DeleteDuplicates @ Flatten[
      Table[First @ Sort[{{r, s}, {m - r, m + 2 - s}}],
        {r, 1, m - 1}, {s, 1, m + 1}], 1];
   SortBy[
     <|"Label" -> #, "Sector" -> If[EvenQ[#[[1]] - #[[2]]], "NS", "R"],
       "h" -> SVirWeight[m, #]|> & /@ labels,
     {#["Sector"] &, #["h"] &}]];

(* ---------- integer series in x = q^(1/2), length kx+1 ---------- *)

truncConv[a_, b_, n_] := Take[ListConvolve[a, b, {1, -1}, 0], n + 1];

(* Prod_{n>=1} (1 + x^(2n-1)) / (1 - x^(2n)) *)
nsPref[kx_Integer] := nsPref[kx] = Module[
   {d = ConstantArray[0, kx + 1], pe = ConstantArray[0, kx + 1]},
   d[[1]] = 1;
   Do[d = d + Join[ConstantArray[0, j], Take[d, kx + 1 - j]], {j, 1, kx, 2}];
   pe[[1 ;; ;; 2]] = PartitionsP[Range[0, Floor[kx/2]]];
   truncConv[d, pe, kx]];

(* Prod_{n>=1} (1 + x^(2n)) / (1 - x^(2n))  =  (PartitionsQ * PartitionsP)(q) *)
rPref[kx_Integer] := rPref[kx] = Module[
   {n = Floor[kx/2], ce = ConstantArray[0, kx + 1]},
   ce[[1 ;; ;; 2]] =
     truncConv[PartitionsQ[Range[0, n]], PartitionsP[Range[0, n]], n];
   ce];

(* Sum_k [x^(2(a_k(r,s)-a_0)) - x^(2(a_k(r,-s)-a_0))] *)
theta[m_, {r_, s_}, kx_] := Module[
   {mm = m (m + 2), a0, kmax, th = ConstantArray[0, kx + 1], e},
   a0 = ((m + 2) r - m s)^2;
   kmax = Ceiling[Sqrt[kx/mm]] + 3;
   Do[
     e = ((2 mm k + (m + 2) r - m s)^2 - a0)/(4 mm);
     If[0 <= e <= kx, th[[e + 1]] += 1];
     e = ((2 mm k + (m + 2) r + m s)^2 - a0)/(4 mm);
     If[0 <= e <= kx, th[[e + 1]] -= 1],
     {k, -kmax, kmax}];
   th];

(* ---------- output formatting ---------- *)

formatOut[coeffs_, off_, "List", meta_, _] :=
  Join[meta, <|"Offset" -> off, "Step" -> 1/2, "Coefficients" -> coeffs|>];

formatOut[coeffs_, off_, "Expression", _, var_] :=
  coeffs . (var^(off + Range[0, Length[coeffs] - 1]/2));

formatOut[coeffs_, off_, "Series", _, var_] := Module[
   {den = LCM[2, Denominator[off]], step, a, kx = Length[coeffs] - 1, full},
   step = den/2;
   a = off den;
   full = ConstantArray[0, kx step + 1];
   full[[1 ;; ;; step]] = coeffs;
   SeriesData[var, 0, full, a, a + kx step + 1, den]];

formatOut[_, _, other_, _, _] := (Message[SVirCharacter::badoutput, other]; $Failed);

(* ---------- main ---------- *)

Options[SVirCharacter] = {"Output" -> "Series", "Grading" -> "Full",
   "Variable" -> Automatic};

SVirCharacter["k" -> k_Integer, rest___] := SVirCharacter[k + 2, rest];

SVirCharacter[m_Integer /; m >= 3, {r_Integer, s_Integer},
   nmax_Integer?NonNegative, opts : OptionsPattern[]] := Module[
   {kx = 2 nmax, h, c, sector, base, delta, fixedQ, grading, coeffs, var},
   grading = OptionValue["Grading"];
   If[grading === "Both",
     Return[{SVirCharacter[m, {r, s}, nmax, "Grading" -> "Even", opts],
             SVirCharacter[m, {r, s}, nmax, "Grading" -> "Odd", opts]}]];
   If[! (1 <= r <= m - 1 && 1 <= s <= m + 1),
     Message[SVirCharacter::badlabel, {r, s}, m]; Return[$Failed]];
   h = SVirWeight[m, {r, s}];
   c = SVirCentralCharge[m];
   sector = If[EvenQ[r - s], "NS", "R"];
   base = truncConv[If[sector === "NS", nsPref[kx], rPref[kx]],
     theta[m, {r, s}, kx], kx];
   fixedQ = sector === "R" && h == c/24;
   delta = UnitVector[kx + 1, 1];
   coeffs = Switch[{sector, grading},
     {"NS", "Full" | "Ungraded"}, base,
     {"NS", "Even"}, base Table[Boole[EvenQ[j]], {j, 0, kx}],
     {"NS", "Odd"}, base Table[Boole[OddQ[j]], {j, 0, kx}],
     {"R", "Ungraded"}, base,
     {"R", "Full"}, If[fixedQ, base, 2 base],
     {"R", "Even"}, If[fixedQ, (base + delta)/2, base],
     {"R", "Odd"}, If[fixedQ, (base - delta)/2, base],
     _, Message[SVirCharacter::badgrading, grading]; Return[$Failed]];
   var = OptionValue["Variable"] /. Automatic -> Global`q;
   formatOut[coeffs, h - c/24, OptionValue["Output"],
     <|"Label" -> {r, s}, "Sector" -> sector, "h" -> h, "c" -> c,
       "Grading" -> grading|>, var]];

Options[SVirCharacters] = Join[Options[SVirCharacter], {"Sector" -> All}];

SVirCharacters["k" -> k_Integer, rest___] := SVirCharacters[k + 2, rest];

SVirCharacters[m_Integer /; m >= 3, nmax_Integer?NonNegative,
   opts : OptionsPattern[]] := Module[{sec = OptionValue["Sector"], reps},
   reps = If[sec === All, SVirKacTable[m],
     Select[SVirKacTable[m], #["Sector"] === sec &]];
   Association[
     (#["Label"] -> SVirCharacter[m, #["Label"], nmax,
         FilterRules[{opts}, Options[SVirCharacter]]]) & /@ reps]];

(* ---------- products ---------- *)
(* Each factor is exact through x^(2 nmax) above its own ground state, so the
   truncated convolution is exact through x^(2 nmax) above the product's. *)

SVirProduct::badgrading =
  "Grading must be one value or a list with one value per factor, and \"Both\" is not allowed.";

Options[SVirProduct] = Options[SVirCharacter];

SVirProduct["k" -> k_Integer, rest___] := SVirProduct[k + 2, rest];

SVirProduct[m_Integer /; m >= 3, labels : {{_Integer, _Integer} ..},
   nmax_Integer?NonNegative, opts : OptionsPattern[]] :=
  SVirProduct[{m, #} & /@ labels, nmax, opts];

SVirProduct[factors : {{_Integer, {_Integer, _Integer}} ..},
   nmax_Integer?NonNegative, opts : OptionsPattern[]] := Module[
   {kx = 2 nmax, gr = OptionValue["Grading"], grs, data, coeffs, var},
   grs = If[ListQ[gr], gr, ConstantArray[gr, Length[factors]]];
   If[Length[grs] =!= Length[factors] || MemberQ[grs, "Both"],
     Message[SVirProduct::badgrading]; Return[$Failed]];
   data = MapThread[
     SVirCharacter[#1[[1]], #1[[2]], nmax, "Output" -> "List",
        "Grading" -> #2] &, {factors, grs}];
   If[MemberQ[data, $Failed], Return[$Failed]];
   coeffs = Fold[truncConv[#1, #2, kx] &,
     First[data]["Coefficients"], #["Coefficients"] & /@ Rest[data]];
   var = OptionValue["Variable"] /. Automatic -> Global`q;
   formatOut[coeffs, Total[#["Offset"] & /@ data], OptionValue["Output"],
     <|"Factors" -> factors, "Sectors" -> (#["Sector"] & /@ data),
       "h" -> Total[#["h"] & /@ data], "c" -> Total[#["c"] & /@ data],
       "Grading" -> grs|>, var]];

Options[SVirProductTable] = Join[Options[SVirProduct], {"Sector" -> "Matched"}];

SVirProductTable["k" -> k_Integer, rest___] := SVirProductTable[k + 2, rest];

SVirProductTable[m_Integer /; m >= 3, nmax_Integer?NonNegative,
   opts : OptionsPattern[]] := Module[{reps, pairs},
   reps = SVirKacTable[m];
   pairs = Flatten[Table[{reps[[i]], reps[[j]]},
      {i, Length[reps]}, {j, i, Length[reps]}], 1];
   If[OptionValue["Sector"] === "Matched",
     pairs = Select[pairs, #[[1]]["Sector"] === #[[2]]["Sector"] &]];
   Association[
     ({#[[1]]["Label"], #[[2]]["Label"]} ->
         SVirProduct[m, {#[[1]]["Label"], #[[2]]["Label"]}, nmax,
           FilterRules[{opts}, Options[SVirProduct]]]) & /@ pairs]];

(* ---------- decomposition a = K . b ---------- *)
(* Each b starts with a single power q^(h_n - c/24) and the weights are distinct
   within a sector, so the coefficient-matching system is triangular: K is
   unique whenever it exists. Nonnegativity and integrality are therefore
   checks on the answer, not constraints to optimize over (no LP needed). *)

SVirDecompose::input =
  "Entries must be \"List\" associations from SVirCharacter/SVirProduct or SeriesData objects.";
SVirDecompose::cmismatch =
  "The characters have different central charges `1`; the q^(-c/24) offsets will not line up.";

(* {exponent -> coefficient, highest exponent known exactly, c (if recorded)} *)
termsOf[a_Association] := Module[
   {cs = a["Coefficients"], off = a["Offset"], st = a["Step"]},
   {AssociationThread[off + st Range[0, Length[cs] - 1], cs],
    off + st (Length[cs] - 1), Lookup[a, "c", Missing[]]}];
termsOf[s_SeriesData] := Module[
   {cs = s[[3]], nmin = s[[4]], nmax = s[[5]], den = s[[6]]},
   {AssociationThread[(nmin + Range[0, Length[cs] - 1])/den, cs],
    (nmax - 1)/den, Missing[]}];
termsOf[_] := $Failed;

leadingExp[assoc_] := Min[Keys[Select[assoc, # != 0 &]]];

decomposeOne[{ta_, topA_, _}, bt_, bLead_] := Module[
   {nb = Length[bt], cand, effTop, cols, exps, mat, rhs, sol, rank = 0,
    row, status = {}},
   (* only b's whose leading power is visible in a can be tested *)
   cand = Flatten @ Position[bLead, l_ /; l <= topA, {1}];
   effTop = Min[Prepend[bt[[cand, 2]], topA]];
   If[effTop < topA,
     AppendTo[status, "some b known to lower order than a: raise nmax for b"]];
   cols = Flatten @ Position[bLead, l_ /; l <= effTop, {1}];
   exps = Select[Union @@ Join[{Keys[ta]}, (Keys[bt[[#, 1]]] & /@ cols)],
     # <= effTop &];
   rhs = Lookup[ta, exps, 0];
   row = ConstantArray[0, nb];
   If[cols === {},
     If[AnyTrue[rhs, # != 0 &],
       AppendTo[status, "no b can match the leading term"]],
     mat = Table[Lookup[bt[[c, 1]], e, 0], {e, exps}, {c, cols}];
     rank = MatrixRank[mat];
     If[rank < Length[cols],
       AppendTo[status, "b characters dependent at this order: K not unique"]];
     sol = Quiet @ Check[LinearSolve[mat, rhs], $Failed];
     If[sol === $Failed,
       AppendTo[status, "no exact solution: a is not a combination of b"],
       row[[cols]] = sol;
       If[! AllTrue[sol, IntegerQ], AppendTo[status, "non-integer multiplicity"]];
       If[AnyTrue[sol, # < 0 &], AppendTo[status, "negative multiplicity"]]]];
   <|"Row" -> row, "Status" -> If[status === {}, "OK", status],
     "Rank" -> rank, "KnownThrough" -> effTop,
     "Untested" -> Complement[Range[nb], cols]|>];

SVirDecompose[as_List, bs_List] := Module[{at, bt, cs, bLead, rows},
   at = termsOf /@ as;
   bt = termsOf /@ bs;
   If[MemberQ[Join[at, bt], $Failed], Message[SVirDecompose::input]; Return[$Failed]];
   cs = Union @ DeleteMissing[Join[at[[All, 3]], bt[[All, 3]]]];
   If[Length[cs] > 1, Message[SVirDecompose::cmismatch, cs]];
   bLead = leadingExp[#[[1]]] & /@ bt;
   rows = decomposeOne[#, bt, bLead] & /@ at;
   <|"K" -> rows[[All, "Row"]], "Rows" -> rows|>];

(* ---------- self test ---------- *)

(* Rocha-Caridi character of Virasoro M(p,p'), integer levels 0..n *)
virChar[p_, pp_, {r_, s_}, n_] := Module[
   {th = ConstantArray[0, n + 1], a0 = (p r - pp s)^2, kmax, e},
   kmax = Ceiling[Sqrt[n/(p pp)]] + 3;
   Do[
     e = ((2 p pp k + p r - pp s)^2 - a0)/(4 p pp);
     If[0 <= e <= n, th[[e + 1]] += 1];
     e = ((2 p pp k + p r + pp s)^2 - a0)/(4 p pp);
     If[0 <= e <= n, th[[e + 1]] -= 1],
     {k, -kmax, kmax}];
   truncConv[PartitionsP[Range[0, n]], th, n]];

SVirSelfTest[n_Integer : 20] := Module[{vir, sv},
   (* M(5,4) labels: h_{r,s} = ((5r-4s)^2-1)/80 *)
   vir[rs_] := virChar[5, 4, rs, n];
   sv[rs_, g_] := SVirCharacter[3, rs, n, "Output" -> "List",
      "Grading" -> g]["Coefficients"];
   <|
    "NS (1,1) even = Vir h=0"     -> (sv[{1, 1}, "Even"][[1 ;; ;; 2]] === vir[{1, 1}]),
    (* odd part of the vacuum starts at level 3/2, i.e. j = 3 *)
    "NS (1,1) odd  = Vir h=3/2"   -> (sv[{1, 1}, "Odd"][[4 ;; ;; 2]] === Take[vir[{3, 1}], n - 1]),
    "NS (1,3) even = Vir h=1/10"  -> (sv[{1, 3}, "Even"][[1 ;; ;; 2]] === vir[{1, 2}]),
    "NS (1,3) odd  = Vir h=3/5"   -> (sv[{1, 3}, "Odd"][[2 ;; ;; 2]] === Most[vir[{3, 2}]]),
    "R (1,2) ungraded = Vir h=3/80" -> (sv[{1, 2}, "Ungraded"][[1 ;; ;; 2]] === vir[{2, 2}]),
    "R (1,4) ungraded = Vir h=7/16" -> (sv[{1, 4}, "Ungraded"][[1 ;; ;; 2]] === vir[{2, 1}]),
    "m=10: 25 NS + 25 R irreps"   -> (Counts[#["Sector"] & /@ SVirKacTable[10]] ===
                                       <|"NS" -> 25, "R" -> 25|>),
    "m=10 R fixed point (5,6): h = c/24, parity split integral" ->
      (SVirWeight[10, {5, 6}] == SVirCentralCharge[10]/24 &&
       AllTrue[SVirCharacter[10, {5, 6}, n, "Output" -> "List",
          "Grading" -> "Odd"]["Coefficients"], IntegerQ]),
    (* independent path: Mathematica's own SeriesData multiplication *)
    "SVirProduct = SeriesData product (NS x NS, Full)" ->
      (Normal[SVirProduct[3, {{1, 1}, {1, 3}}, n] -
          SVirCharacter[3, {1, 1}, n] SVirCharacter[3, {1, 3}, n]] === 0),
    "SVirProduct = SeriesData product (R x R, mixed grading)" ->
      (Normal[SVirProduct[3, {{1, 2}, {1, 4}}, n, "Grading" -> {"Ungraded", "Full"}] -
          SVirCharacter[3, {1, 2}, n, "Grading" -> "Ungraded"] *
          SVirCharacter[3, {1, 4}, n]] === 0),
    (* synthetic decomposition with a known answer: a = 2 (1,1) + (1,5) at m = 10 *)
    "SVirDecompose recovers a known K" ->
      Module[{bs, a, res},
        bs = SVirCharacter[10, #, n, "Output" -> "List"] & /@ {{1, 1}, {1, 3}, {1, 5}};
        a = 2 SVirCharacter[10, {1, 1}, n] + SVirCharacter[10, {1, 5}, n];
        res = SVirDecompose[{a}, bs];
        res["K"] === {{2, 0, 1}} && res["Rows"][[1, "Status"]] === "OK"]
   |>];

End[];
EndPackage[];
