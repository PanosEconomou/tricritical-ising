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
     "Tilde"    : trace with (-1)^F inserted (Even - Odd). Zero for R unless
                  h == c/24, where it is the constant 1 (Witten index).
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
\"Odd\" | \"Tilde\" | \"Both\"; \"Variable\" -> q. SVirCharacter[\"k\" -> k, ...] uses m = k+2.";
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
SVirModularData::usage =
  "SVirModularData[m, opts] gives the modular data of a fermionic super-Virasoro \
minimal model. Options: \"Invariant\" -> \"AA\" (default) | \"DE\" | {\"D\",\"E\"} | \
{Mr, Ms}, where the first letter / matrix Mr is an SU(2)_(m-2) invariant acting on r \
and the second Ms an SU(2)_m invariant acting on s (see SVirSU2Invariant). \
\"Sector\" -> All (default) | \"NS,NS\" | \"NS,R\" (spatial, temporal boundary \
conditions of the fermion). All gives the full object: \"irreps\" (labels \
{\"NS\"|\"NS~\"|\"R\", {r,s}}), \"S\", \"T\", \"Z\" (one coupling matrix per spin \
structure \"NS\", \"NS~\", \"R\"), \"modules\", \"fixed\", \"c\", \"theory\", \
\"algebra\", \"conventions\". \"NS,NS\" and \"NS,R\" give objects ready for \
ModularBootstrapTraces: a single \"S\" and \"Z\" matrix and \"modules\" with \
\"positions\". \"NS,NS\": traces Tr_NS(L q^..), positivity in the NS defect space. \
\"NS,R\": traces Tr_NS(L (-1)^F q^..), positivity in the R defect space, written in \
ungraded R characters. \"Joint\" imposes both on the same traces, with Ramond \
multiplicities in units of 2 (1 when both labels are fixed points); it needs the \
\"units\" key, so ModularBootstrapTraces must divide by it. This is the one to use \
for the bootstrap of a fermionic theory. SVirModularData[\"k\" -> k, opts] uses m = k+2.";
SVirSU2Invariant::usage =
  "SVirSU2Invariant[type, k] gives the A, D or E modular invariant of SU(2)_k as a \
(k+1) x (k+1) matrix in labels 1..k+1 (label = 2j+1). D needs even k >= 4, E needs \
k = 10, 16 or 28.";
SVirSelfTest::usage =
  "SVirSelfTest[nmax:20] checks m = 3 against independent Virasoro M(5,4) \
characters and the m = 10 Ramond fixed point. All entries should be True.";

Begin["`Private`"];

SVirCharacter::badlabel =
  "`1` is not a Kac label for m = `2`: need 1 <= r <= m-1, 1 <= s <= m+1.";
SVirCharacter::badgrading =
  "Unknown grading `1`. Use \"Full\", \"Ungraded\", \"Even\", \"Odd\", \"Tilde\" or \"Both\".";
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
     (* trace with (-1)^F inserted, ground state even *)
     {"NS", "Tilde"}, base Table[(-1)^j, {j, 0, kx}],
     {"R", "Tilde"}, If[fixedQ, delta, 0 base],
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

(* ---------- SU(2)_k modular invariants (Cappelli-Itzykson-Zuber) ---------- *)

SVirSU2Invariant::none =
  "SU(2)_`2` has no `1`-type invariant (D needs even k >= 4, E needs k = 10, 16 or 28).";

SVirSU2Invariant["A", k_Integer?Positive] := IdentityMatrix[k + 1];

SVirSU2Invariant["D", k_Integer /; EvenQ[k] && k >= 4] := Table[
   If[Divisible[k, 4],
     (* D_even: |chi_r + chi_(k+2-r)|^2 over odd r, middle label doubled *)
     Boole[OddQ[r] && OddQ[rp]] (Boole[rp == r] + Boole[rp == k + 2 - r]),
     (* D_odd: diagonal on odd r, r -> k+2-r on even r *)
     If[OddQ[r], Boole[rp == r], Boole[rp == k + 2 - r]]],
   {r, k + 1}, {rp, k + 1}];

SVirSU2Invariant["E", k : (10 | 16 | 28)] := Module[
   {e = UnitVector[k + 1, #] &, blocks, extra = 0},
   blocks = Switch[k,
     10, {e[1] + e[7], e[4] + e[8], e[5] + e[11]},                      (* E6 *)
     16, {e[1] + e[17], e[5] + e[13], e[7] + e[11], e[9]},              (* E7 *)
     28, {e[1] + e[11] + e[19] + e[29], e[7] + e[13] + e[17] + e[23]}]; (* E8 *)
   (* E7 also has (chi_3 + chi_15) chi_9bar + c.c. *)
   If[k == 16, extra = Outer[Times, e[3] + e[15], e[9]]; extra += Transpose[extra]];
   Total[Outer[Times, #, #] & /@ blocks] + extra];

SVirSU2Invariant[t_, k_] := (Message[SVirSU2Invariant::none, t, k]; $Failed);

(* "DE" or {"D","E"} or {Mr, Ms} -> {Mr, Ms}; Mr acts on r (SU(2)_(m-2)), Ms on s (SU(2)_m) *)
invariantMatrices[m_, s_String /; StringLength[s] == 2] := invariantMatrices[m, Characters[s]];
invariantMatrices[m_, {t1_String, t2_String}] :=
  {SVirSU2Invariant[t1, m - 2], SVirSU2Invariant[t2, m]};
invariantMatrices[m_, {mr_?MatrixQ, ms_?MatrixQ}] /;
   Dimensions[mr] == {m - 1, m - 1} && Dimensions[ms] == {m + 1, m + 1} := {mr, ms};
invariantMatrices[_, _] := (Message[SVirModularData::badinv]; $Failed);

invariantName[s_String] := s;
invariantName[{t1_String, t2_String}] := t1 <> t2;
invariantName[_] := "custom";

theoryName[m_, "AA"] := "SM(" <> ToString[m] <> ")";
theoryName[m_, name_] := "SM(" <> ToString[m] <> ")_" <> name;

normSector[All | "All"] := All;
normSector[{a_String, b_String}] := normSector[a <> "," <> b];
normSector[s_String] := With[{t = StringReplace[ToUpperCase[s], " " -> ""]},
   Which[MemberQ[{"NS,NS", "NS,R"}, t], t,
     MemberQ[{"JOINT", "NS,NS+NS,R"}, t], "Joint",
     True, $Failed]];
normSector[_] := $Failed;

(* ---------- modular data ---------- *)
(* Basis: NS_a = Tr_NS q^..., NS~_a = Tr_NS (-1)^F q^... (ground state even),
   R_rho = full graded Ramond irrep ("Full").  With chi(-1/tau) = S.chi(tau):
     NS  -> NS  : 4 Shat
     NS~ -> R   : 2 Sqrt[2] eps_a Shat_{a rho}
     R   -> NS~ : 2 Sqrt[2] d_rho eps_a Shat_{a rho},   d = 2 (h != c/24), 1 (fixed point)
   Shat_{ab} = sin(pi r r'/m) sin(pi s s'/(m+2)) / Sqrt[m(m+2)], and eps_a = +-1 records
   which coset branching function holds the NS ground state.  The NS~/R blocks come
   from the GKO coset; all three were checked against the q-series at m = 3,...,6, 10.

   Modular invariants.  From the coset SU(2)_(m-2) x SU(2)_2 / SU(2)_m with an invariant
   pair (Mr, Ms) and the Majorana spin structures on the SU(2)_2 factor, summing over
   both members of each identification orbit [r,s] = {(r,s), (m-r, m+2-s)}:
     Z_NS[a,b]        = 1/2 Sum_{x in a, y in b} Mr[x1,y1] Ms[x2,y2]       (NS orbits have size 2)
     Z_NS~[a,b]       = (-1)^(2(h_a - h_b)) Z_NS[a,b]                      (= T^T Z_NS Conjugate[T])
     Z_R^ungr[rho,sg] =     Sum_{x in rho, y in sg} Mr[x1,y1] Ms[x2,y2]    (ungraded R characters)
     Z_R^Full         = Z_R^ungr / (d_rho d_sg)
   Z_R agrees with S^T Z_NS~ S for every A/D/E pair at m = 3..18, fixed points included.
   This keeps every coupling matrix exact: no simplifying products of sines. *)

SVirModularData::badinv =
  "\"Invariant\" must be a two-letter string such as \"AA\" or \"DE\", a pair such as \
{\"D\",\"E\"}, or a pair of matrices {Mr, Ms} of sizes m-1 and m+1.";
SVirModularData::badsector =
  "Unknown sector `1`. Use All, \"Joint\", \"NS,NS\" or \"NS,R\". (\"R,NS\" has no positivity: its \
S-image is Tr_NS (-1)^F, with signs. \"R,R\" only sees the Witten index.)";
SVirModularData::notinv =
  "The invariant `1` fails consistency checks: `2`. Returning the data anyway.";

modularCore[m_, {mr_, ms_}] := modularCore[m, {mr, ms}] = Module[
   {c = SVirCentralCharge[m], tab = SVirKacTable[m], ns, rr, n, nr, hN, hR,
    sh, eps, d, orb, coupling, sNS, sNR, zNS, zNSt, zRu, problems = {}},
   ns = Select[tab, #["Sector"] === "NS" &][[All, "Label"]];
   rr = Select[tab, #["Sector"] === "R" &][[All, "Label"]];
   n = Length[ns]; nr = Length[rr];
   hN = SVirWeight[m, #] & /@ ns;
   hR = SVirWeight[m, #] & /@ rr;

   sh[a_, b_] := Sin[Pi a[[1]] b[[1]]/m] Sin[Pi a[[2]] b[[2]]/(m + 2)]/Sqrt[m (m + 2)];
   eps = If[IntegerQ[SVirWeight[m, #] -
         ((#[[1]] - 1) (#[[1]] + 1)/(4 m) - (#[[2]] - 1) (#[[2]] + 1)/(4 (m + 2)))],
       1, -1] & /@ ns;
   d = If[# == c/24, 1, 2] & /@ hR;
   sNS = 4 Outer[sh, ns, ns, 1];                           (* NS  -> NS       *)
   sNR = 2 Sqrt[2] eps Outer[sh, ns, rr, 1];               (* NS~ -> R (Full) *)

   orb[{r_, s_}] := DeleteDuplicates[{{r, s}, {m - r, m + 2 - s}}];
   coupling[a_, b_] := Sum[mr[[x[[1]], y[[1]]]] ms[[x[[2]], y[[2]]]], {x, orb[a]}, {y, orb[b]}];
   zNS  = Table[coupling[a, b], {a, ns}, {b, ns}]/2;
   zNSt = zNS Outer[If[IntegerQ[#1 - #2], 1, -1] &, hN, hN];
   zRu  = Table[coupling[a, b], {a, rr}, {b, rr}];

   (* consistency checks on the invariant *)
   If[zNS[[1, 1]] =!= 1,
     AppendTo[problems, "the vacuum appears " <> ToString[zNS[[1, 1]]] <> " times"]];
   If[! AllTrue[Flatten[{zNS, zRu}], IntegerQ[#] && # >= 0 &],
     AppendTo[problems, "multiplicities are not nonnegative integers"]];
   If[! AllTrue[Flatten @ Table[zNS[[i, j]] == 0 || IntegerQ[2 (hN[[i]] - hN[[j]])],
        {i, n}, {j, n}], TrueQ],
     AppendTo[problems, "spins outside Z/2 (not T^2 invariant)"]];
   If[Max[Abs[N[Transpose[sNS] . zNS . sNS - zNS]]] > 10^-8,
     AppendTo[problems, "Z_NS is not S invariant"]];

   <|"c" -> c, "NS" -> ns, "R" -> rr, "hNS" -> hN, "hR" -> hR, "d" -> d,
     "SNS" -> sNS, "SNR" -> sNR, "ZNS" -> zNS, "ZNSt" -> zNSt, "ZRu" -> zRu,
     "Problems" -> problems|>];

(* one module per nonzero entry of z; "positions" index z, as ModularBootstrapTraces wants *)
modulesOf[z_, labs_, h_, mult_] := Module[{pos, lab},
   pos = Select[Tuples[Range /@ Dimensions[z]], Extract[z, #] != 0 &];
   lab = {labs[[#[[1]]]], labs[[#[[2]]]]} & /@ pos;
   <|"labels"         -> lab,
     "index"          -> AssociationThread[lab -> Range[Length[pos]]],
     "positions"      -> pos,
     "weights"        -> ({h[[#[[1]]]], h[[#[[2]]]]} & /@ pos),
     "multiplicities" -> Extract[mult, pos]|>];

irrepBlock[tag_, labs_, h_] := <|
   "labels"  -> ({tag, #} & /@ labs),
   "index"   -> AssociationThread[({tag, #} & /@ labs) -> Range[Length[labs]]],
   "weights" -> h,
   "sectors" -> ConstantArray[tag, Length[labs]]|>;

(* all three spin structures in one block basis, as before *)
fullData[m_, core_, name_] := Module[
   {ns = core["NS"], rr = core["R"], hN = core["hNS"], d = core["d"], c = core["c"],
    n, nr, labels, weights, z, sMat, ph, tMat, zMats, modules},
   n = Length[ns]; nr = Length[rr];
   labels  = Join[{"NS", #} & /@ ns, {"NS~", #} & /@ ns, {"R", #} & /@ rr];
   weights = Join[hN, hN, core["hR"]];
   z[a_, b_] := ConstantArray[0, {a, b}];

   sMat = ArrayFlatten[{
      {core["SNS"], z[n, n], z[n, nr]},
      {z[n, n], z[n, n], core["SNR"]},
      {z[nr, n], d Transpose[core["SNR"]], z[nr, nr]}}];

   ph = Exp[2 Pi I (# - c/24)] & /@ weights;
   tMat = ArrayFlatten[{
      {z[n, n], DiagonalMatrix[ph[[1 ;; n]]], z[n, nr]},
      {DiagonalMatrix[ph[[1 ;; n]]], z[n, n], z[n, nr]},
      {z[nr, n], z[nr, n], DiagonalMatrix[ph[[2 n + 1 ;;]]]}}];

   (* Z[sigma] = chi . Z[sigma] . Conjugate[chi] *)
   zMats = <|
     "NS"  -> ArrayFlatten[{{core["ZNS"], z[n, n], z[n, nr]},
                {z[n, n], z[n, n], z[n, nr]}, {z[nr, n], z[nr, n], z[nr, nr]}}],
     "NS~" -> ArrayFlatten[{{z[n, n], z[n, n], z[n, nr]},
                {z[n, n], core["ZNSt"], z[n, nr]}, {z[nr, n], z[nr, n], z[nr, nr]}}],
     "R"   -> ArrayFlatten[{{z[n, n], z[n, n], z[n, nr]},
                {z[n, n], z[n, n], z[n, nr]},
                {z[nr, n], z[nr, n], core["ZRu"]/Outer[Times, d, d]}}]|>;

   (* local operators = NS sector; positions index the NS block of the full matrix *)
   modules = modulesOf[core["ZNS"], ns, hN, core["ZNS"]];

   <|
    "irreps" -> <|
       "labels"  -> labels,
       "index"   -> AssociationThread[labels -> Range[Length[labels]]],
       "weights" -> weights,
       "sectors" -> (First /@ labels)|>,
    "S" -> sMat,
    "T" -> tMat,
    "Z" -> zMats,
    "modules"   -> modules,
    "fixed"     -> Select[rr, SVirWeight[m, #] == c/24 &],
    "c"         -> c,
    "theory"    -> theoryName[m, name],
    "algebra"   -> "SVir(" <> ToString[InputForm[c]] <> ")",
    "invariant" -> name,
    "sector"    -> All,
    "conventions" -> "chi_i(-1/tau) = S_ij chi_j(tau), chi_i(tau+1) = T_ij chi_j(tau). \
NS = Tr_NS, NS~ = Tr_NS (-1)^F with the ground state even, R = full graded Ramond \
irrep (\"Full\"). Z[sigma] = chi.Z[sigma].Conjugate[chi]."
   |>];

(* (NS,NS): Tr_NS L q^.. ; S-closed, dual channel is the NS defect Hilbert space *)
nsnsData[m_, core_, name_] := With[{ns = core["NS"], hN = core["hNS"], c = core["c"]},
   <|
    "irreps"    -> irrepBlock["NS", ns, hN],
    "S"         -> core["SNS"],
    "T2"        -> DiagonalMatrix[Exp[4 Pi I (hN - c/24)]],
    "Z"         -> core["ZNS"],
    "modules"   -> modulesOf[core["ZNS"], ns, hN, core["ZNS"]],
    "fixed"     -> {},
    "c"         -> c,
    "theory"    -> theoryName[m, name],
    "algebra"   -> "SVir(" <> ToString[InputForm[c]] <> ")",
    "invariant" -> name,
    "sector"    -> "NS,NS",
    "conventions" -> "Spin structure (space, time) = (NS, NS). chi = NS characters \
(Tr_NS, no (-1)^F); chi(-1/tau) = S.chi(tau) with S = 4 Shat, real symmetric. \
Z_L = chi.Z_L.Conjugate[chi] with Z_L the trace of L on each module; \
Transpose[S].Z_L.S gives the NS defect Hilbert space. T maps NS to NS~, so only T^2 \
acts within this sector: chi(tau+2) = T2.chi(tau)."
   |>];

(* (NS,R): Tr_NS L (-1)^F q^.. ; S maps it to Tr over the R defect Hilbert space.
   S is n x nr, from NS~ characters to UNGRADED R characters, so the identity line
   gives integer dual multiplicities (Full R characters would give halves). *)
nsrData[m_, core_, name_] := Module[
   {ns = core["NS"], rr = core["R"], hN = core["hNS"], d = core["d"], c = core["c"], mods},
   mods = modulesOf[core["ZNSt"], ns, hN, core["ZNS"]];
   mods["parities"] = Extract[Sign[core["ZNSt"]], mods["positions"]];
   <|
    "irreps"     -> irrepBlock["NS~", ns, hN],
    "dualIrreps" -> Append[irrepBlock["R", rr, core["hR"]], "grading" -> "Ungraded"],
    "S"          -> core["SNR"] . DiagonalMatrix[d],
    "Z"          -> core["ZNSt"],
    "dualZ"      -> core["ZRu"],
    "modules"    -> mods,
    "fixed"      -> Select[rr, SVirWeight[m, #] == c/24 &],
    "c"          -> c,
    "theory"     -> theoryName[m, name],
    "algebra"    -> "SVir(" <> ToString[InputForm[c]] <> ")",
    "invariant"  -> name,
    "sector"     -> "NS,R",
    "conventions" -> "Spin structure (space, time) = (NS, R). chi = NS~ characters \
(Tr_NS (-1)^F, ground state even); chi(-1/tau) = S.chi_R(tau) with chi_R the UNGRADED \
Ramond characters, so S is rectangular (NS x R). Z_L = chi.Z_L.Conjugate[chi]; \
Transpose[S].Z_L.S gives the R defect Hilbert space in ungraded characters, and \
Transpose[S].Z.S = dualZ for the identity line. \"modules\" -> \"multiplicities\" are \
the NS multiplicities and \"parities\" the (-1)^F sign of each module's ground state, \
so the identity line has trace multiplicity*parity = Z. Integrality in the ungraded \
basis is a necessary condition only: a non-fixed R x R module of an ordinary line \
comes in units of 2."
   |>];

(* Joint: (NS,NS) and (NS,R) on the SAME traces X, which is what a line must satisfy.
   Block-diagonal S and Z; every module sits at two positions, one per block.
   The NS~ sign (-1)^(2(h_a - h_b)) is moved into S as phases: with
   s = DiagonalMatrix[Exp[-2 Pi I h]] . S_R, ConjugateTranspose[s] . Z . s gives back the
   signed Z~ automatically, so Z itself carries no signs.
   "units" divides the dual multiplicities entrywise: 1 on NS, and on R x R
   2 (ordinary lines: G_0, Gbar_0 act as a Clifford algebra on the ground states),
   or 1 when both labels are Ramond fixed points (G_0 = 0 there). q-type lines are
   therefore NOT allowed by this sector. *)
jointData[m_, core_, name_] := Module[
   {ns = core["NS"], rr = core["R"], hN = core["hNS"], d = core["d"], c = core["c"],
    n, nr, z, sR, sMat, zMat, units, mods},
   n = Length[ns]; nr = Length[rr];
   z[a_, b_] := ConstantArray[0, {a, b}];
   sR = DiagonalMatrix[Exp[-2 Pi I hN]] . core["SNR"] . DiagonalMatrix[d];
   sMat = ArrayFlatten[{{core["SNS"], z[n, nr]}, {z[n, n], sR}}];         (* 2n x (n+nr) *)
   zMat = ArrayFlatten[{{core["ZNS"], z[n, n]}, {z[n, n], core["ZNS"]}}];  (* 2n x 2n *)
   units = ArrayFlatten[{
      {ConstantArray[1, {n, n}], ConstantArray[1, {n, nr}]},
      {ConstantArray[1, {nr, n}], Outer[If[#1 == 1 && #2 == 1, 1, 2] &, d, d]}}];
   mods = modulesOf[core["ZNS"], ns, hN, core["ZNS"]];
   mods["positions"] = ({#, # + n} & /@ mods["positions"]);
   mods["parities"] = Extract[Sign[core["ZNSt"]], First /@ mods["positions"]];
   <|
    "irreps"     -> With[{lab = Join[({"NS", #} & /@ ns), ({"NS~", #} & /@ ns)]},
                      <|"labels" -> lab, "index" -> AssociationThread[lab -> Range[2 n]],
                        "weights" -> Join[hN, hN], "sectors" -> (First /@ lab)|>],
    "dualIrreps" -> <|"labels" -> Join[({"NS", #} & /@ ns), ({"R", #} & /@ rr)],
                      "weights" -> Join[hN, core["hR"]], "grading" -> "R: Ungraded"|>,
    "S"          -> sMat,
    "Z"          -> zMat,
    "units"      -> units,
    "dualZ"      -> ArrayFlatten[{{core["ZNS"], z[n, nr]},
                      {z[nr, n], core["ZRu"]/units[[n + 1 ;;, n + 1 ;;]]}}],
    "modules"    -> mods,
    "fixed"      -> Select[rr, SVirWeight[m, #] == c/24 &],
    "c"          -> c,
    "theory"     -> theoryName[m, name],
    "algebra"    -> "SVir(" <> ToString[InputForm[c]] <> ")",
    "invariant"  -> name,
    "sector"     -> "Joint",
    "conventions" -> "(NS,NS) and (NS,R) imposed on the same traces. Each module has two \
positions (NS block, NS~ block) carrying the same variable. (ConjugateTranspose[S] . Z_L . S) / \
units must be a nonnegative integer matrix: its NS x NS block is the NS defect Hilbert space, \
its R x R block the R defect Hilbert space in ungraded characters divided by 2 (by 1 for \
fixed point x fixed point). Off-diagonal blocks vanish identically. For the identity line \
it equals dualZ."
   |>];

Options[SVirModularData] = {"Invariant" -> "AA", "Sector" -> All};

SVirModularData["k" -> k_Integer, rest___] := SVirModularData[k + 2, rest];

SVirModularData[m_Integer /; m >= 3, opts : OptionsPattern[]] := Module[
   {inv = OptionValue["Invariant"], sec = normSector[OptionValue["Sector"]], mats, core, name},
   If[sec === $Failed,
     Message[SVirModularData::badsector, OptionValue["Sector"]]; Return[$Failed]];
   mats = invariantMatrices[m, inv];
   If[mats === $Failed || MemberQ[mats, $Failed], Return[$Failed]];
   name = invariantName[inv];
   core = modularCore[m, mats];
   If[core["Problems"] =!= {}, Message[SVirModularData::notinv, name, core["Problems"]]];
   Switch[sec,
     All,     fullData[m, core, name],
     "NS,NS", nsnsData[m, core, name],
     "NS,R",  nsrData[m, core, name],
     "Joint", jointData[m, core, name]]];

(* character values at tau = i t in the SVirModularData basis, for checks *)
charVector[m_, t_?NumericQ, nmax_] := Module[{q = Exp[-2. Pi t]},
   With[{a = SVirCharacter[m, Last[#], nmax, "Output" -> "List",
         "Grading" -> If[First[#] === "NS~", "Tilde", "Full"]]},
      a["Coefficients"] . (q^(a["Offset"] + Range[0, Length[a["Coefficients"]] - 1]/2))
     ] & /@ SVirModularData[m]["irreps"]["labels"]];

(* S- and T-invariance of an SU(2)_k invariant, for the self test *)
su2InvariantOK[t_, k_] := Module[{mm = SVirSU2Invariant[t, k], s, h},
   s = Table[N[Sqrt[2/(k + 2)] Sin[Pi a b/(k + 2)]], {a, k + 1}, {b, k + 1}];
   h = Table[(a^2 - 1)/(4 (k + 2)), {a, k + 1}];
   Max[Abs[s . mm . s - mm]] < 10^-10 &&
    AllTrue[Select[Tuples[Range[k + 1], 2], mm[[#[[1]], #[[2]]]] != 0 &],
      IntegerQ[h[[#[[1]]]] - h[[#[[2]]]]] &]];

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
        res["K"] === {{2, 0, 1}} && res["Rows"][[1, "Status"]] === "OK"],
    (* modular data, checked against the q-series themselves *)
    "S chi(t) = chi(1/t) numerically (m = 3, 10)" ->
      AllTrue[{3, 10}, Max[Abs[charVector[#, 1/1.3, 40] -
           N[SVirModularData[#]["S"]] . charVector[#, 1.3, 40]]] < 10^-8 &],
    "S.S = 1 (m = 3, 10)" ->
      AllTrue[{3, 10}, With[{s = N[SVirModularData[#]["S"]]},
          Max[Abs[s . s - IdentityMatrix[Length[s]]]] < 10^-10] &],
    "S^T Z[NS~] S = Z[R],  T^T Z[NS] T* = Z[NS~]  (m = 10)" ->
      With[{md = SVirModularData[10]},
        Max[Abs[N[Transpose[md["S"]] . md["Z"]["NS~"] . md["S"] - md["Z"]["R"]]]] < 10^-10 &&
        Max[Abs[N[Transpose[md["T"]] . md["Z"]["NS"] . Conjugate[md["T"]] - md["Z"]["NS~"]]]] < 10^-10],
    (* modular invariants *)
    "SU(2)_k A/D/E invariants are S and T invariant (k = 4..28)" ->
      AllTrue[{{"A", 8}, {"D", 4}, {"D", 6}, {"D", 8}, {"D", 10}, {"E", 10}, {"E", 16}, {"E", 28}},
        su2InvariantOK @@ # &],
    "AA default unchanged: Z[R] = diag(1/d) (m = 10)" ->
      With[{md = SVirModularData[10]},
        md["Z"]["R"][[51 ;;, 51 ;;]] ===
          DiagonalMatrix[If[SVirWeight[10, #] == SVirCentralCharge[10]/24, 1, 1/2] & /@
            md["irreps"]["labels"][[51 ;;, 2]]]],
    "m = 10 DE (NS,NS): one vacuum, S^T Z S = Z" ->
      With[{md = SVirModularData[10, "Invariant" -> "DE", "Sector" -> "NS,NS"]},
        md["Z"][[1, 1]] === 1 && md["modules"]["weights"][[1]] === {0, 0} &&
        Max[Abs[N[Transpose[md["S"]] . md["Z"] . md["S"] - md["Z"]]]] < 10^-10],
    "m = 10 DE (NS,R): S^T Z S = dualZ, integer and >= 0" ->
      With[{md = SVirModularData[10, "Invariant" -> "DE", "Sector" -> "NS,R"]},
        Max[Abs[N[Transpose[md["S"]] . md["Z"] . md["S"] - md["dualZ"]]]] < 10^-10 &&
        AllTrue[Flatten[md["dualZ"]], IntegerQ[#] && # >= 0 &]],
    "m = 10 DE (full): S^T Z[NS~] S = Z[R],  T^T Z[NS] T* = Z[NS~]" ->
      With[{md = SVirModularData[10, "Invariant" -> "DE"]},
        Max[Abs[N[Transpose[md["S"]] . md["Z"]["NS~"] . md["S"] - md["Z"]["R"]]]] < 10^-10 &&
        Max[Abs[N[Transpose[md["T"]] . md["Z"]["NS"] . Conjugate[md["T"]] - md["Z"]["NS~"]]]] < 10^-10],
    "m = 10 DE (Joint): identity line gives dualZ, integer and >= 0" ->
      With[{md = SVirModularData[10, "Invariant" -> "DE", "Sector" -> "Joint"]},
        Max[Abs[N[(ConjugateTranspose[md["S"]] . md["Z"] . md["S"])/md["units"] - md["dualZ"]]]] < 10^-10 &&
        AllTrue[Flatten[md["dualZ"]], IntegerQ[#] && # >= 0 &]]
   |>];

End[];
EndPackage[];
