# (line wrapping off: GAP otherwise inserts backslash continuations)
SetPrintFormattingStatus("*stdout*", false);
##############################################################################
# candidates.g -- which finite groups are compatible with a set of module
# characters coming from the g = 1 families of a modular bootstrap.
#
# The invertible topological lines of a theory form a finite group G. The
# g = 1 families are its CONJUGACY CLASSES, not its elements, and the trace
# vector of a class restricted to a module of multiplicity m is the character
# of an m-dimensional representation of G.
#
# So: the number of g = 1 families is the number of classes, every module
# multiplicity bounds an irrep degree, and every module character must
# decompose into Irr(G) with non-negative integer multiplicities. This
# searches SmallGroup for everything satisfying all of that.
#
# Written by Mathematica; output is Mathematica syntax on stdout.
##############################################################################

MMANumber := function(x)
  if IsInt(x) then return String(x); fi;
  if IsRat(x) then return Concatenation(String(NumeratorRat(x)), "/",
                                        String(DenominatorRat(x))); fi;
  # Cyclotomics need a real translation; flag rather than emit something wrong.
  return Concatenation("$Failed(*", String(x), "*)");
end;

MMAList := function(l, f)
  local s, i;
  s := "{";
  for i in [1..Length(l)] do
    if i > 1 then Append(s, ","); fi;
    Append(s, f(l[i]));
  od;
  Append(s, "}");
  return s;
end;

MMAMatrix := function(m) return MMAList(m, r -> MMAList(r, MMANumber)); end;

# chars: module characters, one list per module, indexed by OUR class order
#        with the identity class first.
CandidateGroups := function(nclasses, maxdeg, chars, maxOrder)
  local n, i, g, cc, irr, idpos, rest, p, ok, m, vec, sol, out, dims, ordering,
        decomp, reps, mats, e;
  out := [];
  for n in [nclasses..maxOrder] do
    if NrSmallGroups(n) = 0 then continue; fi;
    for i in [1..NrSmallGroups(n)] do
      g   := SmallGroup(n, i);
      cc  := ConjugacyClasses(g);
      if Length(cc) <> nclasses then continue; fi;
      irr  := Irr(g);
      dims := List(irr, Degree);
      if Maximum(dims) > maxdeg then continue; fi;
      idpos := PositionProperty(cc, c -> Representative(c) = One(g));
      rest  := Difference([1..nclasses], [idpos]);
      for p in PermutationsList(rest) do
        ordering := Concatenation([idpos], p);   # our class j  ->  GAP class ordering[j]
        ok := true; decomp := [];
        for m in chars do
          vec := List([1..nclasses], j -> m[Position(ordering, j)]);
          sol := SolutionMat(List(irr, x -> List([1..nclasses], j -> x[j])), vec);
          if sol = fail or ForAny(sol, x -> not IsInt(x) or x < 0) then
            ok := false; break;
          fi;
          Add(decomp, sol);
        od;
        if ok then
          Add(out, rec(order := n, id := i, name := StructureDescription(g),
                       classOrder := ordering, degrees := dims, decomp := decomp));
        fi;
      od;
    od;
  od;
  return out;
end;

# Explicit matrices for every irrep, evaluated on every group element, together
# with which class each element belongs to. This is what lets Mathematica
# assemble the blocks: a module's block is the direct sum of the irreps in its
# decomposition.
ReportGroup := function(n, id)
  local g, cc, reps, elts, r;
  g    := SmallGroup(n, id);
  cc   := ConjugacyClasses(g);
  elts := Elements(g);
  reps := IrreducibleRepresentations(g);
  Print("<|\"order\"->", n, ",\"id\"->", id,
        ",\"name\"->\"", StructureDescription(g), "\"",
        ",\"elementClass\"->", MMAList(elts, e ->
            String(PositionProperty(cc, c -> e in c))),
        ",\"degrees\"->", MMAList(reps, r ->
            String(DimensionOfMatrixGroup(Image(r)))),
        ",\"matrices\"->", MMAList(reps, r ->
            MMAList(elts, e -> MMAMatrix(Image(r, e)))),
        "|>");
end;
    

##############################################################################
# Writers. All the quoting lives here rather than in generated Mathematica
# strings, and each opens its own stream with formatting off -- GAP otherwise
# wraps at screen width and inserts backslash continuations mid-token, which
# is a syntax error on the Mathematica side.
#
# Every name assigned inside a GAP function must appear in its `local` list.
# A missing one is only a read-time WARNING, and the definition then silently
# fails to land -- you find out at call time with "must have a value".
##############################################################################

MMAOpen := function(path)
  local o;
  o := OutputTextFile(path, false);
  SetPrintFormattingStatus(o, false);
  return o;
end;

WriteCandidates := function(path, res)
  local o, i;
  o := MMAOpen(path);
  PrintTo(o, "{");
  for i in [1..Length(res)] do
    if i > 1 then AppendTo(o, ","); fi;
    AppendTo(o, "<|\"name\"->\"", res[i].name, "\",\"order\"->", res[i].order,
             ",\"id\"->", res[i].id,
             ",\"classOrder\"->", MMAList(res[i].classOrder, String),
             ",\"degrees\"->", MMAList(res[i].degrees, String),
             ",\"decomposition\"->", MMAList(res[i].decomp, d -> MMAList(d, MMANumber)),
             "|>");
  od;
  AppendTo(o, "}");
  CloseStream(o);
end;

WriteGroupReport := function(path, n, id)
  local o, g, cc, reps, elts;
  g    := SmallGroup(n, id);
  cc   := ConjugacyClasses(g);
  elts := Elements(g);
  reps := IrreducibleRepresentations(g);
  o    := MMAOpen(path);
  PrintTo(o, "<|\"order\"->", n, ",\"id\"->", id,
             ",\"name\"->\"", StructureDescription(g), "\"",
             ",\"elementClass\"->", MMAList(elts, e ->
                 String(PositionProperty(cc, c -> e in c))),
             ",\"degrees\"->", MMAList(reps, r ->
                 String(DimensionOfMatrixGroup(Image(r)))),
             ",\"matrices\"->", MMAList(reps, r ->
                 MMAList(elts, e -> MMAMatrix(Image(r, e)))),
             "|>");
  CloseStream(o);
end;
