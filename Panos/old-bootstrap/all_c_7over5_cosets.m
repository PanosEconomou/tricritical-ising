(* ::Package:: *)

(* ::Package:: *)

(* Conformal embeddings g_k / (subalgebra factors), converted from the
   original Python/NumPy dump into Wolfram Language.

   conformalEmbeddings is a List (NOT an Association) of 57 associations.
   A List is required because three labels occur twice with different
   projection matrices, i.e. inequivalent embeddings with the same name:
     B7_1 / (A1_10 x A1_3 x A1_3)
     D8_1 / (A1_2 x A1_2 x A1_3 x A1_3)
     E8_1 / (A1_2 x A1_2 x A1_3 x A1_3)
   Keying by label would silently discard one of each pair.

   Each entry:
     "Label"       canonical string, e.g. "A7_2 / (D4_4)"
     "Ambient"     <|"Type" -> "A7", "Level" -> 2|>
     "U1Factors"   number of u(1) directions (these carry no projection matrix)
     "Subalgebras" list of <|"Type", "Level", "Projection"|>

   "Projection" is an integer matrix with Rank[Type] rows and
   Rank[Ambient] columns; every entry in the file satisfies this.

   Usage:
     Select[conformalEmbeddings, #Ambient["Type"] == "E8" &]
     Cases[conformalEmbeddings, KeyValuePattern["U1Factors" -> 0]]
     Lookup[#, "Projection"] & /@ conformalEmbeddings[[1, "Subalgebras"]]
*)

conformalEmbeddings = {
  <|
    "Label" -> "A7_2 / (D4_4)",
    "Ambient" -> <|"Type" -> "A7", "Level" -> 2|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "D4", "Level" -> 4, "Projection" -> {
        {1, 0, 0, 0, 0, 0, 1},
        {0, 1, 0, 0, 0, 1, 0},
        {0, 0, 1, 0, 1, 0, 0},
        {0, 0, 1, 2, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "B6_1 / (A1_2 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "B6", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 2, 3, 4, 4, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 1, 0, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B7_1 / (A1_10 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "B7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 10, "Projection" -> {2, 4, 6, 8, 10, 10, 5}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 2, 2, 1, 0, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 0, 0, -1, 0, 1, 0}|>
    }
  |>,
  <|
    "Label" -> "B7_1 / (A1_10 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "B7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 10, "Projection" -> {4, 6, 6, 6, 6, 6, 3}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 2, 3, 4, 4, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 1, 0, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B7_1 / (B2_1 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "B7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "B2", "Level" -> 1, "Projection" -> {
        {1, 0, 0, 0, 0, 0, 0},
        {0, 2, 2, 2, 2, 2, 1}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 2, 3, 4, 4, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 1, 0, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x A1_1 x A1_1 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 0, 0, 0, 0, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 4, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x A1_28 x A1_28)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 6, 10, 12, 12, 12, 12, 6}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 0, 6, 10, 12, 6}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x A1_28 x A1_3 x A1_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 6, 10, 12, 12, 12, 12, 6}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 1, 2, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x A1_28 x G2_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 6, 10, 12, 12, 12, 12, 6}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 0, 1, 0, 2, 1},
        {0, 0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x A1_3 x A1_1 x A1_28)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 0, 6, 10, 12, 6}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x A1_3 x A1_3 x A1_1 x A1_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 1, 2, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x A1_3 x A1_1 x G2_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 0, 1, 0, 2, 1},
        {0, 0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x G2_1 x A1_28)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 1, 0, 2, 2, 2, 2, 1},
        {0, 0, 1, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 0, 6, 10, 12, 6}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x G2_1 x A1_3 x A1_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 1, 0, 2, 2, 2, 2, 1},
        {0, 0, 1, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 1, 2, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_2 x G2_1 x G2_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 1, 0, 2, 2, 2, 2, 1},
        {0, 0, 1, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 0, 1, 0, 2, 1},
        {0, 0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_1 x A1_1 x A1_2 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 0, 0, 0, 0, 0, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 4, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_28 x A1_28 x A1_2)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {6, 10, 12, 12, 12, 12, 12, 6}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 0, 6, 10, 12, 6}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_28 x A1_2 x A1_3 x A1_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {6, 10, 12, 12, 12, 12, 12, 6}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 1, 2, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_28 x A1_2 x G2_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {6, 10, 12, 12, 12, 12, 12, 6}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 0, 1, 0, 2, 1},
        {0, 0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_3 x A1_1 x A1_2 x A1_28)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 0, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 0, 6, 10, 12, 6}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_3 x A1_3 x A1_1 x A1_1 x A1_2)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 0, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 0, 1, 2, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (A1_3 x A1_1 x A1_2 x G2_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 0, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 2, 2, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 0, 1, 0, 2, 1},
        {0, 0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (G2_1 x A1_2 x A1_28)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {1, 0, 2, 2, 2, 2, 2, 1},
        {0, 1, 0, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 0, 6, 10, 12, 6}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (G2_1 x A1_2 x A1_3 x A1_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {1, 0, 2, 2, 2, 2, 2, 1},
        {0, 1, 0, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 1, 2, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (G2_1 x A1_2 x G2_1)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {1, 0, 2, 2, 2, 2, 2, 1},
        {0, 1, 0, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 0, 0, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 0, 1, 0, 2, 1},
        {0, 0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "B8_1 / (B3_1 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "B8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "B3", "Level" -> 1, "Projection" -> {
        {1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 2, 2, 2, 2, 2, 1}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 4, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 2, 1}|>
    }
  |>,
  <|
    "Label" -> "C2_3 / (A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "C2", "Level" -> 3|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1}|>
    }
  |>,
  <|
    "Label" -> "C8_1 / (D4_4 x U(1))",
    "Ambient" -> <|"Type" -> "C8", "Level" -> 1|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "D4", "Level" -> 4, "Projection" -> {
        {1, 0, 0, 0, 0, 0, 1, 0},
        {0, 1, 0, 0, 0, 1, 0, 0},
        {0, 0, 1, 0, 1, 0, 0, 0},
        {0, 0, 1, 2, 1, 0, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "D5_1 / (A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "D5", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {2, 3, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 0, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (A1_1 x A1_1 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 0, 0, 0, 0, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 2, 3, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 1, 0, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (A1_28 x A1_28)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {6, 10, 12, 12, 12, 6, 6}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 6, 10, 6, 6}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (A1_28 x A1_3 x A1_1)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {6, 10, 12, 12, 12, 6, 6}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 1, 0, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 1, 2, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (A1_28 x G2_1)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {6, 10, 12, 12, 12, 6, 6}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 1, 0, 1, 1},
        {0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (A1_3 x A1_1 x A1_28)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 0, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 6, 10, 6, 6}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (A1_3 x A1_3 x A1_1 x A1_1)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 0, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 1, 0, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 1, 2, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (A1_3 x A1_1 x G2_1)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {1, 0, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {1, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 1, 0, 1, 1},
        {0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (G2_1 x A1_28)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {1, 0, 2, 2, 2, 1, 1},
        {0, 1, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 28, "Projection" -> {0, 0, 0, 6, 10, 6, 6}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (G2_1 x A1_3 x A1_1)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {1, 0, 2, 2, 2, 1, 1},
        {0, 1, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 1, 0, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 1, 2, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D7_1 / (G2_1 x G2_1)",
    "Ambient" -> <|"Type" -> "D7", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {1, 0, 2, 2, 2, 1, 1},
        {0, 1, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "G2", "Level" -> 1, "Projection" -> {
        {0, 0, 0, 1, 0, 1, 1},
        {0, 0, 0, 0, 1, 0, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "D8_2 / (A7_2 x U(1))",
    "Ambient" -> <|"Type" -> "D8", "Level" -> 2|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "A7", "Level" -> 2, "Projection" -> {
        {1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0}
      }|>
    }
  |>,
  <|
    "Label" -> "D8_1 / (A1_1 x A1_3 x A1_3 x U(1)^2)",
    "Ambient" -> <|"Type" -> "D8", "Level" -> 1|>,
    "U1Factors" -> 2,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D8_1 / (A2_1 x A1_3 x A1_3 x U(1))",
    "Ambient" -> <|"Type" -> "D8", "Level" -> 1|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "A2", "Level" -> 1, "Projection" -> {
        {0, 1, 2, 2, 2, 2, 1, 1},
        {1, 0, 0, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D8_1 / (A1_2 x A1_2 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "D8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 2, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D8_1 / (A1_1 x A1_1 x A1_3 x A1_3 x U(1))",
    "Ambient" -> <|"Type" -> "D8", "Level" -> 1|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 0, 0, 0, 0, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D8_1 / (A3_1 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "D8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A3", "Level" -> 1, "Projection" -> {
        {0, 1, 2, 2, 2, 2, 1, 1},
        {1, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 0, 0, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 2, 3, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 0, 1, 0, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "D8_1 / (A1_2 x A1_2 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "D8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {2, 2, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 2, 2, 2, 2, 2, 1, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 2, 3, 4, 4, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 0, 0, 1, 0, 2, 1, 1}|>
    }
  |>,
  <|
    "Label" -> "E6_1 / (A1_3 x A1_3 x U(1))",
    "Ambient" -> <|"Type" -> "E6", "Level" -> 1|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -2, -2, -4, -3, -2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -1, -1, 0, -1, 0}|>
    }
  |>,
  <|
    "Label" -> "E7_1 / (A1_3 x A1_3 x U(1)^2)",
    "Ambient" -> <|"Type" -> "E7", "Level" -> 1|>,
    "U1Factors" -> 2,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -2, -2, -4, -3, -2, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -1, -1, 0, -1, 0, 0}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A1_3 x A1_3 x A1_1 x U(1)^2)",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 2,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -2, -2, -4, -3, -2, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -1, -1, 0, -1, 0, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {-2, -3, -4, -6, -5, -4, -3, -2}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A1_3 x A1_3 x A2_1 x U(1))",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -2, -2, -4, -3, -2, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, -1, -1, 0, -1, 0, 0, 0}|>,
      <|"Type" -> "A2", "Level" -> 1, "Projection" -> {
        {-2, -3, -4, -6, -5, -4, -3, -2},
        {0, 0, 0, 0, 0, 0, 0, 1}
      }|>
    }
  |>,
  <|
    "Label" -> "E8_2 / (A7_2 x A1_2)",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 2|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A7", "Level" -> 2, "Projection" -> {
        {-2, -2, -3, -4, -3, -2, -1, 0},
        {1, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 1, 0}
      }|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {-2, -3, -4, -6, -5, -4, -3, -2}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A1_1 x A1_3 x A1_3 x U(1)^2)",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 2,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 1, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 2, 2, 4, 3, 2, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 1, 0, 1, 0, 0, 0}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A2_1 x A1_3 x A1_3 x U(1))",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "A2", "Level" -> 1, "Projection" -> {
        {0, 1, 1, 2, 2, 2, 2, 1},
        {-2, -3, -4, -6, -5, -4, -3, -2}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 2, 2, 4, 3, 2, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 1, 0, 1, 0, 0, 0}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A1_2 x A1_2 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {-4, -5, -7, -10, -8, -6, -4, -2}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 1, 1, 2, 2, 2, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 2, 2, 4, 3, 2, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 1, 0, 1, 0, 0, 0}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A1_1 x A1_1 x A1_3 x A1_3 x U(1))",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 1,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 1, 1, 2, 2, 2, 2, 1}|>,
      <|"Type" -> "A1", "Level" -> 1, "Projection" -> {0, 0, 0, 0, 0, 0, 0, 1}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 2, 2, 4, 3, 2, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 1, 0, 1, 0, 0, 0}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A3_1 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A3", "Level" -> 1, "Projection" -> {
        {0, 1, 1, 2, 2, 2, 2, 1},
        {-2, -3, -4, -6, -5, -4, -3, -2},
        {0, 0, 0, 0, 0, 0, 0, 1}
      }|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 2, 2, 4, 3, 2, 0, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 1, 0, 1, 0, 0, 0}|>
    }
  |>,
  <|
    "Label" -> "E8_1 / (A1_2 x A1_2 x A1_3 x A1_3)",
    "Ambient" -> <|"Type" -> "E8", "Level" -> 1|>,
    "U1Factors" -> 0,
    "Subalgebras" -> {
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {-4, -5, -7, -10, -8, -6, -4, -2}|>,
      <|"Type" -> "A1", "Level" -> 2, "Projection" -> {0, 1, 1, 2, 2, 2, 2, 2}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 2, 2, 4, 4, 3, 2, 0}|>,
      <|"Type" -> "A1", "Level" -> 3, "Projection" -> {0, 1, 1, 2, 0, 1, 0, 0}|>
    }
  |>
};
