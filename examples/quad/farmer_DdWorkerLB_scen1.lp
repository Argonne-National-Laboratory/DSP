\ENCODING=ISO-8859-1
\Problem name: OsiDefaultName

Minimize
 obj1: 49.9999995 x1 + 76.6666659 x2 + 86.6666658 x3 + 79.33333254 x4
       + 69.9999993 x5 - 56.6666661 x6 - 49.9999995 x7 - 11.99999988 x8
       - 3.3333333 x9 + 0 x10
Subject To
 c1: x1 + x2 + x3 <= 500.5
 c2: 2.5 x1 + x4 - x6 >= 200
 c3: 3 x2 + x5 - x7 >= 240
 c4: - 20 x3 + x8 + x9 <= 0
Bounds
 0 <= x8 <= 6000
      x10 Free
End
