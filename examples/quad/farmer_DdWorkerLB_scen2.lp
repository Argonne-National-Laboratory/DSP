\ENCODING=ISO-8859-1
\Problem name: OsiDefaultName

Minimize
 obj1: 50.000001 x1 + 76.6666682 x2 + 86.6666684 x3 + 79.33333492 x4
       + 70.0000014 x5 - 56.6666678 x6 - 50.000001 x7 - 12.00000024 x8
       - 3.3333334 x9 + 0 x10
Subject To
 c1: x1 + x2 + x3 <= 500.5
 c2: 2 x1 + x4 - x6 >= 200
 c3: 2.4 x2 + x5 - x7 >= 240
 c4: - 16 x3 + x8 + x9 <= 0
Bounds
 0 <= x8 <= 6000
      x10 Free
End
