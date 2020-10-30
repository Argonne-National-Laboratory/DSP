\ENCODING=ISO-8859-1
\Problem name: OsiDefaultName

Minimize
 obj1: 79.33333254 x1 + 69.9999993 x2 - 56.6666661 x3 - 49.9999995 x4
       - 11.99999988 x5 - 3.3333333 x6
Subject To
 c1: x1 - x3 >= 200
 c2: x2 - x4 >= 240
 c3: x5 + x6 <= 0
 q1: [ x3 ^2 ] <= 1000000
Bounds
 0 <= x5 <= 6000
End
