// Gmsh project created on Fri Aug 14 17:16:09 2026
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Point(5) = {2, 0, 0, 0};
//+
Point(6) = {2, 1, 0, 0};
//+
Line(5) = {3, 6};
//+
Line(6) = {6, 5};
//+
Line(7) = {5, 2};
//+
Curve Loop(2) = {7, 2, 5, 6};
//+
Plane Surface(2) = {2};
