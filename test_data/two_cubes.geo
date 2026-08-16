// Gmsh project created on Sun Aug 16 15:54:30 2026
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Box(2) = {1, 0, 0, 1, 1, 1};
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Extrude {0, 1, 0} {
  Surface{4}; 
}
//+
Recursive Delete {
  Curve{20}; Surface{9}; 
}
//+
Recursive Delete {
  Surface{11}; 
}
//+
Recursive Delete {
  Surface{11}; 
}
//+
Recursive Delete {
  Volume{2}; 
}
