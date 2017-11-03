basename = "ch4_save_density";
cell = Import[basename <> ".h5", {"Datasets", "cell"}]
Dimensions[data = Import[basename <> ".h5", {"Datasets", "data"}]]
dataRot = 
  RotateRight[data, Dimensions[data]/2]; (*center at zero*)
p = 
 ListContourPlot3D[dataRot, 
  BoxRatios -> {cell[[1, 1]], cell[[2, 2]], cell[[3, 3]]},
  DataRange -> {{0, cell[[1, 1]]}, {0, cell[[2, 2]]}, {0, 
     cell[[3, 3]]}}, Mesh -> False,
  (*RegionFunction->Function[{x,y,z},y>cell[[2,2]]/2 && z>cell[[3,3]]/
  2 ],*)
  BoundaryStyle -> Red,
  Contours -> {0.007, 0.028},
  AxesLabel -> {x, y, z},
  ContourStyle -> Directive[Blue, Specularity[White, 20], Opacity[0.5]]
  ]
Export[basename <> ".png", p]
