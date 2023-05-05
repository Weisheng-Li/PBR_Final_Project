# PBR Final Project - Weisheng Li

See the PDF report for more details about implementation and correctness validation. Check the release for images and test scenes in the report.

### Code I wrote and stuffs I added
***********************

Anisotropic volume rendering using SGGX:

* `medium.h, medium.cpp` - Add a HomogeneousAnisotropicMedium medium class and a SGGXPhase phase function class. These two needs to work together.
* `parser` - Parse the new medium class
* `vec.h` - Add a few functions that end with 3x3. They still take 4x4 matrices as input and output, but only process the top-left 3x3 submatrix. Other entries are all set to 0.

Layered material with statistical operator:

* `material.h, material.cpp` - Create new material class called Layered. Add a few functions for sampling and evaluating surface normal from Beckmann distribution
* `layered.h` - Implement FGD class that reads precomputed FGD texture from `data/FGD.bin` and supports interpolated texture lookup
* `parser.cpp` - Parse the new material class
* `onb.h` - Now allows converting from world coordinate to local coordinate using ONB
