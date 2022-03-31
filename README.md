# EllieOctree
EllieOctree is a C++ header file (plus tests file) that builds an octree when given objects and a coordinate function.

## Requirements
This code requires C++20 to compile and run.

## Usage
See tests.cpp for examples on using EllieOctree.hpp.

## Other notes
This octree has a number of unique features:
- The octree will prioritise functions depending on whether it can move and/or copy data. This allows, e.g. ```unique_ptr``` to be added as an object into the octree.
  - To achieve this, concepts from C++20 have been utilised.
- A node contains "center of data" variables, for quick access in, e.g. N-body gravity sims.
- The octree has an "updateTree()" function.
  - This function is quicker than rebuilding for when a couple of objects move.
  - Retrieving the objects and destroying and rebuilding the octree is quicker for when all particles move.
  - See tests.cpp for examples

## Future planned updates
- Generalise Node into a template class to allow for custom Node declarations.

## License
[MIT](https://choosealicense.com/licenses/mit/)
