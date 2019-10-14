# NearRepeat

Python code to perform analysis of space-time clustering within event data, primarily for application to crime data. Currently limited to an implementation of the Knox test for spatio-temporal association, as implemented in Jerry Ratcliffe's Near Repeat Calculator.

A number of [notebooks](notebooks) are provided to demonstrate the implementation:

- [Usage of Knox test](notebooks/Usage of Knox test.ipynb) provides a sample application of the test
- [Inconsistency with Near Repeat Calculator](notebooks/Inconsistency with Near Repeat Calculator.ipynb) examines the inconsistency between results found using this implementation and those obtained using the Near Repeat Calculator

In parallel with the development of the code here, Wouter Steenbeek has developed an alternative implementation in R, which can be found [here](https://github.com/wsteenbeek/NearRepeat). We have collaborated on some of the accompanying materials and reviewed each other's code, and as far as we are aware the two implementations should give identical results.