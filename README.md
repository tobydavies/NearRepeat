# NearRepeat

Python code to perform analysis of space-time clustering within event data, primarily for application to crime data. Currently limited to an implementation of the Knox test for spatio-temporal association, equivalent to that implemented in Jerry Ratcliffe's [Near Repeat Calculator](http://www.cla.temple.edu/center-for-security-and-crime-science/projects/#near-repeat-calculator).

One of the advantages of this implementation is that it is fairly fast - in applications to benchmark crime datasets, runtimes are lower than alternative implementations by some margin. Examples of this are included in the notebook below.

A number of [notebooks](notebooks) are provided to demonstrate the implementation:

- [Usage of Knox test](notebooks/Usage%20of%20Knox%20test.ipynb) provides a sample application of the test
- [Inconsistency with Near Repeat Calculator](notebooks/Inconsistency%20with%20Near%20Repeat%20Calculator.ipynb) examines the inconsistency between results found using this implementation and those obtained using the Near Repeat Calculator
- [Runtime benchmarking](notebooks/Runtime%20benchmarking.ipynb) tests the speed of the implementation, and its scaling with parameter settings

In parallel with the development of the code here, Wouter Steenbeek has developed an alternative implementation in R, which can be found [here](https://github.com/wsteenbeek/NearRepeat). We have collaborated on some of the accompanying materials and reviewed each other's code, and as far as we are aware the two implementations should give identical results.