# PhD Thesis

This repository contains all the programming tools used during my PhD Thesis. The files include:
- Mathematica file to generate the functional equations.
- Folder with C++ code to solve Dyson-Schwinger and Bethe-Salpeter equations.
- Folder with C++ code to solve flow equations.
- Folder with numerical routines required during the run.


## Mathematica file FunExp.nb

- Use of external packages:

  - DoFun installed and loaded. Related functions used are described.
  - FormTracer installed and loaded. Related definition of tensors and groups described and implemented.
  
- Template for the Dyson-Schwinger equations:
 
  - Definition of fields and action.
  - Definition of vertices and projections.
  - Derivation of the DSEs in a symbolic and algebraic form.
  - Exportation of the DSE equation to be used in the C++ codes. 
  
- Template for the Renormalisation Group equations:
 
  - Definition of fields and action.
  - Definition of vertices, regulators and projections.
  - Derivation of the RGEs in a symbolic and algebraic form.
  - Exportation performed as in the DSE section.
  
