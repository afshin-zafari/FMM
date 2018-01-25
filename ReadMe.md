# Compiling the Program

 `OMP_TASKS`	if defined, the OMP task is used for computations.
 `FMM_3D`	if defined, the FMM will be run for 3D boxes.



# Running the Program
```
app N L Q P T Tree ops flags 
```

 `N`    Problem size

 `L`    number of levels

 `Q`    number of equivalent points

 `P`    number of minimum points per box

 `T`    number of threads

 `Tree`  path to Tree input file

 `ops`   path to operator files

 `flags` can contain thtese characters with corresponding meanings

   `f`    run fmm or not

   `t`    use task based 

   `s`    run sequentially

   `S`    use SuperGlue

   `a`	   use 'add' for data access in SuperGlue tasks

   `w`    use 'write' for data access in superGlue tasks

   `l`    first keep tasks in list, then submit them at once

   `F`   computes Far Fields

   `N`   computes Near Fields 

   `M`   use OpenMP tasks

# Implementation Overview
`gemv` is the main computation routine in the algorithm. Different instances of `gemv_<nnnn>` are created for computations for leaves, upward (child-to-parent) and downward (parent-to-children) the tree and for translation operations within levels. Corresponding names appear in trace output of SuperGlue.
