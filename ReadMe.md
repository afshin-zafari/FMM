# Compiling the Program

 `OMP_TASKS`	if defined, the OMP task is used for computations.


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
