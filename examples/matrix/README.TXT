LinApi can be compiled for 2 paradigms:

Shared memory (threaded): Lapack
Distributed memory: Scalapack

for each paradigm static and dynamic (shared) linkage libraries are available.

This demo program can be compiled with make and 
executed with run script, for each paradigm.

Individual compilation or execution can be specified in make and run script by:

lapst linApiEx01Lapackst (Shared memory LAPACK dynamic linkage)
lapsh linApiEx01Lapacksh (Shared memory LAPACK dynamic linkage)

scast linApiEx01Scalapackst (distributed memory SCALAPACK dynamic linkage)
scash linApiEx01Scalapacksh (distributed memory SCALAPACK dynamic linkage)

