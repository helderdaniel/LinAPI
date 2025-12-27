# LinAPI

Virtualization layer to abstract code from parallel/distributed backend.
Suported backends are: Lapack[1] and Scalapack[2].

The same C++ code can be linked with the suitable backend, allowing it to run in symmetric parallel CPU cores or distributed computing nodes.

[1] Anderson, E. and Bai, Z. and Bischof, C. and Blackford, S. and Demmel, J. and Dongarra, J. and Du Croz, J. and Greenbaum, A. and Hammarling, S. and McKenney, A. and Sorensen, D. (1999). LAPACK Users' Guide 3rd Ed.. Society for Industrial and Applied Mathematics,Philadelphia, PA. ISBN = 0-89871-447-8

[2] Blackford, L. S. and Choi, J. and Cleary, A. and D'Azevedo, E. and Demmel, J. and Dhillon, I. and Dongarra, J. and Hammarling, S. and Henry, G. and Petitet, A. and Stanley, K. and Walker, D. and Whaley, R. C. (1997). ScaLAPACK} Users' Guide. Society for Industrial and Applied Mathematics, Philadelphia, PA. ISBN = 0-89871-397-8
