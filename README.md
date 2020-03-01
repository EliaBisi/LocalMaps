# LocalMaps
Mathematica package accompanying the paper "The geometric Burge correspondence and the partition function of polymer replicas" by Elia Bisi, Neil O'Connell, and Nikos Zygouras - arXiv:2001.09145

The very last equation of the appendix of our paper defines a map C
(denoted by an upper-case calligraphic C) as a composition of several "local maps".
The present package aims at verifying that such a map equals the identity.
For the notation, we refer to the appendix of the paper.

To do so, it suffices to run, for example, the following commands:

```
m = 2;
q = 3;
matIn = Table[Subscript[w,i,j],{i,5},{j,4}];
matOut = CMap[mat,m,q];
```

The parameter m can be chosen to be greater than or equal to 1.
The parameter q can be chosen to be greater than or equal to 3.
The input symbolic matrix matIn needs to be at least of size (m+q)x(q+1)
(in the example, we chose a 5x4 matrix).
The verification consists in checking that the output matrix matOut equals the input matrix mat.
The choice of the parameters (m, q, and the dimensions of mat) is arbitrary,
in the sense that verifying the identity for a fixed given choice of parameters
is equivalent to verifying the identity for any possible choice.
