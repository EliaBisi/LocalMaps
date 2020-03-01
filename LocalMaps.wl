(* ::Package:: *)

(* LocalMaps: Mathematica package accompanying the paper
"The geometric Burge correspondence and the partition function of polymer replicas"
by Elia Bisi, Neil O'Connell, and Nikos Zygouras - arXiv:2001.09145 *)


(* The very last equation of the appendix of our paper defines a map C
(denoted by an upper-case calligraphic C) as a composition of several "local maps".
The present package aims at verifying that such a map equals the identity.
For the notation, we refer to the appendix of the paper.

To do so, it suffices to run, for example, the following commands:

m = 2;
q = 3;
matIn = Table[Subscript[w,i,j],{i,5},{j,4}];
matOut = CMap[mat,m,q];

The parameter m can be chosen to be greater than or equal to 1.
The parameter q can be chosen to be greater than or equal to 3.
The input symbolic matrix matIn needs to be at least of size (m+q)x(q+1)
(in the example, we chose a 5x4 matrix).
The verification consists in checking that the output matrix matOut equals the input matrix mat.
The choice of the parameters (m, q, and the dimensions of mat) is arbitrary,
in the sense that verifying the identity for a fixed given choice of parameters
is equivalent to verifying the identity for any possible choice. *)



(* Local map a_{i,j} *)
aMap[matrix_,i_,j_] := Module[{x=matrix,a,b1,b2,c1,c2},
a = x[[i,j]];
{b1, b2} = If[i==1, If[j==1, {1/2,1/2}, {x[[i,j-1]],0}], If[j==1, {0,x[[i-1,j]]}, {x[[i,j-1]],x[[i-1,j]]}]];
{c1, c2} = {x[[i+1,j]],x[[i,j+1]]};
x[[i,j]] = (b1+b2) * a^(-1) * (c1^(-1)+c2^(-1))^(-1);
x];


(* Local map tilde-a_{i,j} *)
aTildeMap[matrix_,i_,j_] := Module[{x=matrix,a,b,c1,c2},
a = x[[i,j]];
b = If[j==1, 1, x[[i,j-1]]];
{c1, c2} = {x[[i+1,j]],x[[i,j+1]]};
x[[i,j]] = b * a^(-1) * (c1^(-1)+c2^(-1))^(-1);
x];


(* Local map b_{i,j} *)
bMap[matrix_,i_,j_] := Module[{x=matrix,a,b1,b2,c},
a = x[[i,j]];
{b1, b2} = If[i==1, If[j==1, {1/2,1/2}, {x[[i,j-1]],0}], If[j==1, {0,x[[i-1,j]]}, {x[[i,j-1]],x[[i-1,j]]}]];
c = x[[i,j+1]];
x[[i,j]] = (b1+b2) * a^(-1) * c;
x];


(* Local map c_{i,j} *)
cMap[matrix_,i_,j_] := Module[{x=matrix,a,b1,b2},
a = x[[i,j]];
{b1, b2} = If[i==1, If[j==1, {1/2,1/2}, {x[[i,j-1]],0}], If[j==1, {0,x[[i-1,j]]}, {x[[i,j-1]],x[[i-1,j]]}]];
x[[i,j]] = (b1+b2) * a;
x];


(* Local map d_{i,j}^{k,l} *)
dMap[matrix_,i_,j_, k_,l_] := Module[{x=matrix,a,b1,b2,c1,c2,d},
a = x[[i,j]];
{b1, b2} = If[i==1, If[j==1, {1/2,1/2}, {x[[i,j-1]],0}], If[j==1, {0,x[[i-1,j]]}, {x[[i,j-1]],x[[i-1,j]]}]];
{c1, c2} = {x[[i+1,j]],x[[i,j+1]]};
d = x[[k,l]];
x[[i,j]] = (a^(-1) + (d*(b1+b2))^(-1))^(-1);
x[[k,l]] = (d*(b1+b2)*a^(-2) + a^(-1)) * (c1^(-1)+c2^(-1))^(-1);
x];


(* Local map tilde-d_{i,j}^{k,l} *)
dTildeMap[matrix_,i_,j_, k_,l_] := Module[{x=matrix,a,b,c1,c2,d},
a = x[[i,j]];
b = If[j==1, 1, x[[i,j-1]]];
{c1, c2} = {x[[i+1,j]],x[[i,j+1]]};
d = x[[k,l]];
x[[i,j]] = (a^(-1) + (d*b)^(-1))^(-1);
x[[k,l]] = (d*b*a^(-2) + a^(-1)) * (c1^(-1)+c2^(-1))^(-1);
x];


(* Inverse of local map d_{i,j}^{k,l} *)
dInvMap[matrix_,i_,j_, k_,l_] := Module[{x=matrix,a,b1,b2,c1,c2,d},
a = x[[i,j]];
{b1, b2} = If[i==1, If[j==1, {1/2,1/2}, {x[[i,j-1]],0}], If[j==1, {0,x[[i-1,j]]}, {x[[i,j-1]],x[[i-1,j]]}]];
{c1, c2} = {x[[i+1,j]],x[[i,j+1]]};
d = x[[k,l]];
x[[i,j]] = a + d^(-1) * (c1^(-1) + c2^(-1))^(-1);
x[[k,l]] = (a + a^2 * d * (c1^(-1) + c2^(-1))) / (b1+b2);
x];


(* Inverse of local map tilde-d_{i,j}^{k,l} *)
dTildeInvMap[matrix_,i_,j_, k_,l_] := Module[{x=matrix,a,b,c1,c2,d},
a = x[[i,j]];
b = If[j==1, 1, x[[i,j-1]]];
{c1, c2} = {x[[i+1,j]],x[[i,j+1]]};
d = x[[k,l]];
x[[i,j]] = a + d^(-1) * (c1^(-1) + c2^(-1))^(-1);
x[[k,l]] = (a + a^2 * d * (c1^(-1) + c2^(-1))) / b;
x];


(* Map A (composition of local maps) *)
AMap[matrix_,m_] := Module[{x=matrix},
x = aTildeMap[x,m+1,1];
x = aMap[x,m+2,2];
x = aTildeMap[x,m+1,2];
x = aMap[x,m+1,2];
x = aMap[x,m,1];
x = aMap[x,m+2,2];
x = aMap[x,m+1,1];
x = Simplify[x];
x];


(* Map B (composition of local maps) *)
BMap[matrix_,m_,q_] := Module[{x=matrix},
x = dMap[x,m+1,1,m+q,q+1];
x = dTildeInvMap[x,m+1,1,m+q,q+1];
x];


(* Map C (composition of local maps) *)
CMap[matrix_,m_,q_] := Module[{x=matrix},
x = aMap[x,m+1,1];
x = aMap[x,m+2,2];
x = aMap[x,m,1];
x = aMap[x,m+1,2];
x = BMap[x,m,q];
x = aTildeMap[x,m+1,2];
x = aMap[x,m+2,2];
x = aTildeMap[x,m+1,1];
x = dTildeMap[x,m+1,2,m+q,q+1];
x = dMap[x,m+2,3,m+q,q+1];
x = AMap[x,m];
x = dInvMap[x,m+2,3,m+q,q+1];
x = dInvMap[x,m+1,2,m+q,q+1];
x = dInvMap[x,m,1,m+q,q+1];
x = Simplify[x];
x];
