(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



t1=AbsoluteTime[];

M=1;

r0Int  =  ToExpression[($CommandLine)[[-2]]];
r0Frac =  ToExpression[($CommandLine)[[-1]]];

r0=r0Int + r0Frac/10;

lmin=0;
lmax=50;


p=r0;
e=0;


E0=Sqrt[((p-2-2e)(p-2+2e))/(p(p-3-e^2))];
L0=Sqrt[(p^2 M^2)/(p-3-e^2)];
f0=1-2M/r0;
ut=Sqrt[r0/(r0-3M)];
u\[Phi]=\[CapitalOmega]\[Phi] ut;

\[CapitalOmega]\[Phi]=N[Sqrt[M/r0^3],300];

rs[r_]=r+2M Log[r/(2M)-1];


\[Lambda][l_]:=(l+2)(l-1);
f[r_]:=1-(2M)/r;
V[l_,m_,\[Omega]_]:=1/4 ((2M f[r])/r^3+(l(l+1) f[r])/r^2-\[Omega]^2);
ode[1][l_,m_,\[Omega]_]:=f[r]^2 D[R[1][r],r,r]+(2M f[r])/r^2 D[R[1][r],r]-( 4 V [l,m,\[Omega]]R[1][r]+(4M)/r^2 f[r]^2 D[R[3][r],r]+(2 f[r])/r^2 (1-(4M)/r)(R[1][r]-R[5][r]-f[r] R[3][r])-(2f[r]^2)/r^2 (1-(6M)/r)R[6][r])==0;
ode[3][l_,m_,\[Omega]_]:=f[r]^2 D[R[3][r],r,r]+(2M f[r])/r^2 D[R[3][r],r]-(4 V [l,m,\[Omega]]R[3][r]-(2f[r])/r^2 (R[1][r]-R[5][r]-(1-(4M)/r)(R[3][r]+R[6][r])))==0;
ode[5][l_,m_,\[Omega]_]:=f[r]^2 D[R[5][r],r,r]+(2M f[r])/r^2 D[R[5][r],r]-(4V [l,m,\[Omega]]R[5][r]+(4f[r])/r^2 ((1-(9M)/(2r))R[5][r]-(l(l+1))/2 (R[1][r]-f[r] R[3][r])+1/2 (1-(3M)/r)(l(l+1) R[6][r]-R[7][r])))==0;
ode[6][l_,m_,\[Omega]_]:=f[r]^2 D[R[6][r],r,r]+(2M f[r])/r^2 D[R[6][r],r]-(4V[l,m,\[Omega]] R[6][r]-(2f[r])/r^2 (R[1][r]-R[5][r]-(1-(4M)/r)(R[3][r]+R[6][r])))==0;
ode[7][l_,m_,\[Omega]_]:=f[r]^2 D[R[7][r],r,r]+(2M f[r])/r^2 D[R[7][r],r]-(4V[l,m,\[Omega]] R[7][r]-(2f[r])/r^2 (R[7][r]+\[Lambda][l] R[5][r]))==0;
ode[9][l_,m_,\[Omega]_]:=f[r]^2 D[R[9][r],r,r]+(2M f[r])/r^2 D[R[9][r],r]-4( V [l,m,\[Omega]]+f[r]/r^2 (1-(9M)/(2r)))R[9][r]+(2f[r])/r^2 (1-(3M)/r)R[10][r]==0
ode[10][l_,m_,\[Omega]_]:=f[r]^2 D[R[10][r],r,r]+(2M f[r])/r^2 D[R[10][r],r]-4( V [l,m,\[Omega]]-f[r]/(2r^2))R[10][r]+(2f[r] \[Lambda][l])/r^2 R[9][r]==0;


gauge[1][l_,m_,\[Omega]_]:=I \[Omega] R[1][r]+f[r](I \[Omega] R[3][r]+D[R[2][r],r]+R[2][r]/r-R[4][r]/r) ;
gauge[2][l_,m_,\[Omega]_]:=-I \[Omega] R[2][r]-f[r] D[R[1][r],r]+f[r]^2 D[R[3][r],r]-f[r]/r (R[1][r]-R[5][r]-f[r]R[3][r]-2f[r]R[6][r]);
gauge[3][l_,m_,\[Omega]_]:=-I \[Omega] R[4][r]-f[r]/r (r D[R[5][r],r]+2R[5][r]+l(l+1)R[6][r]-R[7][r]);
gauge[4][l_,m_,\[Omega]_]:=-I \[Omega] R[8][r]-f[r]/r (r D[R[9][r],r]+2 R[9][r]-R[10][r]);


outerBCsEven[l1_,m1_,\[Omega]1_,rout1_,Nmax1_,a01_,fields_]:=Block[{l=l1,m=m1,\[Omega]=\[Omega]1,rout=rout1,Nmax=Nmax1,a0=a01,a,k,L,C1,D1,D3,E1,E3,F3,C5,D5,E5,C7,D7,E7,A1,A3,A5,A6,A7,r,rs,res,result,infCoeffs},
L=l(l+1);

C1[k_]=k(k+1)+4M I \[Omega] k-2-L;
D1[k_]=2M(5+L-2k^2-3k);
D3[k_]=2M(2k-8+4M I \[Omega]);
E1[k_]=4M^2 (k^2+2k-3);
E3[k_]=8M^2 (5-2k);
F3[k_]=16M^3 (k-2);
C5[k_]=C1[k]-2;
D5[k_]=2M(12-2k^2-3k+L);
E5[k_]=4M^2 (k^2+2k-8);
C7[k_]=k(k+1)+4M I \[Omega] k -L+2;
D7[k_]=2M(L-3-2k^2-3k);
E7[k_]=4M^2 (k^2+2k+1);

a[1][0]=a0[[1]];
a[3][0]=a0[[2]];
a[5][0]=a0[[3]];
a[6][0]=a0[[4]];
a[7][0]=a0[[5]];
Table[a[j1][j2]=0,{j1,fields},{j2,-3,-1}];


A1[k_]=(C1[k-1] a[1][k-1]+(2-4 M I \[Omega]) a[3][k-1]+2 a[5][k-1]+2 a[6][k-1]+D1[k-2] a[1][k-2]+D3[k-2] a[3][k-2]-12 M a[5][k-2]-20 M a[6][k-2]+E1[k-3] a[1][k-3]+E3[k-3] a[3][k-3]+16 M^2 a[5][k-3]+56 M^2 a[6][k-3]+F3[k-4] a[3][k-4]-48 M^3 a[6][k-4])/(2 I \[Omega] k);
A3[k_]=(C1[k-1] a[3][k-1]+2 a[1][k-1]-2 a[5][k-1]-2 a[6][k-1]+D1[k-2] a[3][k-2]+4 M (-a[1][k-2]+a[5][k-2]+3 a[6][k-2]) +E1[k-3] a[3][k-3]-16 M^2 a[6][k-3])/(2 I \[Omega] k);

A5[k_]=(C5[k-1]a[5][k-1]+2L(a[1][k-1]-a[3][k-1]-a[6][k-1])+2a[7][k-1]+D5[k-2]a[5][k-2]-10M a[7][k-2]+2M L(-2a[1][k-2]+4a[3][k-2]+5a[6][k-2])+E5[k-3]a[5][k-3]+4M^2 (-2L a[3][k-3]-3L a[6][k-3]+3a[7][k-3]))/(2I \[Omega] k);


A6[k_]=(C1[k-1] a[6][k-1]+2 a[1][k-1]-2 a[5][k-1]-2 a[3][k-1]+D1[k-2] a[6][k-2]+4 M (-a[1][k-2]+a[5][k-2]+3 a[3][k-2]) +E1[k-3] a[6][k-3]-16 M^2 a[3][k-3])/(2 I \[Omega] k);
A7[k_]=(C7[k-1]a[7][k-1]+2\[Lambda][l] a[5][k-1]+D7[k-2]a[7][k-2]-4M \[Lambda][l]a[5][k-2]+E7[k-3]a[7][k-3])/(2 I \[Omega] k);

For[k=1,k<=Nmax,k++,a[1][k]=A1[k]; a[3][k]=A3[k];a[5][k]=A5[k]; a[6][k]=A6[k];a[7][k]=A7[k]];
rs=r+2M Log[r/(2M)-1];

infCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
Table[infCoeffs[[i]]=Table[a[i][j],{j,0,7}],{i,fields}];

Table[res[j]=Exp[I \[Omega] rs] Sum[a[j][i]/r^i,{i,0,Nmax}],{j,fields}];
result=Table[ode[j][l,m,\[Omega]][[1]],{j,fields}]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rout};
Join[Table[res[j]/.r->rout,{j,fields}],Table[D[res[j],r]/.r->rout,{j,fields}],result,infCoeffs]
]


outerBCsOdd[l1_,m1_,\[Omega]1_,rout1_,Nmax1_,a01_,fields_]:=Block[{l=l1,m=m1,\[Omega]=\[Omega]1,rout=rout1,Nmax=Nmax1,a0=a01,a,k,L,C1,I1,D1,J1,E1,K1,A9,A10,r,rs,res,result,infCoeffs},
L=l(l+1);

C1[k_]=4M I \[Omega] k + k(k+1)-L-4;
I1[k_]=4M I \[Omega] k+k(k+1)-L+2;
D1[k_]=-6M k-4M k^2+24M+2M L;
J1[k_]=-6M k-4M k^2-6M+2M L;
E1[k_]=4M^2 (k^2+2k-8);
K1[k_]=4M^2 (k^2+2k+1);


a[9][0]=a0[[1]];
a[10][0]=a0[[2]];
Table[a[j1][j2]=0,{j1,fields},{j2,-3,-1}];


A9[k_]=(C1[k-1]a[9][k-1]+D1[k-2]a[9][k-2]+E1[k-3]a[9][k-3]+2a[10][k-1]-10M a[10][k-2]+12M^2 a[10][k-3])/(2I \[Omega] k);
A10[k_]=(I1[k-1]a[10][k-1]+J1[k-2]a[10][k-2]+K1[k-3]a[10][k-3]+2\[Lambda][l] a[9][k-1]-4M \[Lambda][l] a[9][k-2])/(2 I \[Omega] k);

For[k=1,k<=Nmax,k++,a[9][k]=A9[k]; a[10][k]=A10[k];];
rs=r+2M Log[r/(2M)-1];

infCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
Table[infCoeffs[[i]]=Table[a[i][j],{j,0,7}],{i,fields}];

Table[res[j]=Exp[I \[Omega] rs] Sum[a[j][i]/r^i,{i,0,Nmax}],{j,fields}];
result=Table[ode[j][l,m,\[Omega]][[1]],{j,fields}]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rout};
Join[Table[res[j]/.r->rout,{j,fields}],Table[D[res[j],r]/.r->rout,{j,fields}],result,infCoeffs]
]


outerBCsDipole[l1_,m1_,\[Omega]1_,rout1_,Nmax1_,a01_,fields_]:=Block[{l=l1,m=m1,\[Omega]=\[Omega]1,rout=rout1,Nmax=Nmax1,a0=a01,a,k,L,C1,D1,D3,E1,E3,F3,C5,D5,E5,A1,A3,A5,A6,r,rs,res,result,infCoeffs},
L=l(l+1);

C1[k_]=k(k+1)+4M I \[Omega] k-2-L;
D1[k_]=2M(5+L-2k^2-3k);
D3[k_]=2M(2k-8+4M I \[Omega]);
E1[k_]=4M^2 (k^2+2k-3);
E3[k_]=8M^2 (5-2k);
F3[k_]=16M^3 (k-2);
C5[k_]=C1[k]-2;
D5[k_]=2M(12-2k^2-3k+L);
E5[k_]=4M^2 (k^2+2k-8);

a[1][0]=a0[[1]];
a[3][0]=a0[[2]];
a[5][0]=a0[[3]];
a[6][0]=a0[[4]];
Table[a[j1][j2]=0,{j1,fields},{j2,-3,-1}];


A1[k_]=(C1[k-1] a[1][k-1]+(2-4 M I \[Omega]) a[3][k-1]+2 a[5][k-1]+2 a[6][k-1]+D1[k-2] a[1][k-2]+D3[k-2] a[3][k-2]-12 M a[5][k-2]-20 M a[6][k-2]+E1[k-3] a[1][k-3]+E3[k-3] a[3][k-3]+16 M^2 a[5][k-3]+56 M^2 a[6][k-3]+F3[k-4] a[3][k-4]-48 M^3 a[6][k-4])/(2 I \[Omega] k);
A3[k_]=(C1[k-1] a[3][k-1]+2 a[1][k-1]-2 a[5][k-1]-2 a[6][k-1]+D1[k-2] a[3][k-2]+4 M (-a[1][k-2]+a[5][k-2]+3 a[6][k-2]) +E1[k-3] a[3][k-3]-16 M^2 a[6][k-3])/(2 I \[Omega] k);
A5[k_]=(C5[k-1] a[5][k-1]+2 L a[1][k-1]-2 L a[3][k-1]-2 L a[6][k-1]+D5[k-2] a[5][k-2]+2 M L (-2 a[1][k-2]+4 a[3][k-2]+5 a[6][k-2]) +E5[k-3] a[5][k-3]+4 M^2 (-2 L a[3][k-3]-3 L a[6][k-3]))/(2 I \[Omega] k);
A6[k_]=(C1[k-1] a[6][k-1]+2 a[1][k-1]-2 a[5][k-1]-2 a[3][k-1]+D1[k-2] a[6][k-2]+4 M (-a[1][k-2]+a[5][k-2]+3 a[3][k-2]) +E1[k-3] a[6][k-3]-16 M^2 a[3][k-3])/(2 I \[Omega] k);

For[k=1,k<=Nmax,k++,a[1][k]=A1[k]; a[3][k]=A3[k];a[5][k]=A5[k]; a[6][k]=A6[k]];
rs=r+2M Log[r/(2M)-1];

infCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
Table[infCoeffs[[i]]=Table[a[i][j],{j,0,7}],{i,fields}];

Table[res[j]=Exp[I \[Omega] rs] Sum[a[j][i]/r^i,{i,0,Nmax}],{j,fields}];
R[7][r]=0;
result=Table[ode[j][l,m,\[Omega]][[1]],{j,fields}]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rout};
Clear[R];
Join[Table[res[j]/.r->rout,{j,fields}],Table[D[res[j],r]/.r->rout,{j,fields}],result,infCoeffs]
]


outerBCsStaticEven[l1_,rout1_,Nmax1_,a01_,fields_]:=Block[{l=l1,rout=rout1,Nmax=Nmax1,a0=a01,a,ab,k,L,C1,C3,C5,D1,D3,D5,I3,J3,E1,E3,E5,H3,F1,F3,G1,G3,r,rs,res,result,c1,c2,c3,c4,r6replace,r7replace,Mat,AR,BR,CR,DR,ER,FR,fl,infCoeffs,infLogCoeffs},
L=l(l+1);

C1[k_]=L+1-k^2;
D1[k_]=k(k-1);
D3[k_]=2(k+1);
E1[k_]=1-2k;
F3[k_]=k+1;

G3[k_]=-(2k^2-2-L);
G1[k_]=-2k;
H3[k_]=4k;
I3[k_]=k^2-1;
J3[k_]=-2k;

C5[k_]=L+k(1-k);
D5[k_]=2k-1;
E5[k_]=k(1-k)+2;


a[3][l]=a01[[1]];
a[5][l]=a01[[2]];
a[5][l+2]=a01[[3]];

Table[{a[j1][j2]=0,ab[j1][j2]=0},{j1,fields},{j2,0(*l-4*),l-1}];

a[1][l]=a[3][l]+a[5][l]/(l+1);
a[1][l+1]=(l-1)/l (l a[3][l]+a[5][l]);
a[3][l+1]=(l+1)a[3][l]+a[5][l]/l;
a[5][l+1]=(2l(l+1)a[3][l]+a[5][l](l^3-l^2-2))/(l(l+1));
ab[1][l]=0;
ab[3][l]=0;
ab[5][l]=0;
ab[1][l+1]=0;
ab[3][l+1]=0;
ab[5][l+1]=0;

fl=(l+1)(1+2l)(3+2l);
ab[1][l+2]=-((2(l+2))/fl)(2a[3][l](l+1)-a[5][l](l-3));
ab[3][l+2]=-ab[1][l+2];
ab[5][l+2]=2l ab[3][l+2];

a[1][l+2]=a[3][l]/(l fl) (-6-2l+12l^2+9l^3+4l^4+5l^5+2l^6)+a[5][l]/(2l^2 fl) (12+8l+6l^2-15l^3-8l^4+15l^5+6l^6)-a[5][l+2]/(2l);
a[3][l+2]=a[5][l]/(2l^2 fl) (-12-8l+6l^2+51l^3+42l^4+7l^5-2l^6)+a[3][l]/(l fl) (6+6l+6l^2+23l^3+28l^4+13l^5+2l^6)+a[5][l+2]/(2l);


Mat[k_]={
{C1[k],2k,-(k+1),1,-1,0},
{0,C1[k],0,-(k+1),0,-1},
{-(k+1),1,C1[k],2k,1,0},
{0,-(k+1),0,C1[k],0,1},
{-2L,0,2L,0,C5[k],D5[k]},
{0,-2L,0,2L,0,C5[k]}
};


AR[k_]=-2M(D1[k-1]a[1][k-1]+D3[k-1]a[3][k-1]+a[5][k-1]+E1[k-1]ab[1][k-1]-2ab[3][k-1])+4M^2 (F3[k-2]a[3][k-2]-ab[3][k-2]);
BR[k_]=-2M(D1[k-1]ab[1][k-1]+D3[k-1]ab[3][k-1]+ab[5][k-1])+4M^2 F3[k-2]ab[3][k-2];
CR[k_]=2M(G3[k-1]a[3][k-1]+G1[k-1]a[1][k-1]+H3[k-1]ab[3][k-1]+2ab[1][k-1])+4M^2 (I3[k-2]a[3][k-2]+J3[k-2]ab[3][k-2]);
DR[k_]=2M(G3[k-1]ab[3][k-1]+G1[k-1]ab[1][k-1])+4M^2 I3[k-2]ab[3][k-2];
ER[k_]=2M(E5[k-1]a[5][k-1]+2L a[3][k-1]+D5[k-1]ab[5][k-1]);
FR[k_]=2M(E5[k-1]ab[5][k-1]+2L ab[3][k-1]);



For[k=l+3,k<=Nmax,k++,{a[1][k],ab[1][k],a[3][k],ab[3][k],a[5][k],ab[5][k]}=Inverse[Mat[k]] . {AR[k],BR[k],CR[k],DR[k],ER[k],FR[k]}];

rs=r+2M Log[r/(2M)-1];

r6replace=Solve[gauge[2][l,0,0]==0,R[6][r]][[1,1]];
r7replace=Solve[gauge[3][l,0,0]==0,R[7][r]][[1,1]];

infCoeffs=ConstantArray[{0,0,0,0,0,0,0,0,0,0},10];
infLogCoeffs=ConstantArray[{0,0,0,0,0,0,0,0,0,0},10];
Table[infCoeffs[[i]]=Table[a[i][j],{j,0,9}],{i,fields}];
Table[infLogCoeffs[[i]]=Table[ab[i][j],{j,0,9}],{i,fields}];

Table[res[j]=Sum[(a[j][i]+ab[j][i]Log[r])/r^i,{i,l,Nmax-1}],{j,{1,3,5}}];
result=Table[ode[j][l,0,0][[1]]/.r7replace/.r6replace,{j,fields}]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rout};
Join[Table[res[j]/.r->rout,{j,fields}],Table[D[res[j],r]/.r->rout,{j,fields}],result,infCoeffs,infLogCoeffs]
]


innerBCsEven[l1_,m1_,\[Omega]1_,rin1_,Nmax1_,b01_,fields_]:=Block[{l=l1,m=m1,\[Omega]=\[Omega]1,rin=rin1,Nmax=Nmax1,b0=b01,b,k,L,C1,C3,D1,D3,E1,E3,F1,F3,G3,H3,I5,J5,C5,D5,E5,B1,B3,B5,B6,B7,K7,r,rs,res,result,c1,c2,c3,c4,horizCoeffs},
L=l(l+1);

C1[k_]=4M^2 (1+3k^2-L-16M I \[Omega] k -k);
C3[k_]=8M^2 (2M I \[Omega]-k);
D1[k_]=2M(3k^2-2k-2L-1-24 M I \[Omega] k);
D3[k_]=4M(4M I \[Omega]-k-1);
E1[k_]=k(k-1)-16M I \[Omega] k - 2-L;
F1[k_]=2 I \[Omega] k;
G3[k_]=2M(2k^2-k-L+1-12M I \[Omega] k);
H3[k_]=k(k-1)-L-2-12 I \[Omega] k;
I5[k_]=2M(2k^2-k-L+4-12M I \[Omega] k);
J5[k_]=k(k-1)-L-4-12M I \[Omega] k;
K7[k_]=k(k-1)-L+2-12M I \[Omega] k;

b[1][0]=b0[[1]];
b[3][0]=b0[[2]];
b[5][0]=b0[[3]];
b[6][0]=b0[[4]];
b[7][0]=b0[[5]];
Table[b[j1][j2]=0,{j1,fields},{j2,-4,-1}];

B1[k_]=(C1[k-1]b[1][k-1]+C3[k-1]b[3][k-1]-8M^2 b[5][k-1]+D1[k-2]b[1][k-2]+D3[k-2]b[3][k-2]-8M b[6][k-2]+E1[k-3]b[1][k-3]+2(1+2M I \[Omega])b[3][k-3]+2b[5][k-3]+2b[6][k-3]-F1[k-4]b[1][k-4])/(8M^3 k(4M I \[Omega]-k));
B3[k_]=(G3[k-1]b[3][k-1]+4M(b[1][k-1]-b[5][k-1]+b[6][k-1])+H3[k-2]b[3][k-2]+2(b[1][k-2]-b[5][k-2]-b[6][k-2])-F1[k-3]b[3][k-3])/(4M^2 k(4M I \[Omega]-k));
B5[k_]=(I5[k-1]b[5][k-1]+2M L(2b[1][k-1]+b[6][k-1])+-2M b[7][k-1]+J5[k-2]b[5][k-2]+2L(b[1][k-2]-b[3][k-2]-b[6][k-2])+2b[7][k-2]-F1[k-3]b[5][k-3])/(4M^2 k(4M I \[Omega]-k));
B6[k_]=(G3[k-1]b[6][k-1]+4M(b[1][k-1]-b[5][k-1]+b[3][k-1])+H3[k-2]b[6][k-2]+2(b[1][k-2]-b[5][k-2]-b[3][k-2])-F1[k-3]b[6][k-3])/(4M^2 k(4M I \[Omega]-k));
B7[k_]=(G3[k-1]b[7][k-1]+4M \[Lambda][l] b[5][k-1]+K7[k-2]b[7][k-2]+2\[Lambda][l] b[5][k-2]-F1[k-3]b[7][k-3])/(4M^2 k(4M I \[Omega]-k));

For[k=1,k<=Nmax,k++,b[1][k]=B1[k]; b[3][k]=B3[k];b[5][k]=B5[k]; b[6][k]=B6[k];b[7][k]=B7[k]];
rs=r+2M Log[r/(2M)-1];

horizCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
Table[horizCoeffs[[i]]=Table[b[i][j],{j,0,7}],{i,fields}];

Table[res[j]=Exp[-I \[Omega] rs] Sum[b[j][i](r-2M)^i,{i,0,Nmax}],{j,fields}];
result=Abs[Table[ode[j][l,m,\[Omega]][[1]],{j,fields}]]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rin};
Join[Table[res[j]/.r->rin,{j,fields}],Table[D[res[j],r]/.r->rin,{j,fields}],result,horizCoeffs]
]


innerBCsOdd[l1_,m1_,\[Omega]1_,rin1_,Nmax1_,b01_,fields_]:=Block[{l=l1,m=m1,\[Omega]=\[Omega]1,rin=rin1,Nmax=Nmax1,b0=b01,b,k,L,C1,H1,D1,J1,E1,B9,B10,r,rs,res,result,c1,c2,c3,c4,horizCoeffs},
L=l(l+1);

C1[k_]=2M(k+12M I \[Omega] k-2k^2+L-4);
H1[k_]=2M(k+12M I \[Omega] k-2k^2+L-1);
D1[k_]=4+12M I \[Omega] k+L-k(k-1);
J1[k_]=-2+12M I \[Omega] k+L-k(k-1);
E1[k_]=2 I \[Omega] k;

b[9][0]=b0[[1]];
b[10][0]=b0[[2]];
Table[b[j1][j2]=0,{j1,fields},{j2,-4,-1}];

B9[k_]=(C1[k-1]b[9][k-1]+D1[k-2]b[9][k-2]+E1[k-3]b[9][k-3]+2M b[10][k-1]-2b[10][k-2])/(4M^2 k(k-4M I \[Omega]));
B10[k_]=(H1[k-1]b[10][k-1]+J1[k-2]b[10][k-2]+E1[k-3]b[10][k-3]-4M \[Lambda][l] b[9][k-1]-2\[Lambda][l] b[9][k-2])/(4M^2 k(k-4M I \[Omega]));


For[k=1,k<=Nmax,k++,b[9][k]=B9[k]; b[10][k]=B10[k];];
rs=r+2M Log[r/(2M)-1];

horizCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
Table[horizCoeffs[[i]]=Table[b[i][j],{j,0,7}],{i,fields}];

Table[res[j]=Exp[-I \[Omega] rs] Sum[b[j][i](r-2M)^i,{i,0,Nmax}],{j,fields}];
result=Abs[Table[ode[j][l,m,\[Omega]][[1]],{j,fields}]]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rin};
Join[Table[res[j]/.r->rin,{j,fields}],Table[D[res[j],r]/.r->rin,{j,fields}],result,horizCoeffs]
]


innerBCsDipole[l1_,m1_,\[Omega]1_,rin1_,Nmax1_,b01_,fields_]:=Block[{l=l1,m=m1,\[Omega]=\[Omega]1,rin=rin1,Nmax=Nmax1,b0=b01,b,k,L,C1,C3,D1,D3,E1,E3,F1,F3,G3,H3,I5,J5,C5,D5,E5,B1,B3,B5,B6,r,rs,res,result,c1,c2,c3,c4,horizCoeffs},
L=l(l+1);

C1[k_]=4M^2 (1+3k^2-L-16M I \[Omega] k -k);
C3[k_]=8M^2 (2M I \[Omega]-k);
D1[k_]=2M(3k^2-2k-2L-1-24 M I \[Omega] k);
D3[k_]=4M(4M I \[Omega]-k-1);
E1[k_]=k(k-1)-16M I \[Omega] k - 2-L;
F1[k_]=2 I \[Omega] k;
G3[k_]=2M(2k^2-k-L+1-12M I \[Omega] k);
H3[k_]=k(k-1)-L-2-12 I \[Omega] k;
I5[k_]=2M(2k^2-k-L+4-12M I \[Omega] k);
J5[k_]=k(k-1)-L-4-12M I \[Omega] k;

b[1][0]=b0[[1]];
b[3][0]=b0[[2]];
b[5][0]=b0[[3]];
b[6][0]=b0[[4]];
Table[b[j1][j2]=0,{j1,fields},{j2,-4,-1}];

B1[k_]=(C1[k-1]b[1][k-1]+C3[k-1]b[3][k-1]-8M^2 b[5][k-1]+D1[k-2]b[1][k-2]+D3[k-2]b[3][k-2]-8M b[6][k-2]+E1[k-3]b[1][k-3]+2(1+2M I \[Omega])b[3][k-3]+2b[5][k-3]+2b[6][k-3]-F1[k-4]b[1][k-4])/(8M^3 k(4M I \[Omega]-k));
B3[k_]=(G3[k-1]b[3][k-1]+4M(b[1][k-1]-b[5][k-1]+b[6][k-1])+H3[k-2]b[3][k-2]+2(b[1][k-2]-b[5][k-2]-b[6][k-2])-F1[k-3]b[3][k-3])/(4M^2 k(4M I \[Omega]-k));
B5[k_]=(I5[k-1]b[5][k-1]+2M L(2b[1][k-1]+b[6][k-1])+J5[k-2]b[5][k-2]+2L(b[1][k-2]-b[3][k-2]-b[6][k-2])-F1[k-3]b[5][k-3])/(4M^2 k(4M I \[Omega]-k));
B6[k_]=(G3[k-1]b[6][k-1]+4M(b[1][k-1]-b[5][k-1]+b[3][k-1])+H3[k-2]b[6][k-2]+2(b[1][k-2]-b[5][k-2]-b[3][k-2])-F1[k-3]b[6][k-3])/(4M^2 k(4M I \[Omega]-k));

For[k=1,k<=Nmax,k++,b[1][k]=B1[k]; b[3][k]=B3[k];b[5][k]=B5[k]; b[6][k]=B6[k]];
rs=r+2M Log[r/(2M)-1];

horizCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
Table[horizCoeffs[[i]]=Table[b[i][j],{j,0,7}],{i,fields}];

Table[res[j]=Exp[-I \[Omega] rs] Sum[b[j][i](r-2M)^i,{i,0,Nmax}],{j,{1,3,5,6}}];
result=Abs[Table[ode[j][l,m,\[Omega]][[1]],{j,fields}]]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rin};
Join[Table[res[j]/.r->rin,{j,fields}],Table[D[res[j],r]/.r->rin,{j,fields}],result,horizCoeffs]
]


innerBCsStaticEven[l1_,rin1_,Nmax1_,b01_,fields_]:=Block[{l=l1,rin=rin1,Nmax=Nmax1,b0=b01,b,k,L,C1,C3,C5,D1,D3,D5,E1,E3,F1,F3,G1,G3,B1,B3,B5,B6,r,rs,res,result,c1,c2,c3,c4,r6replace,r7replace,horizCoeffs},
L=l(l+1);

C1[k_]=-4 M^2 (k+1);
C3[k_]=4M^2 k(k-1);
C5[k_]=2M L-4M(1+k^2);
D1[k_]=2M(k-2);
D3[k_]=2M(L+k(1-2k));
D5[k_]=L-k(k+1);
E1[k_]=L+1-k^2;
E3[k_]=k-1;
G3[k_]=2M k;
F1[k_]=4M^2 (L+1+4k-3k^2);
G1[k_]=2M(2L+2+2k-3k^2);

b[1][0]=0;
b[1][1]=0;
b[3][0]=b0[[1]];
b[3][1]=b0[[2]];
b[5][0]=0;
b[5][1]=b0[[3]];

Table[b[j1][j2]=0,{j1,fields},{j2,-4,-1}];


B1[k_]=(F1[k-1]b[1][k-1]+G1[k-2]b[1][k-2]+G3[k-2]b[3][k-2]-2M b[5][k-2]+E1[k-3]b[1][k-3]+E3[k-3]b[3][k-3]-b[5][k-3])/(8M^3 k(k-2));
B5[k_]=(C5[k-1]b[5][k-1]-4M L b[1][k-1]+D5[k-2]b[5][k-2]+2L(b[3][k-2]-b[1][k-2]))/(4M k(k-1));
B3[k_]=(C1[k-1]b[1][k-1]-8M^3 k b[1][k]+4M^2 b[5][k-1]+D3[k-2]b[3][k-2]+D1[k-2]b[1][k-2]+4M b[5][k-2]+E1[k-3]b[3][k-3]+E3[k-3]b[1][k-3]+b[5][k-3])/(C3[k-1]);

b[1][2]=(4M^2 b[5][k-1]+D3[k-2] b[3][k-2])/(8k M^3)/.k->2;
b[5][2]=B5[2];

For[k=3,k<=Nmax,k++,b[1][k]=B1[k];b[5][k]=B5[k];  b[3][k-1]=B3[k];];

rs=r+2M Log[r/(2M)-1];

r6replace=Solve[gauge[2][l,0,0]==0,R[6][r]][[1,1]];
r7replace=Solve[gauge[3][l,0,0]==0,R[7][r]][[1,1]];

horizCoeffs=ConstantArray[{0,0,0,0,0,0,0,0,0,0},10];
Table[horizCoeffs[[i]]=Table[b[i][j],{j,0,9}],{i,fields}];

Table[res[j]=Sum[b[j][i](r-2M)^i,{i,0,Nmax-1}],{j,{1,3,5}}];
result=Abs[Table[ode[j][l,0,0][[1]]/.r7replace/.r6replace,{j,fields}]]//.{R[i_][r]:>res[i],Derivative[n_][R[i_]][r]:>D[res[i],{r,n}],r->rin};
Join[Table[res[j]/.r->rin,{j,fields}],Table[D[res[j],r]/.r->rin,{j,fields}],result,horizCoeffs]
]


(*innerBCsStaticEven[N[4,100],2+10^-2,20,{0,1,0},{1,3,5}][[-14;;-11]]//N
outerBCsStaticEven[N[2,100],10,70,{0,1,0},{1,3,5}][[-3;;]]//N*)


(*N[innerBCsStaticEven[N[2,100],2+10^-2,30,{0,0,1},{1,3,5}],16]*)


(*N[innerBCsEven[N[2,100],2,1/10,2+10^-2,10,{1,0,0,0,0},{1,3,5,6,7}],10][[-10;;]]//MatrixForm
N[innerBCsOdd[N[2,100],2,1/10,2+10^-2,10,{1,0},{9,10}],16][[-10;;]]//MatrixForm*)


(*outerBCsEven[2,2,0.1`100,1000,20,{1,0,0,0,0},{1,3,5,6,7}][[-9]]*)


(*outerBCsEven[2,2,0.1`100,1000,10,{1,0,0,0,0},{1,3,5,6,7}][[-5-10;;-11]]//N*)


sourceEvenDipole[l_,m_]:=(-16 Pi E0)/f0^2 SphericalHarmonicY[l,m,Pi/2,0]{0,0,0,0,f0^2/r0,f0/r0,0,r0 \[CapitalOmega]\[Phi]^2};
sourceEven[l_,m_]:=(-16 Pi E0)/f0^2 SphericalHarmonicY[l,m,Pi/2,0]{0,0,0,0,0,f0^2/r0,f0/r0,0,r0 \[CapitalOmega]\[Phi]^2,r0 \[CapitalOmega]\[Phi]^2 (l(l+1)-2m^2)};
sourceOdd[l_,m_]:=(-16 Pi E0)/f0^2 D[SphericalHarmonicY[l,m,\[Theta],0],\[Theta]]{0,0,0,2I m r0 \[CapitalOmega]\[Phi]^2}/.\[Theta]->\[Pi]/2
sourceEvenStatic[l_,m_]:=(-16 Pi E0)/f0^2 SphericalHarmonicY[l,m,Pi/2,0]{0,0,0,f0^2/r0,f0/r0,0};


REq[i_]:=R[i][r]/.Solve[gauge[Switch[i,2,2,4,3,8,4]][l,m,\[Omega]]==0,R[i][r]][[1,1]]


nMax=12;


F[i_][r]:=Sum[a[i][n]/r^n,{n,0,nMax}]


derivF=Simplify[D[Exp[I \[Omega] rs[r]] F[r],r]/Exp[I \[Omega] rs[r]]]


dF[i_][r]:=derivF/.{F[r]->F[i][r],F'[r]->D[F[i][r],r]}


replace={R[i_][r]:>F[i][r],R[i_]'[r]:>dF[i][r]}


Do[GaugeExpansion\[Infinity][i]=Series[{R[i][r],REq[i]}/.replace/.{\[Omega]->\[Omega]1,l->l1},{r,\[Infinity],10}]//Simplify,{i,{2,4,8}}]


Do[GaugeExpansionReplacement\[Infinity][i]=Table[SeriesCoefficient[GaugeExpansion\[Infinity][i][[1]],n]->SeriesCoefficient[GaugeExpansion\[Infinity][i][[2]],n],{n,0,nMax}],{i,{2,4,8}}]


expStaticEven [i_][r_]:=Sum[ (a[i][n]+al[i][n]Log[r])/r^n,{n,0,nMax}];


r6replace=Solve[gauge[2][l,0,0]==0,R[6][r]][[1,1]];
r7replace=Simplify[Solve[gauge[3][l,0,0]==0,R[7][r]][[1]]/.r6replace][[1]];


Do[staticEvenExpansion[j]=Series[{R[j][r],R[j][r]/.r6replace/.r7replace/.l->l1}/.{R[i_][r]:>expStaticEven[i][r],Derivative[n_][R[i_]][r]:>D[expStaticEven[i][r],{r,n}]},{r,\[Infinity],nMax}];,{j,{6,7}}]


Do[GaugeStaticEvenNonLogReplacement\[Infinity][j]=Table[Coefficient[SeriesCoefficient[staticEvenExpansion[j][[1]],n],Log[r],0]->Coefficient[SeriesCoefficient[staticEvenExpansion[j][[2]],n],Log[r],0],{n,0,nMax}];

GaugeStaticEvenLogReplacement\[Infinity][j]=Table[Coefficient[SeriesCoefficient[staticEvenExpansion[j][[1]],n],Log[r],1]->Coefficient[SeriesCoefficient[staticEvenExpansion[j][[2]],n],Log[r],1],{n,0,nMax}];
,{j,{6,7}}]


G[i_][\[Rho]_]:=Sum[b[i][n] \[Rho]^n,{n,0,nMax}]


derivG=Simplify[D[Exp[-I \[Omega] rs[r]] G[\[Rho]]/.r->\[Rho]+2M,\[Rho]]/Exp[-I \[Omega] rs[r]]]/.r->\[Rho]+2M


dG[i_][r]:=derivG/.{G[\[Rho]]->G[i][\[Rho]],G'[\[Rho]]->D[G[i][\[Rho]],\[Rho]]}


replaceH={R[i_][r]:>G[i][\[Rho]],R[i_]'[r]:>dG[i][r]}


Do[GaugeExpansionH[i]=Series[{R[i][r],REq[i]}/.replaceH/.r->\[Rho]+2M/.{\[Omega]->\[Omega]1,l->l1},{\[Rho],0,10}]//Simplify,{i,{2,4,8}}]


Do[GaugeExpansionReplacementH[i]=Table[SeriesCoefficient[GaugeExpansionH[i][[1]],n]->SeriesCoefficient[GaugeExpansionH[i][[2]],n],{n,0,nMax}],{i,{2,4,8}}]


expStaticEvenH [i_][r_]:=Sum[ b[i][n](r-2M)^n,{n,0,nMax}];


Do[staticEvenExpansionH[j]=Series[{R[j][r],R[j][r]/.r6replace/.r7replace/.l->l1}/.{R[i_][r]:>expStaticEvenH[i][r],Derivative[n_][R[i_]][r]:>D[expStaticEvenH[i][r],{r,n}]}/.r->\[Rho]+2M,{\[Rho],0,nMax}];,{j,{6,7}}]


Do[GaugeStaticEvenReplacementH[j]=Table[SeriesCoefficient[staticEvenExpansionH[j][[1]],n]->SeriesCoefficient[staticEvenExpansionH[j][[2]],n],{n,0,nMax}];
,{j,{6,7}}]


kmaxMax=100;
BCThreshold=10^-20;


MostCircLorenzModes[l_,m_,WP_,PG_]:=Module[{nf,fields,\[Omega],rin, rout,kmax,\[CapitalPhi]0,RIn,ROut,outerBCBasis,innerBCBasis,c,RInhom,dR,d2R,f,source,gauge2,greplace,R6rep,R7rep,horizCoeffs,replace, infCoeffs,infLogCoeffs},
Clear[R];
\[Omega]=m \[CapitalOmega]\[Phi];

If[Mod[l+m,2]==0,
If[m==0,
fields={1,3,5};nf=3;,
fields={1,3,5,6,7};nf=5;],
fields={9,10};nf=2;
];
If[l==1&&m==1,fields={1,3,5,6}; nf=4;];

If[m==0,rout=r0+10;kmax=l+40;,rout=Abs[10/\[Omega]];kmax=15;];

While[True,
If[Mod[l+m,2]==0,
If[l==1&&m==1,
Table[outerBCBasis[i]=outerBCsDipole[l,m,\[Omega],rout,kmax,IdentityMatrix[nf][[i]],fields],{i,1,nf}];,
If[m==0,
(*Table[outerBCBasis[i]={rout^l,rout^l,rout^l,rout^l,rout^l,rout^l,1,1,1,1}outerBCsStaticEven[N[l,100],rout,kmax,IdentityMatrix[nf][[i]],fields],{i,1,nf}]*)
Table[outerBCBasis[i]=Join[{rout^l,rout^l,rout^l,rout^l,rout^l,rout^l},ConstantArray[1,23]]outerBCsStaticEven[l,N[rout,100],kmax,IdentityMatrix[nf][[i]],fields],{i,1,nf}];
,
Table[outerBCBasis[i]=outerBCsEven[l,m,\[Omega],rout,kmax,IdentityMatrix[nf][[i]],fields],{i,1,nf}];
]
],
Table[outerBCBasis[i]=outerBCsOdd[l,m,\[Omega],rout,kmax,IdentityMatrix[nf][[i]],fields],{i,1,nf}];
];

(*Print[N[Total[Abs[Table[outerBCBasis[i][[-nf-10;;-11]],{i,nf}]],2]]," ", rout," ",kmax];*)

If[Total[Abs[Table[outerBCBasis[i][[-nf-10+If[m==0,-10,0];;-11+If[m==0,-10,0]]],{i,nf}]],2]<BCThreshold,
Break[],
kmax+=20;
If[kmax>kmaxMax,rout*=1+2/10;kmax=30];
]
];

rin=2+10^-2;

If[Mod[l+m,2]==0,
If[l==1&&m==1,
Table[innerBCBasis[i]=innerBCsDipole[l,m,\[Omega],rin,20,IdentityMatrix[nf][[i]],fields],{i,1,nf}];,
If[m==0,
Table[innerBCBasis[i]=innerBCsStaticEven[N[l,100],rin,20,IdentityMatrix[nf][[i]],fields],{i,1,nf}],
Table[innerBCBasis[i]=innerBCsEven[l,m,\[Omega],rin,20,IdentityMatrix[nf][[i]],fields],{i,1,nf}];
]
],
Table[innerBCBasis[i]=innerBCsOdd[l,m,\[Omega],rin,20,IdentityMatrix[nf][[i]],fields],{i,1,nf}];
];

R6rep=Solve[gauge[2][l,0,0]==0,R[6][r]][[1]];
R7rep=Solve[gauge[3][l,0,0]==0,R[7][r]][[1]]/.R6rep;

If[l==1&&m==1,R[7][r]=0;];
If[Mod[l+m,2]==0 && m==0,R[6][r]=R[6][r]/.R6rep;R[7][r]=R[7][r]/.R7rep;];

\[CapitalPhi]0=ConstantArray[0,{2nf,2nf}];
Do[
RIn=NDSolve[
Flatten[Table[{ode[fields[[j]]][l,m,\[Omega]],R[fields[[j]]][rin]==innerBCBasis[i][[j]],R[fields[[j]]]'[rin]==innerBCBasis[i][[j+nf]]},{j,1,nf}]],
Flatten[Table[{R[i][r0],R[i]'[r0]},{i,fields}]],{r,rin,r0},WorkingPrecision->WP,PrecisionGoal->PG,MaxSteps->Infinity,Method->"StiffnessSwitching"];
ROut=NDSolve[
Flatten[Table[{ode[fields[[j]]][l,m,\[Omega]],R[fields[[j]]][rout]==outerBCBasis[i][[j]],R[fields[[j]]]'[rout]==outerBCBasis[i][[j+nf]]},{j,1,nf}]],
Flatten[Table[{R[i][r0],R[i]'[r0]},{i,fields}]],{r,rout,r0},WorkingPrecision->WP,PrecisionGoal->PG,MaxSteps->Infinity,Method->"StiffnessSwitching"];
Table[\[CapitalPhi]0[[j,i]]=(-R[fields[[j]]][r0]/.RIn)[[1]],{j,1,nf}];
Table[\[CapitalPhi]0[[j+nf,i]]=(-R[fields[[j]]]'[r0]/.RIn)[[1]],{j,1,nf}];

Table[\[CapitalPhi]0[[j,i+nf]]=(R[fields[[j]]][r0]/.ROut)[[1]],{j,1,nf}];
Table[\[CapitalPhi]0[[j+nf,i+nf]]=(R[fields[[j]]]'[r0]/.ROut)[[1]],{j,1,nf}];

,{i,1,nf}
];

If[EvenQ[l+m],If[l==1,source=sourceEvenDipole,If[m==0,source=sourceEvenStatic,source=sourceEven]],source=sourceOdd];

c=Inverse[\[CapitalPhi]0] . source[l,m];

Table[RInhom[-1][i]=0,{i,1,10}];
Table[RInhom[+1][i]=0,{i,1,10}];

Table[dR[-1][i]=0,{i,1,10}];
Table[dR[+1][i]=0,{i,1,10}];

Table[d2R[-1][i]=0,{i,1,10}];
Table[d2R[+1][i]=0,{i,1,10}];

Table[RInhom[-1][fields[[j]]]=-Sum[c[[i]] \[CapitalPhi]0[[j,i]],{i,1,nf}],{j,1,nf}];
Table[RInhom[+1][fields[[j]]]=Sum[c[[i+nf]] \[CapitalPhi]0[[j,i+nf]],{i,1,nf}],{j,1,nf}];

Table[dR[-1][fields[[j]]]=-Sum[c[[i]] \[CapitalPhi]0[[j+nf,i]],{i,1,nf}],{j,1,nf}];
Table[dR[+1][fields[[j]]]=Sum[c[[i+nf]] \[CapitalPhi]0[[j+nf,i+nf]],{i,1,nf}],{j,1,nf}];

Table[d2R[-1][j]=R[j]''[r]/.Solve[ode[j][l,m,\[Omega]],R[j]''[r]][[1]]/.{R[i_][r]->RInhom[-1][i],R[i_]'[r]->dR[-1][i]},{j,fields}];
Table[d2R[+1][j]=R[j]''[r]/.Solve[ode[j][l,m,\[Omega]],R[j]''[r]][[1]]/.{R[i_][r]->RInhom[+1][i],R[i_]'[r]->dR[+1][i]},{j,fields}];

(*f=1-2M/r;*)

Clear[R];
greplace={R[i_][r]->RInhom[io][i],R[i_]'[r]->dR[io][i],r->r0,R[i_]''[r]:>d2R[io][i]};

If[m!=0,
Table[
{
RInhom[io][4]=R[4][r]/.Solve[gauge[3][l,m,\[Omega]]==0,R[4][r]][[1]]//.greplace,
dR[io][4]=R[4]'[r]/.Solve[D[gauge[3][l,m,\[Omega]],r]==0,R[4]'[r]][[1]]//.greplace,

RInhom[io][2]=R[2][r]/.Solve[gauge[2][l,m,\[Omega]]==0,R[2][r]][[1]]//.greplace,
dR[io][2]=R[2]'[r]/.Solve[D[gauge[2][l,m,\[Omega]],r]==0,R[2]'[r]][[1]]//.greplace,

RInhom[io][8]=R[8][r]/.Solve[gauge[4][l,m,\[Omega]]==0,R[8][r]][[1]]//.greplace,
dR[io][8]=R[8]'[r]/.Solve[D[gauge[4][l,m,\[Omega]],r]==0,R[8]'[r]][[1]]//.greplace
},
{io,{-1,+1}}];
,
Table[
{
RInhom[io][6]=R[6][r]/.R6rep/.greplace,
dR[io][6]=R[6]'[r]/.Solve[D[gauge[2][l,0,0],r]==0,R[6]'[r]][[1]]//.greplace,

RInhom[io][7]=R[7][r]/.R7rep/.greplace,
dR[io][7]=R[7]'[r]/.Solve[D[gauge[3][l,0,0],r]==0,R[7]'[r]][[1]]//.greplace
},
{io,{-1,+1}}
];
];

horizCoeffs=Total[Table[c[[i]]innerBCBasis[i][[-10;;]],{i,1,nf}]];

If[EvenQ[l+m]&&m==0,

(*rout^l factor below is because we premultiply by this factor earlier on *)
infCoeffs = rout^l Total[Table[c[[i+nf]]outerBCBasis[i][[-20;;-11]],{i,1,nf}]];
infLogCoeffs = rout^l Total[Table[c[[i+nf]]outerBCBasis[i][[-10;;]],{i,1,nf}]];
,
infCoeffs = Total[Table[c[[i+nf]]outerBCBasis[i][[-10;;]],{i,1,nf}]];
];

replace={\[Omega]1->\[Omega],a[i_][n_]:>infCoeffs[[i,n+1]],al[i_][n_]:>infLogCoeffs[[i,n+1]],b[i_][n_]:>horizCoeffs[[i,n+1]],l1->l};

If[EvenQ[l+m]&&m==0,
Table[
{

horizCoeffs[[6,i+1]]=b[6][i]/.GaugeStaticEvenReplacementH[6]/.replace;,
horizCoeffs[[7,i+1]]=b[7][i]/.GaugeStaticEvenReplacementH[7]/.replace;,


infCoeffs[[6,i+1]]=a[6][i]/.GaugeStaticEvenNonLogReplacement\[Infinity][6]/.replace;,
infCoeffs[[7,i+1]]=a[7][i]/.GaugeStaticEvenNonLogReplacement\[Infinity][7]/.replace;,

infLogCoeffs[[6,i+1]]=al[6][i]/.GaugeStaticEvenLogReplacement\[Infinity][6]/.replace;,
infLogCoeffs[[7,i+1]]=al[7][i]/.GaugeStaticEvenLogReplacement\[Infinity][7]/.replace;

}
,{i,0,7}],
Table[
{horizCoeffs[[2,i+1]]=b[2][i]/.GaugeExpansionReplacementH[2]/.replace;,
horizCoeffs[[4,i+1]]=b[4][i]/.GaugeExpansionReplacementH[4]/.replace;,
horizCoeffs[[8,i+1]]=b[8][i]/.GaugeExpansionReplacementH[8]/.replace;,

infCoeffs[[2,i+1]]=a[2][i]/.GaugeExpansionReplacement\[Infinity][2]/.replace;,
infCoeffs[[4,i+1]]=a[4][i]/.GaugeExpansionReplacement\[Infinity][4]/.replace;,
infCoeffs[[8,i+1]]=a[8][i]/.GaugeExpansionReplacement\[Infinity][8]/.replace;
}
,{i,0,7}]
];

Clear[R];
{Table[{RInhom[-1][i],dR[-1][i],dR[+1][i]},{i,1,10}],horizCoeffs,infCoeffs,infLogCoeffs}
]


res=MostCircLorenzModes[2,2,40,20];


Edot\[Infinity]=(m^2 \[CapitalOmega]\[Phi]^2)/(64 \[Pi] \[Lambda][l]l(l+1)) Abs[res[[-2,7,1]]]^2/.{l->2,m->2}


<<Teukolsky`


orbit=KerrGeoOrbit[0,10.`30,0,1];


modeT=TeukolskyPointParticleMode[-2,2,2,0,0,orbit];


1-modeT["Fluxes"]["FluxInf"]/Edot\[Infinity]


(*N[res[[2]],20]//MatrixForm*)


(*EdotEH=(\[Lambda][l] l(l+1))/(256\[Pi] M^2(1+16M^2m^2\[CapitalOmega]\[Phi]^2))Abs[hbar[1]+(1+4 I M m \[CapitalOmega]\[Phi])/(l(l+1))(hbar[5]-I hbar[9]+2I M m \[CapitalOmega]\[Phi] \[Lambda][l]^-1(hbar[7]-I hbar[10]))]^2;*)


(*2EdotEH/.{l\[Rule]2,m\[Rule]2,hbar[a_]\[RuleDelayed]res[[2,a,1]]}*)


lOddm0HorizSeriesCoeffs[l1_]:=Module[{series,x1},
series=Series[With[{l=l1},(x1 (1+x1)^2)/(l+2) D[(x1/(x1+1))^(1/2) LegendreP[l,1,3,1+2x1 ],{x1,2}]]/.x1->r2/(2M)-1/.r2->\[Rho]+2M,{\[Rho],0,10},Assumptions->r\[Element]Reals];
Table[SeriesCoefficient[series,n],{n,0,10}]
]


lOddm0InfSeriesCoeffs[l1_]:=Module[{series,x1},
series=Series[With[{l=l1},(x1 (1+x1)^2)/(l+2) D[(x1/(x1+1))^(1/2) LegendreQ[l,1,3,1+2x1 ],{x1,2}]]/.x1->r2/(2M)-1,{r2,\[Infinity],10},Assumptions->r\[Element]Reals];
Table[SeriesCoefficient[series,n],{n,0,10}]
]


CircLorenzOddStatic[l1_,r1_]:=Block[{l=l1,r=r1,h8HomIn,h8HomOut,J,\[CapitalPhi],res,Cs,horizCoeffs,A,infCoeffs,r2},

J[l_]:={0,-((32\[Pi] E0)/f0)\[CapitalOmega]\[Phi] D[SphericalHarmonicY[l,0,\[Theta],0],\[Theta]]/.\[Theta]->\[Pi]/2};

If[l>=3,
h8HomIn[l_,r2_]:=I (x1 (1+x1)^2)/(l+2) D[(x1/(x1+1))^(1/2) LegendreP[l,1,1+2x1 ],{x1,2}]/.x1->(r2/(2M)-1);
h8HomOut[l_,r2_]:=I (x1 (1+x1)^2)/(l+2) D[(x1/(x1+1))^(1/2) LegendreQ[l,1,1+2x1 ],{x1,2}]/.x1->(r2/(2M)-1);
,
h8HomIn[l_,r2_]:=r2^2;
h8HomOut[l_,r2_]:=1/r2;
];

\[CapitalPhi][l_,r_]:={
{-h8HomIn[l,r],h8HomOut[l,r]},
{-D[h8HomIn[l,r2],r2],D[h8HomOut[l,r2],r2]}
}/.r2->r;

Cs=Inverse[\[CapitalPhi][l,r0]] . J[l];

A=-(256/3)Sqrt[3\[Pi]] M^4 E0 \[CapitalOmega]\[Phi]/(r0-2M);

res=ConstantArray[0,{10,3}];
horizCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
infCoeffs=ConstantArray[{0,0,0,0,0,0,0,0},10];
If[l>=3,
horizCoeffs[[8]]=Re[Cs[[1]]]lOddm0HorizSeriesCoeffs[l];
infCoeffs[[8]]=Re[Cs[[2]]]lOddm0InfSeriesCoeffs[l];
,
horizCoeffs[[8]]=Re[Cs[[1]]]Table[SeriesCoefficient[Series[r2^2/.r2->\[Rho]+2M,{\[Rho],0,nMax}],n],{n,0,nMax}];
horizCoeffs[[9]]=Re[A Table[SeriesCoefficient[Series[r2^-2/.r2->\[Rho]+2M,{\[Rho],0,nMax}],n],{n,0,nMax}]];


infCoeffs[[8]]=Re[Cs[[2]]]Table[SeriesCoefficient[Series[r2^-1,{r2,\[Infinity],nMax}],n],{n,0,nMax}];
infCoeffs[[9]]=Re[A Table[SeriesCoefficient[Series[r2^-2,{r2,\[Infinity],nMax}],n],{n,0,nMax}]];

];


res[[8,1]]=Re[Cs[[1]] h8HomIn[l,r]];
res[[8,2]]=Re[Cs[[1]]D[h8HomIn[l,r2],r2]/.r2->r];
res[[8,3]]=Re[Cs[[2]]D[h8HomOut[l,r2],r2]/.r2->r];

(*Fix to make the odd dipole regular at the horizon, without changing the angular-momentum*)
If[l==1,
res[[9,1]]=A r^-2;
res[[9,2]]=-2A r^-3;
res[[9,3]]=-2A r^-3;
];

{res,horizCoeffs,infCoeffs}
]


(*MonopoleIrreg[]:=Block[{A=(2E0)/(3M r0 f0)(M-(r0-3M)Log[f0]),f=1-2M/r,P=(r^2+2M r+4M^2),Q=r^3-M r^2-2M^2r+12M^3,httl0In,hrrl0In,h\[Phi]\[Phi]l0In,h1l0In,h6l0In,h3l0In,httl0Out,hrrl0Out,h\[Phi]\[Phi]l0Out,h1l0Out,h6l0Out,h3l0Out,horizCoeffs,res},
httl0In=-((A f M)/r^3)P;
hrrl0In=A/(r^3f)Q;
h\[Phi]\[Phi]l0In=A f P;
h1l0In=2Sqrt[\[Pi]]r(httl0In+f^2hrrl0In);
h6l0In=2Sqrt[\[Pi]]r/f(httl0In-f^2hrrl0In);
h3l0In=4Sqrt[\[Pi]]1/rh\[Phi]\[Phi]l0In;

httl0Out=(2E0)/(3r^4r0 f0)(3r^3(r0-r)+M^2(r0^2-12M r0+8M^2)+(r0-3M)(-r M(r+4M)+r P f Log[f]+8M^3Log[r0/r]));
hrrl0Out=-((2E0)/(3r^4r0 f0 f^2))(-r^3r0-2M r(r0^2-6M r0-10M^2)+3M^2(r0^2-12M r0+8M^2)+(r0-3M)(5M r^2+r/M Q f Log[f]-8M^2(2r-3M)Log[r0/r]));
h\[Phi]\[Phi]l0Out=-((2E0)/(9r r0 f0))(3r0^2M-80M^2r0+156M^3+(r0-3M)(-3r^2-12M r+3r/M P f Log[f]+44M^2+24M^2Log[r0/r]));
h1l0Out=2Sqrt[\[Pi]]r(httl0Out+f^2hrrl0Out);
h6l0Out=2Sqrt[\[Pi]]r/f(httl0Out-f^2hrrl0Out);
h3l0Out=4Sqrt[\[Pi]]1/rh\[Phi]\[Phi]l0Out;

horizCoeffs=Transpose[{
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],0],
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],1],
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],2],
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],3],
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],4],
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],5],
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],6],
Coefficient[Normal[Series[{h1l0In,0,h3l0In,0,0,h6l0In,0,0,0,0}/.r\[Rule]\[Rho]+2M,{\[Rho],0,10}]],\[Rho],7]

}];

res=ConstantArray[{0,0,0},10];
res[[1]]={h1l0In,D[h1l0In,r],D[h1l0Out,r]}/.r\[Rule]r0;
res[[3]]={h3l0In,D[h3l0In,r],D[h3l0Out,r]}/.r\[Rule]r0;
res[[6]]={h6l0In,D[h6l0In,r],D[h6l0Out,r]}/.r\[Rule]r0;
{res,horizCoeffs}
]*)


MonopoleReg[]:=Block[{P,Q,W,K,r,H,f,R,\[CapitalPhi],J1,source,coeffs,hretIn,hretOut,res,horizCoeffs,infCoeffs,hretOutExp,infLogCoeffs,zeros},
f=1-2M/r;
P=r^2+2r M+4M^2;
Q=r^3-r^2 M-2r M^2+12M^3;
W=3r^3-r^2 M-4r M^2-28 M^3/3;
K=r^3 M-5r^2 M^2-20r M^3/3+28 M^4;

H[1]={-f,f^-1,1};
H[2]={-((f M)/r^3)P,f^-1/r^3 Q,f/r^2 P};
H[3]={-(M^4/r^4),(M^3 f^-2 (3M-2r))/r^4,M^3/r^3};
H[4]={
M/r^4 (W+r P f Log[f]-8M^3 Log[r/M]),
f^-2/r^4 (K-r Q f Log[f]-8M^3 (2r-3M)Log[r/M]),
1/r^3 (3r^3-W-r P f Log[f]+8M^3 Log[r/M])
};

Table[Simplify[R[1][i]=2Sqrt[\[Pi]] r/M (H[i][[1]]+f^2 H[i][[2]])],{i,1,4}];
Table[Simplify[R[3][i]=4Sqrt[\[Pi]]r/M H[i][[3]]],{i,1,4}];
Table[Simplify[R[6][i]=2Sqrt[\[Pi]]r/(f M)(H[i][[1]]-f^2 H[i][[2]])],{i,1,4}];

\[CapitalPhi]={
{-R[1][1],-R[1][2],R[1][3],R[1][4]},
{-R[3][1],-R[3][2],R[3][3],R[3][4]},
D[{-R[1][1],-R[1][2],R[1][3],R[1][4]},r],
D[{-R[3][1],-R[3][2],R[3][3],R[3][4]},r]
};

J1=-16 \[Pi] E0/r0 SphericalHarmonicY[0,0,\[Pi]/2,0];
source={0,0,J1,J1/f};

coeffs=Simplify[Inverse[\[CapitalPhi]] . source]/.r->r0;
(*To get the irregular solution use:*)
(*Table[hretIn[i]=(coeffs[[1]]+coeffs[[2]])R[i][2],{i,{1,3,6}}];
Table[hretOut[i]=-coeffs[[1]](R[i][1]-R[i][2])+coeffs[[3]]R[i][3]+coeffs[[4]]R[i][4],{i,{1,3,6}}];*)

(*To get the regular solution use:*)
Table[hretIn[i]=coeffs[[1]]R[i][1]+coeffs[[2]]R[i][2],{i,{1,3,6}}];
Table[hretOut[i]=coeffs[[3]]R[i][3]+coeffs[[4]]R[i][4],{i,{1,3,6}}];

horizCoeffs=Transpose[{
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],0],
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],1],
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],2],
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],3],
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],4],
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],5],
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],6],
Coefficient[Normal[Series[{hretIn[1],0,hretIn[3],0,0,hretIn[6],0,0,0,0}/.r->\[Rho]+2M,{\[Rho],0,10}]],\[Rho],7]
}];

hretOutExp[i_,lg_]:=Table[Coefficient[Coefficient[Normal[Series[hretOut[i],{r,\[Infinity],nMax}]],r,-n],Log[r],lg],{n,0,nMax}];
zeros=ConstantArray[0,10];
infCoeffs={hretOutExp[1,0],zeros,hretOutExp[3,0],zeros,zeros,hretOutExp[6,0],zeros,zeros,zeros,zeros};
infLogCoeffs={hretOutExp[1,1],zeros,hretOutExp[3,1],zeros,zeros,hretOutExp[6,1],zeros,zeros,zeros,zeros};

res=ConstantArray[{0,0,0},10];
res[[1]]={hretIn[1],D[hretIn[1],r],D[hretOut[1],r]}/.r->r0;
res[[3]]={hretIn[3],D[hretIn[3],r],D[hretOut[3],r]}/.r->r0;
res[[6]]={hretIn[6],D[hretIn[6],r],D[hretOut[6],r]}/.r->r0;

{res,horizCoeffs,infCoeffs,infLogCoeffs}
]


MonopoleReg[][[3,2]]


CircLorenz[l_,m_]:=Module[{res},
If[l==0&&m==0,
res=MonopoleReg[],
If[OddQ[l]&&m==0,
res=CircLorenzOddStatic[l,r0],
res=MostCircLorenzModes[l,m,40,20]
]];
res
]


(*test=CircLorenz[4,0];
1-test[[1,1]]/1.50635113857301083936655833685409076332835356030978220631302303633536637001233`40.*)


modes=Flatten[Table[{l,m},{l,lmin,lmax},{m,0,l}],1];
Length[modes]


Clear[res]
count=0;
SetSharedVariable[count]
SetSharedFunction[res]
Monitor[ParallelTable[Block[{$MaxExtraPrecision=100,l=modes[[i1,1]],m=modes[[i1,2]]},res[l,m]=CircLorenz[l,m]];count++;,{i1,1,Length[modes]}],ToString[count]<>"/"<>ToString[Length[modes]]];


output=Table[Block[{l=modes[[j,1]],m=modes[[j,2]]},Table[{l,m,i,Re[res[l,m][[1,i,1]]],Im[res[l,m][[1,i,1]]],Re[res[l,m][[1,i,2]]],Im[res[l,m][[1,i,2]]],Re[res[l,m][[1,i,3]]],Im[res[l,m][[1,i,3]]]},{i,1,10}]],{j,1,Length[modes]}];


outputHorizCoeffs=Table[Block[{l=modes[[k,1]],m=modes[[k,2]]},Table[
Flatten[{l,m,i,Table[N[{Re[res[l,m][[2,i,j]]],Im[res[l,m][[2,i,j]]]},30],{j,1,8}]}],{i,1,10}]],{k,1,Length[modes]}
];


outputInfCoeffs=Table[Block[{l=modes[[k,1]],m=modes[[k,2]]},Table[
Flatten[{l,m,i,Table[N[{Re[res[l,m][[3,i,j]]],Im[res[l,m][[3,i,j]]]},30],{j,1,8}]}],{i,1,10}]],{k,1,Length[modes]}
];


outputInfLogCoeffs=Table[Block[{l=modes[[k,1]],m=modes[[k,2]]},If[EvenQ[l+m]&&m==0,
Table[
Flatten[{l,m,i,Table[N[{Re[res[l,m][[4,i,j]]],Im[res[l,m][[4,i,j]]]},30],{j,1,8}]}],{i,1,10}],
##&[]
]],{k,1,Length[modes]}
];


r0ToString[r_Integer]:=ToString[r]
r0ToString[r_]:=ToString[N[r]]


(*SetDirectory["~/Temp/"];*)
(*Export["hlmi_r"<>r0ToString[r0]<>"_high_acc_l"<>ToString[lmin]<>"-"<>ToString[lmax]<>".dat",Flatten[N[output,30],1]]*)
SetDirectory["~/Code/SecondOrder/FirstOrderData/data/ret_data"];
Export["hlmi_r"<>r0ToString[r0]<>"_high_acc.dat",Flatten[N[output,30],1]]

SetDirectory["~/Code/SecondOrder/FirstOrderData/data/boundary_data/horiz_boundary_data"];
Export["hlmi_r"<>r0ToString[r0]<>"_high_acc_horiz_coeffs_reg_monopole.dat",Flatten[outputHorizCoeffs,1]]

SetDirectory["~/Code/SecondOrder/FirstOrderData/data/boundary_data/inf_boundary_data"];
Export["hlmi_r"<>r0ToString[r0]<>"_high_acc_inf_coeffs_reg_monopole.dat",Flatten[outputInfCoeffs,1]]
Export["hlmi_r"<>r0ToString[r0]<>"_high_acc_inf_log_coeffs_reg_monopole.dat",Flatten[outputInfLogCoeffs,1]]


Print[(AbsoluteTime[]-t1)/60]
