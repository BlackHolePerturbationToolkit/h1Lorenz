(* ::Package:: *)

\[CapitalDelta]r0=If[r0<=20,2,4];


M=1;


V1Boundary = 2+1/2M;
rin = 2+10^-5;


rs[r_]:=r+2 Log[-1+r/2]


rOfrs[rs1_]:=r/.Quiet@FindRoot[rs[r]==rs1,{r,2+10^-10,100},WorkingPrecision->30]


rsValues=Subdivide[rs[rin],rs[V1Boundary],500];


rValues=rOfrs/@rsValues;


V1Region = rValues;
V2Region = Range[V1Boundary,r0-\[CapitalDelta]r0,0.01][[2;;]];
TRegion  = Range[r0-\[CapitalDelta]r0,r0+\[CapitalDelta]r0,0.01][[2;;]];
U1Region = Range[r0+\[CapitalDelta]r0,100,0.1][[2;;]];
U2Region = Range[100,1000,0.5][[2;;]];
U3Region = Range[1000,10^4,4][[2;;]];


grid=Union@N[Join[V1Region,V2Region,TRegion,U1Region,U2Region,U3Region]];


importantIndexes=Flatten[{Nearest[grid->"Index",r0-\[CapitalDelta]r0],Nearest[grid->"Index",r0],Nearest[grid->"Index",r0+\[CapitalDelta]r0]}];


If[Length[importantIndexes]=!=3,Abort[]];


Nr0[r_Integer]:=r
Nr0[r_]:=N[r]


Export[gridFile,{grid,importantIndexes},{"Datasets",{"r","ImportantIndexes"}}]
