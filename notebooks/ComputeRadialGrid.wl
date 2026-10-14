(* ::Package:: *)

(* Radial grid for h1Lorenz: set r0 and gridFile, then Get this script. Optional settings (unset: the defaults in brackets):
     gridDeltaR0  worldtube half-width; spacing 0.01 across r0 +- gridDeltaR0  [2 for r0 <= 20, 4 above]
     gridR1, gridR2  outside the worldtube the spacing is 0.1 up to gridR1, 0.5 up to gridR2, 4 up to 10^4  [100, 1000]
     gridRefine   k: k times finer spacings outside the worldtube  [1]
     gridRmax     extend the grid beyond 10^4 with spacing 8/gridRefine up to gridRmax  [10^4]
   With the defaults (and r0 + gridDeltaR0 < 100) the grid is the original construction, unchanged. *)
\[CapitalDelta]r0=If[NumericQ[gridDeltaR0],gridDeltaR0,If[r0<=20,2,4]];
gridOptR1=If[NumericQ[gridR1],gridR1,100];
gridOptR2=If[NumericQ[gridR2],gridR2,1000];
gridOptK=If[NumericQ[gridRefine],gridRefine,1];
gridOptRmax=If[NumericQ[gridRmax],gridRmax,10^4];


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
(* outside the worldtube. r0 + Delta r0 >= 100: Range[r0 + Delta r0, 100, 0.1] is empty, so the regions are selected by r instead *)
URegions = If[r0+\[CapitalDelta]r0<100,
  {Range[r0+\[CapitalDelta]r0,100,0.1][[2;;]], Range[100,1000,0.5][[2;;]], Range[1000,10^4,4][[2;;]]},
  {Select[Join[Range[r0+\[CapitalDelta]r0,100,0.1],Range[100,1000,0.5],Range[1000,10^4,4]],#>r0+\[CapitalDelta]r0+10^-6&]}];
(* other spacings outside the worldtube: rebuilt with exact steps *)
If[gridOptK!=1||gridOptR1!=100||gridOptR2!=1000,
  URegions = {Select[Join[Range[r0+\[CapitalDelta]r0,gridOptR1,1/(10gridOptK)],Range[gridOptR1,gridOptR2,1/(2gridOptK)],Range[gridOptR2,10^4,4/gridOptK]],N[#]>r0+\[CapitalDelta]r0+10^-6&]}];


grid=Union@N[Join[V1Region,V2Region,TRegion,Sequence@@URegions]];
If[gridOptRmax>10^4,grid=Union[grid,N[Range[10^4,gridOptRmax,8/gridOptK][[2;;]]]]];


importantIndexes=Flatten[{Nearest[grid->"Index",r0-\[CapitalDelta]r0],Nearest[grid->"Index",r0],Nearest[grid->"Index",r0+\[CapitalDelta]r0]}];


If[Length[importantIndexes]=!=3,Abort[]];


Nr0[r_Integer]:=r
Nr0[r_]:=N[r]


Export[gridFile,{grid,importantIndexes},{"Datasets",{"r","ImportantIndexes"}}]
