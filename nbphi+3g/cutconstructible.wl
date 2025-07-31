(* ::Package:: *)

(* ::Title:: *)
(*phi+3g cc*)


SetDirectory["/Users/christianbiello/Desktop/onshellamps"];
<< "SpinorHelicity.m";
<< "MomentumTwistors.m";
SetDirectory["/Users/christianbiello/Desktop/Mathnb/SBfiniteflow"];
<< "InitTwoLoopToolsFF.m";
<< "InitDiagramsFF.m";
<< "setupfiles/Setup_phi4g.m";


(* ::Subsubsection::Closed:: *)
(*Momentum twistor parametrization*)


GetRandomPoint[n_]:=Rule@@@Transpose@{ex/@Range[3*n-10],RandomInteger[{-10^2,10^2},3*n-10]/10^2};
PS5 = GetPSfromMT[GetZ[5]];
PSnum=GetPSfromMT[GetZ[5]]/. GetRandomPoint[5];
ReduceVars = Solve[{
  spA[PS5[[2]],PS5[[5]]]==0,
  spB[PS5[[2]],PS5[[5]]]==0
},{ex[2],ex[3]}][[1]];
ReLabelVars = {ex[4]->ex[2], ex[5]->ex[3]};
PS5 = PS5 /. ReduceVars /. ReLabelVars;
PH = Plus@@(ToMOM4/@{PS5[[4]],PS5[[5]]}) // MOM4SUM // Simplify;
PS4 = Join[PS5[[1;;3]],{PH}] (*3dofs=2+1mass*)
pspt =Join[GetRandomPoint[4],{ex[3]->RandomInteger[{-10^2,10^2},1][[1]]/10^2}];
PS4num=PS4/. pspt
s[PS4num[[4]]]
ToMOM4/@PS4;
MOM4SUM[Plus@@%]//Simplify


(* ::Subsubsection:: *)
(*Analytical expression for scalar integrals*)


I41m[s_,t_,m_]:=2/(s*t)*(((-s)^(-eps)+(-t)^(-eps)-(-m)^(-eps))/eps^2+Li1m[s,t,m])
Li1m[s_,t_,m_]:=-PolyLog[2,1-m/s]-PolyLog[2,1-m/t]-Log[s/t]^2/2-Pi^2/6
I32m[m1_,m2_]:=1/(m1-m2)*((-m1)^(-eps)-(-m2)^(-eps))/eps^2
I2[m_]:=(1/eps+2)(mu^2/(-m))^(eps)


(* ::Subsubsection:: *)
(*Result from cuts*)


(* ::Text:: *)
(*From boxes and triangles we have the following contribition.*)


d[x_,y_,z_]:=-(s[x,4]s[z,4])/2
c[x_,y_,z_]:=s[x,4]-s[4]


h[x_,y_,z_]:=d[x,y,z] I41m[s[x,4],s[z,4],s[4]]+c[x,y,z]I32m[s[x,4],s[4]]//Expand


Permute[p[1,2,3],CyclicGroup[3]]//.p[x_,y_,z_]:>h[x,y,z];
Plus@@%//Simplify;
Finite43=%+Sum[(-s[i,4])^(-eps)/eps^2,{i,1,3}]//Simplify (*subtraction of IR divergences*)


(* ::Text:: *)
(*From bubbles*)


Bub=Series[-(ds-2)/6*I2[s[4]],{eps,0,0}]
Finite2=Bub+(ds-2)/6/eps


(* ::Text:: *)
(*Summing together the contributions from the two sectors*)


ris=(Finite43+Finite2)


Amp1L=-2s[5]^2/spA[1,2]/spA[2,3]/spA[3,4]/spA[4,1];
finalres=ris*%


(* ::Text:: *)
(*"finalres" is the result from cut-constructible part and we can use it for a comparison with the complete computation.*)
