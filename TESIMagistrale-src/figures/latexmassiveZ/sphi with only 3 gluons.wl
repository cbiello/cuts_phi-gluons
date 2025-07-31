(* ::Package:: *)

SetDirectory["/Users/christianbiello/Desktop/onshellamps"];
<< "SpinorHelicity.m";
<< "MomentumTwistors.m";
SetDirectory["/Users/christianbiello/Desktop/Mathnb/SBfiniteflow"];
<< "InitTwoLoopToolsFF.m";
<< "InitDiagramsFF.m";
(*<< "InitFORMpaths.m";*)
<< "setupfiles/Setup_phi4g.m";


(*Momentum twistor*)
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


A1ym[p1_,p2_,p3_,p4_,p5_]:=1/3/(spA[p1,p2]*spA[p2,p3]*spA[p3,p4]*spA[p4,p5]*spA[p5,p1])*(
trm[p1,p2,p3,p4]+trm[p1,p2,p3,p5]+trm[p1,p3,p4,p5]+trm[p1,p2,p4,p5]+trm[p2,p3,p4,p5]
)


A0phi[p1_,p2_,pH_]:=-spA[p1,p2]^2;


ordersp={spA[i_,i_]:>0,spB[i_,i_]:>0,
spA[2,1]:>-spA[1,2],spB[2,1]:>-spB[1,2],
spA[3,1]:>-spA[1,3],spB[3,1]:>-spB[1,3],
spA[4,1]:>-spA[1,4],spB[4,1]:>-spB[1,4],
spA[3,2]:>-spA[2,3],spB[3,2]:>-spB[2,3],
spA[4,2]:>-spA[2,4],spB[4,2]:>-spB[2,4],
spA[4,3]:>-spA[3,4],spB[4,3]:>-spB[3,4],
spA[l1,1]:>-spA[1,l1],spB[l1,1]:>-spB[1,l1],
spA[l1,2]:>-spA[2,l1],spB[l1,2]:>-spB[2,l1],
spA[l1,3]:>-spA[3,l1],spB[l1,3]:>-spB[3,l1],
spA[l1,4]:>-spA[4,l1],spB[l1,4]:>-spB[4,l1],
spA[l2,1]:>-spA[1,l2],spB[l2,1]:>-spB[1,l2],
spA[l2,2]:>-spA[2,l2],spB[l2,2]:>-spB[2,l2],
spA[l2,3]:>-spA[3,l2],spB[l2,3]:>-spB[3,l2],
spA[l2,4]:>-spA[4,l2],spB[l2,4]:>-spB[4,l2], spA[l2,l1]:>-spA[l1,l2]
};


cut = A0phi[l1,-l2,5]*A1ym[1,2,3,l1,-l2];
%/.spA[x_,-l2]:>-spA[x,l2]/.spA[-l2,x_]:>-spA[l2,x];
%/.trm[a_,b_,c_,-l2]:>-trm[a,b,c,l2];
%/.trm[1,2,3,l2]:>trm[1,2,3,l1]+trm[1,2,3,2];
%/.trm[1,3,l1,l2]:>-trm[l1,3,l1,l2]-trm[2,3,l1,l2];
cut1=%/.trexpandmassless/.ordersp//Expand;
%//SymbolForm


cut1/.spB[l2,l1]:>-s[4]/(spA[l1,l2]);
%/.spB[2,l1]:>(-spB[2,1]spA[1,l2]-spB[2,3]spA[3,l2])/spA[l1,l2]//Expand;
cut2=%/.spB[1,l2]spA[3,l2]:>spB[1,l1]spA[3,l1]+spB[1,2]spA[3,2]/.ordersp//Expand;
%//SymbolForm


cut2/.spA[l1,l2]/spA[a_,l1]/spA[b_,l2]:>((spAB[b,4,b]spAB[a,4,a]-s[4] spAB[a,b,a])/(spAB[a,l1,a]spAB[b,l2,b])-spAB[b,4,b]/spAB[b,l2,b]+spAB[a,4,a]/spAB[a,l1,a])/(2spA[a,b])//Expand;
%/.ordersp/.spA[l1,l2]spB[a_,l1]/spA[b_,l2]:>-s[4]*spB[b,a]/spAB[b,4,b]/.ordersp;
cut3=%/.spA[l1,l2]spB[a_,l2]/spA[b_,l1]:>-s[4]*spB[b,a]/spAB[b,4,b]/.ordersp//Expand;
%//SymbolForm


cut3/.x_/(spAB[1,l2,1]spAB[3,l1,3]):>INT[-x,{1,1,1,1}];
%/.x_/spAB[1,l2,1]:>INT[-x,{1,1,0,1}]/.x_/spAB[3,l1,3]:>INT[x,{1,0,1,1}];
%/.x_/spAB[1,4,1]:>INT[x/spAB[1,4,1],{1,0,1,0}];
%/.x_/spAB[3,4,3]:>INT[x/spAB[3,4,3],{1,0,1,0}];
cut4=%//.INT[x_,pow_]+INT[y_,pow_]:>INT[x+y,pow];
%//SymbolForm


bub=cut4/.INT[x_,{1,1,1,1}]->0/.INT[x_,{1,1,0,1}]->0/.INT[x_,{1,0,1,1}]->0/.INT[x_,{1,0,1,0}]->x;
%*spA[1,2]spA[2,3]spA[3,1]*3/s[4]/.ordersp//Expand;
%//.spA[a_,b_]spB[a_,b_]:>-s[a,b];
%/.spAB[1,4,1]:>s[1,4]-s[4];
bubred=%/.spAB[3,4,3]:>s[3,4]-s[4];
%//SymbolForm


(*phi unordered: numerical result*)
PS4numr1={PS4num[[3]],PS4num[[1]],PS4num[[2]],PS4num[[4]]};
PS4numr2={PS4num[[2]],PS4num[[3]],PS4num[[1]],PS4num[[4]]};
bubr[2]=bubred/. x_spB:>EvaluateSpinors[x,PS4numr1]/.x_spA:>EvaluateSpinors[x,PS4numr1]/.x_s:>EvaluateSpinors[x,PS4numr1];
bubr[3]=bubred/. x_spB:>EvaluateSpinors[x,PS4numr2]/.x_spA:>EvaluateSpinors[x,PS4numr2]/.x_s:>EvaluateSpinors[x,PS4numr2];
bubr[1]=bubred/. x_spB:>EvaluateSpinors[x,PS4num]/.x_spA:>EvaluateSpinors[x,PS4num]/.x_s:>EvaluateSpinors[x,PS4num];
Sum[bubr[i],{i,1,3}]
(*phi unordered: analytical result*)
PS4r1={PS4[[3]],PS4[[1]],PS4[[2]],PS4[[4]]};
PS4r2={PS4[[2]],PS4[[3]],PS4[[1]],PS4[[4]]};
bubr[2]=bubred/. x_spB:>EvaluateSpinors[x,PS4r1]/.x_spA:>EvaluateSpinors[x,PS4r1]/.x_s:>EvaluateSpinors[x,PS4r1];
bubr[3]=bubred/. x_spB:>EvaluateSpinors[x,PS4r2]/.x_spA:>EvaluateSpinors[x,PS4r2]/.x_s:>EvaluateSpinors[x,PS4r2];
bubr[1]=bubred/. x_spB:>EvaluateSpinors[x,PS4]/.x_spA:>EvaluateSpinors[x,PS4]/.x_s:>EvaluateSpinors[x,PS4];
Sum[bubr[i],{i,1,3}]//Simplify


2*s[4]/.x_s:>EvaluateSpinors[x,PS4num]
2*s[4]/.x_s:>EvaluateSpinors[x,PS4]



