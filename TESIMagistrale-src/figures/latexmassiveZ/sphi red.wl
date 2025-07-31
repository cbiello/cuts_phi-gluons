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
PS6 = GetPSfromMT[GetZ[6]];
PSnum=GetPSfromMT[GetZ[6]]/. GetRandomPoint[6];
ReduceVars = Solve[{
  spA[PS6[[2]],PS6[[6]]]==0,
  spB[PS6[[2]],PS6[[6]]]==0
},{ex[4],ex[6]}][[1]];
ReLabelVars = {ex[5]->ex[4], ex[7]->ex[5], ex[8]->ex[6]};
PS6 = PS6 /. ReduceVars /. ReLabelVars;
PH = Plus@@(ToMOM4/@{PS6[[5]],PS6[[6]]}) // MOM4SUM // Simplify;
PS5 = Join[PS6[[1;;4]],{PH}]; (*6dofs=5+1mass*)
pspt =Join[GetRandomPoint[5],{ex[6]->RandomInteger[{-10^2,10^2},1][[1]]/10^2}];
PS5num=PS5/. pspt;
mH=Sqrt[s[PS5num[[5]]]]
ToMOM4/@PS5;
MOM4SUM[Plus@@%]//Simplify


A1ym[p1_,p2_,p3_,p4_,p5_,p6_]:=1/3/(spA[p1,p2]*spA[p2,p3]*spA[p3,p4]*spA[p4,p5]*spA[p5,p6]spA[p6,p1])*
(trm[p1,p2,p3,p4]+trm[p1,p2,p3,p5]+trm[p1,p2,p4,p5]+trm[p1,p3,p4,p5]+trm[p2,p3,p4,p5]+
trm[p1,p2,p3,p6]+trm[p1,p2,p4,p6]+trm[p1,p3,p4,p6]+trm[p2,p3,p4,p6]+
trm[p1,p2,p5,p6]+trm[p1,p3,p5,p6]+trm[p2,p3,p5,p6]+
trm[p1,p4,p5,p6]+trm[p2,p4,p5,p6]+
trm[p3,p4,p5,p6]);


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


cut = A0phi[l1,-l2,5]*A1ym[1,2,3,4,l1,-l2];
%/.spA[x_,-l2]:>-spA[x,l2]/.spA[-l2,x_]:>-spA[l2,x];
%/.trm[a_,b_,c_,-l2]:>-trm[a,b,c,l2];
%/.trm[3,4,l1,l2]:>-trm[l1,4,l1,l2]-trm[1,4,l1,l2]-trm[2,4,l1,l2];
%/.trm[1,2,l1,l2]:>trm[1,l2,l1,l2]-trm[1,3,l1,l2]-trm[1,4,l1,l2];
%/.trm[1,4,l1,l2]:>trm[1,4,l1,2]+trm[1,4,l1,3]+trm[1,4,l1,4];
%/.trm[2,3,l1,l2]:>trm[2,3,l1,1]+trm[2,3,l1,3]+trm[2,3,l1,4]//Expand;
cut1=%/.trexpandmassless/.ordersp//Expand;
%//SymbolForm


cut1/.spB[l2,l1]:>-s[5]/spA[l1,l2];
%/.spA[3,l1]:>(spA[1,3]spA[4,l1]+spA[3,4]spA[1,l1])/spA[1,4]//Expand;
%/.spA[3,l2]:>(spA[1,3]spA[4,l2]+spA[3,4]spA[1,l2])/spA[1,4]//Expand;
%/.spA[2,l1]:>(spA[1,2]spA[4,l1]+spA[2,4]spA[1,l1])/spA[1,4]//Expand;
%/.spA[2,l2]:>(spA[1,2]spA[4,l2]+spA[2,4]spA[1,l2])/spA[1,4]//Expand;
%/.spA[4,l2]spB[x_,l2]:>spA[4,l1]spB[x,l1]+spA[4,1]spB[x,1]+spA[4,2]spB[x,2]+spA[4,3]spB[x,3]//Expand;
cut2=%/.spA[1,l1]spB[x_,l1]:>spA[1,l2]spB[x,l2]-spA[1,2]spB[x,2]-spA[1,3]spB[x,3]-spA[1,4]spB[x,4]/.ordersp//Expand;
%//SymbolForm


cut2/.spA[l1,l2]/spA[a_,l1]/spA[b_,l2]:>((spAB[b,5,b]spAB[a,5,a]-s[5] spAB[a,b,a])/(spAB[a,l1,a]spAB[b,l2,b])-spAB[b,5,b]/spAB[b,l2,b]+spAB[a,5,a]/spAB[a,l1,a])/(2spA[a,b])//Expand;
%/.spA[l1,a_]:>-spA[a,l1];
cut3=%/.spA[l2,a_]:>-spA[a,l2];
%//SymbolForm


cut3/.ordersp/.spA[l1,l2]spB[a_,l1]/spA[b_,l2]:>-s[5]*spB[b,a]/spAB[b,5,b]/.ordersp;
%/.spA[l1,l2]spB[a_,l2]/spA[b_,l1]:>-s[5]*spB[b,a]/spAB[b,5,b]/.ordersp//Expand;
%/.x_/(spAB[1,l2,1]spAB[4,l1,4]):>INT[-x,{1,1,1,1}];
%/.x_/spAB[1,l2,1]:>INT[-x,{1,1,0,1}]/.x_/spAB[4,l1,4]:>INT[x,{1,0,1,1}];
%/.x_/spAB[1,5,1]:>INT[x/spAB[1,5,1],{1,0,1,0}];
%/.x_/spAB[4,5,4]:>INT[x/spAB[4,5,4],{1,0,1,0}];
cut4=%//.INT[x_,pow_]+INT[y_,pow_]:>INT[x+y,pow];
%//SymbolForm


cut4/.spAB[x_,5,x_]:>-spA[x,2]spB[2,x]-spA[x,3]spB[3,x]-spA[x,4]spB[4,x]-spA[x,1]spB[1,x];
%/.spAB[x_,1,x_]:>spA[x,1]spB[1,x];
%/.spAB[x_,2,x_]:>spA[x,2]spB[2,x];
%/.spAB[x_,3,x_]:>spA[x,3]spB[3,x];
%/.spAB[x_,4,x_]:>spA[x,4]spB[4,x];
cutf=%/.ordersp//Expand;
Amp/. x_spB:>EvaluateSpinors[x,PS5num];
%/.x_spA:>EvaluateSpinors[x,PS5num];
%/.x_s:>EvaluateSpinors[x,PS5num]


Box1check=INT[1/3*spB[2,3]^2/(spA[1,4]*spA[1,4])*(-s[4,5]s[5,1]+s[5]s[3,2])/2,{1,1,1,1}];
%/. {x_spB:>EvaluateSpinors[x,PS5num],x_spA:>EvaluateSpinors[x,PS5num],x_s:>EvaluateSpinors[x,PS5num]}


(*Bubble part*)
bub=cut4/.INT[x_,{1,1,1,1}]->0/.INT[x_,{1,1,0,1}]->0/.INT[x_,{1,0,1,1}]->0/.INT[x_,{1,0,1,0}]->x;
%*spA[1,2]spA[2,3]spA[3,4]spA[4,1]/s[5]*3/.ordersp//Expand;
%//.spA[a_,b_]spB[a_,b_]:>-s[a,b];
%/.spAB[1,5,1]:>s[1,5]-s[5];
%/.spAB[4,5,4]:>s[4,5]-s[5];
bubred=%/.spA[1,4]^2spB[1,4]^2:>s[1,4]^2;
%//SymbolForm


(*phi unordered: numerical result*)
PS5numr1={PS5num[[4]],PS5num[[1]],PS5num[[2]],PS5num[[3]],PS5num[[5]]};
PS5numr2={PS5num[[3]],PS5num[[4]],PS5num[[1]],PS5num[[2]],PS5num[[5]]};
PS5numr3={PS5num[[2]],PS5num[[3]],PS5num[[4]],PS5num[[1]],PS5num[[5]]};
bubr[2]=bubred/. x_spB:>EvaluateSpinors[x,PS5numr1]/.x_spA:>EvaluateSpinors[x,PS5numr1]/.x_s:>EvaluateSpinors[x,PS5numr1];
bubr[3]=bubred/. x_spB:>EvaluateSpinors[x,PS5numr2]/.x_spA:>EvaluateSpinors[x,PS5numr2]/.x_s:>EvaluateSpinors[x,PS5numr2];
bubr[4]=bubred/. x_spB:>EvaluateSpinors[x,PS5numr3]/.x_spA:>EvaluateSpinors[x,PS5numr3]/.x_s:>EvaluateSpinors[x,PS5numr3];
bubr[1]=bubred/. x_spB:>EvaluateSpinors[x,PS5num]/.x_spA:>EvaluateSpinors[x,PS5num]/.x_s:>EvaluateSpinors[x,PS5num];
Sum[bubr[i],{i,1,4}]
(*phi unordered: analytical result*)
PS5r1={PS5[[4]],PS5[[1]],PS5[[2]],PS5[[3]],PS5[[5]]};
PS5r2={PS5[[3]],PS5[[4]],PS5[[1]],PS5[[2]],PS5[[5]]};
PS5r3={PS5[[2]],PS5[[3]],PS5[[4]],PS5[[1]],PS5[[5]]};
bubr[2]=bubred/. x_spB:>EvaluateSpinors[x,PS5r1]/.x_spA:>EvaluateSpinors[x,PS5r1]/.x_s:>EvaluateSpinors[x,PS5r1];
bubr[3]=bubred/. x_spB:>EvaluateSpinors[x,PS5r2]/.x_spA:>EvaluateSpinors[x,PS5r2]/.x_s:>EvaluateSpinors[x,PS5r2];
bubr[4]=bubred/. x_spB:>EvaluateSpinors[x,PS5r3]/.x_spA:>EvaluateSpinors[x,PS5r3]/.x_s:>EvaluateSpinors[x,PS5r3];
bubr[1]=bubred/. x_spB:>EvaluateSpinors[x,PS5]/.x_spA:>EvaluateSpinors[x,PS5]/.x_s:>EvaluateSpinors[x,PS5];
Sum[bubr[i],{i,1,4}]//Simplify


2s[5]/.x_s:>EvaluateSpinors[x,PS5]


2*s[5]/.x_s:>EvaluateSpinors[x,PS5num]
