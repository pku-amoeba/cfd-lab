(* ::Package:: *)

(*Nonlinear physics course, Spring semester, 2020 *)
(*Plot bifurcation diagram for given 1-D dynamics system *)

(*dx/dt=rx-x/(1+x)*)
(*Phase space*)
Plot[{x/(1+x),0.1*x,x,20*x,-2*x},{x,-4,2}]
(*Bifurcation diagram*)
g1=Plot[-1+(1/r),{r,0.1,1},PlotStyle->Dashed,PlotRange->All];
g2=Plot[-1+(1/r),{r,1,5}];
g3=Plot[-1+(1/r),{r,-2,-0.01}];
g4=Plot[0*r,{r,-2,1},PlotStyle->Thick];
g5=Plot[0*r,{r,1,5},PlotStyle->Dashed];
Show[g1,g2,g3,g4,g5]


(*dx/dt=rx-x/(1+x^2)*)
(*Phase space*)
Plot[{x/(1+x^2),0.2*x,-0.2*x},{x,-4,4}]
(*Bifurcation diagram*)
h1=Plot[Sqrt[-1+1/r],{r,0.1,1},PlotStyle->Dashed,PlotRange->All];
h2=Plot[-Sqrt[-1+1/r],{r,0.1,1},PlotStyle->Dashed];
h3=Plot[0*r,{r,-1,1}];
h4=Plot[0*r,{r,1,2},PlotStyle->Dashed];
Show[h1,h2,h3,h4]


(*dx/dt=rx+x^3/(1+x^2)*)
(*Phase space*)
Plot[{-x^3/(1+x^2),-0.1*x,0.1*x,0*x,-x},{x,-1,1}]
(*Bifurcation diagram*)
i1=Plot[Sqrt[-r/(r+1)],{r,-0.9,0},PlotStyle->Dashed,PlotRange->All];
i2=Plot[-Sqrt[-r/(r+1)],{r,-0.9,0},PlotStyle->Dashed];
i3=Plot[0*r,{r,0,1},PlotStyle->Dashed];
Show[i1,i2,i3]


(*dx/dt=5-re^(-x^2)*)
(*Phase space*)
j1=Plot[Log[r/5],{r,5,6},PlotStyle->Dashed,PlotRange->All];
j2=Plot[-Log[r/5],{r,5,6}];
j3=Plot[0*r,{r,0,5}];
Show[j1,j2]


(*dx/dt=rx+ax^2-x^3*)
Manipulate[Plot[r*x+x^2-x^3,{x,0,1.5}],{r,-5,0}]


Manipulate[Plot[r*x+x^2-x^3,{x,-20,20}],{r,0,200}]


Manipulate[Plot[r*x-x^2-x^3,{x,-1.5,0}],{r,-5,0}]


Manipulate[Plot[r*x-x^2-x^3,{x,-20,20}],{r,0,200}]


(*Plot in parameter space*)
v1=Plot[Sqrt[-4*r],{r,-1,0},PlotRange->All,PlotStyle->Thick];
v2=Plot[-Sqrt[-4*r],{r,-1,0},PlotStyle->Thick];
v3=Plot[0*r,{r,0,1},PlotStyle->Gray];
Show[v1,v2,v3]


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)


(* ::InheritFromParent:: *)
(**)
