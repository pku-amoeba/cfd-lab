(* ::Package:: *)

(*Mach cone visualization*)
(*A disturbance source is moving at the speed of u (Ma).
Frontwave surfaces activated at t=1,2,3 are ploted, respectively.
By manipulating the value of u, a series of relative positions of the source and the frontwaves are determined.
It is shown that the cone is to be outlined by the envelope of the frontwaves whenever u>1 Ma. 

Herein, we suppose the disturbance moves at constant speed 1 in the medium so that 
the frontwaves could be described by a series of circles. In other words, 
we gives the analytic solutions of the wave equations directly. *)

Manipulate[ContourPlot[{(x-2*Ma)^2+y^2==1,(x-1*Ma)^2+y^2==4,(x-0*Ma)^2+y^2==9},{x,-10,10},{y,-10,10},Axes->True],{Ma,0,3}]
