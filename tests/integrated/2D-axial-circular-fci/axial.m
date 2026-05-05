(* ::Package:: *)

BeginPackage["axial`"]

(* To run the package, execute the following line in a notebook
<<axial.m
 *)

(* 
"Provides MMS solution and source for the anisotropic diffusion model in circular geometry
- Only Dirichlet boundary conditions in the perpendicular direction is considered
- Parallel boundaries (Penalisation) are not taken into account 
- MMS solution is firstly is prescribed in cylindrical coordinates 
(r,p,z) and afterwards a transformed to Cartesian coordinates (x,y,phi) used in GRILLIX 
r = sqrt(x^2+y^2);  p = arctan(y/x);  z = phi, t = time
- Paramters of model are: chi_par, chi_perp, safety factor q, and limiting flux surfaces rmin and rmax"
*)


Print["Computing MMS Terms"];

(*
define x,y and parallel derivatives in terms of r,p,z derivative \
and normalised radial coordinate rn
*)
absb[x_] = Sqrt[1 + x^2/q[x]^2];
pgrad[f_, x_, z_, y_, t_] = (D[f[x,z,y,t],y] + 1/q[x]*D[f[x,z,y,t],z])/absb[x];
ddx[f_, r_, p_, z_, t_] = D[f[r, p, z, t], r];
ddy[f_, r_, p_, z_, t_] = D[f[r, p, z, t], r]*Sin[p] + D[f[r, p, z, t], p]*Cos[p]/r;
d2dx2[f_, r_, p_, z_, t_] = D[ddx[f, r, p, z, t], r];
d2dy2[f_, r_, p_, z_, t_] = D[ddy[f, r, p, z, t], r]*Sin[p] +  D[ddy[f, r, p, z, t], p]*Cos[p]/r;


LaplacePerpMmsSol[f_,r_, p_, z_, t_] =  d2dx2[f, r, p, z, t] + d2dy2[f, r, p, z, t];

LaplacePerp[f_,r_,p_,z_,t_] = D[f[r,p,z,t],{r,2}] + D[f[r,p,z,t],r]/r + 1.0/(r*r)*D[f[r,p,z,t],{p,2}];


  
xn[x_] = (x-xmin)/(xmax-xmin);  
    
d2dpar2[f_, x_, z_, y_, t_] = (D[D[f[x, z, y, t], y], y] + 2/q[x]*D[D[f[x, z, y, t], y], z] +  1/q[x]^2*D[D[f[x, z, y, t], z], z])/absb[x]^2;
    


(*
Define normalised rho and MMS solution in terms of mode numbers \
given above
*)
MmsDens[x_, z_, y_, t_] = amp*Sin[2.0*Pi*kx*xn[x]]*Sin[kz*z - phz]*Cos[ky*y- phy]+offset;
MmsUpar[x_, z_, y_, t_]=1;
rhos = 0.0002284697436697996
Bnorm = 1.0
qe = 1.60217663*^-19
Me = 9.1093837*^-31
e0 = 8.85418781*^-12
Omegaci = qe * Bnorm / (1836.0*Me)


pflux[x_, z_, y_, t_]=MmsDens[x, z, y, t];
(*Smms[x_, z_, y_, t_]=D[MmsDens[x,z,y,t],t]-d2dpar2[MmsDens,x,z,y,t]-Laplaceperpe[MmsDens,x,z,y,t];*)
Smms[x_, z_, y_, t_]=D[MmsDens[x,z,y,t],t]-Dperp/(rhos*rhos*Omegaci) * LaplacePerp[MmsDens,x,z,y,t] * (rhos*rhos)-
					Dpar/(rhos*rhos*Omegaci)*d2dpar2[MmsDens,x,z,y,t]* (rhos*rhos);


Print["Finished MMS Terms"];
Print[Omegaci]
(* Set a dummy return value *)
1


EndPackage[ ]
