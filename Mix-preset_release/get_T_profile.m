function T=get_T_profile(u,U2,T2,Pr,M1,gamma)

T=(T2*(1-u)+(u-U2))/(1-U2) + sqrt(Pr)*0.5*(gamma-1)*M1^2*(1-u).*(u- U2);

end