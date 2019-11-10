clear
clc

%variables
sdpvar x0 x1 x2 x3;
x=[x0;x1;x2;x3];

%inequalities 
g1 = x0^2 + x1^2 + x2^2 + x3^2 -2*x0*x2-2*x1*x3-1;
g=[g1];

%equalities
h1 = x0^2 + x1^2 - 1;
h2 = x2^2 + x3^2 - 1;
h=[h1;h2];


k=0; %the order of relaxation 

PolySys18(x,g,h,k)