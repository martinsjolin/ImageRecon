

A_max = [-3.8268 -9.2388];
A_min = [3.8268 -9.2388];
P = [30 -10];

norm(A_min - P)

norm(A_max - P)

frac = norm(A_min - P)^2/norm(A_max - P)^2


%%

M = 5.97*10^24;
G = 6.67*10^(-11);


r = 39e3 + 6.38*10^6;
g_prim = G*M/r^2;
% g_prim = 9.82;
h = 0.5;
D = 0.001;
t = 1e-4/(pi*(D/2)^2*sqrt(2*g_prim*h))


