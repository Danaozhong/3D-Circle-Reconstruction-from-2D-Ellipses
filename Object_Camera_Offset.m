%
% EQUATION SOLVER TO TRANSLATE AN ELLIPSE
% 19 November 2014
%


%
%
% PART I: Set up the equations
%
%

clear all;


syms A B C D E F
syms At Bt Ct Dt Et Ft

syms u v
syms fx fy cx cy
ut = u + cx;
vt = v*fy/fx + cy;
conicEqt1 = A*u^2 + B*u*v + C*v^2 + D*u + E*v + F;
conicEqt2 = At*ut^2 + Bt*ut*vt + Ct*vt^2 + Dt*ut + Et*vt + Ft;

coeffs1 = coeffs(conicEqt1, [u, v]);
coeffs2 = coeffs(conicEqt2, [u, v]);

eqt1 = coeffs1(1) - coeffs2(1);
eqt2 = coeffs1(2) - coeffs2(2);
eqt3 = coeffs1(3) - coeffs2(3);
eqt4 = coeffs1(4) - coeffs2(4);
eqt5 = coeffs1(5) - coeffs2(5);
eqt6 = coeffs1(6) - coeffs2(6);

%solve ([eqt1 eqt2 eqt3 eqt4 eqt5 eqt6], [At, Bt, Ct, Dt, Et, Ft])

At = solve(eqt6, At)
Ct = solve(eqt3, Ct)
Bt = solve(eqt5, Bt)
Dt = solve(eval(eqt4), Dt)
Et = solve(eval(eqt2), Et)
Ft = simple(solve(eval(eqt1), Ft))
%solve(conicEqt1, conicEqt2, [At])