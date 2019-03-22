%
% EQUATION SOLVER TO CALCULATE THE GENERAL CONIC EQUATION FROM THE STANDARD
% ELLIPSE EQUATION
% 9 March 2015
%


clear all;

syms cx cy a b angle;

%x = cx + a * cos(angle) - b * sin(angle);
%y = cy + a * sin(angle) + b * cos(angle);

syms x y;
xcan = (x - cx)*cos(angle) + (y - cy)*sin(angle);
ycan = -(x - cx)*sin(angle) + (y - cy)*cos(angle);

eqt = xcan*xcan/(a*a) + ycan*ycan/(b*b) - 1;

% Part 2 - find coefficients from conic equation.
syms A C D E F;
conicEqt = A*x^2 + C*y^2 + D*x + E*y + F;

ellipseEqt = (x - cx)^2/(a^2) + (y - cy)^2/(b^2) - 1;

% compare coefficients
coeffs1 = coeffs(collect(conicEqt, [x, y]), [x, y]);
coeffs2 = coeffs(collect(ellipseEqt, [x, y]), [x, y]);

eqt1 = coeffs1(1) - coeffs2(1);
eqt2 = coeffs1(2) - coeffs2(2);
eqt3 = coeffs1(3) - coeffs2(3);
eqt4 = coeffs1(4) - coeffs2(4);
eqt5 = coeffs1(5) - coeffs2(5);
