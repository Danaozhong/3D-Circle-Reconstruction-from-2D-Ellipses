%
% EQUATION SOLVER TO UNPROJECT THE TRAJECTORY OF A 2D IMAGE POINT TO A 3D 
% POINT
% This program searches for a 3D point based on its projection on an
% image, assuming that the motion the point does is circular.
% 13 November 2014
%


clear all;



%
%
% PART I: Set up the equations
%
%

% parameters of the intrinsic camera parameters
syms f fx fy cx cy;

% the known extrinsic parameters
syms tx ty tz;
syms r13 r23 r33

% define the properties of the circle we are looking for. Assume the
% circle's plane to be on the xy plane of the local coordinate
% system. As we assume for all allowed circles that their center is at x= 0
% and y = 0, the only parameters we are looking for is the radius and the
% orientation of the z axis of the local coordinate system.
syms R h;


% define the projection matrix. For simplification, the projection
% parameters were shifted and scaled so that f = fx = fy and cx = cy = 0
MIntrinsic = [f 0 0; 0 f 0; 0 0 1];
MIntrinsicOriginal = [fx 0 cx; 0 fy cy; 0 0 1];


% Idea: values of r11 r21 r31 and r12 r22 r32 are not important as the
% object is rotation symmetric.
% Note that the length of the direction vectors is NOT always 1, thus the X
% and Z parameters must be corrected. ´
MExtrinsic = [-r33 -r13*r23 r13 tx; 0 (r13^2+r33^2) r23 ty; r13 -r23*r33 r33 tz];



syms x y;
% This represents the coordinate system of the 3D circle in the local
% object coordinate system
VCircle = [x; y; h; 1];

% Transform it to camera space...
VCircleInCameraCoordinates = MIntrinsic * MExtrinsic * VCircle;

%... and apply the perspective projection to receive image coordinates.
syms u v
eqU = VCircleInCameraCoordinates(1) / VCircleInCameraCoordinates(3) - u;
eqV = VCircleInCameraCoordinates(2) / VCircleInCameraCoordinates(3) - v;

% solve the perspective projection equations to get x and y based on the
% values of u and v
eq = solve(eqU, eqV);

x = eq(1).x;
y = eq(1).y;

% apply the circle condition. 
%eqt = x^2 + y^2 - R^2;
eqtCircle = x^2*(r33^2+r13^2) + y^2*((r13*r23)^2 + (r13^2 + r33^2)^2 + (r33*r23)^2) - R^2;



% transform to a polynomial of u and v by multiplying (ignoring) with the
% denominator. Just assume that it is nonzero.
[nom, den] = numden(eqtCircle);

% rearrange in terms of u and v
eqt = collect(nom, [u,v]);
% eqt now contains all the neccessary information for the projection. Its
% unknowns are u, v (the image point coordinates) and R, h (the 3D circle
% coordinats we are looking for.


% get the coefficients of the conic equation
uvCoefficients = coeffs(eqt, [u, v]);

% Assign the coefficients to the variables (A = u^2, B = uv, C = v^2, D =
% u, E = v, F = rest)
A = uvCoefficients(6);
B = uvCoefficients(5);
C = uvCoefficients(3);
D = uvCoefficients(4);
E = uvCoefficients(2);
F = uvCoefficients(1);

conicCoeffs = [A; B; C; D; E; F];

% calculate the discriminant of the conic equation
%disc = 4*A*C - B*B;


sqrRCoeffs = sym(1:6);
sqrHCoeffs = sym(1:6);
hCoeffs = sym(1:6);
otherCoeffs = sym(1:6);

test = 0;
for i = 1:6
    currentFunc = collect(conicCoeffs(i), [R, h]);
    [currentCoeffs, coeffTerm] = coeffs(currentFunc, [R, h]);
    
    % set all current coefficients to zero
    
    sqrRCoeffs(i) = 0;
    sqrHCoeffs(i) = 0;
    hCoeffs(i) = 0;
    otherCoeffs(i) = 0;

    
    for k = 1:length(coeffTerm)
        if(coeffTerm(k) == R^2)
            sqrRCoeffs(i) = currentCoeffs(k);
        elseif(coeffTerm(k) == h^2)
            sqrHCoeffs(i) = currentCoeffs(k);
        elseif(coeffTerm(k) == h)
            hCoeffs(i) = currentCoeffs(k);
        elseif(coeffTerm(k) == 1)
            otherCoeffs(i) = currentCoeffs(k);
        else
            throw(ImplementationException);
        end      
    end
end


% Part one is now finished. The values of the coefficients are now known
% and we can continue to solve the conic equation by substituting the
% coefficients by variables.
sSqrRCoeffs = sym('sSqrRCoeffs', [6 1]);
sSqrHCoeffs = sym('sSqrHCoeffs', [6 1]);
shCoeffs = sym('shCoeffs', [6 1]);
sOtherCoeffs = sym('sOtherCoeffs', [6 1]);

sConicCoeffs = sym(1:6);

for i = 1:6
    sConicCoeffs(i) = sSqrRCoeffs(i)*R^2 + sSqrHCoeffs(i)*h^2 + shCoeffs(i)*h + sOtherCoeffs(i);
end

sA = sConicCoeffs(1);
sB = sConicCoeffs(2);
sC = sConicCoeffs(3);
sD = sConicCoeffs(4);
sE = sConicCoeffs(5);
sF = sConicCoeffs(6);

sEqt = sA*u^2 + sB*u*v + sC*v^2 + sD*u + sE*v + sF;


syms a1 a2 a3;
sEqt = R^2 + a1*h^2 + a2*h + a3;

lsEqt = sEqt * sEqt;

diffEqts = sym(1:2);

diffEqts(1) = diff(lsEqt, R);
diffEqts(2) = diff(lsEqt, h);

% solve the differential equations
syms sumSqra1 sumSqra2 sumSqra3
syms suma1a2 suma1a3 suma2a3
syms suma1 suma2 suma3 N
for i = 1:length(diffEqts)
    diffEqts(i) = collect(diffEqts(i), [a1 a2 a3]);
    % as usual, check the coefficients and replace them by sums
    [currentCoeffs, coeffTerm] = coeffs(diffEqts(i), [a1 a2 a3]);
    
    currentEquation = 0;
    
    for k = 1:length(coeffTerm)
        if(coeffTerm(k) == a1^2)
            currentEquation = currentEquation + currentCoeffs(k) * sumSqra1;
        elseif(coeffTerm(k) == a2^2)
            currentEquation = currentEquation + currentCoeffs(k) * sumSqra2;
        elseif(coeffTerm(k) == a3^2)
            currentEquation = currentEquation + currentCoeffs(k) * sumSqra3;
        elseif(coeffTerm(k) == a1*a2)
            currentEquation = currentEquation + currentCoeffs(k) * suma1a2;
        elseif(coeffTerm(k) == a1*a3)
            currentEquation = currentEquation + currentCoeffs(k) * suma1a3;
        elseif(coeffTerm(k) == a2*a3)
            currentEquation = currentEquation + currentCoeffs(k) * suma2a3;
        elseif(coeffTerm(k) == a1)
            currentEquation = currentEquation + currentCoeffs(k) * suma1;
        elseif(coeffTerm(k) == a2)
            currentEquation = currentEquation + currentCoeffs(k) * suma2;
        elseif(coeffTerm(k) == a3)
            currentEquation = currentEquation + currentCoeffs(k) * suma3;
        elseif(coeffTerm(k) == 1)
            currentEquation = currentEquation + N * currentCoeffs(k);
        else
            throw(ImplementationException);
        end     
    end
    
    diffEqts(i) = currentEquation;
end

% simplyfy 
diffEqts(1) = diffEqts(1) / (4*R);

break


% divide by the coefficient of R^2 resuls in the equation:
a1 = cSqrRr13r23/cSqrRsqrr13;
a2 = cSqrRr13r33/cSqrRsqrr13;
a3 = cSqrRsqrr23/cSqrRsqrr13;


break

break

%
%
% PART II: Perform tests
%
%

f = 1.25;
fx = 1.25;
fy = 1.7;
cx = 0.51;
cy = 0.49;
tx = 1.2;
ty = 2.1;
tz = 10;

% rotation matrix with Euler anlges (90, 0, 0)
eulerX = 130 / 180*pi;
eulerY = 35 / 180*pi;
eulerZ = 60 / 180*pi;
matrix = makehgtform('zrotate', eulerZ)*makehgtform('yrotate', eulerY)*makehgtform('xrotate', eulerX);

r11 = matrix(1,1);
r21 = matrix(2,1);
r31 = matrix(3,1);

r12 = matrix(1,2);
r22 = matrix(2,2);
r32 = matrix(3,2);

r13 = matrix(1,3);
r23 = matrix(2,3);
r33 = matrix(3,3);

% evaluate all the coefficients for comparison purpose
sqrRCoeffs = eval(sqrRCoeffs);
sqrHCoeffs = eval(sqrHCoeffs);
hCoeffs = eval(hCoeffs);
otherCoeffs = eval(otherCoeffs);

% This can probably done in a better way, but I have no idea how... and I
% am lazy :(
sSqrRCoeffs1 = sqrRCoeffs(1);
sSqrHCoeffs1 = sqrHCoeffs(1);
sHCoeffs1 = hCoeffs(1);
sOtherCoeffs1 = otherCoeffs(1);

sSqrRCoeffs2 = sqrRCoeffs(2);
sSqrHCoeffs2 = sqrHCoeffs(2);
sHCoeffs2 = hCoeffs(2);
sOtherCoeffs2 = otherCoeffs(2);

sSqrRCoeffs3 = sqrRCoeffs(3);
sSqrHCoeffs3 = sqrHCoeffs(3);
sHCoeffs3 = hCoeffs(3);
sOtherCoeffs3 = otherCoeffs(3);

sSqrRCoeffs4 = sqrRCoeffs(4);
sSqrHCoeffs4 = sqrHCoeffs(4);
sHCoeffs4 = hCoeffs(4);
sOtherCoeffs4 = otherCoeffs(4);

sSqrRCoeffs5 = sqrRCoeffs(5);
sSqrHCoeffs5 = sqrHCoeffs(5);
sHCoeffs5 = hCoeffs(5);
sOtherCoeffs5 = otherCoeffs(5);

sSqrRCoeffs6 = sqrRCoeffs(6);
sSqrHCoeffs6 = sqrHCoeffs(6);
sHCoeffs6 = hCoeffs(6);
sOtherCoeffs6 = otherCoeffs(6);


% set some random image coordinates and evaluate the image-coordinate
% related values
u = 0.6;
v = -0.1;

eval(a1);
eval(a2);
eval(a3);

% test the projection settings
sourceR = 2.0;
sourceh = 3.0;

numOfPoints = 60;

% create the 3D data points (homogeneous coordinates)
points3D = cell(1, numOfPoints);
for i = 1:numOfPoints
    currentAngle = 2*pi*(i-1) / numOfPoints;
    X = cos(currentAngle)*sourceR;
    Y = sin(currentAngle)*sourceR;
    points3D{i} = [X; Y; sourceh; 1];
end

% create the 2D data points
MIntrinsicOriginal = eval(MIntrinsicOriginal);
MExtrinsicOriginal = eval(MExtrinsicOriginal);

points2DOriginal = cell(1, numOfPoints);
points2DAlternative = cell(1, numOfPoints);

for i = 1:numOfPoints
    currentPoint = MIntrinsicOriginal * MExtrinsicOriginal * points3D{i};
    currentu = currentPoint(1) / currentPoint(3);
    currentv = currentPoint(2) / currentPoint(3);
    points2DOriginal{i} = [currentu, currentv];
    
    % As the simplified transformation matrix does not have a unit length
    % of one, we must strech the X and Y parameters of the current3D point
    % to match a 3D circle shape
    current3DPoint = points3D{i};
    current3DPoint(1) = current3DPoint(1) / sqrt(r33^2+r13^2);
    current3DPoint(2) = current3DPoint(2) / sqrt((r13*r23)^2 + (r13^2 + r33^2)^2 + (r33*r23)^2);
    
    currentPoint = eval(MIntrinsic) * eval(MExtrinsic) * current3DPoint;
    currentu = currentPoint(1) / currentPoint(3);
    currentv = currentPoint(2) / currentPoint(3);
    points2DAlternative{i} = [currentu, currentv];
end

% recalculate the image points to fit them to the simplified projection
% matrix
points2D = cell(1, numOfPoints);
point2DNonSimplified = cell(1, numOfPoints);
for i = 1:numOfPoints
    currentu = points2DOriginal{i}(1) - cx;
    currentv = (points2DOriginal{i}(2) - cy)*fx/fy;
    points2D{i} = [currentu, currentv];
    
    % direcly calculate the points based on the simplified projection
    % matrix
    currentPoint = eval(MIntrinsic) * MExtrinsicOriginal * points3D{i};
    compareu = currentPoint(1) / currentPoint(3);
    comparev = currentPoint(2) / currentPoint(3);
    point2DNonSimplified{i} = [currentu, currentv];
    assert(abs(currentu - compareu) < 0.001);
    assert(abs(currentv - comparev) < 0.001);
end

%create a plot with all 2D image points
xvalues = zeros(1, numOfPoints);
yvalues = zeros(1, numOfPoints);
xvaluesNonSimplified = zeros(1, numOfPoints);
yvaluesNonSimplified = zeros(1, numOfPoints);
xvaluesAlternative = zeros(1, numOfPoints);
yvaluesAlternative = zeros(1, numOfPoints);

for i = 1:numOfPoints
    xvalues(i) = points2D{i}(1);
    yvalues(i) = points2D{i}(2);    
    
    xvaluesNonSimplified(i) = point2DNonSimplified{i}(1);
    yvaluesNonSimplified(i) = point2DNonSimplified{i}(2);
    
    xvaluesAlternative(i) = points2DAlternative{i}(1);
    yvaluesAlternative(i) = points2DAlternative{i}(2);
end

% create the plot. If everything is okay, the points from the projection
% with the simplified extrinsic and intrinsics should be the same as the
% ones originated from the non-simplified parameters.
figure();
hold on
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
scatter(xvalues, yvalues, [], [1 0 0], 'filled');
scatter(xvaluesNonSimplified, yvaluesNonSimplified, [], [0 0 1]);
scatter(xvaluesAlternative, yvaluesAlternative, [], [0 1 0]);
hold off

% calculate a, b, c and d (see paper for formulas)
sum_a1 = 0;
sum_a2 = 0;
sum_a3 = 0;
sumSqr_a1 = 0;
sumSqr_a2 = 0;
sum_a1a2 = 0;
sum_a1a3 = 0;
sum_a2a3 = 0;

for i = 1:numOfPoints
    u = points2D{i}(1);
    v = points2D{i}(2);
    
    current_a1 = eval(a1); 
    current_a2 = eval(a2);
    current_a3 = eval(a3);
    
    sum_a1 = sum_a1 + current_a1;
    sum_a2 = sum_a2 + current_a2;
    sum_a3 = sum_a3 + current_a3;
    sumSqr_a1 = sumSqr_a1 + current_a1^2;
    sumSqr_a2 = sumSqr_a2 + current_a2^2;
    
    sum_a1a2 = sum_a1a2 + current_a1*current_a2;
    sum_a1a3 = sum_a1a3 + current_a1*current_a3;
    sum_a2a3 = sum_a2a3 + current_a2*current_a3;
end

N = numOfPoints;

a = 2*N*sumSqr_a1 - 2*sum_a1^2;
b = 3*(N*sum_a1a2 - sum_a1*sum_a2);
c = N*(sumSqr_a2 + 2*sum_a1a3) - 2*sum_a1*sum_a3 - sum_a2^2;
d = N*sum_a2a3 - sum_a2*sum_a3;

partDiffEqt = [a b c d];
resultsh = roots(partDiffEqt);

% calculate the radius for each height value
resultsR = zeros(1, 3);
for i = 1:3
    currenth = resultsh(i);
   resultsR(i) = sqrt((-currenth*currenth*sum_a1 - currenth*sum_a2 - sum_a3)/N);  
end

% choose the value that minimizes the cost equation
lowestCost = 999999.999;
bestIndex = -1;

for i = 1:3
    h = resultsh(i);
    R = resultsR(i);
    
    % ensure that this result meets the conic equation to result in an
    % ellipse.
    if 0 >= 4*eval(A)*eval(C) - eval(B)^2
        continue;
    end
    
    cost = 0;
    for k = 1:N
        u = points2D{k}(1);
        v = points2D{k}(2);
    
        current_a1 = eval(a1);
        current_a2 = eval(a2);
        current_a3 = eval(a3);
        cost = cost + (R^2 + current_a1*h^2 + current_a2*h + current_a3)^2; 
    end
   
    if cost < lowestCost
       lowestCost = cost;
       bestIndex = i;
    end
end


if( -1 == bestIndex)
   error('No solution could be found!');
end
% print out the result
R = resultsR(bestIndex)
h = resultsh(bestIndex)


% check that circle equation is fulfilled
assert(abs(eval(eqtCircle)) < 1.0);

% calculate the x and y coordinates based on their projections.
xvalues = zeros(1, N);
yvalues = zeros(1, N);
zvalues = zeros(1, N);

for i = 1:N
    u = points2D{i}(1);
    v = points2D{i}(2);

    xvalues(i) = eval(x) * sqrt(r33^2+r13^2);
    yvalues(i) = eval(y) * sqrt(r33^2+r23^2);
end

% create a 3D plot of the reconstructed circle
figure();
scatter3(xvalues, yvalues, zvalues);