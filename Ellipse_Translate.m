%
% EQUATION SOLVER TO TRANSLATE ANY OBJECT, GIVEN THE CAMERA OFFSET IN X
% AND Y COORDINATES
% 28 November 2014
%


%
%
% PART I: Set up the equations
%
%

clear all;


syms fx fy cx cy
MIntrinsic = [fx 0 cx; 0 fy cy; 0 0 1];

syms u v uT vT 
syms dx dy;
% create the screen points with a random z value
syms z
p = [u*z; v*z; z];
pT = [(u + dx)*z; (v + dy)*z; z];

% unproject to 3D space
pW = MIntrinsic\p;
pTW = MIntrinsic\pT;

% calculate rotation axis
rotAxis = cross(pW, pTW);

%divide by z^2 as no more relevant
rotAxis = simple(rotAxis / (z^2));


%calculate the rotation angle
cosAngle = pW.'*pTW / (sqrt(pW.'*pW) * sqrt(pTW.'*pTW));
cosAngle = simple(cosAngle);


break


%evaluate
fx = 1.69;
fy = 2.259;
cx = 0.5;
cy = 0.5;

u = 0.3;
v = 0.4;
dx = 0.07;
dy = 0.0;
z = 100;

cosAngle = eval(cosAngle);

angle = acos(cosAngle)
rotAxis = eval(rotAxis);
rotAxis = rotAxis / norm(rotAxis)
pW = eval(pW);
pTW = eval(pTW);
% verify code with matlab function
vrrotvec(pW, pTW)


% to validate, calculate the rotation matrix and see it the projected point
% matches the shifted image point

MRotation = vrrotvec2mat([rotAxis(1) rotAxis(2) rotAxis(3) angle]);
imagePoint = eval(MIntrinsic)*MRotation*pW;

% draw a little diagram
uC = imagePoint(1) / imagePoint(3);
vC = imagePoint(2) / imagePoint(3);

