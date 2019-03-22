%
% EQUATION SOLVER TO UNPROJECT A 2D ELLIPSE TO A 3D CIRCLE
% This program searches for a 3D circle and its radius and orientation if
% based on multiple image views
% 5 January 2015
%


%
%
% PART I: Set up equations
%
%


clear all;

% parameters of the intrinsic camera parameters
syms focal_length;

% the unknown transformation from the circle to the world coordinate system
syms ctx cty ctz;
syms cr11 cr12 cr13 cr21 cr22 cr23 cr31 cr32 cr33;
MCircleToWorldTransformation = [cr11 cr12 cr13 ctx; cr21 cr22 cr23 cty; cr31 cr32 cr33 ctz; 0 0 0 1];

% the unknown circle radius
syms cRadius;
    

% define the properties of the circle we are looking for. Six parameters
% are enough to fully describle the 3D circle (C the center point, N the
% circle support plane normal with the length of N equal to the circle
% radius).
syms nx ny nz;
syms cx cy cz;
N = [nx; ny; nz];
C = [cx; cy; cz];


% define the projection matrix. For simplification, the projection
% parameters were shifted and scaled so that f = fx = fy and cx = cy = 0
MIntrinsic = [focal_length 0 0; 0 focal_length 0; 0 0 1];

%
%
% PART II: Generate test data
%
%

cRadius = 3.0;

focal_length = 1.25;
MIntrinsicValues = eval(MIntrinsic);
    
ctx = 2.0;
cty = 1.3;
ctz = 0.3;
% rotation matrix with Euler angles
eulerX = 22 / 180*pi;
eulerY = 73 / 180*pi;
eulerZ = 152 / 180*pi;
matrix = makehgtform('zrotate', eulerZ)*makehgtform('yrotate', eulerY)*makehgtform('xrotate', eulerX);
cr11 = matrix(1,1);
cr12 = matrix(1,2);
cr13 = matrix(1,3);
cr21 = matrix(2,1);
cr22 = matrix(2,2);
cr23 = matrix(2,3);
cr31 = matrix(3,1);
cr32 = matrix(3,2);
cr33 = matrix(3,3);

MCircleToWorldTransformation = eval(MCircleToWorldTransformation);

numOfEllipses = 3;

for i = 1:numOfEllipses
    tx = rand * 2.0;
    ty = rand * 3.0;
    tz = (rand + 1) * 10.0;
    % rotation matrix with Euler angles
    eulerX = rand*360 / 180*pi;
    eulerY = rand*360 / 180*pi;
    eulerZ = rand*360 / 180*pi;
    
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

    % test the projection settings
    numOfPoints = 60;

    % create the 3D data points (homogeneous coordinates)
    points3D = cell(1, numOfPoints);
    for k = 1:numOfPoints
        currentAngle = 2*pi*(k-1) / numOfPoints;
        X = cos(currentAngle)*cRadius;
        Y = sin(currentAngle)*cRadius;
        points3D{k} = [X; Y; 0; 1];
    end

    % create the 2D data points
    MExtrinsicValues = [r11 r12 r13 tx; r21 r22 r23 ty; r31 r32 r33 tz];


    points2DOriginal = cell(1, numOfPoints);

    for k = 1:numOfPoints
        currentPoint = MIntrinsicValues * MExtrinsicValues * MCircleToWorldTransformation * points3D{k};
        currentu = currentPoint(1) / currentPoint(3);
        currentv = currentPoint(2) / currentPoint(3);
        points2DOriginal{k} = [currentu, currentv];
    end

    % find the ellipse parameters
    [a b c d e f] = fitEllipseToPoints(points2DOriginal);
    
        
    % scale all ellipses to ensure that ac - b^2 = 1
    scale = 1/sqrt(a*c - b*b);
    a = a * scale;
    b = b * scale;
    c = c * scale;
    d = d * scale;
    e = e * scale;
    f = f * scale;    
    
    ellipse_data(i).points_original = points2DOriginal;
    ellipse_data(i).a = a;
    ellipse_data(i).b = b;
    ellipse_data(i).c = c;
    ellipse_data(i).d = d;
    ellipse_data(i).e = e;
    ellipse_data(i).f = f;
    
    extrinsic_data{i} = MExtrinsicValues; 
end

%
%
% PART II Validation - Print the generated points of the ellipses and the
% found ellipse using the LS algorithm on the same plot.  If everything is 
% okay, the points from the given ellipse should match the calculated LS
% ellipse.
%
%
figure();
hold on
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);

for i = 1:length(ellipse_data)
    a = ellipse_data(i).a;
    b = ellipse_data(i).b;
    c = ellipse_data(i).c;
    d = ellipse_data(i).d;
    e = ellipse_data(i).e;
    f = ellipse_data(i).f;
    points2DOriginal = ellipse_data(i).points_original;
    
    [cx, cy] = ellipseCenter(a, b, c, d, e, f);
    [angle] = ellipseAngle(a, b, c, d, e, f);
    [majorAxis, minorAxis] = ellipseAxisLength(a, b, c, d, e, f);

    xvalues = zeros(1, length(points2DOriginal));
    yvalues = zeros(1, length(points2DOriginal));
    for k = 1:length(points2DOriginal)
        xvalues(k) = points2DOriginal{k}(1);
        yvalues(k) = points2DOriginal{k}(2);    
    end

    % draw the input point data 
    scatter(xvalues, yvalues, [], [1 0 0], 'filled');

    p = ellipsePoints(cx, cy, majorAxis, minorAxis, angle);
    plot(p(:,1), p(:,2), '.-')
end
hold off




%
%
% PART III: Find the optimum 3D circle
%
%

% Transform the input data. Decompose the intrinsic parameters in a rotation
% matrix and the translation vector
for i = 1:numOfEllipses
    % extrinsic_matrix = R [I3 -S]
    R = extrinsic_data{i}((1:3), (1:3));
    T = extrinsic_data{i}((1:3), 4);
    S = -inv(R)*T;
    s(i).x = S(1);
    s(i).y = S(2);
    s(i).z = S(3);

end


clearvars r11 r12 r13 r21 r22 r23 r31 r32 r33;
syms r11 r12 r13 r21 r22 r23 r31 r32 r33;
syms sx sy sz;

gS = [sx; sy; sz];

gM1 = focal_length*[r11; r12; r13];
gM2 = focal_length*[r21; r22; r23];
gM3 = [r31; r32; r33];
gM = [gM1; gM2; gM3];

% set initial values for N and C
vecValues = [1; 1; 1; 1; 1; 1];

%set the threshold for the finish state
maxIterations = 200;
errorThreshold = 10;
numIterations = 0;

while numIterations < maxIterations %errorThreshold > 0.1
    numIterations = numIterations + 1;
    
    C = [vecValues(1); vecValues(2); vecValues(3)];
    N = [vecValues(4); vecValues(5); vecValues(6)];
    
    % Recalculate the jacobian matrix
    jacobi = zeros(numOfEllipses*5, 6);
    
    r = zeros(5*numOfEllipses, 1);
    
    for i = 1:numOfEllipses
        % read the values from the current ellipse screen coordinates: first,
        % start with the existrinsic parameters
        r11 = extrinsic_data{i}(1, 1);
        r12 = extrinsic_data{i}(1, 2);
        r13 = extrinsic_data{i}(1, 3);

        r21 = extrinsic_data{i}(2, 1);
        r22 = extrinsic_data{i}(2, 2);
        r23 = extrinsic_data{i}(2, 3);

        r31 = extrinsic_data{i}(3, 1);
        r32 = extrinsic_data{i}(3, 2);
        r33 = extrinsic_data{i}(3, 3);

        sx = s(i).x;
        sy = s(i).y;
        sz = s(i).z;

        M{1} = eval(gM1);
        M{2} = eval(gM2);
        M{3} = eval(gM3);

        S = eval(gS);

        % and the ellipse 2D parameters
        a = ellipse_data(i).a;
        b = ellipse_data(i).b;
        c = ellipse_data(i).c;
        d = ellipse_data(i).d;
        e = ellipse_data(i).e;
        f = ellipse_data(i).f;

        % transform the ellipse to their dual quadric 
        eobs11 = c*f - e*e;
        eobs12 = d*e - b*f;
        eobs13 = b*e - c*d;
        eobs22 = a*f - d*d;
        eobs23 = b*d - a*e;
        eobs33 = 1;
        assert(abs(a*c - b*b - 1) < 0.00001);
        
        % matris is symmetric
        Eobs = [eobs11 eobs12 eobs13; eobs12 eobs22 eobs23; eobs13 eobs23 eobs33];

        % anonymous function to calculate the A matrix (required for the
        % differential equations)
        MatrixA = @(i, j) (M{i}*M{j}.' + M{j}*M{i}.');

        % calculate the derivatives that are commonly needed. The
        % derivatives are vectors (size 3)
        dE33dC = (C - S).' * MatrixA(3,3);
        dE33dN = N.' * (MatrixA(3, 3) - (2*M{3}.'*M{3}) * eye(3));

        % This is a scalar
        E33 = dot((C - S), M{3}) * dot((C - S), M{3}) - dot(cross(N, M{3}), cross(N, M{3}));
        
        % each ellipse yields us five equations. Each value represents a
        % different 
        for k = 1:5
            % first calculate the partial derivative. Each derivative is a
            % vector of three elements (x y and z elements for the C or N
            % vector respectively).
            indexRow = 0;
            indexCol = 0;
            if (1 == k)
                % index (1,1)
                indexRow = 1;
                indexCol = 1;
            elseif(2 == k)
                % index (1,2)
                indexRow = 1;
                indexCol = 2;
            elseif(3 == k)
                % index (2,2) 
                indexRow = 2;
                indexCol = 2;
            elseif(4 == k)
                % index (1,3)
                indexRow = 1;
                indexCol = 3;
            else
                % index (2,3)
                indexRow = 2;
                indexCol = 3;
            end

            dEdC = (C - S).' * MatrixA(indexRow,indexCol);
            dEdN = N.' * (MatrixA(indexRow,indexCol) - (2*M{indexRow}.'*M{indexCol}) * eye(3));
            currentEobs = Eobs(indexRow,indexCol); 

            row = (i-1)*5 + k;

            % set the values of the jacobi matrix
            for l = 1:3
                % set up the partial derivatives for dF / dC
                jacobi(row, l) = dEdC(l) - dE33dC(l)*currentEobs;

                % set up the partial derivatives for dF / dN
                jacobi(row, 3 + l) = dEdN(l) - dE33dN(l)*currentEobs;
            end
            
            Ecrowccol = dot(C - S, M{indexRow}) * dot(C - S, M{indexCol}) - dot(cross(N, M{indexRow}), cross(N, M{indexCol}));
            
            r(row) = Ecrowccol - E33*currentEobs;
        end
    end
      
   
    % all pre-calculations done, iterate 
    vecValuesNew = vecValues - inv(jacobi.'*jacobi)*jacobi.'*r;
    
    vecValues = vecValuesNew;
end
