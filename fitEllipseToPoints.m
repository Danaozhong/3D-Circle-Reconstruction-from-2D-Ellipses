%
%
% Function to fit an ellipse into a set of points using a least-squares
% algorithm.
%
%
%
function [a, b, c, d, e, f] = fitEllipseToPoints(vec_points)
    % solve equations:
    % 1) ax^2 + 2bxy + cy^2 + 2dx + 2ey + f = 0
    % 2) ac - b^2 > 0
    %
    % For technical details, see the documentation of the Master's thesis
    % of Clemens Zangl.
    S1 = zeros(3, 3);
    S2 = zeros(3, 3);
    S3 = zeros(3, 3);
    
    for i = 1:length(vec_points)
        x = vec_points{i}(1);
        y = vec_points{i}(2);
        D1 = [x^2, 2*x*y, y^2];
        D2 = [2*x, 2*y, 1];
        
        S1 = S1 + D1.' * D1;
        S2 = S2 + D1.' * D2;
        S3 = S3 + D2.' * D2;
    end
    
    T = -inv(S3) * S2.';
    M = S1 + S2 * T;
    M = [M(3, :) ./ 2; -M(2, :); M(1, :) ./ 2.0];
    
    
    % calculate the eigenvalues
    [V, D] = eig(M);
    
    
    condition = V(1, :) .* V(3, :) - V(2, :) .^2;
    a1 = V(:, find(condition > 0));
    result = [a1; T * a1];
    % the correct ellipse parameters equal to the eigenvector with the
    % greatest eigenvalue. This is because 1/lambda*a = S^-1*C*a with a
    % being the vector of the ellipse parameter that are desired to be
    % found.
    a = result(1);
    b = result(2);
    c = result(3);
    d = result(4);
    e = result(5);
    f = result(6);    
end