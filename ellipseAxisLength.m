function [major, minor] = ellipseAxisLength(a, b, c, d, e, f)
    up = 2*(a*e*e + c*d*d + f*b*b - 2*b*d*e - a*c*f);
    down1 = (b*b - a*c)*((c - a)*sqrt(1 + 4*b*b/((a - c)*(a - c))) - (c + a));
    down2 = (b*b - a*c)*((a - c)*sqrt(1 + 4*b*b/((a - c)*(a - c))) - (c + a));
    
    major = sqrt(up / down1);
    minor = sqrt(up / down2);
end