function [angle] = ellipseAngle(a, b, c, d, e, f)
    if(b == 0)
        if(a < c)    
            angle = 0; 
        else(a > c)
            angle = pi / 2;
        end
    else
       if(a < c)
         angle = 0.5*acot((a - c) / (2*b));    
       else
         % Todo Clemens: Check why + 90° is not necessary.
         angle =  + 0.5*acot((a - c) / (2*b));     
       end
    end
        

end