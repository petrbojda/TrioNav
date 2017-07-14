function [ Cb2n ] = Cb2n( phi,th,psi )
% Calculates Cb2n - Direction cosine matrix
Cb2n =  [  cos(psi)*cos(th),   cos(psi)*sin(phi)*sin(th) - cos(phi)*sin(psi),	sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(th);
            cos(th)*sin(psi),   cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(th),  cos(phi)*sin(psi)*sin(th) - cos(psi)*sin(phi);
            -sin(th),           cos(th)*sin(phi),                               cos(phi)*cos(th)];
    end

