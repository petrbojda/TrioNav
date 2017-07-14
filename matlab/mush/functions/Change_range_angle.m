function angle = Change_range_angle(psi, mode)

switch mode
    case 1 % range (0-2*pi) 
        psi = mod(psi,2*pi);
        if psi<0
            psi = psi + 2*pi;
        elseif psi >= 2*pi
            psi = psi - 2*pi;
        end
end
angle = psi;
       