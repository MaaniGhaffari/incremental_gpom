function [z, sz] = fusionBCM(x, y, sx, sy)
% Bayesian Committee Machine for map fusion
% inputs: 
%       x: first variable
%       y: second variable
%       sx: first variable variance
%       sy: second variable variance
% outputs: 
%       z: fused value
%       sz: fused variance
if sx == 0
    sz = sy;
    z = y;
    return;
elseif sy == 0
    sz = sx;
    z = x;
    return;
elseif sx == 0 && sy == 0
    error('Both sx and sy are zeros.')
else
    sz = 1./(1./sx + 1./sy);
    z = sz .* ( 1./sx .* x + 1./sy .* y);
end
    

    