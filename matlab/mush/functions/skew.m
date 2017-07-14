function [ X ] = skew( x )
%SKEW 
%
%   Performs the cross product or skew symetric form of a 3x1 vector.
%
%   Calling Sequence
%   [X] = skew( x )
%
%   Inputs
%       x  - 3 x 1, a vector 
%   Outputs
%       X  - 3 x 3, a cross product
%
%   Bibliography
%       % [1] Eun-Hwan Shin, Accuracy Improvement of Low Cost INS GPS for
%       Land Applications (3.5)

X = [ 0    -x(3)  x(2);
      x(3)  0    -x(1);
     -x(2)  x(1)  0];

end

