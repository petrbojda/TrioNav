function [xp, Pp, S, in,K] = EKF_correction_update_ekf1_case1( x, P, y,h, H, R, M )
%EKFMEASUREMENTUPDATE
%
%   Performs the EKF measurement update of the state estimate and estimation
%   error covariance matrix
%   For a discrete-time system given as follows:
%
%   x[k] = f[k-1](x[k-1], u[k-1], w[k-1])
%   y[k] = h[k](x[k],v[k])
%
%   The w[k] and v[k] are white, zero-mean, uncorrelated and have known
%   covariance matrices Q and R.
%
%   w[k]   ~ (0,Q)
%   v[k]   ~ (0,R)
%   E[vw'] = 0
%
%   Calling Sequence
%   [xp, Pp] = EKF_correction_update( x, P, y, h, H, R )
%   [xp, Pp, K, in ] = EKF_correction_update( x, P, y, h, H, R, M )
%
%   Inputs
%       x  - n x 1, a priori state estimate
%       P  - n x n, a priori estimation error covariance matrix
%       y  - m x 1, measurement vector
%       h  - measurement equation (as a numeric vector or a function
%            handle)
%       H  - 1 x n, partial derivative of measurement function h with respect to
%            state (as a numeric matrix or a function handle)
%       R  - m x m, measurement noise covariance
%       M  - partial derivative of measurement function h with respect to
%            measurement noise (as a numeric matrix or a function handle)
%   Outputs
%       xp - n x 1, a posteriori state estimate
%       Pp - n x n, a posteriori estimation error covariance matrix
%       K  - Kalman gain
%       IN - inovation  y-h[k](x[k],0)
%       NIS - normalized innovation squared
%
%   Bibliography
%       [1] Simon D.: Optimal state estimation (13.48, 13.49)

% postacuje:
hx = h;
% hx = H*x;

M = eye(size(R,1));

% valid measurement
S = (H*P*H'+M*R*M');
K = (P*H')/S;             %[1] 13.49
% K = ((P*H')/S)*0;             %[1] 13.49

% in = y-hx;
in = y-hx;

if in(4)>345*pi/180 && in(4)<(360)*pi/180
   hx(4)=360*pi/180+hx(4);
elseif in(4)>-(360)*pi/180 && in(4)<-345*pi/180
    y(4)=y(4)+360*pi/180;
else
    y(4)=y(4);
    hx(4)=hx(4);
end
in = y-hx;    
xp = x + K*in;
Pp = (eye(size(P,1)) - K*H)*P;
Pp = 0.5*(Pp + Pp');   % Symmetrization of P
%     NIS = in'*inv(S)*in;
end
% hun = [[x;y;z];...
%     Cb2n*[vbx;vby;vbz];...
%     ([sfx;sfy;sfz]-g*[sin(th);sin(phi)*cos(th);cos(phi)*cos(th)])];



