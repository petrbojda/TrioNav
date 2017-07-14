function [xm, Pm] = EKF_time_update( x, P, f, F, Q, L, u )
%EKF TIME UPDATE 
%    
%   Performs the discrete-time EKF time update of the state estimate and estimation
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
%   [xm, Pm] = ekfTimeUpdate( x, P, f, F, Q)
%   [xm, Pm] = ekfTimeUpdate( x, P, f, F, Q, L, u )
%
%   Inputs
%       x  - n x 1, state estimate of previous step
%       P  - n x n, estimation error covariance matrix of previous step
%       f  - n x 1, system equation (as a numeric vector or a function
%            handle)
%       F  - n x n, partial derivative of system equation f with respect to
%            state (as a numeric matrix or a function handle)
%       Q  - n x n, covariance matrix of the process noise
%       L  - n x n, partial derivative of system equation f with respect to
%            process noise w
%       u  - n x 1, input vector
%   Outputs
%       xm - n x 1, a priori state mean estimate
%       Pm - n x n, a priori estimation error covariance matrix 
%
%   Bibliography
%       [1] Simon D.: Optimal state estimation (13.47)

if nargin < 5
    error('Wrong number of input arguments. See help kfMeasurementUpdate.')
end

% postacuje:
xm = f;
L = eye(size(x,1),size(Q,2));
Pm = F*P*F' + L*Q*L';                           % [1] 13.47
Pm = 0.5*(Pm + Pm');   % Symmetrization of P

end

