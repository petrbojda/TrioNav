function [funk_ekf1,Fk_ekf1,hunk_ekf1,Hk_ekf1,fun_ekf1] = EKF_compute_process_model_EKF1_case1(dt, prev_state,ars)
%% Read Me
% Process Model of EKF1
% Inputs: 3x ARS (x,y,z)
% State: x = [3x attitude; 3x ARS bias]
% Measurments: [3x ACC; Azimuth from GPS] or [3x ACC]

% function inputs are: 1)sampling time, 2)previous state from EKF1,
% 3)Angular rates
% function outputs are: 1) State prediction, 2) Jacobian of the state
% function, 3) Measurement prediction, 4) Jacobian of the measurement function, 5) Increment in the
% state variable

%% -------------------- variable's definition
% Inputs
arx = ars(1); ary = ars(2); arz = ars(3);
% States (3x attitude, 3x ARS bias)
phi = prev_state(1); th =prev_state(2); psi = prev_state(3);
bgx = prev_state(4); bgy = prev_state(5); bgz = prev_state(6);
%--------------------------------------------------------------------------
%% EKF 1 model computation
% Rotation of the angular rates from body frame to navigation frame (NED)
ARS_trans = [1 sin(phi)*tan(th) cos(phi)*tan(th);
             0 cos(phi)         -sin(phi);
             0 sin(phi)*sec(th) cos(phi)*sec(th)];
         
% function definiton of the EKF1
fun_ekf1 = [ARS_trans * ([arx;ary;arz] - [bgx; bgy; bgz]);...
        [0; 0; 0]];
% State prediciton
funk_ekf1 = prev_state + fun_ekf1*dt;
% Jacobian of the state function
Fk_ekf1 =...
[ dt*(cos(phi)*tan(th)*(ary - bgy) - sin(phi)*tan(th)*(arz - bgz)) + 1,         dt*(cos(phi)*(arz - bgz)*(tan(th)^2 + 1) + sin(phi)*(ary - bgy)*(tan(th)^2 + 1)), 0, -dt,   -dt*sin(phi)*tan(th),   -dt*cos(phi)*tan(th);
                     -dt*(cos(phi)*(arz - bgz) + sin(phi)*(ary - bgy)),                                                                                        1, 0,   0,           -dt*cos(phi),            dt*sin(phi);
  dt*((cos(phi)*(ary - bgy))/cos(th) - (sin(phi)*(arz - bgz))/cos(th)), dt*((cos(phi)*sin(th)*(arz - bgz))/cos(th)^2 + (sin(phi)*sin(th)*(ary - bgy))/cos(th)^2), 1,   0, -(dt*sin(phi))/cos(th), -(dt*cos(phi))/cos(th);
                                                                     0,                                                                                        0, 0,   1,                      0,                      0;
                                                                     0,                                                                                        0, 0,   0,                      1,                      0;
                                                                     0,                                                                                        0, 0,   0,                      0,                      1];
%
phi = funk_ekf1(1); th=funk_ekf1(2); psi = funk_ekf1(3);
psi=Change_range_angle(psi, 1);
%% Measurement Model if we use 3x ACC and azimuth from GPS 
% Measurment prediciton
hunk_ekf1=...
   [ sin(th);
 -cos(th)*sin(phi);
 -cos(phi)*cos(th);
    psi];
% Jacobian of the measurement function
Hk_ekf1=...
    [                 0,          cos(th), 0, 0, 0, 0;
  -cos(phi)*cos(th), sin(phi)*sin(th), 0, 0, 0, 0;
   cos(th)*sin(phi), cos(phi)*sin(th), 0, 0, 0, 0;
                  0,                0, 1, 0, 0, 0];
%% Measurement Model if we use 3x ACC only
%%% Measurment prediciton
% hunk=...
%    [ sin(th);
%  -cos(th)*sin(phi);
%  -cos(phi)*cos(th)];
%%% Jacobian of the measurement function
% Hk=...
%     [                 0,          cos(th), 0, 0, 0, 0;
%   -cos(phi)*cos(th), sin(phi)*sin(th), 0, 0, 0, 0;
%    cos(th)*sin(phi), cos(phi)*sin(th), 0, 0, 0, 0];
end

