function [funk_ekf2,Fk_ekf2,hunk_ekf2,Hk_ekf2] = EKF_compute_process_model_EKF2_case1_plus_vel_mag(dt, prev_state,Cb2n,aeb)

%% Read Me
% Process Model of EKF2
% Inputs: [3x Attitude (Euler's angel); 3x body acceleration; 3x
% centripetal force compensation]
% State, x = [3x position (NED); 3x velocity (body-frame); 3x ACC bias ]
% Measurements: [3x GPS position (NED)]

% function inputs are: 1)sampling time, 2)previous state from EKF2,
% 3)Direction cosine matrix, DCM (body frame to navigation frame), 4) body
% acceleration
% function outputs are: 1) State prediction, 2) Jacobian of the state
% function, 3) Measurement prediction, 4) Jacobian of the measurement function
%% ---------variable's definition ------------
% Defining the coefficients of DCM
C11 = Cb2n(1,1);
C12 = Cb2n(1,2);
C13 = Cb2n(1,3);
C21 = Cb2n(2,1);
C22 = Cb2n(2,2);
C23 = Cb2n(2,3);
C31 = Cb2n(3,1);
C32 = Cb2n(3,2);
C33 = Cb2n(3,3);
% body acceleration
aebx = aeb(1); aeby = aeb(2); aebz = aeb(3); 
% States
vbx = prev_state(4);vby=prev_state(5);vbz=prev_state(6);
bax=prev_state(7);bay=prev_state(8);baz=prev_state(9);
 
%% EKF 2 model computation
% function definiton of the EKF2
fun=...
    [C11*vbx + C12*vby + C13*vbz;
 C21*vbx + C22*vby + C23*vbz;
 C31*vbx + C32*vby + C33*vbz;
                  aebx - bax;
                  aeby - bay;
                  aebz - baz;
                           0;
                           0;
                           0];
                       
% State prediciton
funk_ekf2 = prev_state + fun*dt;
% Jacobian of the state function
Fk_ekf2=...
[ 1, 0, 0, C11*dt, C12*dt, C13*dt,   0,   0,   0;
  0, 1, 0, C21*dt, C22*dt, C23*dt,   0,   0,   0;
  0, 0, 1, C31*dt, C32*dt, C33*dt,   0,   0,   0;
  0, 0, 0,      1,      0,      0, -dt,   0,   0;
  0, 0, 0,      0,      1,      0,   0, -dt,   0;
  0, 0, 0,      0,      0,      1,   0,   0, -dt;
  0, 0, 0,      0,      0,      0,   1,   0,   0;
  0, 0, 0,      0,      0,      0,   0,   1,   0;
  0, 0, 0,      0,      0,      0,   0,   0,   1];

%% Measurement Model of EKF2 
% Measurment prediciton
hunk_ekf2 = [funk_ekf2(1);funk_ekf2(2);funk_ekf2(3)];
% Jacobian of the measurement function
Hk_ekf2=...
    [eye(3) zeros(3,6)];
% Hk_ekf2(4,4)=0;
% Hk_ekf2(5,5)=0;

end

