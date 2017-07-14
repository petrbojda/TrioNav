function [states_KF1, meas_KF1,Qc] = EKF_select_model_covariance_edit_MA_EKF1(model_selection)

switch model_selection
   
    case 3
    Q_EA = deg2rad(0.05);
    Q_EA_psi = deg2rad(0.1);
    Q_GYRBias = deg2rad(1e-4);                % ARS bias driving noise when initialized
   Qc =diag([Q_EA^2, Q_EA^2, Q_EA_psi^2,Q_GYRBias^2, Q_GYRBias^2, Q_GYRBias^2]);   

        
%         rx_std = 1.45444e-5; %% rad/s
%         ry_std = 1.45444e-5; %% rad/s
%         rz_std = 1.45444e-5; %% rad/s

        
%     Q_param.sigACCvel =deg2rad(0.4);
%     Q_param.sigACCbias =deg2rad(2e-3);
% 
%         Qc = diag([  (Q_param.sigACCvel)^2,    (Q_param.sigACCvel)^2,     (Q_param.sigACCvel)^2, ...
%                          (Q_param.sigACCbias)^2, (Q_param.sigACCbias)^2, (Q_param.sigACCbias)^2]);

%         Q_KF2 = diag([  (Q_param.sigACCvel)^2,    (Q_param.sigACCvel)^2,     (Q_param.sigACCvel)^2, ...
%                          (Q_param.sigACCbias)^2, (Q_param.sigACCbias)^2, (Q_param.sigACCbias)^2]);
% %                        
%  Q_KF2 =...
% [0.00850009000000000,0,0,0.499950000000000,0,0,0,0,0;0,0.00850008270000000,0,0,0.499950000000000,0,0,0,0;0,0,0.00850008270000000,0,0,0.499950000000000,0,0,0;0.499954601100000,0,0,499.999998800000,0,0,0,0,0;0,0.499954601000000,0,0,499.999998800000,0,0,0,0;0,0,0.499954601000000,0,0,499.999998800000,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0];
%  







%         rx_std = 1.93925e-6; %% rad/s
%         ry_std = 1.93925e-6; %% rad/s
%         rz_std = 1.93925e-6; %% rad/s
% % 0.00174533
% % 
% %         phi_std = 10^-10;
% %         th_std = 10^-10;
% %         psi_std = 10^-10;
%         phi_std = 0.0010;
%         th_std = 0.0011;
%         psi_std = 0.0011;
% 
% 
%         tau = 1*60*60;  %% 1 hour in seconds
        states_KF1 = 6;
        meas_KF1 = 3;
%         
%         Qc =diag([phi_std^2/tau, th_std^2/tau, psi_std^2/tau, (rx_std)^2/tau, (ry_std)^2/tau, (rz_std)^2/tau]);
end

