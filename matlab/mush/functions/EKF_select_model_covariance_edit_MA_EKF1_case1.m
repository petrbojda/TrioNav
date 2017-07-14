function [Qc] = EKF_select_model_covariance_edit_MA_EKF1_case1(model_selection)

switch model_selection
   
    case 3
    Q_EA = deg2rad(10e-2)^2; 
    Q_GYRBias = deg2rad(1e-6)^2;                % ARS bias driving noise when initialized
   Qc =diag([Q_EA/1000, Q_EA/1000, Q_EA,Q_GYRBias, Q_GYRBias, Q_GYRBias]);   


%%
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
% %         rx_std = 1.45444e-5; %% rad/s
% %         ry_std = 1.45444e-5; %% rad/s
% %         rz_std = 1.45444e-5; %% rad/s
%         tau = 1*60*60;  %% 1 hour in seconds
% 
%         
%         Qc =diag([phi_std^2/tau, th_std^2/tau, psi_std^2/tau, (rx_std)^2/tau, (ry_std)^2/tau, (rz_std)^2/tau]);
%         


end

