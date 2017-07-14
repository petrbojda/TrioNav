function [states_KF2, meas_KF2,Qc] = EKF_select_model_covariance_edit_MA_EKF2_v4(model_selection)
switch model_selection
      case 3
        states_KF2 = 9;
        meas_KF2 = 4;        
        cov_ax = 10e-3*9.80;
        cov_ay = 10e-3*9.80;
        cov_az = 1e-6*9.80;
% Q_param.sigPOSx = 2*1.265*1*10e2*1*1;
% Q_param.sigPOSy = 2*1.265*1*10e2*1*1;
% Q_param.sigPOSz = 2*1.265*1*10e2*1*1;
% Q_param.sigACCvelx = 30;          %smerodatna odchylka ACC (m/s^2)
% Q_param.sigACCvely = 30;          %smerodatna odchylka ACC (m/s^2)
% Q_param.sigACCvelz = 30;          %smerodatna odchylka ACC (m/s^2)
Q_param.sigPOSx = 10;
Q_param.sigPOSy = 10;
Q_param.sigPOSz = 1e-6;
Q_param.sigACCvelx = 10;          %smerodatna odchylka ACC (m/s^2)
Q_param.sigACCvely = 15;          %smerodatna odchylka ACC (m/s^2)
Q_param.sigACCvelz = 1;          %smerodatna odchylka ACC (m/s^2)
        tau = 1*1*1;  %% 1 hour in seconds

        Qc = diag([  (Q_param.sigPOSx)^2,    (Q_param.sigPOSy)^2,     (Q_param.sigPOSz)^2, ...
                        (Q_param.sigACCvelx)^2,(Q_param.sigACCvely)^2, (Q_param.sigACCvelz)^2, ...
                        (cov_ax)^2/tau, (cov_ay)^2/tau, (cov_az)^2/tau]);

end

