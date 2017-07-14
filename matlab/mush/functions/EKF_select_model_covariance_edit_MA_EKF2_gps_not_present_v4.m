function [states_KF2, meas_KF2,Qc] = EKF_select_model_covariance_edit_MA_EKF2_gps_not_present_v4(model_selection)
switch model_selection
       case 3
        states_KF2 = 9;
        meas_KF2 = 4;        
cov_ax = 9.5795e-6*9.81;
cov_ay = 9.8039e-8*9.81;
cov_az = 8.2437e-5*9.81;
Q_param.sigPOSx = 2*1.265*10e-6*1*0;
Q_param.sigPOSy = 2*1.265*10e-6*1*0;
Q_param.sigPOSz = 2*1.265*10e-6*1*0;
Q_param.sigACCvelx = 8*0.3162*1;          %smerodatna odchylka ACC (m/s^2)
Q_param.sigACCvely = 8*0.3162*1;          %smerodatna odchylka ACC (m/s^2)
Q_param.sigACCvelz = 8*0.3162*1;          %smerodatna odchylka ACC (m/s^2)
tau = 1*60*60;  %% 1 hour in seconds

        Qc = diag([  (Q_param.sigPOSx)^2,    (Q_param.sigPOSy)^2,     (Q_param.sigPOSz)^2, ...
                        (Q_param.sigACCvelx)^2,(Q_param.sigACCvely)^2/10000, (Q_param.sigACCvelz)^2, ...
                        (cov_ax)^2/tau, (cov_ay)^2/tau, (cov_az)^2/tau])*1000;

end


