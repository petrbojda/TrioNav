clc, close all; clear all;
load('data/full_length_gps_data_saved.mat')
load Data_Outages/20160514_INS_DMU10_Flight_01_MAXM8W_GEO.mat
% load Data_Outages/20160514_INS_DMU10_Flight_01_MAXM8W_GEO_shift_0.mat


warning ('off','all');
addpath('.\functions');
% data= data(:,:);
data= data(1:(400000),:);
Tini = 130; 
YAW_gps = 12;   %%% initial yaw from GPS in degrees
vel_gps =sqrt(data(:,26).^2+data(:,27).^2+data(:,28).^2);
%% Initilization set up
Fs = 200;
Ts = 1/Fs;          % sampling period of the IMU data
M = length(data);    %%% total number of data
t = (0:1:M-1)*Ts;
Mini=Tini*Fs;	% initialization length (seconds)
GPS_INIT = 1:Tini*200;
GPS_lever_arm = [0; 0; 0];               
Ts_GPS = 0.2;       % GPS sampling time
%------------------------ sensoric ans system data
% system data ranges
acc = [data(:,3) data(:,5) data(:,7)];
ars = [data(:,11) data(:,13) data(:,15)];
ars = deg2rad(ars);
%% GPS DATA Setting
GPSr = 20:22; 
GPS = data(:,GPSr);
% PN_gps = NaN(3,M);       % pozice NED na GPS frekvenci (mimo prijata data NaN)
% Geodetic constant
a = 6378137.0;                  % semi-major axis (m)
b = 6356752.3142;               % semi-minor axis (m)
e = sqrt(1 - b^2/a^2);          % ecentricity
%------------ MAG
iniLAT = nanmedian(GPS(GPS_INIT,1));  % Initial latitude (deg)
iniLON = nanmedian(GPS(GPS_INIT,2));  % Initial longitude (deg)
iniALT = nanmedian(GPS(GPS_INIT,3));  % Initial altitude (m)
% Calculation of the gravitational acceleration on the basis of initial position from GPS
gN = comp_gravity([iniLAT iniLON iniALT]');   % Gravitational accleration (m/s^2)
% Radius of the position depending on the initial latitude and longitude
R_N = a*(1 - e^2)/(1 - e^2*(sin(deg2rad(iniLAT))^2))^(3/2);  % Radius of curvature in the prime vertical
R_M = a/(1 - e^2*(sin(deg2rad(iniLAT))^2))^(1/2);            % Meridian radius of curvature
PN_gps_mem=NaN(3,M);
%% EKF Model
% EKF1
% [No_Xk_states_ekf1, No_meas_states_ekf1,Qc_ekf1] = EKF_select_model_covariance_edit_MA_EKF1(3);  %%%% EKF_select_model_covariance file FIXED
[No_Xk_states_ekf1, No_meas_states_ekf1,~] = EKF_select_model_covariance_edit_MA_EKF1(3);  %%%% EKF_select_model_covariance file FIXED
[No_Xk_states_ekf2, No_meas_states_ekf2,Qc_ekf2] = EKF_select_model_covariance_edit_MA_EKF2_v4(3);  %%%% EKF_select_model_covariance file FIXED
% [~,~,Qc_ekf2_gps_not_present] = EKF_select_model_covariance_edit_MA_EKF2_gps_not_present_v4(3);

Qc_ekf1_case1 = EKF_select_model_covariance_edit_MA_EKF1_case1(3);
Qc_ekf1_case2 = EKF_select_model_covariance_edit_MA_EKF1_case2(3);
Qc_ekf1_case3 = EKF_select_model_covariance_edit_MA_EKF1_case3(3);
Qc_ekf1_case4 = EKF_select_model_covariance_edit_MA_EKF1_case4(3);

% EKF1 Variable Sets
Xk_ekf1 = zeros(No_Xk_states_ekf1,1);     % State vector
Pk_ekf1 = zeros(No_Xk_states_ekf1);       % Covariance matrix
Zk_ekf1 = zeros(No_meas_states_ekf1,1);   % Measurement vector
Xp_ekf1 = Xk_ekf1;
Pp_ekf1 = Pk_ekf1;
%EKF2 Variable sets
Xk_ekf2 = zeros(No_Xk_states_ekf2,1);     % State vector
Pk_ekf2 = zeros(No_Xk_states_ekf2);       % Covariance matrix
Zk_ekf2 = zeros(No_meas_states_ekf2,1);   % Measurement vector
Xp_ekf2 = Xk_ekf2;
Pp_ekf2 = Pk_ekf2;

%Data Saving
Acentrifugal=zeros(3,M);
aeb=zeros(3,M);
X_state =zeros(No_Xk_states_ekf1+No_Xk_states_ekf2,M);
%% Motion Detection Variable Definition
motion_detec_a = nan(1,M);
motion_a = nan(1,M);
motion_detec_r = nan(1,M);
motion_r = nan(1,M);
motion_detec_dpsi = nan(1,M);
motion_dpsi = nan(1,M);
motion =nan(1,M);
aq=zeros(1,M);
as=zeros(1,M);
V_body_mem =nan(3,M);
acc_th = zeros(1,M);
ars_th = zeros(1,M);
vel_gps_mag =nan(1,M);
psidot_th =nan(1,M);
delta_euler=zeros(6,M);
V_NED_gps_mem=nan(3,M);
% V_NED_gps_mem(1:3,j)=V_NED_gps;
V_NED_ekf_estimate_mem =nan(3,M);


%% Various Data Saving
error_pos_est_RTK =nan(3,M);
error_pos_RTK_GPS=nan(3,M);
error_pos_est_GPS=nan(3,M);
error_pos_est_GPSvel = nan(3,M);
error_pos_est_GPS_withalltime=nan(3,M);
error_pos_est_GPSvel_withalltime=nan(3,M);
error_pos_est_GPS_RTK_withalltime=nan(3,M);

%% Dynamic Motion Detection
ACC_motion_threshold = 50/1000;         % 50 mg
ARS_motion_threshold = deg2rad(3);    % 3 deg/s converted to rad/s
delta_euler_motion_threshold = deg2rad(5);
normACC = zeros(1,M);           % Total magnitude of the ACC
normARS = zeros(1,M);           % Total magnitude of the ARS

%% FUSION and Parameter Settings
% EKF1
Xk_mem_ekf1 = zeros(No_Xk_states_ekf1, M);
Pk_mem_ekf1 = zeros(No_Xk_states_ekf1,No_Xk_states_ekf1,M);
% Innov_mem_ekf1 = NaN(3,M);
Innov_mem_ekf1 = NaN(4,M);
% Innov_mem_ekf1 = zeros(4,M);
% residual_mem_ekf1 = ones(3,M);
% smoothed_innv_mem =NaN(4,M);
% Zk_mem_ekf1 = NaN(3,M);
Zk_mem_ekf1 = NaN(4,M);
EUL_GYR_mem = zeros(3,M);       % Euler angles from GYROs
EUL_ACC_MAG_mem = zeros(3,M);   % Euler angles from ACC&MAG
% EKF2
C_b2n_mem = zeros(3,3,M);                       % Rotational matrix from body frame to navigation frame
Xk_mem_ekf2 = zeros(No_Xk_states_ekf2, M);
Pk_mem_ekf2 = zeros(No_Xk_states_ekf2,No_Xk_states_ekf2,M);
Innov_mem_ekf2 = NaN(6,M);
% Zk_mem_ekf2 = NaN(No_meas_states_ekf2,M);
Zk_mem_ekf2 = NaN(6,M);
% residual_mem_ekf2 = ones(3,M);

%% === definice parametru
%------------ IMU
ROT = eye(3);   % ROT = [0 1 0; 1 0 0; 0 0 1];  % YX-Z (ENU) do XYZ (NED~body) 
K_acc = [   1    0    0
            0    1    0
            0    0    1];
SF_acc = [  1    0    0
            0    1    0
            0    0    1];
b_fmu =  [  0;   0;   0];
% b_fmu =  [  0.0897;   -0.056;   0.0403];
ACC_deterministic_comp = ROT*SF_acc*K_acc;        
        
K_ars = [   1    0    0
            0    1    0
            0    0    1];
SF_ars = [  1    0    0
            0    1    0
            0    0    1];
ARS_deterministic_comp = SF_ars*K_ars*ROT;         


% design simple lowpass FIR filter
h=fdesign.lowpass('N,Fc', 40, 5, Fs);
Hd=design(h);



ARS_mem(:,1) = filter(Hd.Numerator,1,ars(:,1));
ARS_mem(:,2) = filter(Hd.Numerator,1,ars(:,2));
ARS_mem(:,3) = filter(Hd.Numerator,1,ars(:,3));
ACC_mem(:,1) = filter(Hd.Numerator,1,acc(:,1));
ACC_mem(:,2) = filter(Hd.Numerator,1,acc(:,2));
ACC_mem(:,3) = filter(Hd.Numerator,1,acc(:,3));

% ACC_mem(:,1) = acc(:,1);
% ACC_mem(:,2) = acc(:,2);
% ACC_mem(:,3) = acc(:,3);
% ARS_mem(:,1) = ars(:,1);
% ARS_mem(:,2) = ars(:,2);
% ARS_mem(:,3) = ars(:,3);

% Compute initial conditions for a period Tini
% Variable settings 
bias_a = zeros(3,1) ;   % primary ACC = bias ACC
bias_w = zeros(3,1) ;   % primary ARS = bias ARS

for j=1:Mini
    ACC = (ACC_deterministic_comp * (ACC_mem(j,1:3)' - b_fmu))';
    ARS = (ARS_deterministic_comp * ARS_mem(j,1:3)')'; 
    % primary ACC
    bias_a(1) = bias_a(1) + ACC(1)/Mini;
    bias_a(2) = bias_a(2) + ACC(2)/Mini;
    bias_a(3) = bias_a(3) + ACC(3)/Mini;
    % prumery ARS = bias ARS
    bias_w(1) = bias_w(1) + ARS(1)/Mini;
    bias_w(2) = bias_w(2) + ARS(2)/Mini;
    bias_w(3) = bias_w(3) + ARS(3)/Mini;
end
% bias_a = (bias_a+[0;0;1])*9.81;
bias_a = [0;0;-0]*9.81./1000;

% bias_w(1:3)=0;
% bias_a =[0;0;0];
% -------------------- Euler angles from ACC & MAG
phi = atan2(mean(-ACC_mem(5:Mini,2)), sqrt( mean(ACC_mem(5:Mini,1))^2 + mean(ACC_mem(5:Mini,3))^2 ));
th = atan2(mean(ACC_mem(5:Mini,1)), mean(-ACC_mem(5:Mini,3)));
psi_gps_init = deg2rad(YAW_gps);

% psi =Change_range_angle(psi_gps_init,1);
% psi =psi_gps_init*pi/180;
psi =psi_gps_init;


EUL_ACC_MAG_mem(:,1:Mini) = repmat([phi;th;psi],1,Mini);
Xk_mem_ekf1(4:6,1:Mini) = repmat(bias_w,1,Mini);
psi_gps(1:1:Mini)=psi_gps_init;
EUL_GYR = [phi; th; psi];
C2bn = Cb2n(EUL_GYR(1),EUL_GYR(2),EUL_GYR(3));
EUL_GYR_mem(:,1:Mini) = repmat(EUL_GYR,1,Mini);
Xk_mem_ekf1(1:3,1:Mini) = repmat(EUL_GYR,1,Mini);

Xk_mem_ekf2(7:9,1:Mini) = repmat(bias_a,1,Mini);

%%  EKF initialization
% EKF1
Xp_ekf1(1) = EUL_GYR(1);
Xp_ekf1(2) = EUL_GYR(2);
Xp_ekf1(3) = EUL_GYR(3);
Xp_ekf1(4) = bias_w(1);
Xp_ekf1(5) = bias_w(2);
Xp_ekf1(6) = bias_w(3);

Pp_ekf1 = zeros(No_Xk_states_ekf1,No_Xk_states_ekf1);

Pp_ekf1(1,1) = 4.5e-3;
Pp_ekf1(2,2) = 4e-3;
Pp_ekf1(3,3) = 1e-3;
Pp_ekf1(4,4) = 8e-11;
Pp_ekf1(5,5) = 2.5e-11;
Pp_ekf1(6,6) = 1.5e-10;

%% EKF2
Xp_ekf2(1) = 0;
Xp_ekf2(2) = 0;
Xp_ekf2(3) = 0;
Xp_ekf2(4) = 0;
Xp_ekf2(5) = 0;
Xp_ekf2(6) = 0;
% Xp_ekf2(7) = bias_a(1);
% Xp_ekf2(8) = bias_a(2);
% Xp_ekf2(9) = bias_a(3);
Xp_ekf2(7) = bias_a(1);
Xp_ekf2(8) = bias_a(2);
Xp_ekf2(9) = bias_a(3);

Pp_ekf2 = zeros(No_Xk_states_ekf2,No_Xk_states_ekf2);
% 
Pp_ekf2(1,1) = 0.5*1;
Pp_ekf2(2,2) = 0.5*1;
Pp_ekf2(3,3) = 0.5*2*1;
Pp_ekf2(4,4) = 1*0.1*1;
Pp_ekf2(5,5) = 1*0.1*1;
Pp_ekf2(6,6) = 1*2*0.1*1;
% Pp_ekf2(7,7) = 10*9.81*0/1000;
% Pp_ekf2(8,8) = 10*9.81*0/1000;
% Pp_ekf2(9,9) = 2*10*9.81*0/1000;
Pp_ekf2(7,7) = 5e-2;
Pp_ekf2(8,8) = 5e-2;
Pp_ekf2(9,9) = 4e-2;
%% GPS from RTK using Geodetic X91+GNSS receiver
GPS_RTK=data(:,47:49);
iniLAT_RTK = nanmedian(GPS_RTK(1:Mini,1));
iniLON_RTK = nanmedian(GPS_RTK(1:Mini,2));
iniALT_RTK = nanmedian(GPS_RTK(1:Mini,3));
PN_gps_RTK_mem=NaN(3,M);
%%
PN_gps_x = NaN(1,M);
PN_gps_y = NaN(1,M);
PN_gps_z = NaN(1,M);
V_gps_z = NaN(1,M);
window_size_m_pred= 300;
num_x =zeros(1,M);
num_y =zeros(1,M);
num_z =zeros(1,M);
num_vz =zeros(1,M);
measurement_error =NaN(4,M);
%%
 for j =1:Mini
     %% For Position estimation from GPS
    if all(GPS(j,:) == 0) || any(isnan(GPS(j,:)))
        GPS(j,:) = [NaN NaN NaN];
        PN_gps = [NaN NaN NaN]';
    else
        if ~isnan(GPS(j,:)) 
            PN_gps(1,1) = (GPS(j,1) - iniLAT)/rad2deg(atan(1/R_M));
            PN_gps(2,1) = (GPS(j,2) - iniLON)/rad2deg(atan(1/(R_N*cos(deg2rad(iniLAT)))));
            PN_gps(3,1) = -GPS(j,3) + iniALT;
            PN_gps = PN_gps - (C2bn * GPS_lever_arm);
        end
    end
PN_gps_mem(:,j) = PN_gps;
%% For GPS position from RTK
    if all(GPS_RTK(j,:) == 0) || any(isnan(GPS_RTK(j,:)))
        GPS_RTK(j,:) = [NaN NaN NaN];
        PN_gps_RTK = [NaN NaN NaN]';
    else
        if ~isnan(GPS_RTK(j,:)) 
            PN_gps_RTK(1,1) = (GPS_RTK(j,1) - iniLAT_RTK)/rad2deg(atan(1/R_M));
            PN_gps_RTK(2,1) = (GPS_RTK(j,2) - iniLON_RTK)/rad2deg(atan(1/(R_N*cos(deg2rad(iniLAT_RTK)))));
            PN_gps_RTK(3,1) = -GPS_RTK(j,3) + iniALT_RTK;
            PN_gps_RTK = PN_gps_RTK - (C2bn * GPS_lever_arm);
        end
    end
PN_gps_RTK_mem(:,j) = PN_gps_RTK;
 end
 
%%
K_ekf1_mem = NaN(6,M);
% K_ekf2_mem = NaN(4,M);
K_ekf2_mem = zeros(9,6,M);
%%
% load('ACF_C2bn.mat')
% load('ACC_data_saved.mat');
fprintf('Percentage Completed (%%)\n')
for j=Mini+1:M
    % Cycle progress in percentage
    if mod(j,floor(M/10))==0
        fprintf('%d ', round(j/M*100));
    end
% %% Input vector processing    
%     % ====================================== ACC+ARS
    ACC = (ACC_deterministic_comp *(ACC_mem(j,1:3)'));     % Deterministic error compensation
    ARS = (ARS_deterministic_comp *(ARS_mem(j,1:3)'));     % Deterministic error compensation
%     % vector magnitudes
    normACC(j) = norm(ACC);                 % Total magnitude of the ACC
    normARS(j) = norm(ARS);     % Total magnitude of the ARS
    
    acc_th(j)=abs(normACC(j) - 1);
    ars_th(j)=abs(normARS(j));
%%
    % ==================================== Euler angles from GYROS only
    phi = EUL_GYR(1); th = EUL_GYR(2);
    ARS_trans = [   1 sin(phi)*tan(th) cos(phi)*tan(th);...
                    0 cos(phi)         -sin(phi);...
                    0 sin(phi)*sec(th) cos(phi)*sec(th)];
                d_Euler = ARS_trans * (ARS - bias_w);

    EUL_GYR = EUL_GYR + d_Euler*Ts;
    EUL_GYR(3) = Change_range_angle(EUL_GYR(3), 1);
    psi = EUL_GYR(3);
    EUL_GYR_mem(:,j) = EUL_GYR;
% ==================================== Eulers angles from ACC+MAG
    ax = ACC(1)/normACC(j); ay = ACC(2)/normACC(j); az = ACC(3)/normACC(j);   
    phi = atan2(-ay, sqrt(ax^2 + az^2));
    th = atan2(ax,-az);
    EUL_ACC_MAG_mem(:,j) = [phi;th;psi];
    
    %% For Position estimation from GPS % ====================================== GPS data processing
    if all(GPS(j,:) == 0) || any(isnan(GPS(j,:)))
        GPS(j,:) = [NaN NaN NaN];
        PN_gps = [NaN NaN NaN]';
    else
        if ~isnan(GPS(j,:)) 
            PN_gps(1,1) = (GPS(j,1) - iniLAT)/rad2deg(atan(1/R_M));
            PN_gps(2,1) = (GPS(j,2) - iniLON)/rad2deg(atan(1/(R_N*cos(deg2rad(iniLAT)))));
            PN_gps(3,1) = -GPS(j,3) + iniALT;
            PN_gps = PN_gps - (C2bn * GPS_lever_arm);
        end
    end
PN_gps_mem(:,j) = PN_gps;
%% For GPS position from RTK
    if all(GPS_RTK(j,:) == 0) || any(isnan(GPS_RTK(j,:)))
        GPS_RTK(j,:) = [NaN NaN NaN];
        PN_gps_RTK = [NaN NaN NaN]';
    else
        if ~isnan(GPS_RTK(j,:)) 
            PN_gps_RTK(1,1) = (GPS_RTK(j,1) - iniLAT_RTK)/rad2deg(atan(1/R_M));
            PN_gps_RTK(2,1) = (GPS_RTK(j,2) - iniLON_RTK)/rad2deg(atan(1/(R_N*cos(deg2rad(iniLAT_RTK)))));
            PN_gps_RTK(3,1) = -GPS_RTK(j,3) + iniALT_RTK;
            PN_gps_RTK = PN_gps_RTK - (C2bn * GPS_lever_arm);
        end
    end
PN_gps_RTK_mem(:,j) = PN_gps_RTK;
%% YAW from GPS
    if ~any(isnan(data(j,28))) && ~any(isnan(data(j,29)))
        psi_gps(1,j)=atan2(data(j,29),data(j,28));
        psi_gps(1,j)=Change_range_angle(psi_gps(1,j), 1);
    else
        psi_gps(1,j) =NaN;
    end
%% EKF1 time update EKF 1 
% Control Input
%     ARS_input_ekf1 = ARS - Xk_mem_ekf1(4:6,j-1);
        ARS_input_ekf1 = ARS ;
 % Measurement   
%     accx = ACC(1); 
%     accy = ACC(2); 
%     accz = ACC(3);
    accx = ACC(1)- (Xk_mem_ekf2(7,j-1))./gN; 
    accy = ACC(2)- (Xk_mem_ekf2(8,j-1))./gN; 
    accz = ACC(3)- (Xk_mem_ekf2(9,j-1))./gN;
vel_gps_mag(j) = sqrt(data(j,28).^2+data(j,29).^2+data(j,27).^2);
%% motion detection
[~,~,~,~,fun_ekf1] = EKF_compute_process_model_EKF1_case1( Ts, Xp_ekf1,ARS_input_ekf1');
delta_euler (:,j)=fun_ekf1;
            if (abs(normACC(j) - 1) <= ACC_motion_threshold)
                motion_detec_a(j) = 0;
            else
                motion_detec_a(j) = 1;
            end
            if sum(motion_detec_a(j-150:j))==0
                motion_a(j) = 0;
            else 
                motion_a(j) =1;
            end
            
            if (abs(normARS(j)) <= ARS_motion_threshold)
                motion_detec_r(j) = 0;
            else
                motion_detec_r(j) = 1;
            end
            
            if sum(motion_detec_r(j-150:j))==0
                motion_r(j) = 0;
            else 
                motion_r(j) =1;
            end
            
            if abs(delta_euler (3,j))<=delta_euler_motion_threshold
                motion_detec_dpsi(j) = 0;
            else
                motion_detec_dpsi(j)=1;
            end

            if sum(motion_detec_dpsi(j-150:j))==0
                motion_dpsi(j) = 0;
            else 
                motion_dpsi(j) =1;
            end
            
            if motion_r(j) ==0 && motion_a(j)==0 &&  motion_dpsi(j)==0
                motion(j) = 0;
            else
                motion(j) =1;
            end
 %%           
    if motion(1,j)== 0

%% time update  PREDICTION EKF1
    if ~any(isnan(psi_gps(1,j))) && vel_gps_mag(j)>=3
        as(j) = 1; %%% Case where  NO dynamic motion + Psi gps available (ACC + Psi from GPS)
    [funk_ekf1,Fk_ekf1,hunk_ekf1,Hk_ekf1,fun_ekf1] = EKF_compute_process_model_EKF1_case1( Ts, Xp_ekf1,ARS_input_ekf1');
%     hunk_ekf1_mem(1:4,j)=hunk_ekf1;
        delta_euler (:,j)=fun_ekf1;

            %%% Process Covariance Calculation
                  G = eye(6);                
            Qk_ekf1 = 0.5*(Fk_ekf1*G*Qc_ekf1_case1*G'+G*Qc_ekf1_case1*G'*Fk_ekf1')*Ts;
%            Qk_ekf1_mem(:,:,j)=Qk_ekf1;


    [Xkp_ekf1, Pkp_ekf1] = EKF_time_update(Xp_ekf1, Pp_ekf1, funk_ekf1, Fk_ekf1,Qk_ekf1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cov_ax = 50e-3;
        cov_ay = 50e-3;
        cov_az = 50e-3;
        cov_psi_gps =(deg2rad(100));
        Rd_ekf1 = [cov_ax^2 0 0 0;0 cov_ay^2 0 0;0 0 cov_az^2 0;0 0 0 cov_psi_gps^2]*1e0; 
    psi_gps(1,j) =psi_gps(1,j);
        Zk_ekf1 = [accx;accy;accz];

    Zk_ekf1 = [Zk_ekf1;psi_gps(1,j)];
    
    
    [Xk_ekf1,Pk_ekf1,S_ekf1,innovations_ekf1,K_ekf1] = EKF_correction_update_ekf1_case1(Xkp_ekf1, Pkp_ekf1, Zk_ekf1,hunk_ekf1, Hk_ekf1, Rd_ekf1);
%     in_correction_check_as1(j)=AA;
    Xk_ekf1(3) = Change_range_angle(Xk_ekf1(3), 1);
    
    
%     K_ekf1_mem_case1(:,:,j)=K_ekf1;

    Xp_ekf1 = Xk_ekf1;
    Pp_ekf1 = Pk_ekf1;

    Innov_mem_ekf1(1:4,j) = innovations_ekf1;
    Zk_mem_ekf1(1:4,j) = Zk_ekf1;
    Xk_mem_ekf1(:,j) = Xk_ekf1;
    Pk_mem_ekf1(:,:,j)=Pk_ekf1;
            funk_ekf1(3)=Change_range_angle(funk_ekf1(3), 1);    

    else % if Psi GPS not available
        as(j) = 2; %%% NO dynamic motion and Psi GPS not available (only ACC avialable)

        [funk_ekf1,Fk_ekf1,hunk_ekf1,Hk_ekf1,fun_ekf1] = EKF_compute_process_model_EKF1_case2( Ts, Xp_ekf1,ARS_input_ekf1');
        delta_euler (:,j)=fun_ekf1;
       G = eye(6);
       Qk_ekf1 = 0.5*(Fk_ekf1*G*Qc_ekf1_case2*G'+G*Qc_ekf1_case2*G'*Fk_ekf1')*Ts;
        [Xkp_ekf1, Pkp_ekf1] = EKF_time_update(Xp_ekf1, Pp_ekf1, funk_ekf1, Fk_ekf1,Qk_ekf1);
        
        cov_ax = 50e-3;
        cov_ay = 50e-3;
        cov_az = 50e-3;
        Rd_ekf1 = [cov_ax^2 0 0 ;0 cov_ay^2 0 ;0 0 cov_az^2]*1e0;
         Zk_ekf1 = [accx;accy;accz];
         
        [Xk_ekf1,Pk_ekf1,S_ekf1,innovations_ekf1,K_ekf1] = EKF_correction_update(Xkp_ekf1, Pkp_ekf1, Zk_ekf1,hunk_ekf1, Hk_ekf1, Rd_ekf1);
        Xk_ekf1(3) = Change_range_angle(Xk_ekf1(3), 1);
%     K_ekf1_mem_case2(:,:,j)=K_ekf1;

       
        Xp_ekf1 = Xk_ekf1;
        Pp_ekf1 = Pk_ekf1;
    %     % -------------------------- savings
        Innov_mem_ekf1(1:3,j) = innovations_ekf1;
        Zk_mem_ekf1(1:3,j) = Zk_ekf1;
        Xk_mem_ekf1(:,j) = Xk_ekf1;
        Pk_mem_ekf1(:,:,j)=Pk_ekf1;
        
        funk_ekf1(3)=Change_range_angle(funk_ekf1(3), 1);
    end
    else
        if ~any(isnan(psi_gps(1,j))) && vel_gps_mag(j)>=3
                    as(j) = 3; %%% YES dynamic motion but Psi GPS available (only Psi GPS avialable)

             [funk_ekf1,Fk_ekf1,hunk_ekf1,Hk_ekf1,fun_ekf1] = EKF_compute_process_model_EKF1_case3( Ts, Xp_ekf1,ARS_input_ekf1');
             delta_euler (:,j)=fun_ekf1;
       G = eye(6);
       Qk_ekf1 = 0.5*(Fk_ekf1*G*Qc_ekf1_case3*G'+G*Qc_ekf1_case3*G'*Fk_ekf1')*Ts;
            [Xkp_ekf1, Pkp_ekf1] = EKF_time_update(Xp_ekf1, Pp_ekf1, funk_ekf1, Fk_ekf1,Qk_ekf1);

            cov_psi_gps =deg2rad(100);
            Rd_ekf1 = diag(cov_psi_gps^2)*1e0;
            
            
            psi_gps(1,j) =psi_gps(1,j);
            Zk_ekf1 = psi_gps(1,j);

        [Xk_ekf1,Pk_ekf1,S_ekf1,innovations_ekf1,K_ekf1] = EKF_correction_update_ekf1_case3(Xkp_ekf1, Pkp_ekf1, Zk_ekf1,hunk_ekf1, Hk_ekf1, Rd_ekf1);
         Xk_ekf1(3) = Change_range_angle(Xk_ekf1(3), 1);
        
%         K_ekf1_mem_case3(:,:,j)=K_ekf1;


        Xp_ekf1 = Xk_ekf1;
        Pp_ekf1 = Pk_ekf1;
                
        Innov_mem_ekf1(4,j) = innovations_ekf1;
        Zk_mem_ekf1(4,j) = Zk_ekf1;
        Xk_mem_ekf1(:,j) = Xk_ekf1;
        Pk_mem_ekf1(:,:,j)=Pk_ekf1;
        
        funk_ekf1(3)=Change_range_angle(funk_ekf1(3), 1);

        else
            as(j) =4; %%% YES dynamic motion + NO Psi GPS (Free Integration)

        [funk_ekf1,Fk_ekf1,hunk_ekf1,Hk_ekf1,fun_ekf1] = EKF_compute_process_model_EKF1_case4( Ts, Xp_ekf1,ARS_input_ekf1');
        delta_euler (:,j)=fun_ekf1;

       G = eye(6);
       Qk_ekf1 = 0.5*(Fk_ekf1*G*Qc_ekf1_case4*G'+G*Qc_ekf1_case4*G'*Fk_ekf1')*Ts;
        %%%%%%%%
        [Xkp_ekf1, Pkp_ekf1] = EKF_time_update(Xp_ekf1, Pp_ekf1, funk_ekf1, Fk_ekf1,Qk_ekf1);
         

        Xk_ekf1 = Xkp_ekf1;
        Xk_ekf1(3) = Change_range_angle(Xk_ekf1(3), 1);
        Pk_ekf1 =Pkp_ekf1;
        

        Xp_ekf1 = Xk_ekf1;
        Pp_ekf1 = Pk_ekf1;

% % %     % -------------------------- savings
        Xk_mem_ekf1(:,j) = Xk_ekf1;
        Pk_mem_ekf1(:,:,j)=Pk_ekf1;
        
        funk_ekf1(3)=Change_range_angle(funk_ekf1(3), 1);
       end
    end   
    
%% EKF 2
% Input
    C2bn = Cb2n(Xk_ekf1(1),Xk_ekf1(2),Xk_ekf1(3));
    C_b2n_mem(:,:,j)=C2bn;
%     C2bn =  C_b2n_mem(:,:,j);
%%
    ars_x = ARS_input_ekf1(1);
    ars_y = ARS_input_ekf1(2);
    ars_z = ARS_input_ekf1(3);
    
    vbx = Xk_mem_ekf2(4,j-1);
    vby = Xk_mem_ekf2(5,j-1);
    vbz = Xk_mem_ekf2(6,j-1);
  Acentrifugal(:,j) = skew([vbx;vby;vbz])*([ars_x;ars_y;ars_z] - [Xk_ekf1(4); Xk_ekf1(5); Xk_ekf1(6)]);

        ACC = ACC-((Xk_mem_ekf2(7:9,j-1))./gN);
        
%         ACC_mem_input (1:3,j)=ACC; 
%         ACC = ACC;

%   aeb(:,j) = (ACC*gN)+C2bn'*[0;0;1*gN]+Acentrifugal(:,j);
  aeb(:,j) = (ACC*gN)+C2bn'*[0;0;1*gN]+Acentrifugal(:,j);

% %%  Orginal sampling Measurement
    PN_gps = PN_gps_mem(:,j);
%     
    V_NED_gps =[data(j,28) data(j,29) data(j,27)]';
    V_NED_gps_mem(1:3,j)=V_NED_gps;
    V_body = C2bn'*V_NED_gps;
    V_body_mem(:,j)=V_body;
    %%
% % % % %% EKF state correction UPDATE STAGE %%%%%%%%%%%%%%*******************************%%%%%%%%%%%%%%%%
Zk_ekf2 = [PN_gps;V_body];
    if ~any(isnan(Zk_ekf2(1:6))) && vel_gps_mag(j)>=3
        aq(j) =1; %%% GPS data available        
%% prediction

         
        [funk_ekf2,Fk_ekf2,hunk_ekf2,Hk_ekf2] = EKF_compute_process_model_EKF2_case1_plus_velgps_3axis( Ts, Xp_ekf2,C2bn,aeb(:,j));
 
       G_ekf2 = eye(9);
%        Qk_ekf2 = 0.5*(Fk_ekf2*G_ekf2*Qc_ekf2*G_ekf2'+G_ekf2*Qc_ekf2*G_ekf2'*Fk_ekf2')*Ts;
        Qk_ekf2=Qc_ekf2;

%     
    R_gps_x=3;
    R_gps_y=3;
    R_gps_z=1;
    R_gps_vx=2;
    R_gps_vy=2;
    R_gps_vz=1;
    
%         R_gps_x=0.08;
%     R_gps_y=0.08;
%     R_gps_z=0.08;
%     R_gps_vx=0.2;
%     R_gps_vy=0.08;
%     R_gps_vz=0.4;
    
        [Xkp_ekf2, Pkp_ekf2] = EKF_time_update(Xp_ekf2,Pp_ekf2,funk_ekf2,Fk_ekf2,Qk_ekf2);
        Rd_ekf2 = diag([ (R_gps_x)^2, (R_gps_y)^2, (R_gps_z)^2, (R_gps_vx)^2, (R_gps_vy)^2, (R_gps_vz)^2])*1;
        hk_ekf2 = hunk_ekf2;
%        
        [Xk_ekf2,Pk_ekf2,S_ekf2,innovations_ekf2,K_ekf2] = EKF_correction_update(Xkp_ekf2,Pkp_ekf2,Zk_ekf2,hk_ekf2,Hk_ekf2,Rd_ekf2);
%          K_ekf2_mem(1:4,j)=diag(K_ekf2);
        K_ekf2_mem(:,:,j)=K_ekf2;
     
%         %%  
        Xp_ekf2 = Xk_ekf2;
        Pp_ekf2 = Pk_ekf2;
   
        Innov_mem_ekf2(1:6,j) = innovations_ekf2;
        Zk_mem_ekf2(1:6,j) = Zk_ekf2;
        Xk_mem_ekf2(:,j) = Xk_ekf2;
        Pk_mem_ekf2(:,:,j) = Pk_ekf2;
    elseif ~any(isnan(Zk_ekf2(1:6))) && vel_gps_mag(j)< 3
                aq(j) =2;  %%% GPS DATA available ...GPS mag less than 6 m/s

        Zk_ekf2 = PN_gps;
        
%             R_gps_x=5;
%     R_gps_y=5;
%     R_gps_z=5;
    
               R_gps_x=3;
    R_gps_y=3;
    R_gps_z=1;
        
        [funk_ekf2,Fk_ekf2,hunk_ekf2,Hk_ekf2] = EKF_compute_process_model_EKF2_case1_plus_vel_mag( Ts, Xp_ekf2,C2bn,aeb(:,j));
 
       G_ekf2 = eye(9);
       
        Qk_ekf2=Qc_ekf2;
%        Qk_ekf2 = 0.5*(Fk_ekf2*G_ekf2*Qc_ekf2*G_ekf2'+G_ekf2*Qc_ekf2*G_ekf2'*Fk_ekf2')*Ts;

%     
        [Xkp_ekf2, Pkp_ekf2] = EKF_time_update(Xp_ekf2,Pp_ekf2,funk_ekf2,Fk_ekf2,Qk_ekf2);
        Rd_ekf2 = diag([ (R_gps_x)^2, (R_gps_y)^2, (R_gps_z)^2])*1;
        hk_ekf2 = hunk_ekf2;
%        
        [Xk_ekf2,Pk_ekf2,S_ekf2,innovations_ekf2,K_ekf2] = EKF_correction_update(Xkp_ekf2,Pkp_ekf2,Zk_ekf2,hk_ekf2,Hk_ekf2,Rd_ekf2);
        K_ekf2_mem(:,1:3,j)=K_ekf2;
     
%         %%  
        Xp_ekf2 = Xk_ekf2;
        Pp_ekf2 = Pk_ekf2;
   
        Innov_mem_ekf2(1:3,j) = innovations_ekf2;
        Zk_mem_ekf2(1:3,j) = Zk_ekf2;
        Xk_mem_ekf2(:,j) = Xk_ekf2;
        Pk_mem_ekf2(:,:,j) = Pk_ekf2;
   
    else
%         %%
        aq(j) =3;  %%% No GPS DATA available ...process model may be used for perfect measurment use
       
[funk_ekf2,Fk_ekf2,hunk_ekf2,Hk_ekf2] = EKF_compute_process_model_EKF2_case2_Pmeas( Ts, Xp_ekf2,C2bn,aeb(:,j));
%     % Covariance
       G_ekf2 = eye(9);
%        Qk_ekf2 = 0.5*(Fk_ekf2*G_ekf2*Qc_ekf2_gps_not_present*G_ekf2'+G_ekf2*Qc_ekf2_gps_not_present*G_ekf2'*Fk_ekf2')*Ts;
       Qk_ekf2 = Qc_ekf2*1e0;
%               Qk_ekf2 = 0.5*(Fk_ekf2*G_ekf2*Qc_ekf2*G_ekf2'+G_ekf2*Qc_ekf2*G_ekf2'*Fk_ekf2')*Ts;


% %   % time update  PREDICTION
    [Xkp_ekf2, Pkp_ekf2] = EKF_time_update(Xp_ekf2,Pp_ekf2,funk_ekf2,Fk_ekf2,Qk_ekf2);

  
        Xk_ekf2 = Xkp_ekf2;
        Pk_ekf2 = Pkp_ekf2;
        Xp_ekf2 = Xk_ekf2;
        Pp_ekf2 = Pk_ekf2;

        Pk_mem_ekf2(:,:,j) = Pk_ekf2;
        Xk_mem_ekf2(:,j) = Xk_ekf2;

    end
    V_NED_ekf_estimate_mem(1:3,j) =C2bn *Xk_mem_ekf2(4:6,j); 
    %% error calculation from RTK with estimation
    if ~any(isnan(PN_gps_RTK_mem(:,j)))
            error_pos_est_RTK(1:3,j) = PN_gps_RTK_mem(:,j)-Xk_mem_ekf2(1:3,j);
            error_pos_RTK_GPS(1:3,j) = PN_gps_RTK_mem(:,j)-PN_gps_mem(:,j);
    else
            error_pos_est_RTK(1:3,j) = nan;
            error_pos_RTK_GPS(1:3,j) =nan;
    end
    
    
    if ~any(isnan(PN_gps_mem(:,j)))
            error_pos_est_GPS(1:3,j) = PN_gps_mem(:,j)-Xk_mem_ekf2(1:3,j);
            error_pos_est_GPSvel(1:3,j) = V_NED_gps_mem(:,j)-V_NED_ekf_estimate_mem(1:3,j);
    else
            error_pos_est_GPS(1:3,j) = nan;
            error_pos_est_GPSvel(1:3,j) = nan;
    end
    error_pos_est_GPS_withalltime(1:3,j) = PN_gps_mem_alltime(:,j)-Xk_mem_ekf2(1:3,j);
    error_pos_est_GPSvel_withalltime(1:3,j) = V_NED_gps_mem_alltime(:,j)-V_NED_ekf_estimate_mem(1:3,j);
    error_pos_est_GPS_RTK_withalltime(1:3,j) = PN_gps_RTK_mem_alltime(:,j)-Xk_mem_ekf2(1:3,j);

%%
    
end
max(abs(error_pos_est_GPS_withalltime(1,:)))
max(abs(error_pos_est_GPS_withalltime(2,:)))
max(abs(error_pos_est_GPS_withalltime(3,:)))

max(abs(error_pos_est_GPS_RTK_withalltime(1,:)))
max(abs(error_pos_est_GPS_RTK_withalltime(2,:)))
max(abs(error_pos_est_GPS_RTK_withalltime(3,:)))

% save('MA_shift_0.mat','Xk_mem_ekf2','V_NED_ekf_estimate_mem','intervals')
%
%% 
figure;
set(gcf,'color','w');
x = -4:0.05:4;
subplot(421), plot(t,Innov_mem_ekf1(1,:),'.'), zoom on, grid on; title('Innovation for ACC (g)')
xlabel('Time (s)'); ylabel('ACC_x (g)'); 
% xlim([0 2000])
subplot(422), hist(Innov_mem_ekf1(1,:),x), title('Histogram');zoom on, grid on;
subplot(423), plot(t,Innov_mem_ekf1(2,:),'.'), zoom on, grid on; 
xlabel('Time (s)'); ylabel('ACC_y (g)'); 
% ylim([-0.8 0.8]);xlim([0 2000]);
subplot(424), hist(Innov_mem_ekf1(2,:),x), zoom on, grid on;
subplot(425), plot(t,Innov_mem_ekf1(3,:),'.'), zoom on, grid on;
% ylim([-0.8 0.8]);xlim([0 2000]);
xlabel('Time (s)'); ylabel('ACC_z (g)')
subplot(426), hist(Innov_mem_ekf1(3,:),x), zoom on, grid on;
subplot(427), plot(t,Innov_mem_ekf1(4,:)*180/pi,'.'), zoom on, grid on; title('Innovation for Heading ({\circ})')
xlabel('Time (s)'); ylabel('\psi ({\circ})'); 
% xlim([0 2000]); ylim([-0.25 0.25])
subplot(428), hist(Innov_mem_ekf1(4,:),x), grid on;

%% with dynamic motion included
figure;
set(gcf,'color','w');
subplot(311), plot(t, EUL_ACC_MAG_mem(1,:)*180/pi,'b',t,Xk_mem_ekf1(1,:)*180/pi,'r',t, EUL_GYR_mem (1,:)*180/pi, 'k'),
hold on
plot(t,motion*20,'g');
title('Euler angle (EA)')
ylabel('\phi ({\circ})');
xlabel('Time (s)');
legend('EA from ACC','EKF1 estimate','EA from ARS'), zoom on, grid on;
subplot(312), plot(t, EUL_ACC_MAG_mem(2,:)*180/pi,'b',t,Xk_mem_ekf1(2,:)*180/pi,'r',t, EUL_GYR_mem (2,:)*180/pi, 'k'),
hold on
plot(t,motion*20,'g');
xlabel('Time (s)');
ylabel('\theta ({\circ})');
zoom on, grid on;
subplot(313), plot(t, psi_gps*180/pi,'*b',t,Xk_mem_ekf1(3,:)*180/pi,'r',t, EUL_GYR_mem (3,:)*180/pi, 'k'),
hold on
plot(t,motion*400,'g');
xlabel('Time (s)');
ylabel('\psi ({\circ})');
zoom on, grid on;

%%
figure;
set(gcf,'color','w');
subplot(311), plot(t,Xk_mem_ekf1(4,:)*180/pi,'r'), zoom on, grid on;
hold on
% plot(t,motion*0.2,'g');
title('ARS bias ({\circ}/s)')
ylabel('ARS_x ({\circ}/s)');
xlabel('Time (s)');
subplot(312), plot(t,Xk_mem_ekf1(5,:)*180/pi,'r'), zoom on, grid on;
hold on
% plot(t,motion*0.1,'g');
xlabel('Time (s)');
ylabel('ARS_y ({\circ}/s)');
zoom on, grid on;
subplot(313), plot(t,Xk_mem_ekf1(6,:)*180/pi,'r'), zoom on, grid on;
hold on
% plot(t,motion*0.1,'g');
xlabel('Time (s)');
ylabel('ARS_z ({\circ}/s)');
%%
figure;
set(gcf,'color','w');
subplot(311), plot(t(1:M),PN_gps_RTK_mem_alltime(1,1:M),'*k'), zoom on, grid on;
hold on
plot(t(1:M),PN_gps_mem_alltime(1,1:M),'.b')
plot(t(1:M),PN_gps_mem(1,1:M),'.g')
plot(t(1:M),Xk_mem_ekf2(1,1:M),'r')
legend('RTK All Time','PN-GPS all time','GPS withoutage','Estimated (EKF2)');
title('Position (NED)')
xlabel('Time (s)');ylabel('North (m)')
subplot(312), plot(t(1:M),PN_gps_RTK_mem_alltime(2,1:M),'*k'), zoom on, grid on;
hold on
plot(t,PN_gps_mem_alltime(2,1:M),'.b')
plot(t,PN_gps_mem(2,1:M),'.g')
plot(t,Xk_mem_ekf2(2,1:M),'r')
xlabel('Time (s)');ylabel('East (m)')
subplot(313), plot(t,PN_gps_RTK_mem_alltime(3,1:M),'*k'), zoom on, grid on;
hold on
plot(t(1:M),PN_gps_mem_alltime(3,1:M),'.b')
plot(t(1:M),PN_gps_mem(3,1:M),'.g')
plot(t(1:M),Xk_mem_ekf2(3,1:M),'r')
% legend('PD-GPS','D-est (EKF2)','Reference');
xlabel('Time (s)');ylabel('Down (m)')

% 
figure;
set(gcf,'color','w');
subplot(311), plot(t(1:M),V_body_mem(1,1:M),'.b'), zoom on, grid on; 
hold on
plot(t(1:M),Xk_mem_ekf2(4,:),'r')
title('Velocity, bodyframe')
legend('V_body All time','Velocity body - GPS','Estimated (EKF2)');
xlabel('Time (s)');ylabel('V_x (m/s)')
subplot(312), plot(t,V_body_mem(2,1:M),'.b'), zoom on, grid on;
hold on
plot(t(1:M),Xk_mem_ekf2(5,:),'r')
% legend('Estimated (EKF2)','Reference');
xlabel('Time (s)');ylabel('V_y (m/s)')
subplot(313), plot(t,V_body_mem(3,1:M),'.b'), zoom on, grid on;
hold on
plot(t(1:M),Xk_mem_ekf2(6,:),'r')
% legend('Estimated (EKF2)','Reference');
xlabel('Time (s)');ylabel('V_z (m/s)')

% 
figure;
set(gcf,'color','w');
subplot(311), plot(t,(Xk_mem_ekf2(7,:)./gN)*1000,'r'), zoom on, grid on;
title('ACC bias (mg)');
xlabel('Time (s)');ylabel('ACC_x (mg)')
subplot(312), plot(t,(Xk_mem_ekf2(8,:)./gN)*1000,'r'), zoom on, grid on;
xlabel('Time (s)');ylabel('ACC_y (mg)')
subplot(313), plot(t,(Xk_mem_ekf2(9,:)./gN)*1000,'r'), zoom on, grid on;
xlabel('Time (s)');ylabel('ACC_z (mg)')


figure;
x = -4:0.05:4;
subplot(421), plot(t,Innov_mem_ekf2(1,:),'.'), zoom on, grid on;
% ylim([-0.5 0.5]);xlim([0 1200]);
title('Innovation for Pn Pe Pd (m)'); xlabel('Time (s)'); ylabel('Postion - N (m)')
subplot(422), hist(Innov_mem_ekf2(1,:),x), zoom on, grid on;
% xlim([-5 5]);
title('Histogram')
subplot(423), plot(t,Innov_mem_ekf2(2,:),'.'), zoom on, grid on; 
xlabel('Time (s)'); ylabel('Postion - E (m)')
% ylim([-0.5 0.5]);xlim([0 1200]);
subplot(424), hist(Innov_mem_ekf2(2,:),x), zoom on, grid on;
% xlim([-5 5]);
subplot(425), plot(t,Innov_mem_ekf2(3,:),'.'), zoom on, grid on; 
xlabel('Time (s)'); ylabel('Postion - D (m)')
% ylim([-0.25 0.25]);xlim([0 1200]);
subplot(426), hist(Innov_mem_ekf2(3,:),x), zoom on, grid on;
% xlim([-5 5]);
subplot(427), plot(t,Innov_mem_ekf2(4,:),'.'), zoom on, grid on; 
xlabel('Time (s)'); ylabel('Velocity - D (m/s)')
% ylim([-0.25 0.25]);xlim([0 1200]);
subplot(428), hist(Innov_mem_ekf2(4,:),x), zoom on, grid on;

figure;
set(gcf,'color','w');
subplot(311), plot(t(1:M),V_NED_gps_mem_alltime(1,1:M),'*g',t(1:M),V_NED_ekf_estimate_mem(1,1:M),'r',t(1:M),V_NED_gps_mem(1,1:M),'.b'), zoom on, grid on; legend('GPS NED all time','Estimated (EKF2) NED','Velocity NED - GPS');
title('Velocity, NED')
xlabel('Time (s)');ylabel('V_N (m/s)')
subplot(312), plot(t(1:M),V_NED_gps_mem_alltime(2,1:M),'*g',t(1:M),V_NED_ekf_estimate_mem(2,:),'r',t(1:M),V_NED_gps_mem(2,1:M),'.b'), zoom on, grid on; 
% legend('Estimated (EKF2)','Reference');
xlabel('Time (s)');ylabel('V_E (m/s)')
subplot(313), plot(t(1:M),V_NED_gps_mem_alltime(3,1:M),'*g',t(1:M),V_NED_ekf_estimate_mem(3,1:M),'r',t(1:M),V_NED_gps_mem(3,1:M),'.b'), zoom on, grid on; 
% legend('Estimated (EKF2)','Reference');
xlabel('Time (s)');ylabel('V_D (m/s)')

figure;
x = -4:0.05:4;
subplot(321), plot(t,error_pos_est_GPS_withalltime(1,:),'.r'), zoom on, grid on;
hold on
plot(t,error_pos_est_GPS(1,:),'.')
% ylim([-0.5 0.5]);xlim([0 1200]);
title('Innovation for Pn Pe Pd (m) in post estimation'); xlabel('Time (s)'); ylabel('Postion - N (m)')
subplot(322), hist(error_pos_est_GPS(1,:),x), zoom on, grid on;
% xlim([-5 5]);
title('Histogram')
subplot(323), plot(t,error_pos_est_GPS_withalltime(2,:),'.r'), zoom on, grid on; 
hold on
 plot(t,error_pos_est_GPS(2,:),'.')
xlabel('Time (s)'); ylabel('Postion - E (m)')
% ylim([-0.5 0.5]);xlim([0 1200]);
subplot(324), hist(error_pos_est_GPS(2,:),x), zoom on, grid on;
% xlim([-5 5]);
subplot(325), plot(t,error_pos_est_GPS_withalltime(3,:),'.r'), zoom on, grid on; 
hold on
plot(t,error_pos_est_GPS(3,:),'.')
xlabel('Time (s)'); ylabel('Postion - D (m)')
% ylim([-0.25 0.25]);xlim([0 1200]);
subplot(326), hist(error_pos_est_GPS(3,:),x), zoom on, grid on;


figure;
x = -4:0.05:4;
subplot(321), plot(t,error_pos_est_GPSvel_withalltime(1,:),'.r'), zoom on, grid on;
hold on
plot(t,error_pos_est_GPSvel(1,:),'.')
% ylim([-0.5 0.5]);xlim([0 1200]);
title('Innovation for Vn Ve Vd (m) in vel estimation'); xlabel('Time (s)'); ylabel('Velocity - N (m/s)')
subplot(322), hist(error_pos_est_GPSvel(1,:),x), zoom on, grid on;
% xlim([-5 5]);
title('Histogram')
subplot(323), plot(t,error_pos_est_GPSvel_withalltime(2,:),'.r'), zoom on, grid on;
hold on
plot(t,error_pos_est_GPSvel(2,:),'.')
xlabel('Time (s)'); ylabel('Velocity - E (m/s)')
% ylim([-0.5 0.5]);xlim([0 1200]);
subplot(324), hist(error_pos_est_GPSvel(2,:),x), zoom on, grid on;
% xlim([-5 5]);
subplot(325), plot(t,error_pos_est_GPSvel_withalltime(3,:),'.r'), zoom on, grid on; 
hold on
plot(t,error_pos_est_GPSvel(3,:),'.')
xlabel('Time (s)'); ylabel('Velocity - D (m/s)')
% ylim([-0.25 0.25]);xlim([0 1200]);
subplot(326), hist(error_pos_est_GPSvel(3,:),x), zoom on, grid on;
