% -------------------------------------------------------------------------
%        'Gated Sensor Fusion: A way to Improve the Precision of 
%                Ambulatory Human Body Motion Estimation'
% -------------------------------------------------------------------------
%
% *************************************************************************
% - Author: Alberto Olivares-Vicente.
% - Entity: University of Granada, Spain.
% - Last revision: 01/22/2015.
% *************************************************************************
%
% - DESCRIPTION:
% The following file is the main routine of the experiments reflected in 
% the paper. 
% The main routine is composed of a series of algorithms which compute the
% orientation of a Wagyromag Magnetic Inertial Measurement Unit attached to
% a mechanical device which provides the real value of the orientation 
% angle.
% Gathered data are loaded from .CSV files, calibrated and used as the 
% input of the different attitude/orientation algorithms. An optimization
% process is carried out to minimize the RMSE of the estimated angle with
% respect to the reference angle. The algorithms are included inside a
% Monte Carlo simulation so the average RMSE for each attitude estimation
% algorithm is provided at the end of the execution. 

%% 0) GENERAL CONFIGURATION \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% -------------------------------------------------------------------------
clear all; close all; clc;

% Load Wagyromag's functions library.
wag = wagLibrary;

% Set flag which indicates whether internal plots of the algorithms are to
% be displayed or not. Possible options : 'yes', 'no'.
showInternalPlots = 'no';

% Set flag which indicates whether the final results of the algorithm are
% to be displayed or not. Possible options : 'yes', 'no'.
showPlots = 'yes';

% Set sampling frequency with which the data were gathered by the MIMU.
f = 40;               

% Set value of the magnitude of the gravity vector in the location in which
% data were gathered. (In our case: Granada, Spain, 37°10'4''N 3°36'3''O, 
% 738 meters over sea level). 
g = 9.797024;

% Set RMSE offset. The computation of the RMSE is done from the Xth signal
% to the Nth signal, where X is an initial offset and N is the length of
% the signal. This is done to allow slower filters to reach convergence.
rmse_offset = 300;

%% 1) LOAD AND CALIBRATE RAW DATA \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% -------------------------------------------------------------------------
%  |
%  |_ This section reads the following data from a .CSV file which was
%     previously gathered using the MIMU. Raw data are stored in the 
%     following vectors of size (1 x N), where N is the total number of 
%     samples. 
%      - ax: Acceleration measured along X axis. Raw values are within the 
%        [0-1023] range. 
%      - ay: Acceleration measured along Y axis. Raw values are within the 
%        [0-1023] range.
%      - az: Acceleration measured along Z axis. Raw values are within the 
%        [0-1023] range.
%      - gx: Angular rate measured along X axis. Raw values are within the 
%        [0-1023] range.
%      - gy: Angular rate measured along Y axis. Raw values are within the 
%        [0-1023] range.
%      - gz: Angular rate measured along Z axis. Raw values are within the 
%        [0-1023] range.
%      - hx: Magnetic field measured along X axis. Raw values are within 
%        the [0-1023] range.
%      - hy: Magnetic field measured along Y axis. Raw values are within 
%        the [0-1023] range.
%      - hz: Magnetic field measured along Z axis. Raw values are within 
%        the [0-1023] range.
%      - Angle: Actual orientation. This vector is used as a reference to
%        check the accuracy of the different methods.

% -------------------------------------------------------------------------
% 1.2) Browse for data file and load data.
% -------------------------------------------------------------------------

% Show the user a folder browser so he can select the folder where the data
% are stored. The '\data' folder is opened by default. 
folder_name = uigetdir('data');

% Get the number of files in the folder. The command 'dir' builds a 
% structure containing different information about the folder. The field 
% 'name' of the structure is a string in which the two first elements are 
% '.' and '..' and the remaining ones are the names of the .CSV files. The
% goal here is to build a for loop which automatically loads all the files
% in the folder.
list = dir(folder_name); 

% -------------------------------------------------------------------------
% 1.3) Initialize Monte Carlo Simulation variables.
% -------------------------------------------------------------------------
% The for loop which iterates through all the signals in the folder is
% actually a Monte Carlo Simulation, since we will apply all the algorithms
% over all the signals, store the results and, at the end of the loop, we
% will average the results. Here we initizalize all the vectors where the
% optimal parameters of each algorithm as well as the resulting RMSE are
% stored. 

% QUEST algorithm: It has two parameters: 'a1' and 'a2'.
a1s_opt_QUEST = zeros(1, length(list) - 2);
a2s_opt_QUEST = zeros(1, length(list) - 2);
rmse_QUEST = zeros(1, length(list) - 2);

% Regular KF algorithm: It has two parameters: 'alpha' and 'beta'.
opt_alphas_KF = zeros(1, length(list) - 2);
opt_betas_KF = zeros(1, length(list) - 2);
rmse_KF = zeros(1, length(list) - 2);

% Regular KF algorithm + QUEST: It has two parameters: 'alpha' and 'beta'.
opt_alphas_KF_QUEST = zeros(1, length(list) - 2);
opt_betas_KF_QUEST = zeros(1, length(list) - 2);
rmse_KF_QUEST = zeros(1, length(list) - 2);

% Extended Kalman Filter algorithm: It has two parameters: 'sigma_w' and
% 'sigma_obs'.
opt_sigmas_w_EKF = zeros(1, length(list) - 2);
opt_sigmas_obs_EKF = zeros(1, length(list) - 2);
rmse_EKF = zeros(1, length(list) - 2);

% Extended Kalman Filter algorithm + QUEST: It has two parameters:
% 'sigma_w' and 'sigma_obs'.
opt_sigmas_w_EKF_QUEST = zeros(1, length(list) - 2);
opt_sigmas_obs_EKF_QUEST = zeros(1, length(list) - 2);
rmse_EKF_QUEST = zeros(1, length(list) - 2);

% Gated Kalman Filter algorithm: It has four parameters: 'alpha1',
% 'alpha2', 'beta1' and 'beta2':
opt_alphas1_GKF = zeros(1, length(list) - 2);
opt_alphas2_GKF = zeros(1, length(list) - 2);
opt_betas1_GKF = zeros(1, length(list) - 2);
opt_betas2_GKF = zeros(1, length(list) - 2);
rmse_GKF = zeros(1, length(list) - 2);

% Gated Kalman Filter + QUEST algorithm: It has four parameters: 'alpha1',
% 'alpha2', 'beta1' and 'beta2':
opt_alphas1_GKF_QUEST = zeros(1, length(list) - 2);
opt_alphas2_GKF_QUEST = zeros(1, length(list) - 2);
opt_betas1_GKF_QUEST = zeros(1, length(list) - 2);
opt_betas2_GKF_QUEST = zeros(1, length(list) - 2);
rmse_GKF_QUEST = zeros(1, length(list) - 2);

% Gated Extended Kalman Filter: It has four parameters: 'sigma_w1',
% 'sigma_w2', 'sigma_obs1' and 'sigma_obs2'.
opt_sigmas_w1_GEKF = zeros(1, length(list) - 2);
opt_sigmas_obs1_GEKF = zeros(1, length(list) - 2);
opt_sigmas_w2_GEKF = zeros(1, length(list) - 2);
opt_sigmas_obs2_GEKF = zeros(1, length(list) - 2);
rmse_GEKF = zeros(1, length(list) - 2);

% Gated Extended Kalman Filter + QUEST: It has four parameters: 'sigma_w1',
% 'sigma_w2', 'sigma_obs1' and 'sigma_obs2'.
opt_sigmas_w1_GEKF_QUEST= zeros(1, length(list) - 2);
opt_sigmas_obs1_GEKF_QUEST = zeros(1, length(list) - 2);
opt_sigmas_w2_GEKF_QUEST = zeros(1, length(list) - 2);
opt_sigmas_obs2_GEKF_QUEST = zeros(1, length(list) - 2);
rmse_GEKF_QUEST = zeros(1, length(list) - 2);

% Improvement vectors.
imp_GKF = zeros(1, length(list) - 2);
imp_GKF_QUEST = zeros(1, length(list) - 2);
imp_GEKF = zeros(1, length(list) - 2);
imp_GEKF_QUEST = zeros(1, length(list) - 2);

% -------------------------------------------------------------------------
% 1.4) Monte Carlo Simulation.
% -------------------------------------------------------------------------
% Start the for loop which iterates through all the data files.
for iter = 3:length(list)

    % Add a slash to complete the folder name.
    folder_name = strcat(folder_name, '/');
    
    % Get the file name.
    file_name = list(iter).name;

    % Print the name of the file being loaded.
    fprintf('Loading %s\n', file_name); 

    % Load raw data from selected file.
    [ax, ay, az, gx, gy, gz, hx, hy, hz, dummy1, dummy2, Angle] = ...
        wag.load_data(fullfile(folder_name, file_name)); 

    % ---------------------------------------------------------------------
    % 1.3) Check the name of the Euler rotation angle.
    % ---------------------------------------------------------------------
    % Since the mechanical device to which the MIMU is attached only allows
    % the rotation around a single Euler angle, each one of the files will
    % contain random rotations around just one Euler angle. The name of the
    % .CSV file contains the name of the Euler angle so we can extract it.
    if strfind(file_name, 'pitch')
        selectedAngle = 'pitch';
    elseif strfind(file_name, 'roll')
        selectedAngle = 'roll';
    elseif strfind(filen_ame, 'yaw')
        selectedAngle = 'yaw';
    end

    % Build an iterator which starts by j = 1 and increases with the value 
    % of the variable 'iter'. This is done to avoid using 'iter - 2' 
    % everytime we want to store a value in a vector.
    j = iter - 2;

    % ---------------------------------------------------------------------
    % 1.4) Calibrate raw data. 
    % ---------------------------------------------------------------------
    % Once data are loaded, we proceed to calibrate them, i.e. we transform
    % the raw units [0-1023] into meaningful physical units (g, deg/s and 
    % Gauss). 

    % Calibrate reference angle. The reference angle also needs to be
    % translated from raw units to degrees. To that purpose, we carried out
    % the characterization of the potentiometer included in the mechanical
    % device in order to build a Look Up Table in which raw values are
    % associated to inclination angles in degrees. The 'calibDigAng'
    % function applies such a translation.
    true_ang = wag.calibDigAng(Angle);  

    % Load Acceleromter, Magnetometer and Gyroscope calibration parameters.
    load 'data/calib_wagyromag.mat'   

    % Calibrate ACCELERATION.
    [axC, ayC, azC] = wag.calibrateAcc(ax, ay, az, Ka, Ra, ba); 

    % Truncate noise (hardcore style).
    axC(abs(axC(:)) < 1e-6) = 0;      
    ayC(abs(ayC(:)) < 1e-6) = 0;      
    azC(abs(azC(:)) < 1e-6) = 0;  

    % Transform units (from 'm/s^2' to 'g').
    axC = axC / g;                    
    ayC = ayC / g;                   
    azC = azC / g;                   

    % Axis adjustment (fitting the MIMU's body frame to the inertial 
    % frame). In this case axes X and Z are interchanged. 
    axCaux = azC;       
    azCaux = axC;       
    axC = axCaux;         
    azC = azCaux;

    % Calibrate MAGNETIC FIELD. 
    [hxC, hyC, hzC] = wag.calibrateMag(hx, hy, hz, T, S, b); 

    % Axis adjustment (fitting the MIMU's body frame to the inertial 
    % frame). In this case, the sense of axes Y and Z is inverted. 
    hyC = -hyC;                       
    hzC = -hzC;                       

    % Calibrate ANGULAR RATE. 
    [gxC, gyC, gzC] = wag.calibrateGyro(gx, gy, gz, Kg, Rg, bg); 

    % Compensate remaining deviation (hardcore style). This is done since
    % in our case the MIMU was static and already stable at the beginning
    % of the data gathering process.
    gxC = gxC - gxC(1);          
    gyC = gyC - gyC(1);         
    gzC = gzC - gzC(1);  

    % Transform units from deg/s to rad/s.
    gxC = gxC * pi/180;          
    gyC = gyC * pi/180;         
    gzC = -gzC * pi/180;         

    % Build the acceleration, angular rate and magnetic field matrices.
    Accelerometer = [axC ayC azC]; 
    Gyroscope = [gxC gyC -gzC];    
    Magnetometer = [hxC' hyC' hzC']; 

    % Build the time vector.
    time = zeros(1, length(axC))';  
    for k = 2: length(axC)
        time(k) = time(k - 1) + 1 / f;
    end

    % ---------------------------------------------------------------------
    % 1.5) Plot reference angle.
    % ---------------------------------------------------------------------
    if strcmpi(showPlots, 'yes')
        angles_figure = figure(1);
        wag.create_figures('first', 'Reference', time, true_ang, 'black')
    end
    
    %% 2) ATTITUDE ESTIMATION \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % ---------------------------------------------------------------------
    %  |
    %  |_ Once all data are properly compensated and calibrated, we proceed
    %     to estimate the orientation of the body which is being monitored. 
    %     We will use different methods since the goal is to carry out a 
    %     comparative study among them to determine their performance and 
    %     accuracy under different motion conditions. 

    %% --------------------------------------------------------------------
    % 2.1) Projection of gravity and magnetic field. 
    % ---------------------------------------------------------------------
    % Compute pitch, roll and yaw. Pitch and roll are computed using the
    % gravity projections and yaw is computed using the magnetic field
    % projections. 
    [pitchAcc, rollAcc] = wag.pitch_roll_decomp(axC, ayC, azC, 'rad'); 
    [yawAccM, magXh, magYh] = wag.yaw_decomp(pitchAcc, rollAcc, hxC',...
        hyC', hzC');

    % Compute RMSE with respect to the reference angle and plot signals. 
    % The RMSE is computed from the 300th sample. This is done to allow 
    % slower algorithms reach their convergence region and therefore allow
    % for a fairer comparison between methods. 
    if strcmpi(selectedAngle, 'pitch')
        acc_angle = pitchAcc;
    elseif strcmpi(selectedAngle, 'roll')
        acc_angle = rollAcc;
    elseif strcmp(selectedAngle, 'yaw')
        acc_angle = yawAccM;
    end
    rmse_acc = wag.compute_rmse(true_ang, 180 / pi * acc_angle, ...
        rmse_offset);
    
    % Plot estimated orientation angle. 
    if strcmpi(showPlots, 'yes')
        wag.create_figures('other','Acc+Mag. Projections', time, ...
            180 / pi * acc_angle, '--b'); 
    end
  
    %% --------------------------------------------------------------------
    % 2.2) QUEST.
    % ---------------------------------------------------------------------
    % The QUEST algorithm was proposed by M. D. Shuster and S. D. Oh in the
    % work entitled: "Three-axis attitude determination from vector 
    % observations". Journal of Guidance Control and Dynamics, 4(1):70–77, 
    % 1981.
    
    % The QUEST algorithm has two variable parameters: 'a1' and 'a2'. These
    % parameters control the weight which is given to the measured
    % acceleration and magnetic field, respectively, to the estimation of
    % the orientation angle. Therefore, in order to find their optimal
    % value for each one of the evaluated signals, we carry out an
    % optimization process which aims to find the optimal parameters by
    % minimizing the RMSE with respect to the reference. By doing this we
    % ensure that each algorithm is executed with its optimal parameters.
    p0_QUEST = [0.3 0.9];    
    [xmin, fmin, ct] = wag.optimizeQUEST(axC, ayC, azC, hxC, hyC, hzC, ...
    true_ang, selectedAngle, p0_QUEST, rmse_offset);
    
    fprintf('------------------QUEST OPTIMIZATION----------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf('Optimal parameters are: \n -A1: %0.4f\n -A2: %0.4f\n',...
        xmin(1), xmin(2))
    fprintf('----------------------------------------------------------\n')

    % Extract and store optimal parameters found by the optimization 
    % routine.
    a1_opt_QUEST = xmin(1);
    a2_opt_QUEST = xmin(2);

    a1s_opt_QUEST(j) = a1_opt_QUEST;
    a2s_opt_QUEST(j) = a2_opt_QUEST;

    % Call again the QUEST algorithm to compute the angle estimate using
    % the optimal parameters. 
    [roll_QUEST, pitch_QUEST, yaw_QUEST, q_quest] = wag.quest(axC, ayC,...
        azC, hxC, hyC, hzC, a1_opt_QUEST, a2_opt_QUEST);

    % Compute RMSE and plot signals.
    if strcmp(selectedAngle, 'pitch')
        angle_QUEST = pitch_QUEST;      
    elseif strcmp(selectedAngle, 'roll')
        angle_QUEST = roll_QUEST;  
    elseif strcmp(selectedAngle, 'yaw')
        angle_QUEST = yaw_QUEST;
    end
    rmse_QUEST(j) = wag.compute_rmse(true_ang, 180 / pi * angle_QUEST', ...
        rmse_offset);
    
    if strcmpi(showPlots, 'yes')
        wag.create_figures('other', 'QUEST', time, ...
            180 / pi * angle_QUEST', 'red'); 
    end
    
    %% --------------------------------------------------------------------
    % 2.3) Sensor fusion I: Regular Kalman Filter (projections).
    % ---------------------------------------------------------------------
    % The first sensor fusion algorithm that we will test is the Regular
    % Kalman Filter using the orientation angle computed with the
    % acceleration projections method as the observation. 
    
    % 2.3.1) Definition of the variables of the optimization process. 
    if strcmp(selectedAngle, 'pitch')
        obs_KF = pitchAcc;
        gyro_KF = gyC;
    elseif strcmp(selectedAngle, 'roll')
        obs_KF = rollAcc;
        gyro_KF = gxC;
    elseif strcmp(selectedAngle, 'yaw')
        obs_KF = yawAccM;
        gyro_KF = gzC;
    end

    % 2.3.2) Parameter optimization.
    
    % Set initial value of parameters;
    p0_KF = [50 0.1];
    
    % Call the optimization routine.
    [xmin, fmin, ct] = wag.optimizeKF(obs_KF, gyro_KF, f, true_ang, ...
        p0_KF, rmse_offset);
    
    fprintf('--------------------KF OPTIMIZATION-----------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf('Optimal parameters are: \n -Alpha: %0.4f\n -Beta: %0.4f\n',...
        xmin(1), xmin(2))
    fprintf('----------------------------------------------------------\n')

    % Extract optimal parameters.
    opt_alpha_KF = xmin(1);
    opt_beta_KF = xmin(2);
    
    % Save optimal parameters.
    opt_alphas_KF(j) = opt_alpha_KF;
    opt_betas_KF(j) = opt_beta_KF;
    
    % 2.3.3) Computation of the orientation estimation using optimal 
    %        parameters.
    angle_KF = wag.fusionKF(gyro_KF, obs_KF, f, opt_alpha_KF, opt_beta_KF);

    % 2.3.4) Compute RMSE and plot signals.
    rmse_KF(j) = wag.compute_rmse(true_ang, 180 / pi * angle_KF, ...
        rmse_offset);
    if strcmpi(showPlots, 'yes')
        wag.create_figures('other', 'Regular KF', time, ...
            180 / pi * angle_KF, 'cyan'); 
    end

    %% --------------------------------------------------------------------
    % 2.4) Sensor fusion II: Regular Kalman Filter (QUEST).
    % ---------------------------------------------------------------------
    % The second sensor fusion algorithm that we will test is the Regular
    % Kalman Filter using the orientation angle computed with the QUEST 
    % algorithm as the observation .
    
    % 2.4.1) Definition of the variables of the optimization process. 
    if strcmp(selectedAngle, 'pitch')
        obs_KF_QUEST = real(pitch_QUEST);
        gyro_KF_QUEST = gyC;
    elseif strcmp(selectedAngle, 'roll')
        obs_KF_QUEST = real(roll_QUEST);
        gyro_KF_QUEST = gxC;
    elseif strcmp(selectedAngle, 'yaw')
        obs_KF_QUEST = real(yaw_QUEST);
        gyro_KF_QUEST = gzC;
    end

    % 2.4.2) Parameter optimization.
    
    % Set initial value of parameters;
    p0_KF_QUEST = [100 0.1];
    
    % Call the optimization routine.
        % Call the optimization routine.
    [xmin, fmin, ct] = wag.optimizeKF(obs_KF_QUEST, gyro_KF_QUEST, f, ...
        true_ang, p0_KF_QUEST, rmse_offset);
    
    fprintf('------------------KF-QUEST OPTIMIZATION-------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf('Optimal parameters are: \n -Alpha: %0.4f\n -Beta: %0.4f\n',...
        xmin(1), xmin(2))
    fprintf('----------------------------------------------------------\n')
    
    % Extract optimal parameters.
    opt_alpha_KF_QUEST = xmin(1);
    opt_beta_KF_QUEST = xmin(2);
    
    % Save optimal parameters.
    opt_alphas_KF_QUEST(j)= opt_alpha_KF_QUEST;
    opt_betas_KF_QUEST(j) = opt_beta_KF_QUEST;
        
    % 2.4.3) Computation of the orientation estimation using optimal 
    %        parameters.  
    angle_KF_QUEST = wag.fusionKF(gyro_KF_QUEST, obs_KF_QUEST, f, ...
        opt_alpha_KF_QUEST, opt_beta_KF_QUEST);

    % 2.4.4) Compute RMSE and plot signals.
    rmse_KF_QUEST(j) = wag.compute_rmse(true_ang, 180 / pi * angle_KF, ...
        rmse_offset);
    if strcmpi(showPlots, 'yes')
        wag.create_figures('other', 'Regular KF-QUEST', time, ...
            180 / pi * angle_KF, 'yellow'); 
    end

    %% --------------------------------------------------------------------
    % 2.5) Sensor fusion III: Extended Kalman Filter I.
    % ---------------------------------------------------------------------
    % The third sensor fusion algorithm that we will test is the Extended
    % Kalman Filter using the orientation angle computed with the 
    % accelerometer projections method as the observation.

    % 2.5.1) Parameter optimization.
    
    % Set initial value of parameters.
    p0_EKF = [0.02 100];
    
    % Call optimization routine.
    [xmin, fmin, ct] = wag.optimizeEKF(axC, ayC, azC, gxC, gyC, gzC, ...
    f, true_ang, p0_EKF, rmse_offset, selectedAngle);

    fprintf('---------------------EKF OPTIMIZATION---------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf(['Optimal parameters are: \n -Sigma_w: %0.4f\n', ...
        ' -Sigma_obs: %0.4f\n'], xmin(1), xmin(2))
    fprintf('----------------------------------------------------------\n')
    
    % Extract optimal parameters.
    opt_sigma_w_EKF = xmin(1);
    opt_sigma_obs_EKF = xmin(2);
    
    % Save optimal parameters.
    opt_sigmas_w_EKF(j) = opt_sigma_w_EKF;
    opt_sigmas_obs_EKF(j) = opt_sigma_obs_EKF;

    % 2.5.2) Computation of the orientation estimation using optimal 
    %        parameters.
    a_EKF = [azC ayC axC];
    w_EKF = [gxC gyC gzC];
    ini = [1 0 0 0];
    dt = 1 / f;
    q = wag.fusionEKF(-a_EKF, w_EKF, dt, opt_sigma_w_EKF, ...
        opt_sigma_obs_EKF, ini);

    % Transformation of the orientation quaternion into Euler angles.
    roll_EKF = zeros(1, length(q));
    pitch_EKF = zeros(1, length(q));
    yaw_EKF = zeros(1, length(q));
    for i = 1 : length(q)
         eulerAngles = wag.quat_to_euler(q(i, :));
         roll_EKF(i) = eulerAngles(1);
         pitch_EKF(i) = eulerAngles(2);
         yaw_EKF(i) = eulerAngles(3);
    end

    % 2.5.3) Compute the RMSE and plot signals.
    if strcmp(selectedAngle, 'pitch')
        angle_EKF = pitch_EKF;
    elseif strcmp(selectedAngle, 'roll')
        angle_EKF = roll_EKF;
    elseif strcmp(selectedAngle, 'yaw')
        angle_EKF = yaw_EKF;
    end
    rmse_EKF(j) = wag.compute_rmse(true_ang, 180 / pi * angle_EKF', ...
        rmse_offset);
    if strcmpi(showPlots, 'yes')
        wag.create_figures('other', 'EKF', time, 180 / pi * angle_EKF, ...
            'green'); 
    end
  
    %% --------------------------------------------------------------------
    % 2.6) Sensor fusion IV: Extended Kalman Filter (QUEST).
    % ---------------------------------------------------------------------
    % The fourth sensor fusion algorithm that we will test is the Extended
    % Kalman Filter using the orientation angle computed with the QUEST
    % algorithm as the observation.

    % 2.6.1) Parameter optimization.
    
    % Set initial value of parameters.
    gyro_EKF_QUEST = [gxC gyC gzC];
    p0_EKF_QUEST = [0.1 0.01];
     
    % Call optimization routine.
    [xmin, fmin, ct] = wag.optimizeEKF_QUEST(gyro_EKF_QUEST, q_quest, ...
        p0_EKF_QUEST, true_ang, selectedAngle, f, rmse_offset);

    fprintf('------------------EKF-QUEST OPTIMIZATION------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf(['Optimal parameters are: \n -Sigma_w: %0.4f\n', ...
        ' -Sigma_obs: %0.4f\n'], xmin(1), xmin(2))
    fprintf('----------------------------------------------------------\n')
    
    % Extract optimal parameters.
    opt_sigma_w_EKF_QUEST = xmin(1);
    opt_sigma_obs_EKF_QUEST = xmin(2);

    % Save optimal parameters.
    opt_sigmas_w_EKF_QUEST(j) = opt_sigma_w_EKF_QUEST;
    opt_sigmas_obs_EKF_QUEST(j) = opt_sigma_obs_EKF_QUEST;

    % 2.6.2) Computation of the orientation estimation using optimal 
    %        parameters.
    quat_ini = q_quest(1,:);
    quat_est= wag.fusionEKF_QUEST(q_quest, gyro_EKF_QUEST, dt, ...
        opt_sigma_w_EKF_QUEST, opt_sigma_obs_EKF_QUEST, quat_ini);

    roll_EKF_QUEST = zeros(1, length(quat_est));
    pitch_EKF_QUEST = zeros(1, length(quat_est));
    yaw_EKF_QUEST = zeros(1, length(quat_est));
    for i = 1 : length(quat_est)
         eulerAngles = wag.quat_to_euler(quat_est(i, :));
         roll_EKF_QUEST(i) = eulerAngles(1);
         pitch_EKF_QUEST(i) = -eulerAngles(2);
         yaw_EKF_QUEST(i) = eulerAngles(3);
    end
    
    % 2.6.3) Compute the RMSE and plot signals.
    if strcmp(selectedAngle, 'pitch')
        angle_EKF_QUEST = pitch_EKF_QUEST;
    elseif strcmp(selectedAngle, 'roll')
        angle_EKF_QUEST = roll_EKF_QUEST;
    elseif strcmp(selectedAngle, 'yaw')
        angle_EKF_QUEST = yaw_EKF_QUEST;
    end
    
    rmse_EKF_QUEST(j) = wag.compute_rmse(true_ang, 180 / pi * ...
        angle_EKF_QUEST', rmse_offset);
    if strcmpi(showPlots, 'yes')
        wag.create_figures('other', 'EKF-QUEST', time, 180 / pi * ...
            angle_EKF, '--red'); 
    end

    %% --------------------------------------------------------------------
    % 2.7) Detection of motion intensity (low/high).
    % ---------------------------------------------------------------------
    % Now we will apply a series of motion intensity detectors to determine
    % the degree of intensity of the motion. In this case we will make two
    % distinctions: smooth and intense motion. The detectors build a binary
    % marker which indicates whether the motion is smoot ('0') or intense
    % ('1'). 
     
    % 2.7.1) Initialize parameters of algorithm.
    
    % MBGT (window size and decision threshold).
    lwin_mbgt = 14;       threshold_mbgt = 1.4;
    
    % MBCUSUM (window size and decision threshold).
    lwin_mbcusum = 18;    threshold_mbcusum = 2.8e-6;   
    
    % FSD (window size, decision threshold, overlapping and normalization
    % factor). 
    lwin_fsd = 20;  threshold_fsd = 6;  shift_fsd = 19; lambda = 30;
    
    % LTSD (window size, decision threshold and overlapping).
    lwin_ltsd = 14;       threshold_ltsd = 4;   shift_ltsd = 2;
    
    % AMD (window size and decision threshold).
    lwin_amd = 86;        threshold_amd = 0.0011;
    
    % AMVD (window size and decision threshold).
    lwin_amvd = 16;       threshold_amvd = 0.0173;
    
    % ARED (window size and decision threshold).
    lwin_ared = 8;        threshold_ared = 1;
    
    % SHOD (window size and decision threshold).
    lwin_shod = 19;       threshold_shod = 1;

	% Define input signal.
    input_signal = sqrt(axC .^ 2 + ayC .^ 2 + azC .^ 2);

    % Computation of intensity markers.
    [marker_mbgt, T_mbgt, marker_mbcusum, T_mbcusum, marker_fsd, ...
        T_fsd, T_fsd_expanded, marker_ltsd, T_ltsd, T_ltsd_expanded, ...
        marker_amd, T_amd, marker_amvd, T_amvd, marker_ared, T_ared, ...
        marker_shod, T_shod] = wag.build_int_markers(input_signal, axC, ...
        ayC, azC, gxC, gyC, gzC, threshold_mbgt, threshold_mbcusum, ...
        threshold_fsd, threshold_ltsd, lwin_mbgt, lwin_mbcusum, ...
        lwin_fsd, lwin_ltsd, shift_fsd, shift_ltsd, lambda, lwin_amd, ...
        threshold_amd, lwin_amvd, threshold_amvd, lwin_ared, ...
        threshold_ared, lwin_shod, threshold_shod);
      
    % Plot input signal and markers in order to decide which marker to
    % use. 
    if strcmpi(showPlots, 'yes')
        detectors_figure = figure(2);
        subplot(3, 3, 1)
        plot(T_mbgt)
        hold on
        plot(threshold_mbgt * ones(1, length(T_mbgt)), 'r')
        legend('Detector output (MBGT)', 'Detection threshold')
        subplot(3, 3, 2)
        plot(T_mbcusum)
        hold on
        plot(threshold_mbcusum * ones(1, length(T_mbcusum)), 'r')
        legend('Detector output (MBCUSUM)', 'Detection threshold')
        subplot(3, 3, 3)
        plot(T_fsd_expanded)
        hold on
        plot(threshold_fsd * ones(1, length(T_fsd_expanded)), 'r')
        legend('Detector output (FSD)', 'Detection threshold')
        subplot(3, 3, 4)
        plot(T_ltsd_expanded)
        hold on
        plot(threshold_ltsd * ones(1, length(T_ltsd_expanded)), 'r')
        legend('Detector output (LTSD)', 'Detection threshold')
        subplot(3, 3, 5)
        plot(T_amd)
        hold on
        plot(threshold_amd * ones(1, length(T_amd)), 'r')
        legend('Detector output (AMD)', 'Detection threshold')
        subplot(3, 3, 6)
        plot(T_amvd)
        hold on
        plot(threshold_amvd * ones(1, length(T_amvd)), 'r')
        legend('Detector output (AMVD)', 'Detection threshold')
        subplot(3, 3, 7)
        plot(T_ared)
        hold on
        plot(threshold_ared * ones(1, length(T_ared)), 'r')
        legend('Detector output (ARED)', 'Detection threshold')
        subplot(3, 3, 8)
        plot(T_shod)
        hold on
        plot(threshold_shod * ones(1, length(T_shod)), 'r')
        legend('Detector output (SHOD)', 'Detection threshold')

        markers_figure = figure(3);
        subplot(3, 3, 1)
        plot(input_signal)
        hold on
        plot(marker_mbgt + 1,'r')
        legend('Input signal','MBGT decision')
        subplot(3, 3, 2)
        plot(input_signal)
        hold on
        plot(marker_mbcusum + 1,'r')
        legend('Input signal','MBCUSUM decision')
        subplot(3, 3, 3)
        plot(input_signal)
        hold on
        plot(marker_fsd + 1, 'r')
        legend('Input signal','FSD decision')
        subplot(3, 3, 4)
        plot(input_signal)
        hold on
        plot(marker_ltsd + 1, 'r')
        legend('Input signal','LTSD decision')
        subplot(3, 3, 5)
        plot(input_signal)
        hold on
        plot(marker_ared + 1, 'r')
        legend('Input signal','ARED decision')
        subplot(3, 3, 6)
        plot(input_signal)
        hold on
        plot(marker_amvd + 1,'r')
        legend('Input signal','AMVD decision')
        subplot(3, 3, 7)
        plot(input_signal)
        hold on    
        plot(marker_amd + 1, 'r')
        legend('Input signal','AMD decision')
        subplot(3, 3, 8)
        plot(input_signal)
        hold on       
        plot(marker_shod + 1,'r')
        legend('Input signal','SHOD decision')
    end
    
    % Select motion intensity marker.
    chosen_marker = marker_fsd;

    %% --------------------------------------------------------------------
    % 2.9) Gated Kalman Filter (DKF) I.
    % ---------------------------------------------------------------------
    % After building the motion intensity markers and selecting the most
    % appropriate one, we will start to test the gating strategy on the
    % same sensor fusion algorithms that were tested above. The goal is to
    % check whether the gating strategy allows for an increment in the
    % precision of the orientation estimates. 
    % The first algorithm that we will check is, again, the standard Kalman
    % Filter.
    
    % 2.9.1) Definition of the variables of the optimization process.  
    if strcmp(selectedAngle, 'pitch')
        gyro_GKF = gyC;
        obs_GKF = pitchAcc;
    elseif strcmp(selectedAngle, 'roll')
        gyro_GKF = gxC;
        obs_GKF = rollAcc;
    elseif strcmp(selectedAngle, 'yaw')
        gyro_GKF = gzC;
        obs_GKF = yawAccM;
    end
    
    % 2.9.2) Parameter optimization.
    
    % Definition of the initial value of the parameters. 
    p0_GKF = [100 100 0.1 0.1];
    
    % Call the optimization routine.
    [xmin, fmin, ct] = wag.optimizeGKF(obs_GKF, gyro_GKF, f, true_ang, ...
    p0_GKF, rmse_offset, chosen_marker);

    fprintf('--------------------GKF OPTIMIZATION----------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf(['Optimal parameters are: \n -Alpha1: %0.6f\n',...
        ' -Alpha2: %0.6f\n -Beta1: %0.6f\n -Beta2: %0.6f\n'],...
        xmin(1), xmin(2), xmin(3), xmin(4))
    fprintf('----------------------------------------------------------\n')
    
    % Get the optimal parameters. 
    opt_alpha1_GKF = xmin(1);
    opt_alpha2_GKF = xmin(2);
    opt_beta1_GKF = xmin(3);
    opt_beta2_GKF = xmin(4);

    % Save optimal parameters.
    opt_alphas1_GKF(j) = opt_alpha1_GKF;
    opt_alphas2_GKF(j) = opt_alpha2_GKF;
    opt_betas1_GKF(j) = opt_beta1_GKF;
    opt_betas2_GKF(j) = opt_beta2_GKF;
    
    % 2.9.3) Computation of the orientation estimation using optimal 
    %        parameters.
    angle_GKF = wag.fusionGKF(gyro_GKF, obs_GKF, f, opt_alpha1_GKF, ...
        opt_alpha2_GKF, opt_beta1_GKF, opt_beta2_GKF, chosen_marker);

    % Compute RMSE and plot signals.
    rmse_GKF(j) = wag.compute_rmse(true_ang, 180 / pi * angle_GKF, ...
        rmse_offset);
    
    figure(1)
    if strcmpi(showPlots, 'yes')
        wag.create_figures('other', 'GKF', time, 180 / pi * ...
            angle_GKF, '--cyan'); 
    end

    %% --------------------------------------------------------------------
    % 2.10) Gated Kalman Filter (QUEST).
    % ---------------------------------------------------------------------
    % 2.10.1) Definition of the variables of the optimization process. 
    if strcmp(selectedAngle, 'pitch')
        obs_GKF_QUEST = pitch_QUEST;
        gyro_GKF_QUEST = gyC;
    elseif strcmp(selectedAngle, 'roll')
        obs_GKF_QUEST = roll_QUEST;
        gyro_GKF_QUEST = gxC;
    elseif strcmp(selectedAngle, 'yaw')
        obs_GKF_QUEST = yaw_QUEST;
        gyro_GKF_QUEST = gzC;
    end

    % 2.10.2) Parameter optimization.
    
    % Definition of the initial value of the parameters. 
    p0_GKF_QUEST = [100 100 0.1 0.1];
    
    % Call the optimization routine.
    [xmin, fmin, ct] = wag.optimizeGKF(obs_GKF_QUEST, gyro_GKF_QUEST, ...
        f, true_ang, p0_GKF_QUEST, rmse_offset, chosen_marker);

    fprintf('------------------GKF-QUEST OPTIMIZATION------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf(['Optimal parameters are: \n -Alpha1: %0.9f\n',...
        ' -Alpha2: %0.9f\n -Beta1: %0.9f\n -Beta2: %0.9f\n'],...
        xmin(1), xmin(2), xmin(3), xmin(4))
    fprintf('----------------------------------------------------------\n')
   
    % Extract optimal parameters.
    opt_alpha1_GKF_QUEST = xmin(1);
    opt_alpha2_GKF_QUEST = xmin(2);
    opt_beta1_GKF_QUEST = xmin(3);
    opt_beta2_GKF_QUEST = xmin(4);
    
    % Store optimal parameters.
    opt_alphas1_GKF_QUEST(j) = opt_alpha1_GKF_QUEST;
    opt_alphas2_GKF_QUEST(j) = opt_alpha2_GKF_QUEST;
    opt_betas1_GKF_QUEST(j) = opt_beta1_GKF_QUEST;
    opt_betas2_GKF_QUEST(j) = opt_beta2_GKF_QUEST;

    % 2.10.3) Computation of the orientation estimation using optimal 
    %        parameters.
    angle_GKF_QUEST = wag.fusionGKF(gyro_GKF_QUEST, obs_GKF_QUEST, f, ...
        opt_alpha1_GKF_QUEST, opt_alpha2_GKF_QUEST, opt_beta1_GKF_QUEST,... 
        opt_beta2_GKF_QUEST, chosen_marker);

    % Compute RMSE and plot signals.
    rmse_GKF_QUEST(j) = wag.compute_rmse(true_ang, 180 / pi * ...
        angle_GKF_QUEST, rmse_offset);
    
    if strcmpi(showPlots, 'yes')
        figure(1)
        wag.create_figures('other', 'GKF-QUEST', time, 180 / pi * ...
            angle_GKF_QUEST, '--black'); 
    end

    %% --------------------------------------------------------------------
    % 2.11) Gated Extended Kalman Filter.
    % ---------------------------------------------------------------------
    % 2.11.1) Parameter optimization.
    
    % Set initial value of parameters.
    p0_GEKF = [0.02 3 0.01 5];
    
    % Call optimization routine.
    [xmin, fmin, ct] = wag.optimizeGEKF(axC, ayC, azC, gxC, gyC, gzC, ...
    f, true_ang, p0_GEKF, rmse_offset, selectedAngle, chosen_marker);

    fprintf('--------------------GEKF OPTIMIZATION---------------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf(['Optimal parameters are: \n -Sigma_w1: %0.9f\n',...
        ' -Sigma_w2: %0.9f\n -Sigma_obs1: %0.9f\n -Sigma_obs2:',...
        ' %0.9f\n'], xmin(1), xmin(2), xmin(3), xmin(4))
    fprintf('----------------------------------------------------------\n')
    
    % Extract optimal parameters
    opt_sigma_w1_GEKF = xmin(1);
    opt_sigma_obs1_GEKF = xmin(2);
    opt_sigma_w2_GEKF = xmin(3);
    opt_sigma_obs2_GEKF = xmin(4);
    
    % Store optimal parameters
    opt_sigmas_w1_GEKF(j) = opt_sigma_w1_GEKF;
    opt_sigmas_obs1_GEKF(j) = opt_sigma_obs1_GEKF;
    opt_sigmas_w2_GEKF(j) = opt_sigma_w2_GEKF;
    opt_sigmas_obs2_GEKF(j) = opt_sigma_obs2_GEKF;
    
    % 2.11.2) Computation of the orientation estimation using optimal 
    %         parameters.
    a_GEKF = [azC ayC axC];
    w_GEKF = [gxC gyC gzC];
    ini = [1 0 0 0];
    dt = 1 / f;
    q = wag.fusionGEKF(-a_GEKF, w_GEKF, dt, opt_sigma_w1_GEKF,...
    opt_sigma_w2_GEKF, opt_sigma_obs1_GEKF, opt_sigma_obs2_GEKF, ...
    ini, chosen_marker);

    % Transformation of the orientation quaternion into Euler angles.
    roll_GEKF = zeros(1, length(q));
    pitch_GEKF = zeros(1, length(q));
    yaw_GEKF = zeros(1, length(q));
    for i = 1 : length(q)
        eulerAngles = wag.quat_to_euler(q(i, :));
        roll_GEKF(i) = eulerAngles(1);
        pitch_GEKF(i) = eulerAngles(2);
        yaw_GEKF(i) = eulerAngles(3);
    end

    % 2.11.3) Compute the RMSE and plot signals.
    if strcmp(selectedAngle, 'pitch')
        angle_GEKF = pitch_GEKF;
    elseif strcmp(selectedAngle, 'roll')
        angle_GEKF = roll_GEKF;
    elseif strcmp(selectedAngle, 'yaw')
        angle_GEKF = yaw_GEKF;
    end

    rmse_GEKF(j) = wag.compute_rmse(true_ang, 180 / pi * angle_GKF, ...
        rmse_offset);
       
    if strcmpi(showPlots, 'yes')
        figure(1)
        wag.create_figures('other', 'GEKF', time, 180 / pi * ...
            angle_GEKF, '--green'); 
    end
    %% --------------------------------------------------------------------
    % 2.12) Gated Extended Kalman Filter (QUEST).
    % ---------------------------------------------------------------------

    % 2.12.1) Optimize parameters.
    gyro_GEKF_QUEST = [gxC gyC gzC];
    
    % Set the initial parameters.
    p0_GEKF_QUEST = [1 0.2 1 0.2];
    
    % Call the optimization routine.
    [xmin, fmin, ct] = wag.optimizeGEKF_QUEST(gyro_GEKF_QUEST, q_quest, ...
    p0_GEKF_QUEST, true_ang, selectedAngle, f, rmse_offset, chosen_marker);

    fprintf('------------------GEKF-QUEST OPTIMIZATION-----------------\n')
    fprintf('The optimization process finished in %d iterations.\n', ct)
    fprintf('The minimum RMSE found is: %0.4f\n', fmin);
    fprintf(['Optimal parameters are: \n -Sigma_w1: %0.9f\n',...
        ' -Sigma_w2: %0.9f\n -Sigma_obs1: %0.9f\n -Sigma_obs2:',...
        ' %0.9f\n'], xmin(1), xmin(2), xmin(3), xmin(4))
    fprintf('----------------------------------------------------------\n')
    
    % Extract optimal parameters.
    opt_sigma_w1_GEKF_QUEST = xmin(1);
    opt_sigma_obs1_GEKF_QUEST = xmin(2);
    opt_sigma_w2_GEKF_QUEST = xmin(3);
    opt_sigma_obs2_GEKF_QUEST = xmin(4);

    % Store optimal parameters.
    opt_sigmas_w1_GEKF_QUEST(j) = opt_sigma_w1_GEKF_QUEST;
    opt_sigmas_obs1_GEKF_QUEST(j) = opt_sigma_obs1_GEKF_QUEST;
    opt_sigmas_w2_GEKF_QUEST(j) = opt_sigma_w2_GEKF_QUEST;
    opt_sigmas_obs2_GEKF_QUEST(j) = opt_sigma_obs2_GEKF_QUEST;

    % 2.12.2) Computation of the orientation estimation using optimal 
    %         parameters.
    quat_ini = q_quest(1, :);
    quat_est = wag.fusionGEKF_QUEST(q_quest, gyro_GEKF_QUEST, dt, ...
        opt_sigma_w1_GEKF_QUEST, opt_sigma_w2_GEKF_QUEST, ...
        opt_sigma_obs1_GEKF_QUEST, opt_sigma_obs2_GEKF_QUEST, quat_ini, ...
        chosen_marker);
    
    roll_GEKF_QUEST = zeros(1, length(quat_est));
    pitch_GEKF_QUEST = zeros(1, length(quat_est));
    yaw_GEKF_QUEST = zeros(1, length(quat_est));
    for i = 1 : length(quat_est)
         eulerAngles = wag.quat_to_euler(quat_est(i, :));
         roll_GEKF_QUEST(i) = eulerAngles(1);
         pitch_GEKF_QUEST(i) =- eulerAngles(2);
         yaw_GEKF_QUEST(i) = eulerAngles(3);
    end

    % 2.12.3) Compute the RMSE and plot signals.
    if strcmp(selectedAngle, 'pitch')
        angle_GEKF_QUEST = pitch_GEKF_QUEST;
    elseif strcmp(selectedAngle, 'roll')
        angle_GEKF_QUEST = roll_GEKF_QUEST;
    elseif strcmp(selectedAngle, 'yaw')
        angle_GEKF_QUEST = yaw_GEKF_QUEST;
    end
    
    rmse_GEKF_QUEST(j) = wag.compute_rmse(true_ang, 180 / pi * ...
        angle_GEKF_QUEST', rmse_offset);
           
    if strcmpi(showPlots, 'yes')
        figure(1)
        wag.create_figures('other', 'GEKF-QUEST', time, 180 / pi * ...
            angle_GEKF_QUEST, '--yellow'); 
    end

    %% --------------------------------------------------------------------
    % 2.14) Final adjustments and computation of improvements.
    % ---------------------------------------------------------------------
    % We need to remove the QUEST signals so they are recomputed in the
    % next iteration.
    clear q_quest;
    clear roll_QUEST;
    clear pitch_QUEST;
    clear yaw_QUEST;

    % Computation of improvements.
    imp_GKF(j) = wag.improvement_perc(rmse_KF(j), rmse_GKF(j));
    imp_GKF_QUEST(j) = wag.improvement_perc(rmse_KF_QUEST(j), ...
        rmse_GKF_QUEST(j));
    imp_GEKF(j) = wag.improvement_perc(rmse_EKF(j), rmse_GEKF(j));
    imp_GEKF_QUEST(j) = wag.improvement_perc(rmse_EKF_QUEST(j), ...
        rmse_GEKF_QUEST(j));
    
    % Store figure and close it.
    if strcmpi(showPlots, 'yes')
        savefig(angles_figure, strcat('figures/angles_', ...
            file_name(1 : end - 4)));
        savefig(detectors_figure, strcat('figures/detectors_', ...
            file_name(1 : end - 4)));
        savefig(markers_figure, strcat('figures/markers_', ...
            file_name(1 : end - 4)));
        close all
    end
end

%% ------------------------------------------------------------------------
% 3) Compute mean RMSE and mean Parameters.
% -------------------------------------------------------------------------
% At the end of the Monte Carlo simulation we need to compute the mean and
% standard deviation of all the parameters and RMSEs.

% 3.1) QUEST.
mean_a1s_opt_QUEST = mean(a1s_opt_QUEST);
std_a1s_opt_QUEST = std(a1s_opt_QUEST);
mean_a2s_opt_QUEST = mean(a2s_opt_QUEST);
std_a2s_opt_QUEST = std(a2s_opt_QUEST);
mean_rmse_QUEST = mean(rmse_QUEST);
std_rmse_QUEST = std(rmse_QUEST);

% 3.2) Regular Kalman filter.
mean_opt_alphas_KF = mean(opt_alphas_KF);
std_opt_alphas_KF = std(opt_alphas_KF);
mean_opt_betas_KF = mean(opt_betas_KF);
std_opt_betas_KF = std(opt_betas_KF);
mean_rmse_KF = mean(rmse_KF);
std_rmse_KF = std(rmse_KF);

% 3.3) Regular Kalman filter + QUEST.
mean_opt_alphas_KF_QUEST = mean(opt_alphas_KF_QUEST);
std_opt_alphas_KF_QUEST  = std(opt_alphas_KF_QUEST);
mean_opt_betas_KF_QUEST  = mean(opt_betas_KF_QUEST);
std_opt_betas_KF_QUEST  = std(opt_betas_KF_QUEST);
mean_rmse_KF_QUEST  = mean(rmse_KF_QUEST );
std_rmse_KF_QUEST  = std(rmse_KF_QUEST );

% 3.4) Extended Kalman filter.
mean_opt_sigmas_w_EKF = mean(opt_sigmas_w_EKF);
std_opt_sigmas_w_EKF = std(opt_sigmas_w_EKF);
mean_opt_sigmas_obs_EKF = mean(opt_sigmas_obs_EKF);
std_opt_sigmas_obs_EKF = std(opt_sigmas_obs_EKF);
mean_rmse_EKF = mean(rmse_EKF);
std_rmse_EKF = std(rmse_EKF);

% 3.5) Extended Kalman filter + QUEST.
mean_opt_sigmas_w_EKF_QUEST = mean(opt_sigmas_w_EKF_QUEST);
std_opt_sigmas_w_EKF_QUEST = std(opt_sigmas_w_EKF_QUEST);
mean_opt_sigmas_obs_EKF_QUEST = mean(opt_sigmas_obs_EKF_QUEST);
std_opt_sigmas_obs_EKF_QUEST = std(opt_sigmas_obs_EKF_QUEST);
mean_rmse_EKF_QUEST = mean(rmse_EKF_QUEST);
std_rmse_EKF_QUEST = std(rmse_EKF_QUEST);

% 3.6) Gated Kalman Filter.
mean_opt_alphas1_GKF = mean(opt_alphas1_GKF);
std_opt_alphas1_GKF = std(opt_alphas1_GKF);
mean_opt_alphas2_GKF = mean(opt_alphas2_GKF);
std_opt_alphas2_GKF = std(opt_alphas2_GKF);
mean_opt_betas1_GKF = mean(opt_betas1_GKF);
std_opt_betas1_GKF = std(opt_betas1_GKF);
mean_opt_betas2_GKF = mean(opt_betas2_GKF);
std_opt_betas2_GKF = std(opt_betas2_GKF);
mean_rmse_GKF = mean(rmse_GKF);
std_rmse_GKF = std(rmse_GKF);

mean_imp_GKF = mean(imp_GKF);
std_imp_GKF = std(imp_GKF);

% 3.7) Gated Kalman Filter + QUEST.
mean_opt_alphas1_GKF_QUEST = mean(opt_alphas1_GKF_QUEST);
std_opt_alphas1_GKF_QUEST = std(opt_alphas1_GKF_QUEST);
mean_opt_alphas2_GKF_QUEST = mean(opt_alphas2_GKF_QUEST);
std_opt_alphas2_GKF_QUEST = std(opt_alphas2_GKF_QUEST);
mean_opt_betas1_GKF_QUEST = mean(opt_betas1_GKF_QUEST);
std_opt_betas1_GKF_QUEST = std(opt_betas1_GKF_QUEST);
mean_opt_betas2_GKF_QUEST = mean(opt_betas2_GKF_QUEST);
std_opt_betas2_GKF_QUEST = std(opt_betas2_GKF_QUEST);
mean_rmse_GKF_QUEST = mean(rmse_GKF_QUEST);
std_rmse_GKF_QUEST = std(rmse_GKF_QUEST);

mean_imp_GKF_QUEST = mean(imp_GKF_QUEST );
std_imp_GKF_QUEST = std(imp_GKF_QUEST );

% 3.8) Gated Extended Kalman Filter.
mean_opt_sigmas_w1_GEKF = mean(opt_sigmas_w1_GEKF);
std_opt_sigmas_w1_GEKF = std(opt_sigmas_w1_GEKF);
mean_opt_sigmas_obs1_GEKF = mean(opt_sigmas_obs1_GEKF);
std_opt_sigmas_obs1_GEKF = std(opt_sigmas_obs1_GEKF);
mean_opt_sigmas_w2_GEKF = mean(opt_sigmas_w2_GEKF);
std_opt_sigmas_w2_GEKF = std(opt_sigmas_w2_GEKF);
mean_opt_sigmas_obs2_GEKF = mean(opt_sigmas_obs2_GEKF);
std_opt_sigmas_obs2_GEKF = std(opt_sigmas_obs2_GEKF);
mean_rmse_GEKF = mean(rmse_GEKF);
std_rmse_GEKF = std(rmse_GEKF);

mean_imp_GEKF = mean(imp_GEKF);
std_imp_GEKF = std(imp_GEKF);

% 3.9) Gated Extended Kalman Filter + QUEST.
mean_opt_sigmas_w1_GEKF_QUEST = mean(opt_sigmas_w1_GEKF_QUEST);
std_opt_sigmas_w1_GEKF_QUEST = std(opt_sigmas_w1_GEKF_QUEST);
mean_opt_sigmas_obs1_GEKF_QUEST = mean(opt_sigmas_obs1_GEKF_QUEST);
std_opt_sigmas_obs1_GEKF_QUEST = std(opt_sigmas_obs1_GEKF_QUEST);
mean_opt_sigmas_w2_GEKF_QUEST = mean(opt_sigmas_w2_GEKF_QUEST);
std_opt_sigmas_w2_GEKF_QUEST = std(opt_sigmas_w2_GEKF_QUEST);
mean_opt_sigmas_obs2_GEKF_QUEST = mean(opt_sigmas_obs2_GEKF_QUEST);
std_opt_sigmas_obs2_GEKF_QUEST = std(opt_sigmas_obs2_GEKF_QUEST);
mean_rmse_GEKF_QUEST = mean(rmse_GEKF_QUEST);
std_rmse_GEKF_QUEST = std(rmse_GEKF_QUEST);

mean_imp_GEKF_QUEST = mean(imp_GEKF_QUEST);
std_imp_GEKF_QUEST = std(imp_GEKF_QUEST);

%% ------------------------------------------------------------------------
% 4) Store results in Excel tables.
% -------------------------------------------------------------------------
% Once the experiments are finished, we store the results in an excel file.
report_name = 'report.xls';

headers = {'','KF','KF (QUEST)','EKF','EKF (QUEST)'};
print_RMSE_KF = strcat(num2str(mean_rmse_KF), '±', num2str(std_rmse_KF));
print_RMSE_KF_QUEST = strcat(num2str(mean_rmse_KF_QUEST), '±', ...
    num2str(std_rmse_KF_QUEST));
print_RMSE_EKF = strcat(num2str(mean_rmse_EKF), '±', ...
    num2str(std_rmse_EKF));
print_RMSE_EKF_QUEST = strcat(num2str(mean_rmse_EKF_QUEST), '±', ...
    num2str(std_rmse_EKF_QUEST));

print_rmses = {'RMSE', print_RMSE_KF, print_RMSE_KF_QUEST, ...
    print_RMSE_EKF, print_RMSE_EKF_QUEST};

xlswrite(strcat('reports/', report_name), headers, 1, 'A1')
xlswrite(strcat('reports/', report_name), print_rmses, 1, 'A2')

headers = {'', 'GKF', 'GKF (QUEST)', 'GEKF', 'GEKF (QUEST)'};
print_RMSE_GKF = strcat(num2str(mean_rmse_GKF), '±', ...
    num2str(std_rmse_GKF));
print_RMSE_GKF_QUEST = strcat(num2str(mean_rmse_GKF_QUEST), '±', ...
    num2str(std_rmse_GKF_QUEST));
print_RMSE_GEKF = strcat(num2str(mean_rmse_GEKF), '±', ...
    num2str(std_rmse_GEKF));
print_RMSE_GEKF_QUEST = strcat(num2str(mean_rmse_GEKF_QUEST), '±', ...
    num2str(std_rmse_GEKF_QUEST));

print_rmses = {'RMSE', print_RMSE_GKF, print_RMSE_GKF_QUEST, ...
    print_RMSE_GEKF, print_RMSE_GEKF_QUEST};
xlswrite(strcat('reports/', report_name), headers, 1, 'A3')
xlswrite(strcat('reports/', report_name), print_rmses, 1, 'A4')

print_imp_GKF = strcat(num2str(mean_imp_GKF), '±', num2str(std_imp_GKF));
print_imp_GKF_QUEST = strcat(num2str(mean_imp_GKF_QUEST), '±', ...
    num2str(std_imp_GKF_QUEST));
print_imp_GEKF = strcat(num2str(mean_imp_GEKF), '±', ...
    num2str(std_imp_GEKF));
print_imp_GEKF_QUEST = strcat(num2str(mean_imp_GEKF_QUEST), '±', ...
    num2str(std_imp_GEKF_QUEST));
print_imp = {'Improvement (%)', print_imp_GKF, print_imp_GKF_QUEST, ...
    print_imp_GEKF, print_imp_GEKF_QUEST};
xlswrite(strcat('reports/', report_name), print_imp, 1, 'A5')

fprintf('Simulation ended successfully :)\n');

% -------------------------------------------------------------------------
% |||||||||||||||||||||||||| END OF MAINLOOP FILE |||||||||||||||||||||||||
% -------------------------------------------------------------------------