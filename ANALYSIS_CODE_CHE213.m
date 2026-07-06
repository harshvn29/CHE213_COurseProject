%% 1. SYSTEM CLEANUP
clear; clc; close all;

% CONFIGURATION
videoFile = 'Hot_water.mp4'; 
v = VideoReader(videoFile);
realDyeRadiusCm = 5; % Update per video spread

% Smart Timestamps (25%, 50%, 75% of video)
maxTime = v.Duration * 0.95;
times = linspace(maxTime/4, maxTime, 3); 

%% 2. GLOBAL CALIBRATION (True C0 Reference)
% Capturing the dark intensity at the earliest point to set as C0
v.CurrentTime = times(1);
firstFrame = double(rgb2gray(readFrame(v)));
globalDarkest = min(firstFrame(:)); 
globalLightest = max(firstFrame(:));

%% 3. STORAGE & FIGURES
all_D = zeros(length(times), 1);
fig1 = figure('Name', 'Concentration Profile', 'Position', [50 450 1400 450]);
fig2 = figure('Name', 'Linearized Proof (lnC vs r^2)', 'Position', [50 50 1400 450]);

%% 4. ANALYSIS LOOP
for i = 1:length(times)
    t_target = times(i);
    v.CurrentTime = t_target;
    frame = readFrame(v);
    grayFrame = double(rgb2gray(frame)); 
    
    % --- PER-FRAME CALIBRATION (Accounting for Drift) ---
    f_cal = figure; imshow(uint8(grayFrame));
    title(sprintf('t = %.1fs | 1. Click CURRENT Center | 2. Click CLEAR Edge', t_target));
    [cols, rows] = ginput(2);
    close(f_cal);
    
    % Extraction
    xLine = linspace(cols(1), cols(2), 200);
    yLine = linspace(rows(1), rows(2), 200);
    I = interp2(grayFrame, xLine, yLine);
    I = I(:);
    
    % GLOBAL NORMALIZATION (Peak is allowed to drop) 
    C_exp = (globalLightest - I) ./ (globalLightest - globalDarkest);
    C_exp(C_exp < 1e-4) = 1e-4; C_exp(C_exp > 1.2) = 1.2; 
    
    r_exp = (linspace(0, realDyeRadiusCm, 200))'; 
    rSq = r_exp.^2;
    lnC = log(C_exp);
    
    %  REGRESSION 
    valid = (C_exp > 0.1 & C_exp < 0.9); % Focus on the gradient zone
    
    if sum(valid) > 10
        p = polyfit(rSq(valid), lnC(valid), 1);
        all_D(i) = abs(-1 / (4 * p(1) * t_target));
    else
        all_D(i) = NaN;
    end
    
    % PLOTTING FIG 1: C vs r
    figure(fig1); subplot(1, 3, i);
    plot(r_exp, C_exp, 'k.', 'MarkerSize', 8); hold on;
    if ~isnan(all_D(i))
        peak_at_r0 = C_exp(1); % Use actual center concentration for this t
        plot(r_exp, peak_at_r0 * exp(-(r_exp.^2)/(4*all_D(i)*t_target)), 'r-', 'LineWidth', 2);
    end
    title(sprintf('t = %.1fs', t_target));
    xlabel('r (cm)'); ylabel('C / C_{initial}'); grid on;
    axis([0 realDyeRadiusCm 0 1.2]);
    
    % PLOT FIG 2: ln(C) vs r^2 
    figure(fig2); subplot(1, 3, i);
    if ~isnan(all_D(i))
        plot(rSq(valid), lnC(valid), 'k.', 'MarkerSize', 10); hold on;
        plot(rSq(valid), polyval(p, rSq(valid)), 'b-', 'LineWidth', 2);
        title(sprintf('Slope Fit: D = %.2e', all_D(i)));
        xlabel('r^2 (cm^2)'); ylabel('ln(C/C_{initial})'); grid on;
    end
end

%% 5. FINAL RESULTS SUMMARY
avg_D = mean(all_D, 'omitnan');
std_D = std(all_D, 'omitnan');

fprintf('\nRESULTS SUMMARY: %s\n', videoFile);
fprintf('Time (s)  |  D (cm^2/s)\n');

for i = 1:length(times)
    fprintf('%8.1f  |  %.4e\n', times(i), all_D(i));
end

fprintf('AVERAGE D : %.4e cm^2/s\n', avg_D);
fprintf('ST. DEV   : %.4e\n', std_D);
