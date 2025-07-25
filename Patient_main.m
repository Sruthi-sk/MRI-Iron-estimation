clc
clear all
close all

%% Load the data
% paths = {'Patient Datasets -20240311/1_2/','Patient Datasets -20240311/3_8/', 'Patient Datasets -20240311/2_7/', 'Patient Datasets -20240311/4_15/'};
paths = {'Patient Dataset/1_2/','Patient Dataset/2_7/', 'Patient Dataset/3_8/', 'Patient Dataset/4_15/'};

% Initialize arrays to store results for each patient
numPatients = numel(paths);
for j = 1:numPatients
    fileList = dir(fullfile(paths{j}, '*.IMA'));
    files = numel(fileList);
end

%% Set the ROI
TE_values = zeros(numPatients, files);
meanIntensityValues = zeros(numPatients, files);
Patient = cell(numPatients, files);
ROI_area = 200;

for j = 1:numPatients
    % Initialize patient_combined for each patient
    fileList = dir(fullfile(paths{j}, '*.IMA'));
    filePath = fullfile(paths{j}, fileList(1).name); %get the path of the first image
    
    image = double(dicomread(filePath)); %save the first image
    disp(size(image))
    [files, ~] = size(fileList); % Determine the number of files
    patient_combined = zeros(files, size(image,1), size(image,2)); %initialize variable with the dimensions of the first image
    
    for i = 1:files  % Loop through each file
        % Get the current image for the current patient
        filePath = fullfile(paths{j}, fileList(i).name);

        % Read DICOM
        info = dicominfo(filePath);
        TE_values(j,i) = info.EchoTime;
        
        image = double(dicomread(filePath));
        Patient{j, i} = image;

        % Assign the current image to the appropriate slice of patient_combined
        patient_combined(i, :, :) = image;
    end
    
    % Visualize combined image
    figure(1);
    imagesc(squeeze(patient_combined(i, :, :)))
    axis image
    axis on
    pt = ginput(1); % Get a point from user for ROI

    % Calculate circular ROI
    radius = sqrt(ROI_area);
    points = bbox2points([pt(1) pt(2) radius radius]);
    xunit = points(:, 1);
    yunit = points(:, 2);

   % Create circular mask
    center = pt; % Assuming pt is the center of the circle for your square
    [rows, cols] = size(patient_combined(1, :, :)); % Assuming patient_combined is your image
    [x, y] = meshgrid(1:cols, 1:rows);
    mask = (x - center(1)).^2 + (y - center(2)).^2 <= (radius/2)^2;
    
    % Plot circle
    theta = linspace(0, 2*pi, 100);
    x_circle = center(1) + (radius/2) * cos(theta);
    y_circle = center(2) + (radius/2) * sin(theta);
    hold on;
    plot(x_circle, y_circle, 'k', 'LineWidth', 2);
    hold off;
    
    % Use mask to select ROI
    Mask = roipoly(squeeze(patient_combined(1, :, :)), xunit, yunit);
    idx1 = find(Mask == 1);
    Count_Mask_Ones = sum(sum(Mask));

    % Calculate mean intensity within the ROI for each image
    for i = 1:files 
        image = Patient{j, i};
        meanIntensityValues(j,i) = mean(image(Mask));
    end
end
fprintf('TE values')
disp(TE_values)

%% PLOT SCANS
figure();
for i = 1:files
    subplot(1, 8, i); % Create subplot at appropriate position
    imagesc(Patient{1, i}); % Show image with original aspect ratio
    axis image
    axis off
    colormap hot
end

%% Task 8
figure;
plot(TE_values(1,:), meanIntensityValues(1,:), 'o-', 'LineWidth', 2);
xlabel('Time to Echo (TE)');
ylabel('Mean Signal Intensity in ROI');
title('TE vs. Mean Signal Intensity in ROI');
grid on;

% hold on
% plot(TE_values(2,:), meanIntensityValues(2,:), 'o--', 'LineWidth', 2);
hold on
plot(TE_values(3,:), meanIntensityValues(3,:), 'o--', 'LineWidth', 2);
% hold on
% plot(TE_values(4,:), meanIntensityValues(4,:), 'o--', 'LineWidth', 2);
legend('Patient1_2', 'Patient3_8'); %, 'Patient2_7','Patient4_15');


%% Task 9 - Fit the models
lambdaA = 1e-5;
lambdaR = 1e-5;
radius = round(radius); %14

% Initialize arrays to store results for each patient
a = cell(numPatients, 1);
r = cell(numPatients, 1);
a_lm = cell(numPatients, 1);
r_lm = cell(numPatients, 1);
YY = cell(numPatients, 1);

a_patient = cell(numPatients, 1);
r_patient = cell(numPatients, 1);
resnorm_patient = cell(numPatients,1);
x_patient = [100,0];

tic
for j = 1:numPatients
    for ii = 1:length(TE_values(j,:)')
        patient_combined(ii, :, :) = Patient{j, ii};

        ROI = Mask .* squeeze(patient_combined(ii, :, :));
        ROI_ = ROI(:);
        Y_(:, :, ii) = reshape(ROI_(idx1), radius, radius);
        YY{j} = Y_;
    end

    Y = Y_; 

    [Nrow, Ncol, bands] = size(Y);
    yReshaped = reshape(Y, Nrow*Ncol, bands)';

    [a{j}, r{j}] = relaxationEst(yReshaped, TE_values(j,:)', Nrow, Ncol, lambdaA, lambdaR); 
    params_initial = [200, 1/mean(TE_values(j,:)')];
    [a_lm{j}, r_lm{j}, ~] = fitmodel_lm(yReshaped, TE_values(j,:)', params_initial);
 
end
toc


%% Plot images (SEXP MODEL)
mean_T2Star = cell(numPatients,1);
mean_a = cell(numPatients,1);
for j = 1:numPatients
    figure('Name', ['Patient ' num2str(j)])    
    subplot(131)
    imagesc(reshape(YY{j}(:, :, 1), Nrow, Ncol))
    axis image
    axis off
    colormap hot
    colorbar
    c = colorbar;
    % w = c.FontSize;
    % c.FontSize = 22;
    % set(gca, 'ColorScale', 'log')
    title('Observations')
    
    subplot(132)
    imagesc(reshape(a{j}, Nrow, Ncol))
    axis image
    axis off
    colormap hot
    colorbar
    c = colorbar;
    % w = c.FontSize;
    % c.FontSize = 22;
    % set(gca, 'ColorScale', 'log')
    title('Intensities Map')


    % Plot T2*
    subplot(133)
    imagesc(reshape(1./r{j}, Nrow, Ncol))
    axis image
    axis off
    colormap hot
    colorbar
    c = colorbar;
    % w = c.FontSize;
    % c.FontSize = 22;
    % set(gca, 'ColorScale', 'log')
    title('T2^* Map')
    
    mean_T2Star{j} = [mean(1./r{j})];
    mean_a{j} = mean(a{j});
    fprintf('ADMM')
    fprintf('Mean T2* value for Patient %d: %f\n', j, mean_T2Star{j});
    fprintf('Mean a value for Patient %d: %f\n', j, mean_a{j});
end


%% Plot images (SEXP MODEL) - Levenberg Marquardt
mean_T2StarLM = cell(numPatients,1);
mean_aLM = cell(numPatients,1);
for j = 1:numPatients
    figure('Name', ['Patient ' num2str(j)])    
    subplot(131)
    imagesc(reshape(YY{j}(:, :, 1), Nrow, Ncol))
    axis image
    axis off
    colormap hot
    colorbar
    c = colorbar;
    % w = c.FontSize;
    % c.FontSize = 22;
    % set(gca, 'ColorScale', 'log')
    title('Observations')
    
    subplot(132)
    imagesc(reshape(a_lm{j}, Nrow, Ncol))
    axis image
    axis off
    colormap hot
    colorbar
    c = colorbar;
    % w = c.FontSize;
    % c.FontSize = 22;
    % set(gca, 'ColorScale', 'log')
    title('Intensities Map LM')


    % Plot T2*
    subplot(133)
    imagesc(reshape(1./r_lm{j}, Nrow, Ncol))
    axis image
    axis off
    colormap hot
    colorbar
    c = colorbar;
    % w = c.FontSize;
    % c.FontSize = 22;
    % set(gca, 'ColorScale', 'log')
    title('T2^* Map LM')
    
    mean_T2StarLM{j} = [mean(1./r_lm{j})];
    mean_aLM{j} = mean(a_lm{j});
    fprintf('LM')
    
    fprintf('Mean T2* value for Patient %d: %f\n', j, mean_T2StarLM{j});
    fprintf('Mean s value for Patient %d: %f\n', j, mean_aLM{j});

    sgtitle(['Patient ' num2str(j) ' Analysis with Levenberg-Marquadt algorithm'])
end