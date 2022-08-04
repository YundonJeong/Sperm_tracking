% Sperm tracking analysis
close all
clear all
clc
%% Import sperm tracking data
read_csv = readtable('spot export.csv'); % Read csv file.

% for matlab version upto R2019
data_micron = read_csv(4:end,[3 5 6 8]); % Extract data to use from data table. Here are [TRACK ID, POSITION X, POSITION Y, FRAME] used.
data_micron = str2double(data_micron{:,:}); % Transfer datatype from table to numeric double martrix.

% % for matlab version R2020 and later
% data_micron = read_csv(1:end,[3 5 6 8]); % Extract data to use from data table. Here are [TRACK ID, POSITION X, POSITION Y, FRAME] used.
% data_micron = double(data_micron{:,:}); % Transfer datatype from table to numeric double martrix.

%% Wall coordinates
wall_definition;

% Fitting wall into linear function
fit_wall = fit(wall_coord(:,1),wall_coord(:,2),'poly1'); % Fit wall coordinates into linear function.
x_fit = [1:0.1:img_size(1)]'; % Create x domain to represent the fitting function 'fit_wall'.
y_fit = fit_wall(x_fit); % Y values given by 'fit_wall'.
index_fit = (y_fit>=1)&(y_fit<=512);
wall_fit_coord_px = [x_fit(index_fit), y_fit(index_fit)]; % Eliminate coordinates that exceed the area of the image.

% Plot the original wall and the fitting line in scatter graph
figure(101), scatter(wall_coord(:,1),wall_coord(:,2),2,'o');
hold on
scatter(wall_fit_coord_px(:,1),wall_fit_coord_px(:,2),1,'o');
axis image, axis ij, grid on;
xlim([0 img_size(1)]), ylim([0 img_size(2)]);
title('Wall fitting');
legend('Data','Fitting');

%% Acquisition settings and field-of-view (FOV) scale
% Acquisition settings
scan_size = 776.72; % The size in microns of the scanning range, usually the square, at final detection.
pixel_number = 512; % The number of pixels that consists one dimension of the acquired image.
mag = 5; % System magnification in [x.] (times).
fr2sec = 1/30; % Duration time between frames of image in the unit of seconds.
fprintf('\nAcquisition parameter settings :\n      scan size = %.2f [um] \n      scanning pixel number = %u [px.] \n      magnification = %.2f [x.] \n      frame rate = %.2f [Hz] \n\n',scan_size, pixel_number, mag, 1/fr2sec);

% FOV
scale = scan_size / mag / pixel_number; % The size of the pixel in the FOV in microns.
fov_size_micron = scale * img_size(1); % The size (one dimension) of the FOV.

%% Unit transfer : pixel to micron
wall_coord_micron = wall_coord * scale; % Wall coordinates defined in the unit of pixel are transferred into microns by multiplying scaling factor 'scale'.
wall_fit_coord_micron = wall_fit_coord_px * scale; % Coordinates of the fitting line defined in the unit of pixel are transferred into microns by multiplying scaling factor 'scale'.

%% Plot entire track sets
data = data_micron;
id_list = [unique(data(:,1))]; % List of 'TRACK ID' of the sperms.

% Plot tracking data
for count = 1 : numel(id_list)
    id = id_list(count);
    data_id = data((id==data(:,1)),:);
    figure(102)
    hold on
    scatter(wall_coord_micron(:,1), wall_coord_micron(:,2), 4, [0.2 0.2 0.2],'filled')
    scatter(wall_fit_coord_micron(:,1),wall_fit_coord_micron(:,2), 4, [0.7 0.7 0.7],'filled')
    scatter(data_id(:,2),data_id(:,3), 4,'filled')
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis image, grid on;
    xlim([0 fov_size_micron]), ylim([0 fov_size_micron]);
    xlabel('x [um]'), ylabel('y [um]');
    title(strcat('ID = ',num2str(id),' (',num2str(count),'th)'))
    pause(0.02)
end
hold off

%% Calculate values
sperm_num = numel(id_list); % The number of the sperms.
calculated_values = zeros([1 14]); % Pre-define a zero-matrix to contain calculated values in the following loop.
% Calculate physical values based on tracking data.
for sperm = 1 : sperm_num % Select an individual sperm. It is repeated for the entire sperms.
id = id_list(sperm); % Bring 'TRACK ID' of the sperm.
data_id = data((id==data(:,1)),:); % Extract data of the sperm with corresponding id.
data_id = sortrows(data_id,4); % Sort data by order of frames.
frame_id = data_id(:,4); % List of the numbered frames that contains tracking data.
frame_num = size(frame_id,1); % The amount of the numbered frames
wall_num = size(wall_coord_micron,1); % The amount of the coordinates that consist the wall defined from the image.
wall_fit_num = size(wall_fit_coord_micron,1); % The amount of the coordinates that consist the fitting line of the wall.

% Pre-define zero-matrices
distance_id_fit = zeros([1 frame_num]);
distance_id_wall = zeros([1 frame_num]);
speed_id = zeros([1 frame_num]);
interlink_vec = zeros([2 frame_num]);
calculated_values_id = zeros([frame_num 14]); % There are 14 outcomes to be induced by each sperm data.

for frame = 1 : frame_num % Here the variable 'frame' works as a loop control variable. Involved contents in this loop will be operated by each frame.
    
    % Minimal distance between the fitting line and the tracked point
    distance_fit = zeros(size([1 wall_fit_num]));
    for coord = 1 : wall_fit_num
        distance_fit(coord) = norm(wall_fit_coord_micron(coord,[1 2])-data_id(frame,[2 3])); % Get distance to every fitting line coordinates from the tracked point.
    end
    distance_id_fit(frame) = min(distance_fit); % Find the minimal distance among 'distance_fit'.
    
    % Minimal distance from between the wall and the tracked point
    % Similarly done with the 'distance_id_fit' case.
    distance_wall = zeros(size([1 wall_num]));
    for coord = 1 : wall_num
        distance_wall(coord) = norm(wall_coord_micron(coord,[1 2])-data_id(frame,[2 3]));
    end
    distance_id_wall(frame) = min(distance_wall);
    
    % Differential speed
    % The differential speed is defined by dividing the length of the position change into the time duration between current frame and the next available frame.
    % The unit of 'speed_id' is [micron/frame number] here. This will be transferred to [micron/second] later using 'fr2sec'.
    if frame == frame_num
        speed_id(frame) = NaN; % Since there is no subsequent frame, the speed at the last frame is not defined.
    else
        speed_id(frame) = norm(data_id(frame+1,[2 3]) - data_id(frame,[2 3]))/(frame_id(frame+1)-frame_id(frame));
    end
    
    % Total distance
    % The total distance represents the total length that the sperm travelled during the acquisition.
    % It is defined by adding up the length of the position change over the entire frames.  
    distance_total = 0; % Initially set the variable to null.
    for c = 1 : frame_num-1
        distance_total = distance_total + norm(data_id(c+1,[2 3])-data_id(c,[2 3])); % Stack up the distance of all frames.
    end
    
    % Total displacement
    % The total displancement is the difference of the position between the initial frame and the final frame.
    % It is defined by subtracting the position of the initial frame from the position of the final frame.
    displacement = data_id(frame_num,[2 3]) - data_id(1,[2 3]);

    % Sperm approaching angle to the wall
    % The approaching angle is defined by the angle between the fitting line of the wall and the displancement vector.
    % The angles are only considered the absolute amount and represented in the unit of [radian].
    wall_tan = fit_wall.p1; % Gradient of the fitting line.
    wall_tan_vec = [1 wall_tan]; % Direction vector of fitting line.
    if norm(displacement) == 0
    angle_rad = NaN; % Return NaN when the size of the displacement is zero.
    else
    angle_cos = displacement * wall_tan_vec' / norm(displacement) / norm(wall_tan_vec); % Cosine of the angle is induced from inner product of two vectors.
    angle_rad = abs(acos(angle_cos)); % The angle is obtained by taking the arc-cosine of 'angle_cos'.
    end
    
    % Midpoint of displacement
    % Midpoint of the displancement is defined by adding the position of the initial frame and that of the final frame then dividing it into half. 
    position_mid = 0.5*(data_id(frame_num,[2 3]) + data_id(1,[2 3]));

    % Minimal distance to the midpoint from the fitting line
    % Similarly done with the 'distance_id_fit' case.
    distance_mid_fit = zeros(size([1 wall_fit_num])); 
    for coord = 1 : wall_fit_num
        distance_mid_fit(coord) = norm(wall_fit_coord_micron(coord,[1 2]) - position_mid);
    end
    distance_id_mid_fit = min(distance_mid_fit);
    
    % Minimal distance to the midpoint from the wall
    % Also similarly done with the 'distance_id_fit' case.
    distance_mid_wall = zeros(size([1 wall_num])); 
    for coord = 1 : wall_num
        distance_mid_wall(coord) = norm(wall_coord_micron(coord,[1 2]) - position_mid);
    end
    distance_id_mid_wall = min(distance_mid_wall);

    % Interlink vector
    % The interlink vector is the displacement vector between two subsequent available frames.
    if frame == frame_num
        interlink_vec(:,frame) = NaN; % Since there is no subsequent frame, the interlink vector at the last frame is not defined.
    else
        interlink_vec(:,frame) = data_id(frame+1,[2 3]) - data_id(frame,[2 3]);
    end

    % Straight line to width ratio
    % The straight line to width ratio (SWR) is the ratio of the length of total displacement to the width that the tracked path is spread transversely.
    % The width is defined by the steps below :
        % 1. Define the line that connects the initial and the final point. The line is described as y=ax+b here. Constant a and b are specifically determined.
        % 2. Calculate the distance to each tracked point from the line.
        % 3. Group the tracked points which are in the left side together from the line, and the right side together.
        % 4. Find the maximum from each group and sum the their absolute values.
    x = data_id(:,2); % X coordinate of the tracked points.
    y = data_id(:,3); % Y coordinate of the tracked points.
    a = (y(frame_num)-y(1)) / (x(frame_num)-x(1)); % Gradient of of the connecting line. 
    b = -a*x(1)+y(1); % Intercept on y axis of the line.
    d = zeros([1 frame_num]); % Distance to each tracked point from the line.
    sgn = zeros([1 frame_num]);
    for tp = 1 : frame_num
        A = data_id(tp,2);
        B = data_id(tp,3);
        d(tp) = sqrt(a^2/(1+a^2)*(A+(b-B)/a)^2); % Distance from the point (xp,yp) to the line kx+ly+m=0 can be expressed as d=abs(k*xp+l*yp+m)/sqrt(k^2+l^2).
        sgn(tp) = -sign((A-x(1))*a-(B-y(1))); % The sign is evaluated as +1 if the point is on the left side from the line, and -1 if it is on the right side. 
    end
    d_sgn = d.*sgn;
    width = max(d_sgn(d_sgn>=0)) + max(abs(d_sgn(d_sgn<=0)));

% Put aside calculate values in pre-defined matrix (except for the last element).
calculated_values_id(frame,:) = [id frame frame_id(frame) distance_id_fit(frame) distance_id_wall(frame) speed_id(frame)/fr2sec frame_num angle_rad*180/pi distance_total norm(displacement) distance_id_mid_fit distance_id_mid_wall width 0];
end

% Interlink angle
% The interlink angle is the angle between two subsequent interlink vectors.
% There will be (frame_num-2) outcomes that are available out of (frame_num) tracking data.
interlink_angle = zeros([1 frame_num]);
for frame = 1 : frame_num
    if frame == frame_num
        interlink_angle(frame) = NaN; % Since there is no subsequent frame, the interlink angle at the last frame is not defined.
    else
        if norm(interlink_vec(:,frame+1))*norm(interlink_vec(:,frame)) == 0
            interlink_angle(frame) = 0; % Return 0 when the size of the interlink vector is zero.
        else
            interlink_cos = interlink_vec(:,frame+1)'*interlink_vec(:,frame)/norm(interlink_vec(:,frame+1))/norm(interlink_vec(:,frame)); % Similarly done with the 'angle_cos' case.
            interlink_angle(frame) = abs(acos(interlink_cos)); % Take only the absolute value. This angle are ranged from 0 to pi.
        end
    end
end

% Mean of interlink angle
% Sum of the interlink angles are divided into (frame_num). 
interlink_angle_mean = sum(interlink_angle(1:end-2))/(frame_num-2);
% Fill the last element of 'calculated_value_id' with 'interlink_angle_mean', transferred to the unit of degree from radian.
calculated_values_id(:,end) = interlink_angle_mean*180/pi;

% Do these processes for the entire sperms and save it in the matrix 'calculated_values'.
calculated_values = [calculated_values; calculated_values_id];
end
calculated_values_table = calculated_values(2:end,:);


%% Write data matrix into csv file
calculated_values_table = array2table(calculated_values_table); % Transform variable datatype from double matrix to table
calculated_values_table.Properties.VariableNames(1:14) = {'Track ID','Order','Frame number','Distance from the wall [um]','Distance from fitting line [um]','Speed [um/s]','Total number of frames','Approaching angle to wall [deg.]','Total travelling distance [um]','Total displacement length [um]','Distance to the midpoint from fitting line [um]','Distance to the midpoint from wall [um]','Width [um]','Interlink angle [deg.]'}; % Designate header.
filename_1 = strcat('tracks_mag',num2str(mag),'x.csv');
writetable(calculated_values_table, filename_1);
disp(strcat('Calculated values are exported as "',filename_1,'"'))


%% Create pivot table
id_num = numel(id_list); % The number of sperm tracking data
pivot = zeros([id_num 13]);

pivot(:,1) = id_list; % List of 'TRACK ID'.
for id = 1 : id_num % Select individual sperm.
    extract_id_matrix = calculated_values(:,1)==id_list(id); % Logical matrix to extract data of the selected sperm.
    
    frame_number = calculated_values(extract_id_matrix,3); % Import frame numbers
    pivot(id,2) = max(frame_number) - min(frame_number); % The number of how many frames are consisting current data.

    pivot(id,3) = pivot(id,2) * fr2sec; % Total time duration in seconds of the current tracking data.

    temp = calculated_values(extract_id_matrix,8); % Import approaching angle to wall.
    pivot(id,4) = temp(1); % All elements of 'temp' have an identical value.  

    temp = calculated_values(extract_id_matrix,9); % Import total distance.
    pivot(id,5) = temp(1);

    pivot(id,6) = pivot(id,5)./pivot(id,3); % Calculate total travel speed (shortened as TTS).

    temp = calculated_values(extract_id_matrix,10); % Import total displacemnet length.
    pivot(id,7) = temp(1);
    
    pivot(id,8) = pivot(id,7)./pivot(id,3); % Calculate straight line speed (shortened as SLS).
    
    pivot(id,9) = pivot(id,8)./pivot(id,6); % Calculate LIN = SLS / TTS
    
    temp = calculated_values(extract_id_matrix,11); % Distance to the midpoint from fitting line.
    pivot(id,10) = temp(1);

    temp = calculated_values(extract_id_matrix,12); % Distance to the midpoint from wall.
    pivot(id,11) = temp(1);
    
    temp = calculated_values(extract_id_matrix,14); % Mean of interlink angle.
    pivot(id,12) = temp(1);
    
    temp = calculated_values(extract_id_matrix,13); % Import width.
    pivot(id,13) = pivot(id,7)/temp(1); % Calculate straight line to width ratio by dividing length of total displacement into width.
end

table_pivot = array2table(pivot);
table_pivot.Properties.VariableNames(1:13) = {'track_id','f_final-f_initial','track_duration','angle_from_wall','total_distance','tts','displacement','sls','lin','distance_to_midpoint_from_fitting_line','distance_to_midpoint_from_wall','mean_of_interlink_angle','straight_line_width_ratio'};
filename_2 = strcat('summary_tracks_mag',num2str(mag),'x.csv');
writetable(table_pivot, filename_2);
disp(strcat('Pivot table are exported as "',filename_2,'"'))




