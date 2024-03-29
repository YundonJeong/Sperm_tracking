% sperm tracking analysis

%% import coordinate data
read_csv = readtable('spot export.csv'); % read csv file

% % for matlab version upto R2019
% data_micron = read_csv(4:end,[3 5 6 8]); % extract data to use
% num = size(data_micron,1);
% data_micron = str2double(data_micron{:,:}); % header = [ID, position_x, position_y, frame];, transfer cell into double

% % for matlab version R2020 and later
data_micron = read_csv(1:end,[3 5 6 8]); % extract data to use
data_micron = double(data_micron{:,:}); % header = [ID, position_x, position_y, frame];, transfer cell into double
num = size(data_micron,1); % header = [ID, position_x, position_y, frame];

%% wall coordinates
wall_definition;

% wall linear regression
fit_wall = fit(wall_coord(:,1),wall_coord(:,2),'poly1');
x_fit = [1:0.1:img_size(1)]';
y_fit = fit_wall(x_fit);
index_fit = (y_fit>=1)&(y_fit<=512);
wall_fit_coord_px = [x_fit(index_fit), y_fit(index_fit)];
% scatter plot
figure(99), scatter(wall_coord(:,1),wall_coord(:,2),3), axis ij, xlim([0 img_size(1)]), ylim([0 img_size(2)]);
hold on
scatter(wall_fit_coord_px(:,1),wall_fit_coord_px(:,2),'.'), xlim([0 img_size(1)]), ylim([0 img_size(2)]);
hold off

%% unit transfer : pixel(wall) to micron(tracker)
% Variables used in this section are initially defined in 'wall_definition.m'.
wall_coord_micron = wall_coord * scale;
wall_fit_coord_micron = wall_fit_coord_px * scale;
fr2sec = 1/30; % [sec] image acquisition rate (frame/second) e.g., 30hz = 1/30

%% plot entire track sets
data = data_micron;
id_list = [unique(data(:,1))]; % sperm id number list

for count = 1 : numel(id_list)
    id = id_list(count);
    data_id = data((id==data(:,1)),:); % extract corresponding id data
    figure(4152)
    hold on
    scatter(wall_coord_micron(:,1),wall_coord_micron(:,2), 4, [0.1 0.2 0.3],'filled')
    scatter(wall_fit_coord_micron(:,1),wall_fit_coord_micron(:,2), 4, [0.6 0.7 0.8],'filled')
    scatter(data_id(:,2),data_id(:,3), 4,'filled')
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    axis image, xlim([0 x_max_micron]), ylim([0 x_max_micron]);
    title(strcat('ID = ',num2str(id),' (',num2str(count),'th)'))
    pause(0.01)
end
hold off

%% calculate parameters

id_list = [unique(data(:,1))]; % list of sperm id number
sperm_num = numel(id_list);
excel_csv = zeros([1 14]);
for sperm = 1 : sperm_num % sperm corresponds to 'Order' in 'spot exprot.csv'
    
% select sperm
id = id_list(sperm); % id corresponds to 'Track ID' in 'spot exprot.csv'
data_id = data((id==data(:,1)),:); % extract data of corresponding id
data_id = sortrows(data_id,4); % sort data by frame number
frame_id = data_id(:,4);
frame_num = size(frame_id,1); % frame_num corresponds to 'Number of frames' in 'spot exprot.csv'
wall_num = size(wall_coord_micron,1); % total size of the coordinates consisting of wall
wall_fit_num = size(wall_fit_coord_micron,1); % total size of the coordinates consisting of fitting line of the wall

% pre-call variable matrices
distance_id_str = zeros([1 frame_num]);
distance_id_cur = zeros([1 frame_num]);
speed_id = zeros([1 frame_num]);
interlink_vec = zeros([2 frame_num]);
excel_csv_frame = zeros([frame_num 14]);

for frame = 1 : frame_num % frame corresponds to 'Frame number' in 'spot exprot.csv'
    
    % minimal distance from the fitted wall(straight line) to sperm point at each frame
    distance_str = zeros(size([1 wall_num]));
    for coord = 1 : wall_fit_num
        distance_str(coord) = norm(wall_fit_coord_micron(coord,[1 2])-data_id(frame,[2 3]));
    end
    distance_id_str(frame) = min(distance_str);
    
    % minimal distance from the curved wall(pixel-like line) to sperm point at each frame
    distance_cur = zeros(size([1 wall_num]));
    for coord = 1 : wall_num
        distance_cur(coord) = norm(wall_coord_micron(coord,[1 2])-data_id(frame,[2 3]));
    end
    distance_id_cur(frame) = min(distance_cur);
    
    % speed (= total travel length / time interval between the frames) at each frame
    if frame == frame_num
        speed_id(frame) = NaN;
    else
        speed_id(frame) = norm(data_id(frame+1,[2 3]) - data_id(frame,[2 3]))/(frame_id(frame+1)-frame_id(frame));
    end
    
    % total distance of a sperm traveled during measurement
    distance_total = 0;
    for c = 1 : frame_num-1
        distance_total = distance_total + norm(data_id(c+1,[2 3])-data_id(c,[2 3]));
    end
    
    % displacement from the start point to the final point of the tracking points
    displacement = data_id(frame_num,[2 3]) - data_id(1,[2 3]);

    % angle between the wall and track (approaching angle of the sperm to the fitted wall)
    wall_tan = fit_wall.p1; % tangent of the fitted wall
    wall_unit = [1 fit_wall.p1]; % unit vector of the fitted wall
    if norm(displacement) == 0
    diff_rad = NaN; % return Nan when the displacement is zero
    else
    diff_cos = displacement*wall_unit'/norm(displacement)/norm(wall_unit);
    diff_rad = acos(diff_cos);
    end
    if diff_rad > pi/2
        diff_rad = pi-diff_rad;
    end
    
    % mid-point of displacement
    displacement_mid = 0.5*(data_id(frame_num,[2 3]) + data_id(1,[2 3]));

    % minimal distance from the fitted wall to the mid-point
    distance_mid_str = zeros(size([1 wall_fit_num])); 
    for coord = 1 : size(wall_fit_coord_micron,1)
        distance_mid_str(coord) = norm(wall_fit_coord_micron(coord,[1 2]) - displacement_mid);
    end
    distance_mid_str = min(distance_mid_str);
    
    % minimal distance from the curved wall to the mid-point
    distance_mid_cur = zeros(size([1 wall_num])); 
    for coord = 1 : size(wall_coord_micron,1)
        distance_mid_cur(coord) = norm(wall_coord_micron(coord,[1 2]) - displacement_mid);
    end
    distance_mid_cur = min(distance_mid_cur);

    % interlink vector (displacement vector at each change of frame)
    if frame == frame_num
        interlink_vec(:,frame) = NaN;
    else
        interlink_vec(:,frame) = data_id(frame+1,[2 3]) - data_id(frame,[2 3]);
    end
    fr2sec = 1/30; % [sec]

    % straight line-to-sideward ratio
        % assuming the line from the start point to the final point as y=ax+b
        x = data_id(:,2);
        y = data_id(:,3);
        a = (y(frame_num)-y(1)) / (x(frame_num)-x(1));
        b = -x(1)*(y(frame_num)-y(1))/(x(frame_num)-x(1))+y(1);
        d = zeros([1 frame_num]); % distance of the tracking points from y=ax+b 
        sgn = zeros([1 frame_num]); % evaluate +1 if the tracking point is on the left side from y=ax+b, or -1 if it is on the right side
    for i = 1 : frame_num
        A = data_id(i,2); B = data_id(i,3);
        d(i) = sqrt(a^2/(1+a^2)*(A+(b-B)/a)^2);
        sgn(i) = -sign((A-x(1))*a-(B-y(1)));
    end
    d_sgn = d.*sgn;
    width = max(d_sgn(d_sgn>=0))+max(abs(d_sgn(d_sgn<=0))); % sum of the longest distance(absolute value) on the left side and the longest distance on the right side

% save above data
excel_csv_frame(frame,:) = [id frame frame_id(frame) distance_id_str(frame) distance_id_cur(frame) speed_id(frame)/fr2sec frame_num diff_rad*180/pi distance_total norm(displacement) distance_mid_str distance_mid_cur width 0]; % save data to csv (interlink rad is temporarily set to zero) 
end

% interlink angle (angle change of interlink vectors)
interlink_rad = zeros([1 frame_num]);
for frame = 1 : frame_num
    if frame == frame_num
        interlink_rad(frame) = NaN;
    else
        if norm(interlink_vec(:,frame+1))*norm(interlink_vec(:,frame)) == 0
            interlink_rad(frame) = 0;
        else
            angle = interlink_vec(:,frame+1)'*interlink_vec(:,frame)/norm(interlink_vec(:,frame+1))/norm(interlink_vec(:,frame));
            interlink_rad(frame) = abs(acos(angle));
        end
    end
end

% average of the interlink angle for each sperm 
interlink_rad_sum = sum(interlink_rad(1:end-2)); % Elements at the last two frames are excluded because they are invalid (physcially not defined).
% save data
excel_csv_frame(:,end) = interlink_rad_sum*180/pi/(frame_num-2);

% collect calculated data for each sperm
excel_csv = [excel_csv; excel_csv_frame];
end
calculated_parameter_csv = excel_csv(2:end,:);

%% write the matrix of the calculated results into csv file
table_calculated_parameter_csv = array2table(calculated_parameter_csv);
table_calculated_parameter_csv.Properties.VariableNames(1:14) = {'Track ID','Order','Frame number','Distance from straight wall [um]','Distance from curved wall [um]','Speed per frame [um/s]','Number of frames','Radian from wall [deg.]','Total distance of target sperm [um]','Displacement of target sperm [um]','Distance from straight wall to mid-point [um]','Distance from curved wall to mid-point [um]','Width [um]','Interlink radian [deg.]'};
filename_1 = strcat('tracks_mag',num2str(mag),'x.csv');
writetable(table_calculated_parameter_csv, filename_1)
disp(strcat('Calculated parameters are exported as "',filename_1,'"'))

%% pivot table
% Pivot table is created at the request of further calculation between the parameters or selected export.
id_num = numel(id_list);
pivot = zeros([id_num 13]);

pivot(:,1) = id_list; % track_id
for id = 1 : id_num
    extract_id_matrix = calculated_parameter_csv(:,1)==id_list(id);
    
    frame_number = calculated_parameter_csv(extract_id_matrix,3);
    pivot(id,2) = max(frame_number) - min(frame_number); % frame_final - frame_initial

    pivot(id,3) = pivot(id,2) * fr2sec; % track_duration [sec]

    temp = calculated_parameter_csv(extract_id_matrix,8);
    pivot(id,4) = temp(1); % angle_from_wall

    temp = calculated_parameter_csv(extract_id_matrix,9);
    pivot(id,5) = temp(1); % total_distance

    pivot(id,6) = pivot(id,5)./pivot(id,3); % vcl

    temp = calculated_parameter_csv(extract_id_matrix,10);
    pivot(id,7) = temp(1); % displacement
    
    pivot(id,8) = pivot(id,7)./pivot(id,3); % vsl
    
    pivot(id,9) = pivot(id,8)./pivot(id,6); % lin
    
    temp = calculated_parameter_csv(extract_id_matrix,11);
    pivot(id,10) = temp(1); % distance_from_fitted_wall_to_mid_point

    temp = calculated_parameter_csv(extract_id_matrix,12);
    pivot(id,11) = temp(1); % distance_from_curved_wall_to_mid_point
    
    temp = calculated_parameter_csv(extract_id_matrix,14);
    pivot(id,12) = temp(1); % mean_of_interlink_angle
    
    temp = calculated_parameter_csv(extract_id_matrix,13);
    pivot(id,13) = pivot(id,7)/temp(1); % straight line/width
end

table_pivot = array2table(pivot);
table_pivot.Properties.VariableNames(1:13) = {'track_id','f_final-f_initial','track_duration','angle_from_wall','total_distance','vcl','displacement','vsl','lin','distance_from_fitted_wall_to_mid_point','distance_from_curved_wall_to_mid_point','mean_of_interlink_angle','straight_width_ratio'};
filename_2 = strcat('summary_tracks_mag',num2str(mag),'x.csv');
writetable(table_pivot, filename_2)
disp(strcat('Pivot table are exported as "',filename_2,'"'))




