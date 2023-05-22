% define wall coordinates from tiff image
% Assumptions
% 1. The tiff image consists of RGB channels with integer-type elements.
% 2. The uterus wall (where does not include sperms) must be represented in
% white color, or identically, (R=255,G=255,B=255).

filename = 'wall.tif'; % input name of the tiff file to define wall
wall_tif = imread(filename); % read tiff file of the averaged sperm tracking images

wall_sum = double(wall_tif(:,:,1))+double(wall_tif(:,:,2))+double(wall_tif(:,:,3));
wall_sum = max(wall_sum(:))-wall_sum;
% figure,imagesc(wall_sum)
wall = double(wall_sum>=1);
% figure,imagesc(wall)

wall_line_x = double(abs(circshift(wall,1,2) - wall)>=1); wall_line_x(:,1) = 0;
% figure,imagesc(wall_line_x)
wall_line_y = double(abs(circshift(wall,1,1) - wall)>=1); wall_line_y(1,:) = 0;
% figure,imagesc(wall_line_y)
wall_line = double(abs(wall_line_x + wall_line_y)>=1);
% figure,imagesc(wall_line)

img_size = size(wall_line);
wall_coord = [0 0;];
for i = 1 : img_size(1)
    for j = 1 : img_size(2)
        if wall_line(i,j) == 1
            wall_coord = [wall_coord; j i;];
        end
    end
end
wall_coord = wall_coord(2:end,:);


%% export coordinate data to csv file

% reverse y-axis (for plotting in R)
wall_coord_yrev = [wall_coord(:,1), img_size(2)-wall_coord(:,2)+1];

% scaling in micron
scan_size = 776.72; % [um]
scanning_pixel = 512; % acquired image pixel size(square), [px.]
mag = 5; % x
fprintf('\nThe parameters are set as :\n      scan size = %.2f\n      scanning pixel = %u\n      magnification = %.2f\n',scan_size, scanning_pixel, mag);
disp("Move to 'wall_definition.m' to modify.")
px_num = img_size(1); % square
x_max_micron = scan_size / mag * img_size(1)/scanning_pixel;
scale = scan_size / mag /scanning_pixel ; % [um/px.] %size of one pixel()
wall_coord_yrev_micron = wall_coord_yrev * scale;


% write coordinate data into csv file
table_wall_coord_csv = array2table([wall_coord_yrev,wall_coord_yrev_micron]);
table_wall_coord_csv.Properties.VariableNames(1:4) = {'x','y_reversed','x_micron','y_reversed_micron'};
writetable(table_wall_coord_csv,strcat(filename(1:end-4),'_coordinates.csv'))
disp(strcat('Wall coordinates are exported as "', filename(1:end-4),'_coordinates.csv"'))

