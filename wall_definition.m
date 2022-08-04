%% Define wall coordinates from image
% Assumptions
% 1. The tiff image consists of RGB channels with integer-type elements.
% 2. The uterus wall (where does not include sperms) is represented in white color, or identically, (R=255,G=255,B=255).

%% Import wall image
filename = 'wall.tif'; % Input name of the image file to define wall.
wall_img = imread(filename); % Read the image file.
% figure, imagesc(wall_1color_inv), axis image, colorbar; title('')

%% Process wall image
wall_1color = double(wall_img(:,:,1)) + double(wall_img(:,:,2)) + double(wall_img(:,:,3)); % Create 1-dimensional color image by summing RGB values of each pixel.
wall_1color_inv = max(wall_1color(:)) - wall_1color; % Invert color of the image
% figure, imagesc(wall_1color_inv), axis image, colorbar; title(sprintf("' %s ', monochrome, inverted",filename));
wall = double(wall_1color_inv >= 1); % Binarize values of the image : pixels with the value of 0 become to 0, and the larger values than 1 become to 1.
% figure, imagesc(wall), axis image, colorbar; title(sprintf("' %s ', monochrome, inverted, binarized",filename));

% Define wall from the processed image
wall_line_x = double(abs(imtranslate(wall,[1 0]) - wall) >= 1); % Evaluate 1 to the pixels whose value change when horizontally shifted by 1 pixel, and evaluate 0 for the other pixels.
wall_line_x(:,1) = 0; % Make the first vertical line 0 because it is meaninglessly evaluated.
% figure, imagesc(wall_line_x), axis image, title('wall line x');
wall_line_y = double(abs(imtranslate(wall,[0 1]) - wall) >= 1); % Evaluate 1 to the pixels whose value change when vertically shifted by 1 pixel, and evaluate 0 for the other pixels.
wall_line_y(1,:) = 0; % Make the first horizontal line 0 because it is meaninglessly evaluated.
% figure, imagesc(wall_line_y), axis image, title('wall line y');
wall_line = double(abs(wall_line_x + wall_line_y) >= 1); % Define 'wall_line' by summing 'wall_line_x' and 'wall_line_y' and binarzing.
% figure, imagesc(wall_line), axis image, title('wall line');

%% Convert binary wall image to coordinates
img_size = size(wall_line);
wall_coord = [0 0;];
for i = 1 : img_size(1)
    for j = 1 : img_size(2)
        if wall_line(i,j) == 1
            wall_coord = [wall_coord; j i;]; % The axis of the image [i j] is reversed to [j i] as converted into coordinates in order to follow the axis of track data. 
        end
    end
end
wall_coord = wall_coord(2:end,:); % 'Wall_coord' contains coordinates of the wall defined in the 2D image 'wall_line' (but the axis is reversed).


