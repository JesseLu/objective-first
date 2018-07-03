clc; clf; close all;

input_mode = 1;
output_mode = 2;

load_epsilon_filename = sprintf('adjoint_eps_bin_%d_%d.csv', input_mode, output_mode);

% Read in a binary epsilon structure in order to start.
eps = csvread(load_epsilon_filename);

[eps_rows, eps_cols] = size(eps);
num_eps_pixels = eps_rows * eps_cols;

max_eps_value = 12.25;
min_eps_value = 1.0;
eps_midpt = 0.5 * (max_eps_value + min_eps_value);

% convert to binary mask
eps_bin = (eps > eps_midpt);

% Length scale of each epsilon pixel
scale = 4.0;

% Generate DeScribe output or Python ouptut
output_type = 'python';

assert((strcmp(output_type, 'python') == 1) || ...
    (strcmp(output_type, 'describe') == 1));

file_name = './output/waveguide_turtle.py';

if (strcmp(output_type, 'describe'))
    file_name = './output/waveguide_describe.txt';
end

file_id = fopen(file_name, 'w');

if (strcmp(output_type, 'python'))
    fprintf(file_id, 'from turtle import Turtle\n');
    fprintf(file_id, 't = Turtle()\n');
    fprintf(file_id, 't.speed(0)\n');
end

call_draw_line = @(pt0, pt1) draw_line(pt0, pt1, scale, file_id, output_type);

% get x-y points

x_pts = zeros(1, num_eps_pixels);
y_pts = zeros(1, num_eps_pixels);
pt_idx = 1;

for y = 1 : 1 : eps_rows
    for x = 1 : 1 : eps_cols
        if eps_bin(y, x)
            x_pts(pt_idx) = x - 1;
            y_pts(pt_idx) = y - 1;
            pt_idx = pt_idx + 1;
        end
    end
end

num_pts = pt_idx - 1;

x_pts_trim = x_pts(1:num_pts);
y_pts_trim = y_pts(1:num_pts);

% Do a triangulation of the structure
triangulation = delaunay(x_pts_trim, y_pts_trim);
[num_triangles, triangle_length] = size(triangulation);

% Remove the holes in the epsilon structure by deleting segments whose
% centroids fall in a gap in the structure;
remove_holes = zeros(num_triangles * 3, 4);
success_idx = 1;
for n = 1 : 1 : num_triangles
    pt0_x = x_pts_trim(triangulation(n, 1));
    pt0_y = y_pts_trim(triangulation(n, 1));

    pt1_x = x_pts_trim(triangulation(n, 2));
    pt1_y = y_pts_trim(triangulation(n, 2));
    
    pt2_x = x_pts_trim(triangulation(n, 3));
    pt2_y = y_pts_trim(triangulation(n, 3));
    
    centroid_0_x = 0.5 * (pt0_x + pt1_x);
    centroid_0_y = 0.5 * (pt0_y + pt1_y);

    centroid_1_x = 0.5 * (pt0_x + pt2_x);
    centroid_1_y = 0.5 * (pt0_y + pt2_y);
    
    centroid_2_x = 0.5 * (pt1_x + pt2_x);
    centroid_2_y = 0.5 * (pt1_y + pt2_y);
    
    snap_0_x = round(centroid_0_x);
    snap_0_y = round(centroid_0_y);
    
    snap_1_x = round(centroid_1_x);
    snap_1_y = round(centroid_1_y);
    
    snap_2_x = round(centroid_2_x);
    snap_2_y = round(centroid_2_y);
    
    check_0 = eps_bin(snap_0_y + 1, snap_0_x + 1);
    check_1 = eps_bin(snap_1_y + 1, snap_1_x + 1);
    check_2 = eps_bin(snap_2_y + 1, snap_2_x + 1);

    if (check_0)
        remove_holes(success_idx, :) = ...
            [pt0_x, pt0_y, pt1_x, pt1_y];
        success_idx = success_idx + 1;
    end
    
    if (check_1)
        remove_holes(success_idx, :) = ...
            [pt0_x, pt0_y, pt2_x, pt2_y];
        success_idx = success_idx + 1;
    end
    
    if (check_2)
        remove_holes(success_idx, :) = ...
            [pt1_x, pt1_y, pt2_x, pt2_y];  
        success_idx = success_idx + 1;
    end
end

num_remaining_lines = success_idx - 1;
remove_holes = remove_holes(1:num_remaining_lines, :);

% figure; triplot(remove_holes, x_pts_trim, y_pts_trim);
% figure; triplot(triangulation, x_pts, y_pts);

filtered_line_segments = zeros(size(remove_holes));
num_filtered_lines = 0;

% Filter for duplicate lines
for n = 1 : 1 : num_remaining_lines   
    get_line_segment = remove_holes(n, :);
    
    already_exists = 0;
    
    for k = 1 : 1 : num_filtered_lines
        filtered_line = filtered_line_segments(k, :);
        
        if ( (get_line_segment(1) == filtered_line(1)) && ...
             (get_line_segment(2) == filtered_line(2)) && ...
             (get_line_segment(3) == filtered_line(3)) && ...
             (get_line_segment(4) == filtered_line(4)) )
            already_exists = 1;
            break;
        end
    end
        
    if (already_exists == 0)
        filtered_line_segments(num_filtered_lines + 1, :) = get_line_segment;
        num_filtered_lines = num_filtered_lines + 1;
    end
end

filtered_line_segments = filtered_line_segments(1:num_filtered_lines,:);

% Sort by x coordinate and then by y
filtered_line_segments = sortrows(filtered_line_segments, 1);
filtered_line_segments = sortrows(filtered_line_segments, 2);
did_change = 1;
while (did_change)
    [filtered_line_segments, did_change] = simplify_triangulation(filtered_line_segments);
end

reduced_dims = size(filtered_line_segments);
num_reduced_segments = reduced_dims(1);

% Re-sort in case things moved around a lot in the reduction
filtered_line_segments = sortrows(filtered_line_segments, 1);
filtered_line_segments = sortrows(filtered_line_segments, 2);

for n = 1 : 1 : num_reduced_segments
    filtered_line = filtered_line_segments(n, :);
    
    call_draw_line([filtered_line(1), filtered_line(2), 0], ...
        [filtered_line(3), filtered_line(4), 0]);
end
    
if (strcmp(output_type, 'python'))
    fprintf(file_id, 't.screen.exitonclick()\n');
end

