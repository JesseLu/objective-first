function [reduced_line_segments, did_change] = simplify_triangulation(line_segments)
    did_change = 0;
    
    original_dims = size(line_segments);
    reduced_line_segments = zeros(original_dims);
    num_reduced = 0;
    
    used_line_mask = zeros(original_dims(1), 1);
        
    for n = 1 : 1 : original_dims(1)
        if (used_line_mask(n) == 1)
            continue;
        end
        
        next_line = line_segments(n, :);
        x_start = next_line(1);
        y_start = next_line(2);
        x_end = next_line(3);
        y_end = next_line(4);

        slope = -1;
        inf_slope = 1;

        if (abs(x_end - x_start) > 0)
            slope = (y_end - y_start) / (x_end - x_start);
            inf_slope = 0;
        end
                
        for k = (n + 1) : 1 : original_dims(1)
            if (used_line_mask(k) == 1)
                continue;
            end
                                    
            new_next_line = line_segments(k, :);
            new_x_start = new_next_line(1);
            new_y_start = new_next_line(2);
            new_x_end = new_next_line(3);
            new_y_end = new_next_line(4);

            cur_slope = -1;
            cur_inf_slope = 1;

            if (abs(new_x_end - new_x_start) > 0)
                cur_slope = (new_y_end - new_y_start) / (new_x_end - new_x_start);
                cur_inf_slope = 0;
            end     
            
            if (((slope == cur_slope) && (inf_slope == 0) && (cur_inf_slope == 0)) || ((inf_slope == 1) && (cur_inf_slope == 1)))
                if ( ((new_x_start == x_start) && (new_y_start == y_start)) || ...
                   ((new_x_start == x_end) && (new_y_start == y_end)) || ...
                   ((new_x_end == x_start) && (new_y_end == y_start)) || ...
                   ((new_x_end == x_end) && (new_y_end == y_end)) )
                                           
                    repl_x_start = x_start;
                    repl_x_end = x_end;
                    repl_y_start = y_start;
                    repl_y_end = y_end;

                    % The normal case that we expect where one segment
                    % links up to the other.  The special cases are handled
                    % below.
                    if ((new_x_start == x_start) && (new_y_start == y_start))
                        repl_x_start = new_x_end;
                        repl_y_start = new_y_end;
                    elseif ((new_x_start == x_end) && (new_y_start == y_end))
                        repl_x_end = new_x_end;
                        repl_y_end = new_y_end;
                    elseif ((new_x_end == x_start) && (new_y_end == y_start))
                        repl_x_start = new_x_start;
                        repl_y_start = new_y_start;
                    elseif ((new_x_end == x_end) && (new_y_end == y_end))
                        repl_x_end = new_x_start;
                        repl_y_end = new_y_start;
                    end
                        
                    % Check to see if first line segment consumes the
                    % second.
                    if ( (new_x_start <= max(x_start, x_end)) && (new_x_start >= min(x_start, x_end)) && ...
                         (new_x_end <= max(x_start, x_end)) && (new_x_end >= min(x_start, x_end)) && ...
                         (new_y_start <= max(y_start, y_end)) && (new_y_start >= min(y_start, y_end)) && ...
                         (new_y_end <= max(y_start, y_end)) && (new_y_end >= min(y_start, y_end)) )
                        repl_x_start = x_start;
                        repl_x_end = x_end;
                        repl_y_start = y_start;
                        repl_y_end = y_end;
                    % Check to see if second line segment consumes the
                    % first.
                    elseif ( (x_start <= max(new_x_start, new_x_end)) && (x_start >= min(new_x_start, new_x_end)) && ...
                         (x_end <= max(new_x_start, new_x_end)) && (x_end >= min(new_x_start, new_x_end)) && ...
                         (y_start <= max(new_y_start, new_y_end)) && (y_start >= min(new_y_start, new_y_end)) && ...
                         (y_end <= max(new_y_start, new_y_end)) && (y_end >= min(new_y_start, new_y_end)) )
                        repl_x_start = new_x_start;
                        repl_x_end = new_x_end;
                        repl_y_start = new_y_start;
                        repl_y_end = new_y_end;
                    % Check if both the lines we had were the same.
                    elseif ((repl_x_start == repl_x_end) && (repl_y_start == repl_y_end))
                        repl_x_start = x_start;
                        repl_x_end = x_end;
                        repl_y_start = y_start;
                        repl_y_end = y_end;
                    end
                    
                    x_start = repl_x_start;
                    x_end = repl_x_end;
                    y_start = repl_y_start;
                    y_end = repl_y_end;
                    
                    used_line_mask(k) = 1;
                    break;
                end
            end
        end
        
        used_line_mask(n) = 1;
        
        reduced_line_segments(num_reduced + 1, :) = [x_start, y_start, x_end, y_end];
        num_reduced = num_reduced + 1;
        
    end
    
    reduced_line_segments = reduced_line_segments(1:num_reduced,:);

    if (num_reduced < original_dims(1))
        did_change = 1;
    end

end

