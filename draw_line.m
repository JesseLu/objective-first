function draw_line(pt0, pt1, scale, file_id, output_type)
    pt0 = scale * pt0;
    pt1 = scale * pt1;
    
    command = sprintf('t.penup()\nt.setpos(%f, %f)\nt.pendown()\nt.goto(%f, %f)\n', ...
        pt0(1), pt0(2), pt1(1), pt1(2));    
    if (strcmp(output_type, 'describe'))
        command = sprintf('%f %f %f\n%f %f %f\nwrite\n', pt0(1), pt0(2), pt0(3), pt1(1), pt1(2), pt1(3));
    end

    fprintf(file_id, command);
end

