function fmt = cell2fmt(cell_array)
    fmt = string(cell_array(1));
    for i = 2:length(cell_array)
        fmt = [fmt, strcat('t',string(cell_array(i)))];
    end
    fmt = [fmt, 'n'];
%     if fmt(1) == ''
%        fmt(1) = 'empty';
%     end
    
    fmt = strjoin(fmt, '\');

end