% PARSE_G2O

function [v e] = parse_g2o(filename)
    disp('==================');
    fprintf('parsing g2o 3D file %s...', filename);
    v = [];
    e = [];
    
    f = fopen(filename);
    
    % read line-wise
    
    while ~feof(f)
        line = fgetl(f);
        if all(line(1:15) == 'VERTEX_SE3:QUAT')
            data = line(17:end);
            vtx = textscan(data, '%f', 8);
            v = [v; vtx{1}'];
        elseif all(line(1:13) == 'EDGE_SE3:QUAT')
            data = line(15:end);
            edg = textscan(data, '%f', 30);
            e = [e; edg{1}'];
        else
            error('unhandled tag');
        end
    end
    
    fprintf(' done.\n');
    fprintf('number of poses:   %d \n',size(v,1));
    fprintf('number of edges:   %d \n\n',size(e,1));
end
