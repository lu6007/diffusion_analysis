% function mark = mark_photobleach_region(intial_solution)
% Mark all the nodes with initial solution value < 10000

% Copyright: Shaoying Lu and Yingxiao Wang 2013
function mark = mark_photobleach_region(initial_solution, varargin)
parameter_name = {'bound','file'};
default_value = {10000,''};
[bound file] = parse_parameter(parameter_name, default_value, varargin);

if(isempty(file)||~exist(file, 'file')),
    num_nodes = size(initial_solution,1);
    mark = zeros(num_nodes,1);
    for i=1:num_nodes,
        if initial_solution(i)<bound,
            mark(i) = 1.;
        end;
    end;
    if ~isempty(file),
        save(file, 'mark');
    end;
else % the file exist 
    res = load(file);
    mark = res.mark;
end;

return;