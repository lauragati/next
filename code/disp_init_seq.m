function [init_name, init_title] = disp_init_seq(initialization, n_inputs)

% Give names
if  initialization == 1
    init_name = 'init_TR';
    init_title = ' Initialized at Taylor rule sequence(s) '
elseif initialization == 2
    init_name = 'init_rand';
    init_title = ' Initialized at random sequence(s) '
end

if n_inputs == 3
    disp('Feeding in pi, x, i')
elseif n_inputs==2
    disp('Feeding in x, i')
elseif n_inputs ==1
    disp('Feeding in i')
end