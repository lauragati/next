function [PLM_name, gain_name, gain_title] = give_names(PLM, gain)

% Give names
if  PLM == 1
    PLM_name = 'constant_only'
elseif PLM == -1
    PLM_name = 'mean_only_PLM';
elseif PLM == 2
    PLM_name = 'slope_and_constant'
end

if  gain == 1
    gain_name = 'dgain'
    gain_title = 'Decreasing gain';
elseif gain == 3
    gain_name = 'cgain'
    gain_title = 'Constant gain';
elseif gain == 21
    gain_name = 'again_critCEMP'
    gain_title = 'Endogenous gain, CEMP''s criterion';
elseif gain == 22
    gain_name = 'again_critCUSUM'
    gain_title = 'Endogenous gain, CUSUM criterion';
end