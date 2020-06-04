function [fs, lw, textinterpreter, mathinterpreter, ...
    fs_pres,lw_pres, ...
    grey, silver, maroon, light_coral, light_salmon,dark_green, green, light_green, light_sky_blue, ...
    teal, purple, saddle_brown, light_brown] = plot_configs

fs=40; % fontsize
lw=4; % linewidth
fs_pres = 80; % versions for beamer presi
lw_pres = 6;
textinterpreter = 'latex';
mathinterpreter = 'tex';

% Some color spectra
% grey color (divide by 255)
grey = [128,128,128]/255;
silver = [192,192,192]/255;
maroon = [138 0 0]/255;
light_coral = [240 128 128]/255;
light_salmon = [255,160,122]/255;
dark_green = [0 100 0]/255;
green = [0 128 0]/255;
light_green = [144 238 144]/255;
light_sky_blue = [135 206 250]/255;

teal = [0,128,128]/255;
purple = [128,0,128]/255;
saddle_brown = [139,69,19]/255;
light_brown = [181,101,29]/255;