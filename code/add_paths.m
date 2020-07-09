function [current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath,tryouts_path] = add_paths


current_dir = pwd;
% determine if you're local or on the server
if contains(current_dir, 'BC_Research') % local
    
    % Add all the relevant paths on local
    cd ../ % go up 1 levels
    basepath = pwd;
    cd .. % go up another level to BC_Research
    BC_researchpath = pwd;
    toolpath = [BC_researchpath '/matlab_toolbox'];
    export_figpath = [toolpath '/Export_Fig'];
    figpath = [basepath '/figures'];
    tablepath = [basepath '/tables'];
    datapath = [basepath '/data'];
    tryouts_path = [toolpath '/tryouts'];
    % inputsRyan_path = [current_dir '/input_sequences_for_Ryan'];
    % RyanPS6_path = [current_dir '/Ryan_PS6_replicate'];
    
    cd(current_dir)
    
    addpath(basepath)
    addpath(toolpath)
    addpath(export_figpath)
    addpath(figpath)
    addpath(datapath)
    addpath(tryouts_path)
    % addpath(inputsRyan_path)
    % addpath(RyanPS6_path)
    
elseif contains(current_dir, 'gsfs0') % sirius server
    % This is going to be something like /gsfs0/data/gati/next_server/code
    
    % Add all the relevant paths on server
    cd ../ % go up 1 levels
    basepath = pwd;
    cd .. % go up another level to /gsfs0/data/gati, the server home
    server_researchpath = pwd;
    toolpath = [server_researchpath '/matlab_toolbox'];
    addonspath = [server_researchpath '/addons2017a']; % write your versions of missing functions in Matlab 2017a
    export_figpath = [toolpath '/Export_Fig'];
    figpath = [basepath '/figures'];
    tablepath = [basepath '/tables'];
    datapath = [basepath '/data'];
    
    cd(current_dir)

    
    addpath(basepath)
    addpath(toolpath)
    addpath(addonspath)
    addpath(export_figpath)
    addpath(figpath)
    addpath(datapath)
    
    BC_researchpath = server_researchpath; % for outputs
end

% Note for future projects use:
% Use "genpath" to Generate a path that includes myfolder and all folders below it.
% p = genpath('myfolder')