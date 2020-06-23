function [current_dir, basepath, BC_researchpath,toolpath,export_figpath,figpath,tablepath,datapath,tryouts_path] = add_paths

% Add all the relevant paths
current_dir = pwd;
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

% Note for future projects use:
% Use "genpath" to Generate a path that includes myfolder and all folders below it.
% p = genpath('myfolder')