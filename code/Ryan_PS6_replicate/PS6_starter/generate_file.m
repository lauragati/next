%Load Parameters
[param,set] = parameters;

%Generate a mod object
[modl,param,set] = model(param,set);

%Get dimensions of the model
nx = length(modl.X);
ny = length(modl.Y);
neps = size(modl.shck,2);
neq  = nx+ny;
np   = length(struct2cell(param));

%Generate model file
model_func(modl);

%Save model object and idx variables
load v_idx
save  model_object modl *_idx  ny neps nx np
