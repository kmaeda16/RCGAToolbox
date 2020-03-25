% Script InitTransition
% This script is called by RCGA_Main

Results.Transition.time = [];
Results.Transition.neval = [];
Results.Transition.generation = [];
Results.Transition.f = [];
Results.Transition.x = [];
if n_constraint > 0
    Results.Transition.phi = [];
    Results.Transition.g = [];
end
