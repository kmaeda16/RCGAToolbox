
function [f g r] = b5_obj(par)

r=[];

load b5_data

if(length(par)==86)
    inputs.model.n_par = 120;
    par = [par,1,0,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1];
end

for iexp=1:inputs.exps.n_exp

        if(~isfield(inputs.exps,'pend') || isempty(inputs.exps.pend{iexp}))
            inputs.exps.pend{iexp}=zeros(4,length(privstruct.t_con{iexp})-1);
        end
        
        if(size(inputs.model.J,1)>0)
            jacobian=1;
        else
            jacobian=0;
        end

[x iflag]=feval('cvodesg_b5',...inputs.ivpsol.ivpmex,...
  inputs.model.n_par,...Number of parameters
  0,...activate sensitivity
  par,...parameters
  ones(1,inputs.model.n_par),...isoptpar
  privstruct.t_in{iexp},...t initial
  privstruct.t_f{iexp},...t final
  length(privstruct.t_s{iexp}),...n_times
  privstruct.t_s{iexp},...t
  inputs.model.n_st,...n_states
  privstruct.exp_y0{iexp},...Intiial values for state variables
  inputs.model.n_stimulus,...  inputs.model.n_stimulus
  length(privstruct.t_con{iexp}),...Number of controls changes(handle discontinuities)
  privstruct.t_con{iexp},...Times of such discontinuities
  privstruct.u{iexp}',... Values of stimuli
  inputs.exps.pend{iexp},... Slope of the line
  inputs.ivpsol.rtol,...reltol
  inputs.ivpsol.atol,...atol
  inputs.ivpsol.max_step_size,...max_step_size
  inputs.ivpsol.ivp_maxnumsteps,...max_num_steps
  500,...%max_error_test_fails
  0,...Sensitivity analysis=false
  jacobian,...%use of the jacobian
  iexp-1);%  iexp-1); %experiment number

  r=[r; (x(:,inputs.exps.index_observables{iexp})-inputs.exps.exp_data{iexp})./0.05];
  
end
 
r=r(:);
f=sum(r.^2);
g=0;

end
