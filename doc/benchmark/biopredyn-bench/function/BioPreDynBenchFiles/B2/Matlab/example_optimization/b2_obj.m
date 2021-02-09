
function [objective,constraints,residuals] = b2_obj(par)

global p
p = par;

%% EXPERIMENTAL DATA

% PEP data (x8):
x8   =[
[0, 2.67]   ;
[.15, 1.99] ;
[.3, 2.1]   ;
[.45, 2.09] ;
[.6, 1.84]  ;
[.8, 2.31]  ;
[5.5, 2.76] ;
[12, 3.05]  ;
[21.5, 2.42];
[31.5, 2.23];
[61, 2.52]  ;
[90, 2.81]  ;
[120.5, 2.71];
[180.5, 2.71];
[300.5, 2.7] ];

% Glucose data (x18):
x18 = [
0,    2;               % 0.5277777778e-1 ; <- this would be the value at 0-, not 0+ (that is, before injecting a pulse)
5.5,  1.255555556     ;
13.5, 1.311111111     ;
31,   1.283333333     ;
61,   0.8611111111    ;
91,   0.5972222223    ;
151,  0.9611111112e-1 ;
181,  0.4333333334e-1 ;
212.5,0.5055555556e-1 ;
241,  0.4777777778e-1 ;
270.5,0.4777777778e-1 ;
301,  0.6000000000e-1 ];

% G6p data (x6):
x6   =[
0,     3.48;
0.15,  4.39;
0.3,   4.76;
0.45,  4.86;
0.6,   4.65;
0.8,   4.75;
5.5,   5.52;
12,    5.86;
21.5,  4.39;
31.5,  3.6 ;
61,    3.83;
90,    4.3 ;
120.5, 4.05;
180.5, 3.27;
300.5, 3.38];

% PYR data (x13):
x13   =[
[0, 2.67]   ;
[.15, 4.07] ;
[.3, 3.71]  ;
[.45, 3.19] ;
[.6, 3.57]   ;
[.8, 3.14]   ;
[5.5, 2.38]   ;
[12, 3.71]   ;
[21.5, 3.19]   ;
[31.5, 5.24]   ;
[61, 4.47]   ;
[90, 3.62]   ;
[120.5, 3.62]   ;
[180.5, 2.86]   ;
[300.5, 2.4]   ];

% F6P data (x3):
x3   =[
0,     0.6;
0.15,  0.62;
0.3,   0.66;
0.45,  0.74;
0.6,   0.62;
0.8,   0.75;
5.5,   0.92;
12,    1.15;
21.5,  0.57;
31.5,  0.46;
61,    0.57;
90,    0.57;
120.5, 0.69;
180.5, 0.46;
300.5, 0.46];

% G1P data (x5):
x5   =[
[0, .6525]   ;
[2, 1.35]   ;
[16, .83]   ;
[19, .83]   ;
[31, .78]   ;
[57, .84]   ;
[91.5, .64]   ;
[150.5, .74]   ;
[299, .7]   ];

% 6PG data (x9):
x9 =[
[0, .81];
[3.5, 1.01];
[4, .92];
[12, 1.15];
[12.25, 1.19];
[21, 1.06];
[25.75, 1.1];
[30, 1.05];
[32.25, 1.08];
[58.5, .97];
[59, 1.01];
[119.75, .92];
[124, .89] ;
[178, .74];
[180, .88];
[209, .8] ];

% FDP data (x4):
x4 =[
[0, .27]   ;
[4.5, .19]   ;
[11, .56]   ;
[20, 1]   ;
[30, 2.83]   ;
[60, 1.5]   ;
[90, 2.26]   ;
[119.5, 2.4]   ;
[180, 1.25]   ;
[239.5, .7e-1]   ;
[300, .2e-1]   ];

% GAP data (x7 -- x10 in the original model):
x7   =[
[0, .22] ;
[4.5, .28];
[11, .32];
[20, .31];
[30, .24];
[60, .30];
[90, .18];
[119.5, .22];
[180, .21]  ;
[239.5, .22];
[300, .20]   ];

%=======================================
% EXPERIMENTAL DATA 1: cpep, cg6p, cpyr, cf6p
%=======================================
n_obs{1} = 4; 
exp_data{1}			= [ ...
									1.99,	4.39,	4.07,	0.62;
									2.1,	4.76,	3.71,	0.66;
									2.09,	4.86,	3.19,	0.74;
									1.84,	4.65,	3.57,	0.62;
									2.31,	4.75,	3.14,	0.75;
									2.76,	5.52,	2.38,	0.92;
									3.05,	5.86,	3.71,	1.15;
									2.42,	4.39,	3.19,	0.57;
									2.23,	3.6,	5.24,	0.46;
									2.52,	3.83,	4.47,	0.57;
									2.81,	4.3,	3.62,	0.57;
									2.71,	4.05,	3.62,	0.69;
									2.71,	3.27,	2.86,	0.46;
									2.7,	3.38,	2.4,	0.46;
									];
std_dev{1} = 0.15*ones(1,n_obs{1});

%=======================================
% EXPERIMENTAL DATA 2: cglcex
%=======================================
n_obs{2} = 1; 
exp_data{2}			= [ ...
									1.255555556;
									1.311111111;
									1.283333333;
									0.8611111111;
									0.5972222223;
									0.09611111112;
									0.04333333334;
									0.05055555556;
									0.04777777778;
									0.04777777778;
									0.06;
									];
std_dev{2} = 0.15*ones(1,n_obs{2});

%=======================================
% EXPERIMENTAL DATA 3: cg1p
%=======================================
n_obs{3} = 1; 
exp_data{3}			= [ ...
									1.35;
									0.83;
									0.83;
									0.78;
									0.84;
									0.64;
									0.74;
									0.70;
									];
std_dev{3} = 0.15*ones(1,n_obs{3});

%=======================================
% EXPERIMENTAL DATA 4: cpg
%=======================================
n_obs{4} = 1; 
exp_data{4}			= [ ...
									1.01;
									0.92;
									1.15;
									1.19;
									1.06;
									1.10;
									1.05;
									1.08;
									0.97;
									1.01;
									0.92;
									0.89;
									0.74;
									0.88;
									0.80;
									];
std_dev{4} = 0.15*ones(1,n_obs{4});

%=======================================
% EXPERIMENTAL DATA 5: cfdp, cgap
%=======================================
n_obs{5} = 2; 
exp_data{5}			= [ ...
									0.19,	0.28;
									0.56,	0.32;
									1,	0.31;
									2.83,	0.24;
									1.5,	0.30;
									2.26,	0.18;
									2.4,	0.22;
									1.25,	0.21;
									0.07,	0.22;
									0.02,	0.20;
									];
std_dev{5} = 0.15*ones(1,n_obs{5});


%% INTEGRATION

% Initial concentrations of the 18 states:
cdhap_init						= 0.167;
ce4p_init						= 0.098;
cpg2_init						= 0.399;
cpg3_init						= 2.13;
cpgp_init						= 0.008;
crib5p_init						= 0.398;
cribu5p_init					= 0.111;
csed7p_init						= 0.276;
cxyl5p_init						= 0.138;
cf6p_init						= 0.6;
cfdp_init						= 0.272;
cg1p_init						= 0.653;
cg6p_init						= 3.48;
cgap_init						= 0.218;
cglcex_init						= 2;
cpep_init						= 2.67;
cpg_init						= 0.808;
cpyr_init						= 2.67;

x0								= zeros(1, 18);
x0(1)							= cdhap_init;
x0(2)							= ce4p_init;
x0(3)							= cf6p_init;
x0(4)							= cfdp_init;
x0(5)							= cg1p_init;
x0(6)							= cg6p_init;
x0(7)							= cgap_init;
x0(8)							= cpep_init;
x0(9)							= cpg_init;
x0(10)							= cpg2_init;
x0(11)							= cpg3_init;
x0(12)							= cpgp_init;
x0(13)							= cpyr_init;
x0(14)							= crib5p_init;
x0(15)							= cribu5p_init;
x0(16)							= csed7p_init;
x0(17)							= cxyl5p_init;
x0(18)							= cglcex_init;
   
% Integration timepoints:
tspan   = 0:0.15:301;

% Nominal values of the 5 known parameters:
cfeed							= 110.96;
Dil								= 2.78e-05;
mu								= 2.78e-05;
cytosol							= 1;
extracellular					= 1;

% Complete parameter vector (both unknown and known parameters):
params = [p,cfeed,Dil,mu,cytosol,extracellular];

% Integrate the model with the radau solver:
[x iflag] = radau5g_b2(...
18, ...        % inputs.model.n_st
0,  ...        % inputs.exps.t_in{iexp}
tspan, ...     % privstruct.vtout{iexp}
2007,...       % n_vtout
tspan, ...     % privstruct.t_int{iexp}
2007, ...      % n_int 
[0 301], ...   % privstruct.t_con{iexp}
1,   ...       % ncon{iexp}
x0,...         % initial conditions
1e-4, ...      % inputs.ivpsol.rtol
1e-6, ...      % inputs.ivpsol.atol
params, ...    % par
121,  ...      % size(par,2)
0,     ...     % ipar
1,    ...      % 1
1,    ...      % privstruct.u{iexp}
0,    ...      % inputs.exps.pend{iexp}
1,    ...      % n_u
1    ...       % privstruct.iflag
);


%% OBJECTIVE FUNCTION

% indices of the simulated x's that correspond to existing measurements:
x8_t  = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
x6_t  = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
x13_t = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
x3_t  = [2 3 4 5 6 38 81 144 211 408 601 804 1204 2004]';
x18_t = [38 91 208 408 608 1008 1208 1418 1608 1804 2007]';
x5_t  = [14 108 128 208 381 611 1004 1994]';
x9_t  = [24 28 81 83 141 173 201 216 391 394 799 828 1188 1201 1394]';
x4_t  = [31 74 134 201 401 601 798 1201 1598 2001]';
x7_t  = [31 74 134 201 401 601 798 1201 1598 2001]';

if size(x,1)~=2007
    x = 1e6*ones(2007,18);
end

ms{1} = [x(x8_t,8),x(x6_t,6),x(x13_t,13),x(x3_t,3)];     
ms{2} = x(x18_t,18);
ms{3} = x(x5_t,5);
ms{4} = x(x9_t,9);
ms{5} = [x(x4_t,4),x(x7_t,7)];

% objective function from AMIGO (llk, homo_var):
f        = 0.0;
g        = [];
n_exp    = 5;
atol     = 1e-6;

for iexp=1:n_exp
    error_matrix=[];
    error_data = 0.15*exp_data{iexp};
    error_data(error_data <= atol) = atol;    %To avoid /0
    for i_obs=1:n_obs{iexp}
        error_matrix=[error_matrix; (ms{iexp}(:,i_obs)-exp_data{iexp}(:,i_obs))./(error_data(:,i_obs))];
    end
    g=[g error_matrix'];
    f=sum(g.^2);
    if isreal(f)==0 || isnan(f)==1
        f=1e20;
    end
end

residuals = g';
constraints = 0;
objective = f;



