function [diam_i, total_flux, total_mass, nspec_i, masspec_i, fluxspec_i, p] = coag_driver(nsections, fd, stickiness)
%
% COAG_DRIVER is the main command-line driver for calculating particle size
% spectra. 
%
% USAGE:
%   Type coag_driver at the command line and follow the instructions
%
% HISTORY:
%   23-09-21: First Cut
%
% Adrian Burd, University of Georgia
%
%

%close all
%clear all

%% Set up and get user options
%

[p, opt] = SetUpCoag(nsections, fd);

p.n_sections = nsections;
p.alpha      = stickiness;


%% Calculate the sectionally integrated coagulation kernels
%

disp('Calculating kernels')

b_brown = CalcBetas(p);
b_brown.b1 = b_brown.b1*p.conBr*p.day_to_sec*p.alpha;
b_brown.b2 = b_brown.b2*p.conBr*p.day_to_sec*p.alpha;
b_brown.b3 = b_brown.b3*p.conBr*p.day_to_sec*p.alpha;
b_brown.b4 = b_brown.b4*p.conBr*p.day_to_sec*p.alpha;
b_brown.b5 = b_brown.b5*p.conBr*p.day_to_sec*p.alpha;
 
p.kernel='KernelCurSh';
b_shear = CalcBetas(p);
b_shear.b1 = b_shear.b1*p.gamma*p.day_to_sec*p.alpha;
b_shear.b2 = b_shear.b2*p.gamma*p.day_to_sec*p.alpha;
b_shear.b3 = b_shear.b3*p.gamma*p.day_to_sec*p.alpha;
b_shear.b4 = b_shear.b4*p.gamma*p.day_to_sec*p.alpha;
b_shear.b5 = b_shear.b5*p.gamma*p.day_to_sec*p.alpha;
b_shear.b25 = b_shear.b25*p.gamma*p.day_to_sec*p.alpha;
 
p.kernel='KernelCurDS';
b_ds    = CalcBetas(p);
b_ds.b1 = b_ds.b1*p.setcon*p.day_to_sec*p.alpha;
b_ds.b2 = b_ds.b2*p.setcon*p.day_to_sec*p.alpha;
b_ds.b3 = b_ds.b3*p.setcon*p.day_to_sec*p.alpha;
b_ds.b4 = b_ds.b4*p.setcon*p.day_to_sec*p.alpha;
b_ds.b5 = b_ds.b5*p.setcon*p.day_to_sec*p.alpha;
b_ds.b25 = b_ds.b25*p.setcon*p.day_to_sec*p.alpha;

% Pack up the betas and store them in a new structure that will get passed
% to the derivative and jacobian calculation routines

p2.b1 =  b_brown.b1 + b_shear.b1 + b_ds.b1;
p2.b2 =  b_brown.b2 + b_shear.b2 + b_ds.b2;
p2.b3 =  b_brown.b3 + b_shear.b3 + b_ds.b3;
p2.b4 =  b_brown.b4 + b_shear.b4 + b_ds.b4;
p2.b5 =  b_brown.b5 + b_shear.b5 + b_ds.b5;

p2.b25 = p2.b2 - p2.b3 - p2.b4 - p2.b5;

%% NOTE WE HAVE NOTE MULTIPLIED BY THE STICKINESS!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate the linear terms in the population balance equation
%

p2.growth    = CalcGrowth(p);
p2.sink_loss = CalcSinkingLoss(p);

p2.linear   = p2.growth - p2.sink_loss;

%% Calculate disaggregation terms: Set them to zero here if you don't want them.

p2.disagg_minus = p.c3*diag(p.c4.^(1 : p.n_sections));
p2.disagg_plus  = p.c3*diag(p.c4.^(2:p.n_sections),-1);

%p2.disagg_minus = zeros(size(p2.disagg_minus));
%p2.disagg_plus  = zeros(size(p2.disagg_plus));

%% Initial Size Spectrum
%  Caclculate the initial size spectrum used for both estimation of the
%  steady state or evolving solution - put both in later versions

spec_init = CalcInitialSpec(p, p2);

%% Integrate Coagulation Equations
%  Set up for integrating over time

disp('Solving ODEs')

calcomp = 1:p.n_sections;
abs_tol = 1.0e-18;        % Absolute Tolerance baseline
rel_tol = 3.0e-14;        % Relative tolerance

at = (abs_tol * 1.5 .^ (-(calcomp-1)));

t_span = p.t_init : p.delta_t : p.t_final;

non_neg_indx = [1 : p.n_sections];

ode_options = odeset('RelTol', rel_tol, 'Refine', 0, 'AbsTol', at, 'Jacobian', @CalcCoagJac, 'NonNegative', non_neg_indx);
%ode_options = odeset('Jacobian', @CalcCoagJac);

[t_out, y] = ode15s(@CalcCoagDeriv, t_span, spec_init, ode_options, p2); 

%% Output
    
[diam_i, total_flux, total_mass, nspec_i, masspec_i, fluxspec_i] = CoagOutput(p, p2, t_out, y);

% Numerical spectrum

delta_vol = (p.v_upper - p.v_lower)';
delta_vol_mat = delta_vol(ones(length(t_out),1), :);

av_vol = p.av_vol';
av_vol_mat = av_vol(ones(length(t_out),1), :);

nv_numerical = y./av_vol_mat./delta_vol_mat;

