function spec_init = CalcInitialSpec(p, p2)
%
% CalcInitialSpec calculates the initial spectrum. At the present, this is
% based on George's code and calculates a spectrum with equal volume in
% each section
%

% Initial spectrum for an exponential

term1 = p.vol0*exp(-p.v_lower/p.vol0);
term2 = p.v_lower.*exp(-p.v_lower/p.vol0);
term3 = p.vol0*exp(-2.0*p.v_lower/p.vol0);
term4 = 2.0*p.v_lower.*exp(-2.0*p.v_lower/p.vol0);

spec_init = p.N0*(term1 + term2 - term3 - term4);


% Numerical calculate things

spec_init_quad = zeros(p.n_sections, 1);

for i_sec = 1 : p.n_sections
    
    spec_init_quad(i_sec) = quadl(@(x) specfun(x, p), p.v_lower(i_sec), p.v_upper(i_sec));
        
end


spec_init = p.av_vol(1)*ones(p.n_sections, 1);

tfactor = 10.^(0 : p.n_sections-1)';

spec_init = max(spec_init ./ tfactor, 1.0e-30);

spec_init = spec_init * p.num_1;

spec_init = spec_init*10;


% % Try a steady state solution as G does
% 
% Vcon        = ones(p.n_sections,1)*p.av_vol(1);
% Vcon(2:end) = Vcon(2:end)/10;        % decrease the initial concentration of others
% Vcon        = p.num_1*Vcon;          % use the initialization data
% 
% 
% Vtry = fsolve(@(x) finitial(x, Vcon(1), p2), Vcon(2:end)*1e7);
% 
% Vcon(2:end)=Vtry/1e7;
% 
% spec_init = Vcon;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = finitial(x, s1conc, p)
%
% Function to optimize to find the steady state spectrum
%


xx=[s1conc;x/1e7];      % have to de scale, add first section conc
FF=CalcCoagDeriv(0,xx,p);
F=FF(2:end)*1e7;        % rescale


function f = specfun(v, p)

f = v .* exp( -v/p.vol0);
f = (p.N0/p.vol0)*f;

