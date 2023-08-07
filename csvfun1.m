function f = csvfun1(v, t, p)

vprime = -v./p.vol0;
tau    = p.N0*p.beta_zero*t;

f = exp(2*vprime/(tau+2));
f = (4*p.N0/(p.vol0*(tau+2)*(tau+2)))*f;

