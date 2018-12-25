%%%GZ function

function [gz]= grav_eff_point(x,xm,m,G)
r_dot=x-xm;
r_dot=sqrt(dot(r_dot,r_dot));
r=r_dot;
z=x(3);
zm=xm(3);

gz=G*m*(z-zm)/r^3;
end