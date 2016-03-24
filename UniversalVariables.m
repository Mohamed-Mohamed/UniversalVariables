function [ r,v ] = UniversalVariables ( ro,vo,dt,muo,tol )
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% Inputs:
% ro          : initial position
% vo          : initial velocity
% dt          : time interval
% muo      :
% tol         : tolerance of newton's method
%% Outputs:
% r           : position at time dt
% v           : velocity at time dt
%% Function body
mag_ro=norm(ro)
mag_vo=norm(vo)
vro=dot(ro,vo)/mag_ro
alpha=2/mag_ro-mag_vo^2/muo
syms x n m;
if alpha>0 % ellipse
    s(n)=(sqrt(n)-sin(sqrt(n)))/(sqrt(n))^3
    c(m)=(1-cos(sqrt(m)))/m
elseif alpha<0 % hyperbola
    s(n)=(sinh(sqrt(-n))-sqrt(-n))/(sqrt(-n))^3
    c(m)=(cosh(sqrt(-m))-1)/-m
elseif alpha==0 % parabola
    s(n)=1/6
    c(m)=1/2
end
F(x)=mag_ro*vro/sqrt(muo)*x^2*c(alpha*x^2)+(1-alpha*mag_ro)*x^3*s(alpha*x^2)+mag_ro*x-sqrt(muo)*dt;
F0=sqrt(muo)*abs(alpha)*dt;
cay = NewtonMethod( F(x),F0,tol )
f=1-cay^2/mag_ro*double(c(alpha*cay^2))
g=dt-cay^3/sqrt(muo)*double(s(alpha*cay^2))
r=f*ro+g*vo
f_dot=sqrt(muo)/mag_ro/norm(r)*(alpha*cay^3*double(s(alpha*cay^2))-cay)
g_dot=1-cay^2/norm(r)*double(c(alpha*cay^2))
v=f_dot*ro+g_dot*vo
end