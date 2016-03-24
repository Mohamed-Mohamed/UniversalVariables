function [ sol ] = NewtonMethod( f,f0,tol )
% this function is about Newton's numerical method
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs:
% f : the function to be ittrated as a function of x
% f0 : the initial value of itteration
% tol : solution tolerance
%% outputs:
% sol : the itteration solution 
%% Function body
syms x;
f(x) = matlabFunction(f);
f_dash(x)=diff(f(x));
sol=f0;
ratio=double(((f(f0))/(f_dash(f0))));
while double(abs(ratio))>= tol;
    sol=double(sol-ratio);
    ratio=double((f(sol))/(f_dash(sol)));
end
end

