function [xmax,fmax,alpha,f0] = line_search (L,P1,P2,P3,x0,grad)
% Find maximum of function L in search direction grad
% FORMAT [xmax,fmax,alpha,f0] = line_search (L,P1,P2,P3,x0,grad)
%
% L         Name of function
% P1,P2,P3  Arguments to function
% x0        Expansion point
% grad      gradient
%
% xmax      maximum
% alpha     step size

fold = feval(L,P1,P2,P3,x0);
f0 = fold;

alpha = 0.2;
x = x0 + alpha*grad;
fnew = feval(L,P1,P2,P3,x);

last_move=1;
for i = 1:8,
    if fnew > fold
        new_move=last_move;
        xmax=x;
        fmax=fnew;
        alpha_max=alpha;
    else
        new_move=-1*last_move;
    end
    if new_move
        % increase
        alpha=alpha*(1+rand(1));
    else 
        % decrease
        alpha=alpha*rand(1);
    end
        
    x = x0 + alpha*grad;
    fold = fnew;
    fnew = feval(L,P1,P2,P3,x);
end

%%% Get rid of randomisation
%%% Golden search ?????


