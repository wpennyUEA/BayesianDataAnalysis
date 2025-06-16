function [] = vbmfx_plot_precisions (hier)
% Plot priors and posteriors on group precision
% FORMAT [] = vbmfx_plot_precisions (hier)


M(1).a=hier.a0;
M(1).b=hier.b0;
M(1).name='Priors';
M(1).str='p(\lambda)';

M(2).a=hier.a;
M(2).b=hier.b;
M(2).name='Posteriors';
M(2).str='p(\lambda|Y)';

for i=1:length(M),
    
    a = M(i).a;
    b = M(i).b;
    
    m = a./b;
    v = a./(b.^2);
    s = sqrt(v);
    xmax = max(m)+max(s);
    Nx = 100;
    dx = xmax/Nx;
    x = [dx:dx:xmax];
    
    figure
    hold on
    P = length(a);
    for p=1:P,
        [pr,m(p),v(p)] = mygampdf(x,1/b(p),a(p));
        plot(x,pr);
        disp(sprintf('Mean = %1.3f SD = %1.3f',m(p),sqrt(v(p))));
        
    end
    grid on
    xlabel('\lambda');
   
    ylabel(M(i).str);
    title(M(i).name);
end
