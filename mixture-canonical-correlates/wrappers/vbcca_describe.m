function [cca] = vbcca_describe (cca,reorder)
% Describe VBCCA model
% FORMAT [cca] = vbcca_describe (cca,reorder)
%
% cca       estimated model from vbcca.m
% reorder   (optional)
%               .by = 'state' or 'context' to reorder by a state mean (mu1)
%               or context mean (mu2)
%               .element which state/context variable to reorder by
%

if nargin > 1
    j = reorder.element;
    switch reorder.by
        case 'state'
            disp(' ');
            disp(sprintf('Reordering cluster indices by state variable x1[%d]',reorder.element));
            for k=1:cca.M,
                x(k) = cca.mu2{k}(j);
            end
        case 'context'
            disp(sprintf('Reordering cluster indices by context variable x2[%d]',reorder.element));
            for k=1:cca.M,
                x(k) = cca.mu2{k}(j);
            end
    end
    
    [tmp,ind] = sort(x);
    new_cca = cca;
    for k=1:cca.M,
        j = ind(k);
        new_cca.pi(k) = cca.pi(j);
        new_cca.W1{k} = cca.W1{j};
        new_cca.W2{k} = cca.W2{j};
        new_cca.mu1{k} = cca.mu1{j};
        new_cca.mu2{k} = cca.mu2{j};
        new_cca.psi1{k} = cca.psi1{j};
        new_cca.psi2{k} = cca.psi2{j};
        new_cca.C1{k} = cca.C1{j};
        new_cca.C2{k} = cca.C2{j};
    end
    cca = new_cca;
end



for k=1:cca.M,
    disp(sprintf('Component %d, prior = %1.3f',k,cca.pi(k)));
    disp('State Mean');
    disp(cca.mu1{k}(:)');
    disp('Context Mean');
    disp(cca.mu2{k}(:)');
end

