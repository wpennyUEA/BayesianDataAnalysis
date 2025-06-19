function hinton(W)
    % HINTON Display a Hinton diagram for visualising a weight matrix.

    % clf; hold on;
    [rows, cols] = size(W);
    max_weight = max(abs(W(:)));
    if max_weight == 0
        max_weight = 1; % avoid division by zero
    end

    for r = 1:rows
        for c = 1:cols
            val = W(r,c);
            area = abs(val) / max_weight;
            sz = sqrt(area);
            if val > 0
                colour = [1 1 1]; % white
            elseif val < 0
                colour = [0 0 0]; % black
            else
                continue
            end
            x = c - sz/2;
            y = rows - r + 1 - sz/2; % invert y-axis for visual alignment
            rectangle('Position',[x y sz sz],'FaceColor',colour,'EdgeColor','k'); % black border
        end
    end

    axis equal;
    set(gca,'XTick',1:cols,'YTick',1:rows,'YDir','normal');
    grid on;
    hold off;
end
