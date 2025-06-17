function r2 = compute_r2(true_vals, predicted_vals)
SS_res = sum((true_vals - predicted_vals).^2, 1);
    SS_tot = sum((true_vals - mean(true_vals)).^2, 1);
    r2 = 1 - SS_res ./ SS_tot;
end