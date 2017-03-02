function qp_serum = qp_serum(mu_set)
    
    if mu_set == 0
        qp_serum = 0;
    else
        A = 0.182518;
        B = -0.039406;
        C = 0.003172;
        D = -0.000027;
        qp_serum = A*mu_set^3 + B*mu_set^2 + C*mu_set + D;
    end
end