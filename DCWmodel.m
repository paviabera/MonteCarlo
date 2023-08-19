function [sse, fittedcurve] =  DCWmodel(params, rho, ydata)
    % This function computes the homogeneous semi-infinite solution to
    % Continuaous Wave PDE
    % with instrument response function
    % fittedcurve     -- Estimated Fluence rate
    % rho             -- source detector separation 
   
    
    mu_eff = params(1);

    % define constants
   rho_o=rho(1);
    % Continuous Wave solution
    phi_cw= (rho_o*(exp(-mu_eff*rho)))./(rho*exp(-rho_o*mu_eff));

    %Residuals
    fittedcurve = phi_cw;
    sse = sum((fittedcurve-ydata).^2);
    
end