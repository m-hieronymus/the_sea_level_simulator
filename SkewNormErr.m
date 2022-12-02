function errval = SkewNormErr85(parms)
    % This is the utility function that fminsearch will use to optimize parameters.
    % parms is a vector of the three parameters of the SkewNormal distribution.
    % this function will compute an error score summarizing the deviation
    % of the values of that distribution from the known Pct05, mean, and Pct95
    
    % Example of known target values.
    % You want to search for Skewnormal parameter values
    % giving a distribution with these properties.
    % Adjust these to whatever values you want.

    global target
   
    
    skn = SkewNor(parms(1),parms(2),parms(3));
    
        errval = ( target(1)-skn.InverseCDF(.05) )^2 + ( target(2)-skn.InverseCDF(.17) )^2 + ( target(3)-skn.InverseCDF(.50) )^2 + ( target(4)-skn.InverseCDF(.83) )^2 + ( target(5)-skn.InverseCDF(.95) )^2;
    
end
