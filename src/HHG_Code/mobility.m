function [mu0] = mobility(gas)
%return the mobility of ions in its parent gas (m2/Vs)
if gas == 'Xe'
    mu0 = 6e-5; %TJ + other sources (see report)
elseif gas == 'Kr'
    mu0 = 1e-4; %see report for source
elseif gas == 'Ar'
    mu0 = 1.4e-4; %see report for source
else
    mu0 =0;
    error = 'gas should be Xe, Kr or Ar'
end


end

