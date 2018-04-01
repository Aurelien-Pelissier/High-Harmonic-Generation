function X0 = potde(gas)
%return the potential depth of the gas

if gas == 'Xe'
    X0 = 1.625; % potential depth Xe

elseif gas== 'Kr'
    X0 = 1.365; % potential depth Kr
       
elseif gas == 'Ar'
    X0 = 1.2; % potential depth Ar
    
else
    X0 = 0;
    sprintf('error, gas should be Xe, Kr, or Ar')
end


end

