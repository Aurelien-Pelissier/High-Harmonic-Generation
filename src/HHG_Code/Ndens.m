function Natm = Ndens(gas)
%return the number density at 1atm of the gas
%same for all gas

Natm = 2.7e25;

% Na = 6.022e23; %Avogadro number (mol-1)
% 
% if gas == 'Xe'
%     rho = 5.894; %mass density at stp (kg/m3)
%     M = 131.293; %Molar mass (g/mol)
%     
% elseif gas == 'Kr'
%     rho = 3.749; %mass density at stp (kg/m3)
%     M = 83.798; %Molar mass (g/mol)
%     
% else
%     Natm = 2.7e25;
% end
% 
% Natm = Na/M*rho*1e3;
% %about 2.7e25, nearly the same for both gas.

end

