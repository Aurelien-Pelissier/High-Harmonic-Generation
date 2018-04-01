function Ip = Ipot(gas)
%return the first ionization potentiel of the gas in eV
    %https://dept.astro.lsa.umich.edu/~cowley/ionen.htm

if numel(gas) == 2
    if gas == 'Xe'
        Ip = 12.1298;
    elseif gas == 'Kr'
        Ip = 13.9996;
    elseif gas == 'Ar'
        Ip = 15.7596;
    end

    
elseif numel(gas) == 6

    if gas == 'Xestar'
        Ip = 3;
    elseif gas == 'Krstar'
        Ip = 3;
    elseif gas == 'Arstar'
        Ip = 3;
    end
    
else
    sprintf('error, gas should be Xe, Kr or Ar')
end


end

