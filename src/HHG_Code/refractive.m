function [refrac] = refractive(lambda,gas)
%Calculate the refractive index of the given gas at 1atm using the Sellmeier equation

if gas == 'Xe'
    B1 = 0.00322869;
    B2 = 0.00355393;
    B3 = 0.0606764;
    C1 = 46.301;
    C2 = 50.578;
    C3 = 112.74;
    %http://refractiveindex.info/?shelf=main&book=Xe&page=Bideau-Mehu
    
elseif gas == 'Kr'
    B1 = 0.00253637;
    B2 = 0.00273649;
    B3 = 0.0620802;
    C1 = 65.4742;
    C2 = 73.698;
    C3 = 181.08;
    %http://refractiveindex.info/?shelf=main&book=Kr&page=Bideau-Mehu
    
elseif gas == 'Ar'
    B1 = 2.5041e-3;
    B2 = 5.00283e-4;
    B3 = 5.22343e-2;
    C1 = 91.012;
    C2 = 87.892;
    C3 = 214.02;
    %http://refractiveindex.info/?shelf=main&book=Ar&page=Bideau-Mehu
    
elseif gas == 'He'
    B1 = 0.01470091;
    B2 = 0;
    B3 = 0;
    C1 = 423.98;
    C2 = 0;
    C3 = 0;
    %https://refractiveindex.info/?shelf=main&book=Kr&page=Bideau-Mehu
    
else
    sprintf('error, gas should be Xe, Kr or Ar')
end

lambdau = lambda*1e6; %Converting to um


refrac = B1/(C1 - lambdau^(-2)) + B2/(C2 - lambdau^(-2)) + B3/(C3 - lambdau^(-2));
end

