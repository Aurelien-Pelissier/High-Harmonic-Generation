function [P,dPz,dPr] = Press(z,r,Pm,lp,profile,znozzle)
%Pm = pressure peak (mbar)
%return Pressure and its derivative at position (z,r) (mbar)%

sigmar = 120e-6;
if profile == 'gauss'
    %function - 2D gaussian
    sigmaz = lp/(2*sqrt(2*log(2)));
    P	= Pm*exp(-(z-znozzle)^2/(2*sigmaz^2)-r^2/(2*sigmar^2));
    dPz	= -(z-znozzle)/(sigmaz^2)*Pm*exp(-(z-znozzle)^2/(2*sigmaz^2)-r^2/(2*sigmar^2));
    dPr	= -r/(sigmar^2)*Pm*exp(-(z-znozzle)^2/(2*sigmaz^2)-r^2/(2*sigmar^2));
    
elseif profile == 'squar'
    %function square * gaussian
    if abs(z-znozzle) > lp/2
        P = 0;
        dPz = 0;
        dPr = 0;
        
    else
        P	= Pm*exp(-r^2/(2*sigmar^2));
        dPz	= 0;
        dPr	= -r/(sigmar^2)*Pm*exp(-r^2/(2*sigmar^2));
    end
    
    
elseif profile == 'squ_g'
    %function square * gaussian smoothed with gaussian
    fa=20;
    sigma = lp/fa;
    if z-znozzle >= lp/2
        P = Pm*exp(-r^2/(2*sigmar^2))*exp(-(z-znozzle-lp/2)^2/sigma^2);
        dPz = (Pm*exp(-r^2/(2*sigmar^2))*exp(-(lp/2 - z + znozzle)^2/sigma^2)*(lp - 2*z + 2*znozzle))/sigma^2;
        dPr = -(Pm*r*exp(-r^2/(2*sigmar^2))*exp(-(lp/2 - z + znozzle)^2/sigma^2))/sigmar^2;
        
    elseif z-znozzle <= -lp/2
        P = Pm*exp(-r^2/(2*sigmar^2))*exp(-(z-znozzle+lp/2)^2/sigma^2);
        dPz = -(Pm*exp(-r^2/(2*sigmar^2))*exp(-(lp/2 + z - znozzle)^2/sigma^2)*(lp + 2*z - 2*znozzle))/sigma^2;
        dPr = -(Pm*r*exp(-r^2/(2*sigmar^2))*exp(-(lp/2 + z - znozzle)^2/sigma^2))/sigmar^2;
        
    else
        P	= Pm*exp(-r^2/(2*sigmar^2));
        dPz	= 0;
        dPr	= -r/(sigmar^2)*Pm*exp(-r^2/(2*sigmar^2));
    end
    
    
else
    P=nan;
    dPz = nan;
    dPr = nan;
    display('error, profile should be gauss or square')
end

end

