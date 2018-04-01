function V = velocity(gas)
%return the sound velocity of the gas (m/s), calculated on the report
if gas == 'Xe'
    V = 153;
elseif gas == 'Kr'
    V = 194;
elseif gas == 'Ar'
    V = 280; 
elseif gas == 'He'
    V = 880;
else
    V =0;
    error = 'gas should be Xe, Kr or Ar'
end


end