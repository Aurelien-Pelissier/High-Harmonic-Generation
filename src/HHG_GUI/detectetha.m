function ethaf = detectetha(I0,tp,R0,lambda1,gas,zmax,rmax,nres,t)
%detecting ionization file to avoid doing the calculation everytime%
%if the file doesnt exist or the resolution is too low, a new one is created%

myfolder = fullfile(sprintf('%.0fnm', lambda1*1e9),'Ioniz',sprintf('%.1fum', R0*1e6),gas,'first_pulse');
if exist (myfolder,'dir')  == 0 %checking if the folder exist
    mkdir(myfolder) %create the folder if doesnt exist
end


myfile = fullfile(sprintf('%.0fnm', lambda1*1e9),'Ioniz',sprintf('%.1fum', R0*1e6),gas,'first_pulse',sprintf('ioniz_I0%.2e_tp%.0ffs_zmax%.1fmm_%.0ffs.txt', I0, tp*1e15, zmax*1e3, t*1e15));
if exist(myfile,'file') ~= 0 %checking if the file exist
    ethafi = dlmread(myfile); %reading the ionization matrix
    if size(ethafi,1)-1 >= nres  %checking if the file resolution is higher than nres
        ethaf = imresize(ethafi, [nres+1 nres+1]); %resizing the ionization matrix
        display(sprintf('using ioniz_I0%.2e_t=%.0ffs.txt with resolution %d', I0, t*1e15, size(ethafi,1)-1))
    else
        ethaf = ethaprz(I0,tp,R0,lambda1,gas,zmax,rmax,nres,myfile,t); %creating a new one
    end
else
    ethaf = ethaprz(I0,tp,R0,lambda1,gas,zmax,rmax,nres,myfile,t); %creating a new one
end
end

function etha = ethaprz(I0,tp,R0,lambda1,gas,zmax,rmax,nres,myfile,t)
%calculate and create a new ionization file%
    c		= 299792458;
    w1		= c*2*pi/lambda1;
    k1		= w1/c;																	
    Zr		= k1*R0^2/2;

    zmin = -zmax;
    rmin = -rmax;
    zi		= [zmin : (zmax-zmin)/nres : zmax];
    ri		= [rmin : (rmax-rmin)/nres : rmax];
    etha = zeros(length(ri),length(zi));
    
    display(sprintf('calculating first pulse ionization file for I0=%.2e_t=%.0ffs, resolution = %d', I0,t*1e15,nres))

    for     l = 1:length(zi)*length(ri)
        i = floor(l/(length(zi)+0.0000000000001))+1;
        j = mod(l,length(ri))+1;
        
        Rzz	= R0*sqrt(1+(zi(i)/Zr)^2);
        Irz   = I0*(R0/Rzz)^2*exp((-2*ri(j)^2)/(Rzz^2)); 
        etha(j,i) = ethap(Irz,tp,lambda1,gas,t);
    end

    %write the resulting ionization fraction on a file
    fid = fopen(myfile,'wt');
    for i=1:size(etha,1)
        fprintf(fid, '%d ', etha(i,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
    
end