function [ethaf,ethafb,ethaft] = detectethan(I0,tp,V,P,lp,profile,znozzle,Te,f,R0,lambda1,gas,zmax,rmax,nres,ethai,t)
%detecting ionization file to avoid doing the calculation everytime%
%if the file doesnt exist or the resolution is too low, a new one is created%

%return:
%ethaft = the ionization at time t of the last pulse, array (r,z)
%ethafb = the ionization just before the last pulse, array (r,z)
%ethaf = the ionization at the end of the last pulse, array (r,z)

myfolder = fullfile(sprintf('%.0fnm', lambda1*1e9),'Ioniz',sprintf('%.1fum', R0*1e6),gas,'final',sprintf('%.0fum', lp*1e6),profile);
if exist (myfolder,'dir')  == 0 %checking if the folder exist
    mkdir(myfolder) %create the folder if doesnt exist
end


% ethafb = zeros(nres+1,nres+1); %zero atm but have to do a file creation for this
% ethafb1 = zeros(nres+1,nres+1);


myfilen = fullfile(sprintf('%.0fnm', lambda1*1e9),'Ioniz',sprintf('%.1fum', R0*1e6),gas,'final',sprintf('%.0fum', lp*1e6),profile,sprintf('ioniz_I0%.2e_tp%.0ffs_V%.0fms_P%.0fmbar_znoz%.2fmm_zmax%.1fmm_ap.txt', I0, tp*1e15, V, P, abs(znozzle)*1e3, zmax*1e3));
myfilenb = strcat(myfilen(1:end-7),'_bp.txt');
if exist(myfilen,'file') ~= 0 %checking if the file exist
    ethafi = dlmread(myfilen); %reading the ionization matrix
    ethafbi = dlmread(myfilenb);
    if size(ethafi,1)-1 >= nres  %checking if the file resolution is higher than nres
        ethaf1 = imresize(ethafi, [nres+1 nres+1]); %resizing the ionization matrix
        ethafb1 = imresize(ethafbi, [nres+1 nres+1]);
        display(sprintf('using ioniz_I0%.2e_V%.0fms_P%.0fmbar_lp%.0fum_znoz%.2fmm.txt with resolution %d', I0, V, P, lp*1e6, znozzle*1e3, size(ethafi,1)-1))
    else
        [ethaf1,ethafb1] = etharz(I0,R0,lambda1,gas,V,P,lp,profile,abs(znozzle),Te,f,zmax,rmax,nres,ethai,myfilen); %creating a new one
    end   
else
    [ethaf1,ethafb1] = etharz(I0,R0,lambda1,gas,V,P,lp,profile,abs(znozzle),Te,f,zmax,rmax,nres,ethai,myfilen); %creating a new one
end


if znozzle < 0 %  taking the symetry to avoid a second calculation
    ethaf = zeros(size(ethaf1,1),size(ethaf1,2));
    ethafb = zeros(size(ethafb1,1),size(ethafb1,2));
    for i = 1: size(ethaf,1) % ethaf and ethafb has the same size
        ethaf(:,i) = ethaf1(:,size(ethaf1,2)-i+1);
        ethafb(:,i) = ethafb1(:,size(ethafb1,2)-i+1);
    end
else
    ethaf = ethaf1;
    ethafb = ethafb1;
end



%%calculating ionization at time t
if t <= -tp
    myfilent = myfilenb;
elseif t >= tp
    myfilent = myfilen;    
else
    myfilent = strcat(myfilen(1:end-7),sprintf('_%.0ffs.txt',t*1e15));
end


if exist(myfilent,'file') ~= 0 %checking if the file exist
    ethafti = dlmread(myfilent); %reading the ionization matrix
    if size(ethafti,1)-1 >= nres  %checking if the file resolution is higher than nres
        ethaft = imresize(ethafti, [nres+1 nres+1]); %resizing the ionization matrix
        display(sprintf('using ioniz_t=%.0ffs', t*1e15))
    else
        ethaft = etharzt(I0,tp,R0,lambda1,gas,zmax,rmax,nres,ethafb,t); %creating a new one
    end   
else
    ethaft = etharzt(I0,tp,R0,lambda1,gas,zmax,rmax,nres,ethafb,t); %creating a new one
end

end

function ethaft = etharzt(I0,tp,R0,lambda1,gas,zmax,rmax,nres,ethafb,t)
%from an ionitial ionization fraction (ethabp), calculate the ionization at time t (in the pulse).

ethapt = detectetha(I0,tp,R0,lambda1,gas,zmax,rmax,nres,t);
ethaft = zeros(size(ethafb,1),size(ethafb,2));

for k = 1:size(ethafb,2)
    for l = 1:size(ethafb,1)
        ethaft(l,k) = ethafb(l,k) + (1-ethafb(l,k))*ethapt(l,k);
    end
end

end
