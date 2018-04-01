function [ethaf,ethafb] = etharz(I0,R0,lambda1,gas,V,Pm,lp,profile,znozzle,Te,f,zmax,rmax,nres,ethai,myfilen)
%calculate the global ionization fraction after np pulses considering gas velocity and ionization recombination/diffusion



%f pulse frequency (60MHz)
tpp = 1/f;         %time between 2 pulses
d = 2.5*R0;        %interaction area

c		= 299792458;
w1		= c*2*pi/lambda1;
k1		= w1/c;																	
Zr		= k1*R0^2/2;


zmin = -zmax;
rmin = -rmax;


zi		= [zmin : (zmax-zmin)/nres : zmax];
ri		= [rmin : (rmax-rmin)/nres : rmax];

if Pm <= 1 % diffusion is really fast and everything diffuse away
        ethaf = ethai;
else



    %-------------- calculating the ionization after n pulse --------------%
    
    np = floor(d*f/V)+1;  %avr number of pulses that atoms will "see"

    % metaZ(:,:,n) = metastable atoms just before the nth pulse
    % ionGas(:,:,n) = ions just before the nth pulse
    % sumGas(:,:,n) = sum of metastable atoms and ions just before the nth pulse 
    
    ethaZ = zeros(length(ri),length(zi),np);
%     metaZ = zeros(length(ri),length(zi),np);
    sumGas = zeros(length(ri),length(zi),np);
%     ionGas = zeros(length(ri),length(zi),np);
    
    
    ethaZ(:,:,1) = ethai;
    
    

    lr = length(ri)*V/f/(rmax-rmin); %average step, have to be integern so we create lstep
    lstep = zeros(np,1);
    for n = 2:np
        display(sprintf('Calculating new ionization file for I0 = %.2e, V = %.0fms,  P = %.0fmbar, lp = %.0fum, znozzle = %.2fmm after %d pulses \nPulse num %d',I0,V,Pm,lp*1e6,znozzle*1e3,np,n))

        [ethaZbp,metaZbp] = ethabp(ethaZ(:,:,n-1),Pm,lp,profile,znozzle,Te,tpp,gas,zmax,rmax,nres);
        
        lstep(n) = floor((n-1)*lr)-sum(lstep(1:n-1));
         %for axample, if lr = 1.8, going in step of 2 by 2 is wrong, 1 by 1 too,
         %so we have to create a variable step to make sure that the average step at the end is 1.8
         
        for k = 1:length(zi)
            for l = 1:length(ri)
                if l < lstep(n) + 1
                    ethaZ(l,k,n) = 0;
%                     metaZ(l,k,n) = 0;
                else
                    sumGas(l,k,n) = metaZbp(l-lstep(n),k) + ethaZbp(l-lstep(n),k);
%                     metaZ(l,k,n) = metaZbp(l-lp,k);
                    ethaZ(l,k,n) = sumGas(l,k,n) + (1-sumGas(l,k,n))*ethaZ(l,k,1);
                end
            end
        end
%         ionGas(:,:,n) = ethaZbp;
    end
    ethafb = sumGas(:,:,end);
    ethaf = ethaZ(:,:,end);
end   
    
    
    %write the resulting ionization fraction on a file
    fid = fopen(myfilen,'wt');
    for i=1:size(ethaf,1)
        fprintf(fid, '%d ', ethaf(i,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
    
    
    myfilenb = strcat(myfilen(1:end-7),'_bp.txt');
    fid = fopen(myfilenb,'wt');
    for i=1:size(ethafb,1)
        fprintf(fid, '%d ', ethafb(i,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
    
    
    % ---------- ploting Xe+, Xe* and sum fraction (not mandatory) ----------- %
%     figure;
%     imagesc(zi.*1e3,ri.*1e6,ethaZ(:,:,np))
%     title(sprintf('ionization after n pulses for V%.0dms and P=%.0dmbar', V, Pm))
%     xlabel('z[mm]'); ylabel('r[um]')
%     caxis([0,1])
%     colormap(jet)
%     colorbar
%     saveas(gcf, fullfile('Ioniz',gas,'final',sprintf('ioniz_I0%.2e_V%.0dms_P=%.0dmbar.png', I0, V, Pm)), 'png')
%     
%     figure('visible','off');
%     imagesc(zi.*1e3,ri.*1e6,metaZ(:,:,np))
%     title('Xe* before')
%     xlabel('z[mm]'); ylabel('r[um]')
%     caxis([0,1])
%     colormap(jet)
%     colorbar
%     saveas(gcf, fullfile('Results',sprintf('meta_I0%.2e_V%dms.png', I0(i), V(j))), 'png')
%     
%     figure('visible','off');
%     imagesc(zi.*1e3,ri.*1e6,ionGas(:,:,np))
%     title('Xe+ before')
%     xlabel('z[mm]'); ylabel('r[um]')
%     caxis([0,1])
%     colormap(jet)
%     colorbar
%     saveas(gcf, fullfile('Results',sprintf('ionXe_I0%.2e_V%dms.png', I0(i), V(j))), 'png')
    
end