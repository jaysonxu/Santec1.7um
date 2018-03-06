%%% OMAG for VSCEL, with phase compensation V10012016 S.Song
clear;clc;close all;

if ~exist('.\outputStru\','dir'),mkdir('.\outputStru\'),end
outputdir=[cd,'\outputStru\'];

% outputdir='G:\exp\VCSEL\Brain\OMAG\';
Thr=5e4; % Threshold to remove low OCT signal
do_PhComp=0;
useref=0;
use_nR1=0;

saveDicom=1;
show_img=1;
do_medianshift=0;
%IMcropRg=51:1050;  nZcrop=numel(IMcropRg); imshowrgZ=1:nZcrop;
distcomp.feature( 'LocalUseMpiexec', false ); % parallel processing

tmp=dir('.\*.oct');
for iFile=1:size(tmp,1)
fprintf(['[%d]\t',tmp(iFile).name(1:end-4),'\n'],iFile);end
%%
tic
RgStr=[65 160];
IMcropRg=1:550;
nZimg=numel(IMcropRg);
for iFile=1:size(tmp,1)
% for iFile=[10:12]
    filename=tmp(iFile).name;
    fid=fopen(filename);
    bob=fread(fid,1,'uint32');
    SPL=fread(fid,1,'double');
    nX=fread(fid,1,'uint32');
    nY=fread(fid,1,'uint32');
    Boffset=fread(fid,1,'uint32');
    Blength=fread(fid,1,'uint32')+1;
    Xcenter=fread(fid,1,'double');
    Xspan=fread(fid,1,'double');
    Ycenter=fread(fid,1,'double');
    Yspan=fread(fid,1,'double');
    nR=fread(fid,1,'uint32');
    if use_nR1,nR=1;end
    n_dataset=fread(fid,1,'uint32');
    ProtMode=fread(fid,1,'uint32');
    fseek(fid,4,'cof');
    sizeBck=fread(fid,1,'uint32');
    Bck=fread(fid,sizeBck,'int16');
    sizeKES=fread(fid,2,'uint32');
    KES=(fread(fid,sizeKES(2),'double'))'*sizeKES(2);
    
    
    IMGheight=floor(Blength/2);
    
    nY=nY/nR;
    
    K=1:Blength;

    use_autoRg=0;
    winG=tukeywin(Blength,0.25);
    %%
    if useref==0%use the first 50k a-lines to calc ref
        fseek(fid,bob,'bof'); %fseek(fid,4,'cof');
        RefBs=fread(fid,[SPL*nX+2 min(50,nY)],'int16');
        RefBs(1:2,:)=[];
        Ref = repmat(mean(reshape(RefBs,SPL,[]),2),1,nX);
    else
        Ref = repmat(Bck,[1,nX]);%Ref = repmat(importdata(refname)',1,linenum);
    end
    Ref = Ref(Boffset+1:Boffset+Blength,:);
    %%    
    if saveDicom
        STR=zeros(nZimg,nX,1,nY);
%         BLD_ed=STR;
    end    

    if nX>1200||IMGheight>2048,warning('off', 'Images:initSize:adjustingMag');end
    %%
    if show_img
        fig_im=figure('position',[10 400 800 800]);

    end
%%
    pb = ProgressBar(nY);
    parfor If=1:nY

%     for If=1:nY
%     for If=147
        Bs=zeros(Blength,nX,nR);
        fid=fopen(filename);
        fseek(fid,bob+(SPL*nX+2)*nR*(If-1)*2,'bof');
        for Ic=1:nR
            fseek(fid,4,'cof');
            B=double(reshape(fread(fid,SPL*nX,'int16'),[SPL,nX]));
            Bc=B(Boffset+1:Boffset+Blength,:);
            Bs(:,:,Ic)=Bc-Ref;
%             Bs(:,:,Ic)=Bc;
%             Bs(:,:,Ic)=bsxfun(@minus,Bc,mean(Bc,2));
        end
        Bs=Bs.*winG; %apply window
        Bd=interp1(KES,Bs,K,'linenumar',0);  % Interpolation + background subtraction
        Bimg=fft(Bd,SPL,1);  %complex value
        IMG=Bimg(IMcropRg,:,:);
        if do_PhComp
            [IMG] = SSOCT_F_PhCompV3(IMG,IMcropRg(1),Thr,do_medianshift);
        end

        %%
        Im=20*log10(abs(IMG));
        Im=mean(Im,3);
        if show_img            
%             RgStrS=[80 130];
            figure(fig_im)
            imagesc(Im,RgStr);
            colormap gray;
            drawnow
        end
        %%
        STR(:,:,1,If)=Im;
%         STRd(:,:,1,If)=Imd;
%         BLD_ed(:,:,If)=20*log10(abs(Im-circshift(Im,[-319 -374])));
        fclose(fid);
        pb.progress;
    end
    pb.stop;
%%
    if saveDicom
        fprintf('calculating range....'); 
        if use_autoRg
            nP=numel(STR);
            RgStr  =  prctile(STR(1:5:nP),[2.5 99.99]);
%             RgStr_d  =  prctile(STRd(1:5:nP),[2.5 99.99]);
%             RgFlow_ed=prctile(BLD_ed(1:37:nP),[2.5 99.99]);
        end
        %%
        fprintf('saving DICOM....\n'); 
        
        STRimg=uint16(65536*(STR-RgStr(1))/(RgStr(2)-RgStr(1)));
     
        
        if do_PhComp
            dicomwrite(STRimg,[outputdir,filename,'-StruC1.dcm']);
%             dicomwrite(BLDimg_ed,[outputdir,filename,'-Flow_edC1.dcm']);
        else
            
            dicomwrite(STRimg,[outputdir,filename,'-StruC0.dcm']);
%             dicomwrite(STRimg_d,[outputdir,filename,'-StruC0_d.dcm']);
%             dicomwrite(BLDimg_ed,[outputdir,filename,'-Flow_edC0.dcm']);
        end
    end
end
toc