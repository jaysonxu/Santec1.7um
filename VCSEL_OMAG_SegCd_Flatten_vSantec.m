%%% OMAG Segmentation + projection V10012016 S.Song
clear;clc
close all;

%%

    tmp=dir('.\*StruC0.dcm');
    for Nfile=1:size(tmp,1), fprintf(['[%d]\t',tmp(Nfile).name(1:end-5),'\n'],Nfile);end
    tic
    %%
    ListSeg={...
        'Flow_edC1'...
%         ,'Flow_omC1'...
        }; 

    do_loadSeg=1;
    do_showSeg=0;
%     do_test_seg=1; if do_test_seg, do_loadSeg=0;end  % Test parameters for segmenting IMG

    do_Blend_Flatten=0;
    do_Blend=0;

    use_autoRg=1;RgFlow=[0.1 0.8];

    num_vol=5;
    nRB=1;
    
    Info.nZshift=0;
    Info.nZseg=350;
    BandLines=[35,45,67,110,205];
    
    Thrsh_ratio=[0.6 0.55];
    
    %%

    StrRg=[0 65535];

    distcomp.feature( 'LocalUseMpiexec', false ); % parallel processing


%     for Nfile=[1]   
    for Nfile=1%2:size(tmp,1)

        filename=tmp(Nfile).name;


        STR=double(squeeze(dicomread(filename)));
        nX=size(STR,2);nY=size(STR,3)/nRB;

        if nRB>1, STR=squeeze(mean(reshape(STR,[],nX,nRB,nY),3));end

        STR(1:50,:,:)=0;
        %%
        if do_loadSeg&&exist([filename,'_SegInfo.mat'],'file')
            load([filename,'_SegInfo.mat']);
        else
            test_iY=400;

            Info.WinSize=70;
            Info.thrsh_ratio=Thrsh_ratio(Nfile);
            Info.medfilt=[50 50];
            Info.gaussfilt=[7 9];

            seg_top=zeros(nX,nY);
            seg_f=seg_top;


                    pb = ProgressBar(nY);
                    
                    for iY=1:nY
%                     parfor iY=1:nY
                        figure(1);
    %                     cimgc=medfilt2(STR(:,:,iY),[Info.WinSize 13]);
                        cimgc=imgaussfilt(STR(:,:,iY),Info.gaussfilt);
                        cimgMax=max(cimgc);
                        for iX=1:nX
                            vol=cimgc(:,iX);
    %                         [~,loc] = findpeaks(vol,'MINPEAKHEIGHT', max(vol)-Info.cut_thrsh);
%                             loc = find(vol>Info.cut_thrsh,1,'first');
                            loc = find(vol>cimgMax(iX)*Info.thrsh_ratio,1,'first');
                            while isempty(loc)
                                cimgMax(iX)=cimgMax(iX)*Info.thrsh_ratio;
                                loc = find(vol>cimgMax(iX)*Info.thrsh_ratio,1,'first');
                            end
                            seg_top(iX,iY)=loc(1);
                        end
                        imagesc(STR(:,:,iY)); colormap gray;
                        hold on,
                        plot(seg_top(:,iY),'-.r');
                        seg_f(:,iY)=medfilt1(seg_top(:,iY),Info.medfilt(1));
                        plot(seg_f(:,iY),'-.g');
                        hold off;
                        drawnow;
                        pb.progress;
                    end
                    pb.stop;
                    fprintf('\n')
                    imagesc(seg_f)
                    save([filename,'_SegInfo.mat'],'seg_top','seg_f','Info')
%                     close all;
                    %% smooth segmentation result

                    end

                    
                
%         end
        %%    
        
         [O,P,Q]=ndgrid(1:Info.nZseg,1:nX,1:nY);
            O=bsxfun(@plus,O,permute(round(seg_f),[3,1,2]));
            O(O>size(STR,1))=size(STR,1);
            O(O<1)=1;
            ind=sub2ind(size(STR),O,P,Q);
            clear O P Q
            
            STR1=uint16(STR(ind));
            STR2=permute(STR1,[1,2,4,3]);
            if ~exist('.\Flatten\','dir'),mkdir('.\Flatten\'),end
            dicomwrite(STR2,['.\Flatten\',filename(1:end-4),'.dcm']);
            %% Projection
    end
    %%
    
    
