%% Dynamical Structure Factor S_XX_wk()
close all
clear
clc

%% Create the Hilbert space and op representations
tic
o1 = makespace_confinedLimitp(3,2,508);
toc

%% Graph - 1 : energy spectrum as a fnc of theta
close
figure(1)
for i = 1
    LTh = 21;
    thetaVector = linspace(0,0.5,LTh);
    eigen_tt = zeros(length(o1.idS),LTh);
    f = figure;
    for i = 1:LTh
        [eigen_t ,eigVec_t] = o1.mbspectrum(thetaVector(i)*pi);
        %         [eigen_tt,sInd]  = sort(real(eigen_t(:).'));
        E0(i) = abs(eigen_tt(2)-eigen_tt(1));

        %         klm = [];
        %         for ii = 1:prod(size(eigen_t))
        %             klm = [klm;sort(real((eigen_t{ii})))];
        %         end
        %         eigen_tt(:,i)  = sort(klm);
        for ii = 1:length(o1.ksec.kvec)
            ei = sort(real(eigen_t{ii}));
            scatter(ii+0*ei(1:100),ei(1:100)-min(eigen_t{6}),38,'b','filled')
            hold on
        end
        hold off
        pause(1)
    end
    %%
    f = figure;
    plot(thetaVector,eigen_tt(1:8,:)-repmat(eigen_tt(1,:),8,1),'k','LineWidth',3)
    xlabel('$\theta/\pi$','Interpreter','latex','FontSize',25)
    ylabel('$E(\theta)$','Interpreter','latex','FontSize',25)
    ylim([0,4.1])
    %
    %     saveas(f,"E_vs_theta.eps")
    %     saveas(f,"E_vs_theta.jpg")

    %Graph - 2 : band-structure in k-space

    close
    klmn = (reshape(cell2mat(o1.ksec.kvec),[2,length(cell2mat(o1.ksec.kvec))/2]));
    subplot(1,2,1)
    scatter([klmn(1,:),klmn(1,:),klmn(1,:)+2,klmn(1,:)+2],[klmn(2,:),klmn(2,:)+2,klmn(2,:),klmn(2,:)+2],60,'filled')
    xlabel('$k_x/\pi$','Interpreter','latex','FontSize',25)
    ylabel('$k_y/\pi$','Interpreter','latex','FontSize',25)

    subplot(1,2,2)
    kIndc = [6,3,12,4,2,6,1,6];
    kIndc = [6,3,12,9,6,7,6,2,4,9,6];
    [eigen_t ,eigVec_t] = o1.mbspectrum(0.4*pi);
    for i = 1:length(kIndc)
        scatter(0*real(eigen_t{kIndc(i)})+i,real(eigen_t{kIndc(i)})-min(real(eigen_t{6})),60,'filled')
        hold on
    end
    %     pause
end
%
% savefig("E_kxky.fig")

%% Graph - 2: fixed (kx,ky) and (E vs theta)
% klmn = (reshape(cell2mat(o1.ksec.kvec),[2,length(cell2mat(o1.ksec.kvec))/2]))
% close
figure(2)

for ikl=1
    % for a 3 x 2 plaquette lattice momentum values [0,0] and [1,1]*pi
    kvecm = cell2mat(o1.ksec.kvec.').';
    kIndc(1) = find(kvecm(1,:)==0 & kvecm(2,:)==0);
    kIndc(2) = find(kvecm(1,:)==1 & kvecm(2,:)==1);
    
    LTh = 11; 
    thetaVector = linspace(0,.5,LTh);
    
    eigen_0 = zeros(length(o1.ksec.proj{kIndc(1)}),LTh);
    eigen_pi = zeros(length(o1.ksec.proj{kIndc(2)}),LTh);

    Xk_weight6  = zeros(length(o1.ksec.proj{kIndc(1) }),LTh);
    Xk_weight10 = zeros(length(o1.ksec.proj{kIndc(2)}),LTh);
    Z4k_weight6 = zeros(length(o1.ksec.proj{kIndc(1) }),LTh);
    Z4k_weight10= zeros(length(o1.ksec.proj{kIndc(2)}),LTh);

    for i = 1:LTh
        [eigen_t ,eigVec_t] = o1.mbspectrum(thetaVector(i)*pi);
        [XkME ,kk] = o1.S_XX_kkp(eigen_t,eigVec_t);
        [Z4kME, ~] = o1.S_Z4Z4_kkp(eigen_t,eigVec_t);
        eigen_tt  = real(cell2mat(eigen_t.'));

        eigen_0(:,i)  = sort(real(eigen_t{kIndc(1) })-min(real(eigen_t{kIndc(1) })));
        eigen_pi(:,i) = sort(real(eigen_t{kIndc(2)})-min(real(eigen_t{kIndc(1)})));

        Xk_weight6(:,i)  = XkME{kIndc(1)} ; Xk_weight10(:,i)  = XkME{kIndc(2)} ;
        Z4k_weight6(:,i) = Z4kME{kIndc(1)}; Z4k_weight10(:,i) = Z4kME{kIndc(2)};

        for ink = 1:length(o1.ksec.kvec)
            abcd1(i,ink)=sum(abs(XkME{ink}'));  
            abcd2(i,ink)=sum(abs(Z4kME{ink}'));
        end
    end
%     mmX(ikl,:)=[max(abs(Xk_weight10(:))),max(abs(Xk_weight6(:)))]
%     mmU(ikl,:)=[max(abs(Z4k_weight10(:))),max(abs(Z4k_weight6(:)))]    %
%
    %%
    SSize = 180;

    e0S = size(eigen_0); e0S = e0S(1);
    epiS = size(eigen_pi); epiS = epiS(1);

    mmax= max([max(abs(Xk_weight10(:))),max(abs(Xk_weight6(:)))]);
    mmaU= max([max(abs(Z4k_weight10(:))),max(abs(Z4k_weight6(:)))]);

    for den = 1
        subplot(2,2,1)
        for ii = 1:e0S
            scatter(thetaVector,eigen_0(ii,:), SSize*abs(Xk_weight6(ii,:))/mmax+1e-2,...
                'b','filled','LineWidth',3)
            hold on
            plot(thetaVector,eigen_0(ii,:),'k:','LineWidth',0.1)
            ax = gca;
            ax.FontSize = 20; 
            xlabel('$\theta/\pi$','Interpreter','latex','FontSize',25)
            ylabel('$E_{\lambda}(\theta;\vec{k}=(0,0))$','Interpreter','latex','FontSize',25)
            ylim([0,12]);xlim([0,0.5]);
        end
        hold off
        title('$\langle \vec{k} = (0,0),\lambda| X_{\vec{k} = (0,0)} | \Psi_{GS}(\theta) \rangle$','FontSize',20,'Interpreter','latex')
        subplot(2,2,2)
        for ii = 1:epiS
            scatter(thetaVector,eigen_pi(ii,:), SSize*abs(Xk_weight10(ii,:))/mmax+1e-2,'b','filled','LineWidth',3)
            hold on
            plot(thetaVector,eigen_pi(ii,:),'k:','LineWidth',0.1)
            ax = gca;
            ax.FontSize = 20; 
            xlabel('$\theta/\pi$','Interpreter','latex','FontSize',25)
            ylabel('$E_{\lambda}(\theta;\vec{k}=(1,1)\pi)$','Interpreter','latex','FontSize',25)
            ylim([0,12]);xlim([0,0.5]);
        end
        hold off
        title('$\langle \vec{k} = (1,1)\pi,\lambda| X_{\vec{k} = (1,1)\pi} | \Psi_{GS}(\theta)\rangle$','FontSize',20,'Interpreter','latex')
        subplot(2,2,3)
        for ii = 1:e0S
            scatter(thetaVector,eigen_0(ii,:), SSize*abs(Z4k_weight6(ii,:))/mmaU+1e-2,'b','filled','LineWidth',3)
            hold on
            plot(thetaVector,eigen_0(ii,:),'k:','LineWidth',0.1)
            ax = gca;
            ax.FontSize = 20; 
            xlabel('$\theta/\pi$','Interpreter','latex','FontSize',25)
            ylabel('$E_{\lambda}(\theta;\vec{k}=(0,0)\pi)$','Interpreter','latex','FontSize',25)
            ylim([0,12]);xlim([0,0.5]);
        end
        hold off
        title('$\langle \vec{k} = (0,0),\lambda| U_{\vec{k} = (0,0)} | \Psi_{GS}(\theta) \rangle$','FontSize',20,'Interpreter','latex')
        subplot(2,2,4)
        for ii = 1:epiS
            scatter(thetaVector,eigen_pi(ii,:), SSize*abs(Z4k_weight10(ii,:))/mmaU+1e-2,'b','filled','LineWidth',3)
            hold on
            plot(thetaVector,eigen_pi(ii,:),'k:','LineWidth',0.1)
            ax = gca;
            ax.FontSize = 20; 
            xlabel('$\theta/\pi$','Interpreter','latex','FontSize',25)
            ylabel('$E_{\lambda}(\theta;\vec{k}=(1,1)\pi)$','Interpreter','latex','FontSize',25)
            ylim([0,12]);xlim([0,0.5]);
        end
        hold off
        title('$\langle \vec{k} = (1,1)\pi,\lambda| U_{\vec{k}= (1,1)\pi} | \Psi_{GS}(\theta)\rangle$','FontSize',20,'Interpreter','latex')
    end
    %%
end
%
% savefig("E_vs_theta_fixed_kxky.fig")

%% Dynamical correlations
% Graph - 3.1:  S_XX(k,w)
figure(3)
for i=1
    str = 'SXX_SUU_Lattice3x2.fig';
    Lm = 30;
    omega = linspace(0,Lm,Lm*1001);
    %
    kvecm = cell2mat(o1.ksec.kvec.').';
    [kxx,om] = meshgrid(kvecm(1,:),omega);
    [kyy, ~] = meshgrid(kvecm(2,:),omega);

    de = 1e-4; % Lorentzian broadening infinitesimal term
    SXX = om*0;
    Sz4z4 = SXX;
    thetaVector = [0.3,0.5];
    % S_XX for 4-theta values
    for i = 1:length(thetaVector)
        [eigen_t ,eigVec_t ] = o1.mbspectrum(thetaVector(i)*pi);

        [XkME ,kk] = o1.S_XX_kkp(eigen_t,eigVec_t);
        eigen_tt  = real(cell2mat(eigen_t.'));

        SXX = om*0;
        kken = 0*real(eigen_t{1}).';
        for i1 = 1:length(eigen_t)
            XkMEi1  = XkME{i1};
            eigen_t2 = real(eigen_t{i1}).';
            for i2 = 1:length(eigen_t2)
                SXX( :,i1) = SXX( :,i1) + abs(XkMEi1( i2))^2./(om(:,i1)-(eigen_t2(i2)-min(real(eigen_t{6})))+1j*de);
            end
        end
        logSXX = log(-(imag(SXX)));
        logSXX(logSXX<-10) = -10;
        subplot(2,2,i+2)
        hold on
        kvecm=reshape(cell2mat(o1.ksec.kvec),[2,length(o1.ksec.kvec)]);
        [kvecm(1,:),indkvec]=sort(kvecm(1,:));
        kvecm(2,:) = kvecm(2,indkvec);
        logSXX = logSXX(:,indkvec);
        kvecm(1,:) = kvecm(1,:)-(kvecm(1,2)+kvecm(1,1))/2;
        
        [kxx,om] = meshgrid(kvecm(1,:),omega);
        [kyy, ~] = meshgrid(kvecm(2,:),omega);

        logSSXX=logSXX; logSSXX(:,end+1) = logSSXX(:,1);   % k-space periodicity, stitched to the end 
        kxxx = kxx; kxxx(:,end+1) = 2+0*kxxx(:,1); % k-space periodicity, stitched to the end 
        omm = om; omm(:,end+1) = om(:,1);          % k-space periodicity, stitched to the end

        mesh(kxxx,omm,logSSXX)
        ax = gca;
        ax.FontSize = 20;
        view(2)
        title(strcat('$S_{XX}(\vec k,\omega,\theta/\pi:',num2str(thetaVector(i)),')$'),'Interpreter','latex','FontSize',25)

        ylabel('$\omega$','Interpreter','latex','FontSize',30);
        xlabel('$k_x/\pi$','Interpreter','latex','FontSize',30);
        xlim([0-(kvecm(1,2)-kvecm(1,1))/2,2-(kvecm(1,2)-kvecm(1,1))/2])
        ylim([0,12])
        %         colorbar('location','north')

        aa{i} = logSSXX;
    end
    % Plotting the first plot
    bottom = min(min(min(aa{1})),min(min(aa{2})));
    top  = max(max(max(aa{1})),max(max(aa{2})));
end
% Graph - 3.2:  S_UU(k,w)
for ii9 = 1
    for i = 1:length(thetaVector)
        [eigen_t ,eigVec_t ] = o1.mbspectrum(thetaVector(i)*pi);
        [Z4kME, ~] = o1.S_Z4Z4_kkp(eigen_t,eigVec_t);
        Sz4z4 = 0*SXX;
        for i1 = 1:length(eigen_t)
            Z4kMEi1 = Z4kME{i1};
            eigen_t2 = real(eigen_t{i1}).';
            for i2 = 1:length(eigen_t2)
                Sz4z4(:,i1) = Sz4z4(:,i1) + abs(Z4kMEi1(i2))^2./(om(:,i1)-(eigen_t2(i2)-min(real(eigen_t{6})))+1j*de);
            end
        end

        logSz4z4 = log(-(imag(Sz4z4)));
        logSz4z4(logSz4z4<-10) = -10;

        logSz4z4 = logSz4z4(:,indkvec);
        logSzz4zz4=logSz4z4; logSzz4zz4(:,end+1) = logSz4z4(:,1);

        subplot(2,2,i)
        hold on
        mesh(kxxx,omm,logSzz4zz4)
        ax = gca;
        ax.FontSize = 20;
        view(2)
        title(strcat('$S_{UU}(\vec k,\omega,\theta/\pi:',num2str(thetaVector(i)),')$'),'Interpreter','latex','FontSize',25)

        ylabel('$\omega$','Interpreter','latex','FontSize',30);
        xlabel('$k_x/\pi$','Interpreter','latex','FontSize',30);
        xlim([0-(kvecm(1,2)-kvecm(1,1))/2,2-(kvecm(1,2)-kvecm(1,1))/2])
        ylim([0,24/2])

        %colorbar('location','north')
        ab{i} = logSzz4zz4; 
    end
    drawnow
    bottom1 = min(bottom,min(min(min(ab{1})),min(min(ab{2}))));
    top1  = max(top,max(max(max(ab{1})),max(max(ab{2}))));
    subplot(2,2,1)
%     shading interp;
    % This sets the limits of the colorbar to manual for the first plot
    caxis manual
    caxis([bottom1 top1]);
    subplot(2,2,2)
%     shading interp;
    % This sets the limits of the colorbar to manual for the first plot
    caxis manual
    caxis([bottom1 top1]);
    subplot(2,2,3)
%     shading interp;
    % This sets the limits of the colorbar to manual for the first plot
    caxis manual
    caxis([bottom1 top1]);
    subplot(2,2,4)
%     shading interp;
    % This sets the limits of the colorbar to manual for the first plot
    caxis manual
    caxis([bottom1 top1]);
end

%% Graph - 4: <W_p> vs theta
figure(4)
for i=1
    Lt = 121;

    thetaVector = linspace(0,0.5,Lt);
    Xi = zeros(length(o1.coordxy(1,:)),length(thetaVector));
    for ii = 1:length(thetaVector)
        [eigen_t,eigVec_t] = o1.mbspectrum(thetaVector(ii)*pi);
        [eigen_t,indc] = sort(real(eigen_t{6}));
        eigVec_t = eigVec_t{6};
        eigVec_t = eigVec_t(:,indc);
        psi0 = (o1.ksec.projMsubs{6}.'*eigVec_t(:,1));
        for pla = 1:length(o1.coordxy(1,:))
            Xi(pla,ii) = psi0'*o1.Z4.Z4mat{pla}*psi0;
        end
    end
    plot(thetaVector,Xi(5,:),'LineWidth',3)
    ylabel('$\langle W_p \rangle$','Interpreter','latex','FontSize',25)
    xlabel('$\theta/\pi$','Interpreter','latex','FontSize',25);
end

%% Graph - 5: < psi | W_p exp(iHt) X_i exp(-iHt) W_p | psi > vs theta - time
figure(5)
for i=1
    Lt   = 4;
    Ltim = 121;
    time        = linspace(0,9,Ltim);
    thetaVector = linspace(0.01,0.5,Lt);

    pla = 5;
    Xi = zeros(o1.numbits,Ltim,Lt);
    for ii = 1:Lt
        [eigen_t,eigVec_t] = o1.mbspectrum(thetaVector(ii)*pi);
        eigVec_tall = [];
        eigen_tall = [];
        for iiV = 1:length(eigVec_t)
            [eigen_tt,indc] = sort(real(eigen_t{iiV}));
            eigVec_tt = eigVec_t{iiV};
            eigVec_tt = eigVec_tt(:,indc);
            eigVec_tall=[eigVec_tall,o1.ksec.projMsubs{iiV}.'*eigVec_tt];
            eigen_tall =[eigen_tall;reshape(eigen_tt,[length(eigen_tt),1])];
        end

        [eigen_t0,eigVec_t0] = o1.mbspectrum(0.495*pi);

        [eigen_tt0,indc0] = sort(real(eigen_t0{6}));
        eigVec_tt0 = eigVec_t0{6};
        eigVec_tt0 = eigVec_tt0(:,indc0);
        psi0 = (o1.ksec.projMsubs{6}.'*eigVec_tt0(:,1));
        psi0 = psi0/sqrt(sum(abs(psi0).^2));
        %         psi0 = eigVec_tt0(:,1);
        psi0_p = o1.Z4.Z4mat{pla}*psi0;
        psi0_p = psi0_p/sqrt(sum(abs(psi0_p).^2));
        %         psi_p = psi0;

        cm_l = eigVec_tall'*psi0_p;
        cm_l = cm_l/sqrt(sum(abs(cm_l).^2));
        ei_l = eigen_tall-min(eigen_tall);

        for it = 1:Ltim
            psi_t = sum(diag(exp(-1j*ei_l*time(it)).*cm_l)*eigVec_tall.').';
            psi_t = psi_t/sqrt(sum(abs(psi_t).^2));
            for indXi = 1:o1.numbits
                Xi(indXi,it,ii) =  psi_t'*o1.x_imat{indXi}*psi_t;
            end
        end
    end
end

%%
flag = 1;
close all
str = 'Z4_timeEvolve_Lattice3x2_v1.gif';
for iii=1
    % Burada zamana bağlı olarak figürler elde edilmeli
    close
    % choose a theta slice
    for indTime = 1:Ltim
        % plot the bond expectation values at each time step with a pause functionality
        % visualise the bonds
        % instead of spy, I would like to use scatter to change the marker size
        % I thereby need the location of each element,
        % y-axis is reversed, for m x n matrix
        % top-most row y - location is m
        % bottom-most row y - location is 1
        % left-most column x - location is 1
        [spylocX,spylocY] = meshgrid(1:2*o1.Lx,2*o1.Ly:-1:1);

        aa = reshape(str2num(o1.idS{1}.'),[2*o1.Ly,2*o1.Lx]);
        [nn,mm] = meshgrid(1:o1.Lx,1:o1.Ly);
        nn = nn(:)*2-0.5; mm = mm(:)*2-1.5;

        % visualise the initial quench
        [inta,~]=find(diag(o1.xmat==4));
        aa = reshape(str2num(o1.Z4.Z4Loc{5}.'),[2*o1.Ly,2*o1.Lx]);
        %     % visualise the time evolution
        for thetaIndex = 1:Lt
            hold off
            subplot(2,2,thetaIndex)
            scatter(spylocX,spylocY,60*5,'ro');
            hold on

            for ip = 1:length(nn)
                poly2.Vertices(:,1) = poly1.Vertices(:,1) + nn(ip);
                poly2.Vertices(:,2) = poly1.Vertices(:,2) + mm(ip);
                plot(poly2);
            end
            axis equal;
            % now, introduce the plaquettes
            abcd=abs(reshape(Xi(:,indTime,thetaIndex)-1,[o1.Ly*2,o1.Lx*2]));
            %             abcd=abcd(end:-1:1,:);
            abcd(abcd<1e-2)=0;

            scatter(spylocX,spylocY,fix(60*abcd)*5+1e-1,'blue','filled')
            title(strcat("$\theta/\pi : ", num2str(thetaVector(thetaIndex)),'$'),'Interpreter','latex')
            xlim([0,1+2*o1.Lx])
            ylim([0,1+2*o1.Ly])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
        end
        drawnow
        pause(0.1)
        if flag == 1
            frame = getframe(gcf);
            img =  frame2im(frame);
            [img,cmap] = rgb2ind(img,256);
            if indTime == 1
                imwrite(img,cmap,str,'gif','LoopCount',Inf,'DelayTime',0.1);
            else
                imwrite(img,cmap,str,'gif','WriteMode','append','DelayTime',0.1);
            end
        end
    end

end

%%

%% Graph - 6: < psi | X_i exp(iHt) W_p exp(-iHt) X_i | psi > vs theta - time
figure(6)
for i=1
    Lt   = 4;
    Ltim = 121;
    time        = linspace(0,9,Ltim);
    thetaVector = linspace(0.005,0.495,Lt);
    polyin = polyshape([0 0 sqrt(2) sqrt(2)],...
        [sqrt(2) 0 0 sqrt(2)]);
    poly1 = rotate(polyin,45);
    poly2 = poly1;
    pla = 9;
    Z4p = zeros(length(o1.Z4.Z4mat),Ltim,Lt);
    for ii = 1:Lt
        [eigen_t,eigVec_t] = o1.mbspectrum(thetaVector(ii)*pi);
        eigVec_tall = [];
        eigen_tall = [];
        for iiV = 1:length(eigVec_t)
            [eigen_tt,indc] = sort(real(eigen_t{iiV}));
            eigVec_tt = eigVec_t{iiV};
            eigVec_tt = eigVec_tt(:,indc);
            eigVec_tall=[eigVec_tall,o1.ksec.projMsubs{iiV}.'*eigVec_tt];
            eigen_tall =[eigen_tall;reshape(eigen_tt,[length(eigen_tt),1])];
        end

        [eigen_t0,eigVec_t0] = o1.mbspectrum(0.01*pi);

        [eigen_tt0,indc0] = sort(real(eigen_t0{6}));
        eigVec_tt0 = eigVec_t0{6};
        eigVec_tt0 = eigVec_tt0(:,indc0);
        psi0 = (o1.ksec.projMsubs{6}.'*eigVec_tt0(:,1));
        psi0 = psi0/sqrt(sum(abs(psi0).^2));
        %         psi0 = eigVec_tt0(:,1);
        psi0_p = o1.x_imat{pla}*psi0;
        psi0_p = psi0_p/sqrt(sum(abs(psi0_p).^2));
        %         psi_p = psi0;

        cm_l = eigVec_tall'*psi0_p;
        cm_l = cm_l/sqrt(sum(abs(cm_l).^2));
        ei_l = eigen_tall-min(eigen_tall);

        for it = 1:Ltim
            psi_t = sum(diag(exp(-1j*ei_l*time(it)).*cm_l)*eigVec_tall.').';
            psi_t = psi_t/sqrt(sum(abs(psi_t).^2));
            for indZp = 1:length(o1.Z4.Z4mat)
                Z4p(indZp,it,ii) =  psi_t'*o1.Z4.Z4mat{indZp}*psi_t;
            end
        end
    end
end

%%
flag = 1;
close all
figure
str = 'Xi_timeEvolve_Lattice3x2_v1.gif';
for iii = 1
    % Burada zamana bağlı olarak figürler elde edilmeli
    close
    % choose a theta slice
    for indTime = 1:Ltim
        % plot the bond expectation values at each time step with a pause functionality
        % visualise the bonds
        % instead of spy, I would like to use scatter to change the marker size
        % I thereby need the location of each element,
        % y-axis is reversed, for m x n matrix
        % top-most row y - location is m
        % bottom-most row y - location is 1
        % left-most column x - location is 1
        [spylocX,spylocY] = meshgrid(1:2*o1.Lx,2*o1.Ly:-1:1);
        [spylocXp,spylocYp] = meshgrid(1:2*o1.Lx,2*o1.Ly-1:-2:1);
        spylocXp = spylocXp+0.5;        spylocYp = spylocYp+0.5;

        spylocYp(:,[2:2:end]) = spylocYp(:,[2:2:end]) +1;
        [nn,mm] = meshgrid(1:o1.Lx,1:o1.Ly);
        nn = nn(:)*2-0.5; mm = mm(:)*2-1.5;

        % visualise the initial quench
        [inta,~]=find(diag(o1.xmat==4));
        %     % visualise the time evolution
        for thetaIndex = 1:Lt
            hold off
            subplot(2,2,thetaIndex)
            scatter(spylocX,spylocY,60*5,'ro');
            hold on

            for ip = 1:length(nn)
                poly2.Vertices(:,1) = poly1.Vertices(:,1) + nn(ip);
                poly2.Vertices(:,2) = poly1.Vertices(:,2) + mm(ip);
                plot(poly2);
            end
            axis equal;
            % now, introduce the plaquettes
            abcd=abs(reshape(Z4p(:,indTime,thetaIndex)-1,[o1.Ly,2*o1.Lx]));
            %             abcd=abcd(end:-1:1,:);
            abcd(abcd<1e-2)=0;

            scatter(spylocXp,spylocYp,fix(60*abcd)*12+1e-1,'blue','filled')
            title(strcat("$\theta/\pi : ", num2str(thetaVector(thetaIndex)),'$'),'Interpreter','latex')
            xlim([0,1+2*o1.Lx])
            ylim([0,1+2*o1.Ly])
            set(gca,'xtick',[])
            set(gca,'ytick',[])
        end
        drawnow
        hold off
        pause(0.1)
        if flag == 1
            frame = getframe(gcf);
            img =  frame2im(frame);
            [img,cmap] = rgb2ind(img,256);
            if indTime == 1
                imwrite(img,cmap,str,'gif','LoopCount',Inf,'DelayTime',0.1);
            else
                imwrite(img,cmap,str,'gif','WriteMode','append','DelayTime',0.1);
            end
        end
    end

end

%% Energy as a func of theta
clear all
close all
figure
for ii = 2:3
    o1 = makespace_confinedLimitp(2,ii,400);
    lTh = 41;
    thetaVector = linspace(0.18,0.25,lTh);
    % figure
    for i=1:lTh
        [eigen_t ,eigVec_t ] = o1.mbspectrum(thetaVector(i)*pi);
        eigen_tt  = sort(real(cell2mat(eigen_t.')));
        scatter(0*real(eigen_tt(1:2))+thetaVector(i),real(eigen_tt(1:2))-min(real(eigen_tt)),10,'b','filled');
        hold on
    end
    hold off
    pause(1)
end

%% Energy gap scaling
clc
clear
close all

N = 11;
LL = 5;
thetaVector = linspace(0,0.5,N);
SXX = zeros(LL,N);
for ii = 1:LL
    ii
    o1 = makespace_confinedLimitp(ii,1,508);
    E0 = 1e5;
    for i = 1:N
        [eigen_t ,eigVec_t] = o1.mbspectrum_noKSector(thetaVector(i)*pi);
        E0 = min(abs(eigen_t(2)-eigen_t(1)),E0);
        %Entanglement entropy
        SXX(ii,i) = o1.entEntropy(thetaVector(i),eigVec_t(:,1));
    end
    EGap(ii) = E0;
end
%% (Energy gap vs theta) and (Entanglement Entropy vs theta)
close all;
subplot(1,2,1);

stem(1./(4*(1:LL)),(EGap),'b','filled',"MarkerSize",20)
ylabel("$\Delta E_{gap}$","Interpreter","latex","FontSize",36)
xlabel("$1/N$","Interpreter","latex","FontSize",26)
xlim([0,1/4])
subplot(1,2,2)
hold off
for ii = 1:LL
    subplot(1,2,2)
    hold on
    plot(thetaVector,real(SXX(ii,:)),'DisplayName',strcat('N = ',num2str(ii)),'LineWidth',2);
    hold off
    xlabel("$\theta/\pi$","Interpreter","latex","FontSize",26)
    ylabel("$\mathcal{S}_{4N}(\theta)$","Interpreter","latex","FontSize",36)
    legend
end
drawnow
