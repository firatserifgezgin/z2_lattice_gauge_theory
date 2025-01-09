classdef makespace_z2_lgt
    % This class is used generate contrained subspace of Z2 basis for a
    % given plaquette lattice size m x n where bond size is 4 m n
    % The basis is generated for the gauge choice G_i = 1.
    % It then generates the permutation cycles and group elements under
    % translations to decompose into k-space subsectors, it then allows for
    % block-diagonal structure

    % below one can find important attiributes and functions
    % object:
    % o1 = makespace_confinedLimit(m,n,StringLength)

    % properties:
    % o1.idS        : generates gauge-invariant basis (all loops of size X)
    % o1.xmat,zmat  : representation X and Z in idS
    % o1.X_k        : representation of X_k for each k-sector
    % o1.Z4         : representation of ZZZZ for each plaquette
    % o1.ksec       : clustering of idS into k-sectors and projectors

    % methods:
    % makespace_z2_lgt(n,m,mxSt): constructor, n x m plaq lattice,
    %                                       mxSt: maxString length 
    % o1.visStateM(idS{1},visualsOn,ifLabelPlaquettes) : visualize idS
    % [eigens,eigVecs] = o1.mbspectrum(theta)       : diag for all k-sec
    % [eigens,eigVecs] = o1.mbzerospectrum(theta)   : diag for k = 0
    % [eigens,eigVecs] = o1.mbspectrumKsector(theta,k) : diag for k = k0
    % [Xij]            = o1.S_XX_ij(obj,eigen,eigVec,i1,i2):S_XX(omega,i,j)
    % [Xkk_abs2,kk]    = o1.S_XX_kkp(eii,eiVec): S_XX(omega,k)
    % [Z4ij]           = o1.S_Z4Z4_ij(obj,eii,eiVec,i1,omega):S_UU(omega,i,j)
    % [Z4kk_abs,kk]    = o1.S_Z4Z4_kkp(eii,eiVec,i1,omega): S_UU(omega,k)

    properties
        id ;
        idS;
        idS_uncut_basis;
        idS_cut;
        idS_cut_horz;
        idS_uncut;
        idS_uncut_horz;
        idS_cut_matrix;
        idS_uncut_matrix;

        idStLen;
        Lx ;
        Ly ;
        st = 0;
        numbits int64;
        oldstatenum;
        newstatenum;
        coordxy;
        coordPartSites;
        coordPartSiteNames;
        coordBondSites;
        coordBondSiteNames;
        newstate;
        newstatePMPM;
        newstateMPMP;
        newstatePPMM;
        newstateMMPP;
        newstateX;
        zpairs;
        zmat;
        xFlipmat;
        xmat;
        x_imat_rep;
        x_imat;
        xtranslatemat;
        ytranslatemat;
        xtranslateperm;
        ytranslateperm;
        staterules1;
        staterules2;
        state;
        ksec struct;
        X_k;
        X_krep;
        U_k;
        Z4;
        X4;
    end

    methods
        function obj = makespace_z2_lgt(n,m,mxSt)
            obj.Lx = n;
            obj.Ly = m;
            obj.numbits = 4*n*m;
            obj.st = 0;
            obj.id(1) = 0;
            % topological sector 1
            obj.idS{1}          = pad('',obj.numbits,'left','0'); % set the first state
            obj.idS_cut{1}      = pad('',obj.numbits/2,'left','0'); % set the first state
            obj.idS_cut_horz{1} = pad('',obj.numbits/2,'left','0'); % set the first state

            obj.idS_uncut{1}        = pad('',obj.numbits/2,'left','0'); % set the first state
            obj.idS_uncut_horz{1}   = pad('',obj.numbits/2,'left','0'); % set the first state

            obj.idS_uncut_basis{1} = pad('',obj.numbits/2,'left','0'); % set the first basis state
            % This section does include the other topological sectors if necessary
            % topological sector 2
            %             if n>m
            %                 idsInt1 = diag(ones(1,2*max(m,n)),1);
            %             elseif m>n
            %                 idsInt1 = diag(ones(1,2*max(m,n)),-1);
            %             else
            %                 idsInt1 = diag(ones(1,2*m));
            %             end
            %             idsInt1= idsInt1(1:2*m,1:2*n);
            %             obj.idS{2} = strrep(num2str(idsInt1(:).'),' ','');
            %             % topological sector 3
            %             idsInt2 = flip(idsInt1);
            %             obj.idS{3} = strrep(num2str(idsInt2(:).'),' ','');
            %             % topological sector 4
            %             idsInt3 = xor(idsInt1,idsInt2);
            %             obj.idS{4} = strrep(num2str(idsInt3(:).'),' ','');
            %

            obj.oldstatenum = length(obj.id);
            obj.newstatenum = length(obj.id)+1;

            coorPartSites = [];
            [cBondx,cBondy] = meshgrid(1:2*obj.Lx,1:2*obj.Ly);
            obj.coordBondSites = [2*cBondx(:),2*cBondy(:)];
            obj.coordBondSiteNames = strrep(join(pad(string(obj.coordBondSites.'),2,'left','0').'),' ','').';
            coorxy = [];
            for i = 1:2*obj.Lx
                if mod(i,2) == 1;j = 1;else;j = 2;end
                for jj = j:2:2*obj.Ly
                    coorxy = [coorxy,[i;jj]];
                    coorPartSites = [coorPartSites,[2*i-1, 2*i+1;2*jj+1, 2*jj-1]];
                end
            end

            obj.coordxy = coorxy;
            strSitesNum = strrep(join(pad(string(coorPartSites),2,'left','0').'),' ','').';
            [~,ia] = unique(strSitesNum);
            obj.coordPartSites = coorPartSites(:,ia);
            obj.coordPartSiteNames = strSitesNum(ia).';

            for indxy = 1:length(coorxy(1,:))
                ix = coorxy(1,indxy); iy = coorxy(2,indxy);
                plaq = cell2mat(obj.zplaqLoc(ix,iy));
                obj.oldstatenum = length(obj.idS);
                for indS = 1:length(obj.idS)
                    idBitForm = reshape(str2num(obj.idS{indS}.'),[2*m,2*n]);
                    obj.newstate = bitxor(idBitForm,plaq);
                    %                     [idBitForm,obj.newstate]

                    if sum(obj.newstate,'all')<=mxSt
                        for ix = 0:2*obj.Lx-1
                            for iy = 0:2*obj.Ly-1
                                nSs = circshift(obj.newstate,[-ix-iy,ix-iy]);
                                intNewSt = strrep(num2str(nSs(:).'),' ','');
                                if sum(strcmp(intNewSt,obj.idS),'all')== 0
                                    obj.idS{end+1} = intNewSt;

                                    nSs_cut_vert = nSs(:,1:end/2);
                                    nSs_uncut_vert = nSs(:,end/2+1:end);
                                    nSs_cut_horz = nSs(1:end/2,:);
                                    nSs_uncut_horz = nSs(end/2+1:end,:);

                                    intNewSt_vert_cut   = strrep(num2str(nSs_cut_vert(:).'),' ','');
                                    intNewSt_vert_uncut = strrep(num2str(nSs_uncut_vert(:).'),' ','');
                                    intNewSt_horz_cut   = strrep(num2str(nSs_cut_horz(:).'),' ','');
                                    intNewSt_horz_uncut = strrep(num2str(nSs_uncut_horz(:).'),' ','');

                                    obj.idS_cut{end+1}        = intNewSt_vert_cut; %vertical default
                                    obj.idS_uncut{end+1}      = intNewSt_vert_uncut; % vertical_default
                                    obj.idS_cut_horz{end+1}   = intNewSt_horz_cut; % horizontal
                                    obj.idS_uncut_horz{end+1} = intNewSt_horz_uncut; % horizontal

                                    obj.idStLen(end+1) = sum(obj.newstate,'all');
                                    if sum(strcmp(intNewSt_vert_uncut,obj.idS_uncut_basis),'all')== 0
                                        obj.idS_uncut_basis{end+1}=intNewSt_vert_uncut;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            L_HS = length(obj.idS);
            zp = [];
            if length(coorxy(1,:))<=2
                coorLength = 1;
            else
                coorLength =length(coorxy(1,:));
            end
                
            for indxy = 1:coorLength
                ix = coorxy(1,indxy); iy = coorxy(2,indxy);
                plaq = cell2mat(obj.zplaqLoc(ix,iy));

                z0 = [];

                for indS = 1:L_HS
                    idBitForm = reshape(str2num(obj.idS{indS}.'),[2*m,2*n]);
                    obj.newstate = bitxor(idBitForm,plaq);
                    intNewSt = strrep(num2str(obj.newstate(:).'),' ','');
                   
                    if sum(strcmp(intNewSt,obj.idS),'all') ~= 0
                        kk = find(strcmp(intNewSt,obj.idS)); lkk = length(kk);
                        z0 = [z0,[indS*ones(1,lkk);kk]];
                        zp = [zp,[indS*ones(1,lkk);kk]];
                    end
                end
                obj.Z4.Z4pairs{indxy} = z0;
                intZmat = sparse(L_HS,L_HS);
                intZmat((z0(2,:)-1)*L_HS+z0(1,:)) = 1;
                obj.Z4.Z4mat{indxy} = intZmat;
                obj.Z4.Z4Loc{indxy} = strrep(num2str(plaq(:).'),' ','');
                obj.Z4.Z4Numbers(indxy) = indxy;
            end
            for indXi = 1:obj.numbits
                for indS = 1:L_HS
                    obj.newstateX = str2num(obj.idS{indS}.');
                    Xi_int(indS)=obj.newstateX(indXi);
                end
                obj.x_imat{indXi} = diag(Xi_int);
                obj.x_imat_rep{indXi} = diag((Xi_int-1/2)*2);
            end

            ss = obj.idS;
            siz = size(ss.');

            obj.zpairs = zp;
            obj.zmat = sparse(zp(1,:),zp(2,:),ones(1,length(zp(1,:))));
            for indSS = 1:prod(siz)
                xmInter(indSS) = sum(double(ss{indSS})-48);
            end
            obj.xFlipmat = sparse(1:L_HS,1:L_HS,xmInter);
            obj.xmat=sparse((2*obj.xFlipmat-double(obj.numbits)*eye(length(obj.xFlipmat))));

            % determine the location of each loop elements
            plaq_vec_indcs = find(diag(obj.xFlipmat) == 4);
            plaq_vec_indcs = plaq_vec_indcs(1:min(length(plaq_vec_indcs),2*obj.Lx*obj.Ly));
            for idP = 1:(length(obj.Z4.Z4mat)-1)
                x_locs= find(double(obj.idS{plaq_vec_indcs(idP)})-48);
%                 plaq_bond_loc{idP} = x_locs;
                obj.X4.X4matSec{idP} = (eye(length(obj.idS))-1j*obj.x_imat_rep{x_locs(1)})*(eye(length(obj.idS))-1j*obj.x_imat_rep{x_locs(2)})*(eye(length(obj.idS))+1j*obj.x_imat_rep{x_locs(3)})*(eye(length(obj.idS))+1j*obj.x_imat_rep{x_locs(4)})/4;
            end
            % note that current X definition is with eigenvalues 0 and 1,
            % normally it is 1 or -1 and SU(2) algebra works for S_i
            % operators

            obj.xtranslateperm = obj.xtranslate();
            obj.xtranslatemat = sparse(obj.xtranslateperm,1:L_HS,ones(1,L_HS));
            obj.ytranslateperm = obj.ytranslate();
            obj.ytranslatemat = sparse(obj.ytranslateperm,1:L_HS,ones(1,L_HS));
            obj.ksec = obj.generateKSectors(obj.xtranslateperm,obj.ytranslateperm);
            obj.X_k = obj.generateX_k;
            for ii1 = 1:length(obj.X_k)
                obj.X_krep{ii1}= obj.X_k{ii1};
%                 obj.X_krep{ii1}=sparse((2*obj.X_k{ii1}-double(obj.numbits)*eye(length(obj.X_k{ii1}))));
            end
            obj.U_k = obj.generateU_k;
            obj.idS_cut_matrix = zeros(L_HS,L_HS);
            obj.idS_uncut_matrix = zeros(L_HS,L_HS);

            for ii = 1:L_HS
                obj.idS_cut_matrix(ii,:)  = double(strcmp(obj.idS_cut{ii},obj.idS_cut));
                obj.idS_uncut_matrix(ii,:) = double(strcmp(obj.idS_uncut{ii},obj.idS_uncut));
            end
            obj.idS_cut_matrix = sparse(obj.idS_cut_matrix);
            obj.idS_cut_matrix = sparse(obj.idS_cut_matrix-eye(L_HS));
            obj.idS_uncut_matrix = sparse(obj.idS_uncut_matrix);
            obj.idS_uncut_matrix = sparse(obj.idS_uncut_matrix-eye(L_HS));
        end

        function outputArg = fromcoord(obj,i,j)
            outputArg = mod(i-1,2*obj.Lx) + mod(j-1,2*obj.Ly)*2*obj.Lx;
        end

        function outputArg = bitvalues(obj,coord)
            outputArg = reshape(str2num(pad(dec2bin(coord),obj.numbits,'left','0').').',[2*obj.Lx,2*obj.Ly]);
        end

        function outputArg = zplaqLoc(obj,ii,jj)
            outputArg = [];
            for iiInd = 1:length(ii)
                i = ii(iiInd);j = jj(iiInd);
                ou = zeros(2*obj.Ly,2*obj.Lx);
                ou(2*obj.Ly-mod([j,j+1,j,j+1]-1,2*obj.Ly),mod([i,i,i+1,i+1]-1,2*obj.Lx)+1) = 1;
                outputArg{iiInd} = ou;
            end
        end

        function outputArg = xtranslate(obj)
            outputArg = [];
            for indS = 1:length(obj.idS)
                aa = reshape(str2num(obj.idS{indS}.'),[2*obj.Ly,2*obj.Lx]);
                intShtX = num2str(circshift(aa,[-1,1]));
                ou1 = strrep(intShtX(:).',' ','');
                abb = find(strcmp(ou1,obj.idS));
                if isempty(abb)
                    for im = 1:2
                        subplot(1,2,im)
                        if im ==2
                            aa = reshape(str2num(ou1.'),[2*obj.Ly,2*obj.Lx]);
                        end
                        polyin = polyshape([0 0 sqrt(2) sqrt(2)],...
                            [sqrt(2) 0 0 sqrt(2)]);
                        poly1 = rotate(polyin,45);
                        poly2 = poly1;
                        [nn,mm] = meshgrid(1:obj.Lx,1:obj.Ly);
                        nn = nn(:)*2-0.5; mm = mm(:)*2-1.5;
                        %                         spy(aa,30);
                        if 1 == 1
                            hold on;
                            spy(ones(2*obj.Ly,2*obj.Lx),'ro',20);
                            for ip = 1:length(nn)
                                poly2.Vertices(:,1) = poly1.Vertices(:,1) + nn(ip);
                                poly2.Vertices(:,2) = poly1.Vertices(:,2) + mm(ip);
                                plot(poly2);
                            end
                            axis equal;
                            [yy,indyy] = sort(obj.coordxy(2,:));
                            text(obj.coordxy(1,indyy)+0.5,obj.Ly*2-yy+0.5,num2str((1:length(obj.coordxy(2,:))')'));

                            hold off;
                        end
                        set(gca,'YTick',[],'XTick',[]);
                    end
                    pause
                end
                outputArg = [outputArg,find(strcmp(ou1,obj.idS))]; % ya 2 satır yukarda transpose almak sorun çıkardı ya da önceki bitvalues fonksiyonu yanlış
            end
        end

        function outputArg = ytranslate(obj)
            outputArg = [];
            for indS = 1:length(obj.idS)
                aa = reshape(str2num(obj.idS{indS}.'),[2*obj.Ly,2*obj.Lx]);
                intShtY = num2str(circshift(aa,[-1,-1]));
                ou1 = strrep(intShtY(:).',' ','');
                outputArg = [outputArg,find(strcmp(ou1,obj.idS))]; % ya 2 satır yukarda transpose almak sorun çıkardı ya da önceki bitvalues fonksiyonu yanlış
            end
        end

        function aa = visStateM(obj,state,vis,plaqInd)
            aa = reshape(str2num(obj.idS{state}.'),[2*obj.Ly,2*obj.Lx]);
            polyin = polyshape([0 0 sqrt(2) sqrt(2)],...
                [sqrt(2) 0 0 sqrt(2)]);
            poly1 = rotate(polyin,45);
            poly2 = poly1;
            [nn,mm] = meshgrid(1:obj.Lx,1:obj.Ly);
            nn = nn(:)*2-0.5; mm = mm(:)*2-1.5;

            spy(aa,30);
            if vis == 1
                hold on;
                spy(ones(2*obj.Ly,2*obj.Lx),'ro',20);
                for ip = 1:length(nn)
                    poly2.Vertices(:,1) = poly1.Vertices(:,1) + nn(ip);
                    poly2.Vertices(:,2) = poly1.Vertices(:,2) + mm(ip);
                    plot(poly2);
                end
                axis equal;
                if plaqInd ==1
                    [~,indyy] = sort(obj.coordxy(1,:));
                    text(obj.coordxy(1,indyy)+0.5,obj.Ly*2-obj.coordxy(2,indyy)+0.5,num2str((1:length(obj.coordxy(2,:))')'));
                end
                hold off;
            end
            set(gca,'YTick',[],'XTick',[]);
        end

        function cycleNotation = permutationCycles(obj,permutation)
            n = length(permutation);
            cycleNotation = {};
            visited = false(1, n);

            for i = 1:n
                if ~visited(i)
                    cycle = [i];
                    j = permutation(i);
                    visited(i) = true;

                    while j ~= i
                        cycle = [cycle, j];
                        visited(j) = true;
                        j = permutation(j);
                    end

                    cycleNotation{end+1} = cycle;
                end
            end
        end

        function permut1 = permutationReplace(obj,permut, elements)
            permut1  = permut;
            n = length(permut);
            for i = 1:n
                permut1{i} = elements(permut{i});
            end
        end

        function [oArgs1,oArgs2] = cyclePerms(obj,xcyc,yp)
            oA = obj.permutationReplace(xcyc,yp);
            oArgs1 = zeros(1,length(xcyc));
            oArgs2 = 0*oArgs1;

            for indx = 1:length(xcyc)
                xcL = length(xcyc{indx});
                if  (xcL-1 > 0)
                    for indy = 1:length(xcyc)
                        %                         sumXcoA = sum(abs(xcyc{indx}-sort(oA{indy})));
                        abb = sort(xcyc{indx});
                        abb=strrep(num2str(abb(:).'),' ','');
                        abc = sort(oA{indy});
                        abc=strrep(num2str(abc(:).'),' ','');
                        %                         if sumXcoA < 1e-5
                        if strcmp(abb,abc) == true
                            oArgs1(indx) = indy;
                            for indz = 0:(xcL-1)
                                if sum(abs(circshift(xcyc{indx},indz)-oA{indy}))< 1e-8
                                    oArgs2(indx) = indz;
                                end
                            end
                            break
                        end
                    end
                else
                    oArgs1(indx) = indx;
                    oArgs2(indx) = 0;
                end
            end
        end

        function ksec1 = generateKSectors(obj,xp,yp)
            xcyc = obj.permutationCycles(xp);
            [ycp,ycs] = obj.cyclePerms(xcyc,yp);
            subspaces = obj.permutationCycles(ycp);
            Lsub = length(subspaces);

            staterules = {};
            statecycles= [];
            kxkypairs = {};
            stateVecs = {};
            ijk = 1;

            for inds = 1:Lsub
                cycind = subspaces{inds};
                cycles = xcyc{cycind};
                shifts = ycs(cycind);
                ashInt = cumsum(shifts);
                ashift = [0,ashInt(1:end-1)];
                xlen = length(cycles);
                ylen = length(cycind);
                totalShift = mod(sum(shifts),xlen);

                if totalShift ==0
                    numLoops = 1;
                else
                    numLoops = (lcm(totalShift,xlen)/totalShift);
                end
                [jj,rr] = meshgrid(1:ylen,1:numLoops);
                jj = jj(:);rr = rr(:);
                [iii,rrr] = meshgrid(1:xlen,rr);
                [~,jjj]   = meshgrid(1:xlen,jj);
                for i1 = 0:(xlen-1)
                    kx = 2*i1/xlen;
                    minKy = mod(kx*totalShift,2)/ylen;
                    kySize = 2/ylen;
                    for j1 = 1:ylen
                        ky = minKy + (j1-1)*kySize;
                        interStateRule = sum(exp(-1i*pi*kx*(iii-1-ashift(jjj) ...
                            -(rrr-1)*totalShift) - 1i*pi*ky*(jjj-1+ylen*(rrr - 1)))) ;
                        staterules{ijk} = interStateRule;
                        statecycles{ijk} = cycles;
                        kxkypairs{ijk} = [kx;mod(ky,2)];
                        ijk = ijk+1;
                    end
                end
            end
            for i = 1:length(staterules)
                a0=(fix(100*(kxkypairs{i}))/100).';
                %                 a0Index(i) = sum(2.^[4*a0(2),2*a0(1)]);
                a0Index{i} = strcat(num2str(a0(2)),num2str(a0(1)));
                a=statecycles{i};
                b=staterules{i};
                c = sparse(ones(1,length(a)),a,b,1,length(obj.idS));
                stateVecs{i} = c/sqrt(sum(abs(c).^2));
            end
            [G,IDd] = findgroups(a0Index);
            for j=1:length(IDd)
                cj = find(strcmp(a0Index,IDd{j}));
                %                 kxkyGroups{j} = cj(:);
                kxkypairsOrd{j} = kxkypairs{cj(1)}.';
                %                 statecyclesOrd{j} = {statecycles{cj(:)}};
                %                 staterulesOrd{j} = {staterules{cj(:)}};
                stateVecsOrd{j} = {stateVecs{cj(:)}};
            end
            ksec1.kvec = kxkypairsOrd;
            ksec1.proj = stateVecsOrd;
            kksecproj = [];
            for i=1:length(ksec1.proj)
                kksecprojsubs{i} = cell2mat(ksec1.proj{i}.');
                kksecproj=[kksecproj;cell2mat(ksec1.proj{i}.')];
            end
            ksec1.projM = kksecproj;
            ksec1.projMsubs = kksecprojsubs;
        end

        function ksecSub = ksecfind(obj,k)
            kvecm = cell2mat(obj.ksec.kvec.').';
            index = find(kvecm(1,:)==k(1) & kvecm(2,:)==k(2));

            ksecSub.kvec = k;
            ksecSub.proj = {obj.ksec.proj{index}};
            ksecSub.projM= cell2mat(obj.ksec.proj{index}.');
        end

        function [eigens,eigVecs] = mbspectrum(obj,theta)
            for i = 1:length(obj.ksec.kvec)
                H{i} = -sparse(cos(theta)*obj.zmat+sin(theta)*obj.xmat);
                basis = obj.ksec.projMsubs{i};
                [V,ei] = eig(full(conj(basis)*H{i}*basis.'));
                [eigens{i},indEi] = sort(real(diag(ei)));
                %                 sort(diag(ei))
                eigVecs{i} = V(:,indEi);
            end
%             H = -sparse(cos(theta)*obj.zmat+sin(theta)*obj.xFlipmat);
            %             basis = obj.ksec.projMsubs{i};
%             [V,ei] = eig(full(H));
%             [eigens{1},indEi] = sort(diag(ei));
%             eigVecs{1} = V(:,indEi);
        end

        function [eigens,eigVecs] = mbspectrum_noKSector(obj,theta)
            H = -sparse(cos(theta)*obj.zmat+sin(theta)*obj.xmat);
            [V,ei] = eig(full(H));
            [eigens,indEi] = sort(real(diag(ei)));
            eigVecs = V(:,indEi);
        end

        function [eigens,eigVecs] = mbzerospectrum(obj,theta)
            H = -sparse(cos(theta)*obj.zmat+sin(theta)*obj.xmat);
            kvecm = cell2mat(obj.ksec.kvec.').';
            index = find(kvecm(1,:)==0 & kvecm(2,:)==0);

            basis = obj.ksec.projMsubs{index};
            [eigVecs,eigens] = eig(full(conj(basis)*H*basis.'));
            [eigens,indEi] = sort(real(diag(eigens)));
            eigVecs = eigVecs(:,indEi);
        end

        function [eigens,eigVecs] = mbspectrumKsector(obj,theta,k)
            H = -sparse(cos(theta)*obj.zmat+sin(theta)*obj.xmat);
            kvecm = cell2mat(obj.ksec.kvec.').';
            index = find(kvecm(1,:)==k(1) & kvecm(2,:)==k(2));

            basis = obj.ksec.projMsubs{index};
            [eigVecs,eigens] = eig(full(conj(basis)*H*basis.'));
            [eigens,indEi] = sort(real(diag(eigens)));
            eigVecs = eigVecs(:,indEi);
        end
   
        function outArgs = generateX_k(obj)
            LL = 1:length(obj.idS);
            kk = obj.ksec.kvec;
            for i1 = 1:length(kk)
                k = kk{i1};
                kx= k(1)*pi;ky = k(2)*pi;
                for i2 = LL
                    [ix_p,iy_p] = meshgrid((0:2*obj.Lx-1)/sqrt(2),(2*obj.Ly-1:-1:0)/sqrt(2));
                    z_prime = 1/sqrt(2)*[1,1;-1,1]*[ix_p(:).';iy_p(:).'];
                    ix = reshape(z_prime(1,:),[2*obj.Ly,2*obj.Lx]);
                    iy = reshape(z_prime(2,:),[2*obj.Ly,2*obj.Lx])-0.5;

                    vecc = 2*(reshape(str2num(obj.idS{i2}.'),[2*obj.Ly,2*obj.Lx])-0.5);
                    cev1 = exp(-1j*(kx*ix+ky*iy))/sqrt(double(obj.numbits));
                    X_ki(i2) = sum(cev1.*vecc,'all');
                end
                outArgs{i1} = sparse(LL,LL,X_ki);
            end
        end

        function [Xij] = S_XX_ij(obj,eigen,eigVec,i1,i2)
            % S_XX_ij creates and excitation at local i bond and measures the
            % propagation of this excitation at an arbitrary j bond. The basis is in X
            % basis, therefore X_i operators can only change the sign of the basis
            % vectors which has a flipped spin, i.e. if X_i crosses with loop-strings
            % within our basis, otherwise it gives a trivial 1.

            U_t = 0*obj.ksec.projM;
            E_t = diag(U_t);
            E_t = E_t-min(E_t);
            k_t = [E_t,E_t];
            ct = 1;
            for ii = 1:length(eigen)
                k0 = obj.ksec.kvec{ii};
                proj = obj.ksec.proj{ii};
                projM= full(obj.ksec.projMsubs{ii});
                vecs = eigVec{ii};
                eigvs = diag(eigen{ii});
                eigvs(abs(eigvs)<1e-9) = 0;
                for ii1 = 1:length(eigvs)
                    k_t(ct,:) = k0;
                    E_t(ct)   = eigvs(ii1);
                    U_t(:,ct) = projM.'*vecs(:,ii1);
                    ct = ct+1;
                end
            end

            % check the states which has a flipped spin at i1
            for j1 = 1:length(obj.idS)
                bbb0 = reshape(str2num(obj.idS{j1}.'),[2*obj.Ly,2*obj.Lx]);
                %                 bbb1 = reshape(str2num(obj.idS{j1}.'),[2*obj.Ly,2*obj.Lx]);
                %                 bbb2 = reshape(str2num(obj.idS{j1}.'),[2*obj.Ly,2*obj.Lx]);

                bb1 = bbb0.*obj.bitvalues(2^obj.fromcoord(i1(2),i1(1)));
                bb2 = bbb0.*obj.bitvalues(2^obj.fromcoord(i2(2),i2(1)));

                overlapi1(j1,1) = -2*sum(bb1(:))+1;
                overlapi2(j1,1) = -2*sum(bb2(:))+1;
            end

            % inject the action of X_i1 and X_i2 on U_t, flip the sign of
            % the vectors which include flipped spins at i1 and i2.
            U_flp1 =  (U_t')*(U_t(:,1).*overlapi1)  ;
            U_flp2 = conj((U_t')*(U_t(:,1).*overlapi2));
            output.matrixElements = U_flp1.*U_flp2;
            output.kvecs = k_t;
            output.energies = E_t-min(E_t);
            % spectral projection to eigenstates
            %             ksecL = length(obj

        end

        function [Xkk_abs,kk] = S_XX_kkp(obj,eii,eiVec)
            kvecm = cell2mat(obj.ksec.kvec.').';
            index = find(kvecm(1,:)==0 & kvecm(2,:)==0);
            eiVecc = eiVec{index};

            vec0 = obj.ksec.projMsubs{index}.'*eiVecc(:,1);
            kk = obj.ksec.kvec{index};

            % calculate the matrix elements
            for i = 1:length(obj.ksec.kvec)
                k0 = obj.ksec.kvec{i};
                X_kmat = obj.X_krep{i};
                eigvecsfull = obj.ksec.projMsubs{i}.'*eiVec{i};
                Xkk_abs{i} = eigvecsfull'*X_kmat*vec0;
                
                eigvecsfull_zero = obj.ksec.projMsubs{index}.'*eiVec{index};
                Xkk_abs_zero{i} = eigvecsfull_zero'*X_kmat*vec0;
            end
            Xkk_abs_zero;
        end

        function U_k = generateU_k(obj)
            kk = obj.ksec.kvec;
            z_p = [1,1;-1,1]*(obj.coordxy-1)/2;

            % locate the positions of single plaq terms for XOR operation
            plaq_vec_indcs = find(diag(obj.xFlipmat) == 4);

            %generate plaq operator representations for each k-sector
            U_k = {};
            for i1 = 1:length(kk)
                k = kk{i1};
                kx= k(1)*pi;ky = k(2)*pi;

                Lp= length(plaq_vec_indcs);
                Lp = length(obj.Z4.Z4mat);
                U_pin = 0*obj.xFlipmat;
                for ii0 = 1:Lp
                    U_pin = U_pin + exp(-1j*(kx*z_p(1,ii0)+ky*z_p(2,ii0)))*obj.Z4.Z4mat{ii0}/sqrt(Lp);
                end
                U_k{i1} = U_pin;
            end
        end

        function [Z4ij] = S_Z4Z4_ij(obj,eii,eiVec,i1,omega)
            [E0,indE0] = min(real(cell2mat(eii.')));
            counter = indE0;
            LL = length(eii);
            for ik = 1:LL
                if counter <= numel(eii{ik})
                    vecc = eiVec{ik};
                    vec0 = obj.ksec.projMsubs{ik}.'*vecc(:,counter);
                    kk = obj.ksec.kvec{ik};
                    break
                else
                    counter = counter-numel(eii{ik});
                end
            end

            % initialize the outputs
            U_t = 0*obj.ksec.projM;
            E_t = diag(U_t);
            E_t = E_t;
            k_t = [E_t,E_t];

            % calculate the full expression together

            Z4mat_i1 = obj.Z4.Z4mat{i1};
            for i2ind= 1:length(obj.coordxy(1,:))
                Z4mat_i2 = obj.Z4.Z4mat{i2ind};
                Z4i = []; Z4j = [];Z4ijw = 0*omega;
                for i = 1:length(obj.ksec.kvec)
                    k0 = obj.ksec.kvec{i};
                    ei3 = eii{i};
                    if size(obj.ksec.projMsubs{i},1)==size(eiVec{i},1)
                        eigvecsfull = obj.ksec.projMsubs{i}.'*eiVec{i};
                        [size(vec0);size(Z4mat_i1);size(eigvecsfull)];
                        Z4i = eigvecsfull'*Z4mat_i1*vec0;
                        Z4j = eigvecsfull'*Z4mat_i2*vec0;
                        Z4ijm = conj(Z4j).*Z4i;
                        for ik = 1: numel(Z4ijm)
                            Z4ijw = Z4ijw + Z4ijm(ik)./(omega-(ei3(ik)-E0));
                        end
                    end
                end
                Z4ij{i2ind} = reshape(Z4ijw,[1,numel(Z4ijw)]);
            end
        end

        function [Z4kk_abs,kk] = S_Z4Z4_kkp(obj,eii,eiVec)
            kvecm = cell2mat(obj.ksec.kvec.').';
            index = find(kvecm(1,:)==0 & kvecm(2,:)==0);
            kk = obj.ksec.kvec{index};
            vv =eiVec{index};
            vec0 = obj.ksec.projMsubs{index}.'*vv(:,1);

            % calculate the matrix elements
            for i = 1:length(obj.ksec.kvec)
                k0 = obj.ksec.kvec{i};
                U_kmat = obj.U_k{i};
                eigvecsfull = obj.ksec.projMsubs{i}.'*eiVec{i};
                Z4kk_abs{i} = eigvecsfull'*U_kmat*vec0;
            end
        end

        function S_A = entEntropy(obj,theta,psi0)
            % for each vector, determine the matching vectors for the spin
            % it is already calculated as obj.idS_cut_matrix
            L_HS = length(obj.idS);
            % define a new basis for the uncut part, each vector are
            % doubled, we acquire the indices
            ss_cut = [];ss_uncut = [];
            ss_cut_full = [];ss_cfull_coeff = [];
            ss_cut_list = {};
            ind = 1;
            for ii = 1:L_HS
                % kesilen kısımları aynı olanlar vektörleri buluyoruz
                ct = find(obj.idS_cut_matrix(ii,:)); 
                ct = ct(ct>ii-1);
                L_ct = length(ct);
                ss_cut       = [ss_cut,       [ii*ones(1,L_ct);ct]];
                ss_cut_full  = [ss_cut_full,  [ii*ones(1,L_ct+1);[ii,ct]]];
                ss_cfull_coeff  = [ss_cfull_coeff,[1/2,ones(1,L_ct)]];

                % burada hangi vektörler bölününce aynı oluyor onu görüyoruz
                unc = find(obj.idS_uncut_matrix(ii,:)); 
                unc = unc(unc>ii-1);
                L_unc = length(unc);
                ss_uncut= [ss_uncut,[ii*ones(1,L_unc);unc]];

                if sum(ss_uncut(2,:)==ii) == 0 
                    ss_uncut_list{ii}=[ii,unc];
                    new_basis_vecs(ind) = ii;
                    ind = ind+1;
                end
            end
            % yekpare vektörler new_basis_vecs vektöründe tutuluyor
            % ss_uncut_list'in de new_basis_vecs elemanlarındaki listelerde
            % aynı cut'ları olan vektörler mevcut

            % şimdi amacımız, ss_cut içindeki tüm elemanları new_basis_vecs
            % içindeki elemanların sayılarıyla değiştirmek
            ss_cut_perm = ss_cut_full;
            L_newB = length(new_basis_vecs);
            for ii = 1:L_newB
                elems = ss_uncut_list{new_basis_vecs(ii)};
                for jj = 1:length(elems)
                    ss_cut_perm(ss_cut_perm==elems(jj)) = ii;
                end
            end

            rhoA = zeros(L_newB,L_newB);
            % determine a basis vectors and matrices for the matches
            % generate density matrix
            for ii = 1:length(ss_cut_full(1,:))
                vec1 = zeros(L_newB,1);   vec1(ss_cut_perm(1,ii),1)  = 1;
                vec2 = zeros(L_newB,1);   vec2(ss_cut_perm(2,ii),1)  = 1;
                coeff = conj(psi0(ss_cut_full(2,ii)))*psi0(ss_cut_full(1,ii))*ss_cfull_coeff(ii); % kesilince eşleşen diğer parçalar
                rhoA = rhoA + coeff*(vec1*vec2');
            end
            rhoA = sparse(rhoA+rhoA');
            
            [eV,eD] = eig(full(rhoA));eD=diag(eD);
            % Entanglement entropy
            S_A = -sum(eD.*log(eD));
        end
    end
end

