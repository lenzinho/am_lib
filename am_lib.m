classdef am_lib

    % 
    % [sc,s2p,p2s] = get_supercell(pc,diag([2,2,2])); sc.u2p=s2p; sc.p2u=p2s; sc.u2i=pc.p2i(s2p);
    % 
    % 
    % nsteps=50; amplitude=1; mode=6; kpt=[.0;.5;.5];
    % [dc_master] = generate_bvk_displacement(bvk,ip,sc,nsteps,kpt,amplitude,mode);
    % 
    % % for i = 1:nsteps
    % %     dc=dc_master; dc.tau=dc.tau(:,:,i);
    % %     save_poscar(get_primitive_cell(dc),sprintf('outfile.POSCAR_%03i',i))
    % % % plot3(dc_master.tau(1,:,i),dc_master.tau(2,:,i),dc_master.tau(3,:,i),'.'); view(2); drawnow; pause(0.1);
    % % end
    % 
    % % % displace phonon mode 
    % % k=[0;0.5;.5]; amp=1; mode=6;
    % % [md] = get_bvk_displacement(bvk,ip,uc,k,amp,mode);
    % 
    % % % run md
    % % dt = 0.1; nsteps = 2000; Q = 5; T = 300;
    % % [md] = run_bvk_md(bvk,ip,uc,dt,nsteps,Q,T);
    %     
    % 
    % 
    % 
    
    properties (Constant)
        tiny      = 1E-3; % precision of atomic coordinates
        eps       = 1E-8; % numerical precision
        units_eV  = 0.06465555; % sqrt( [eV/Ang^2] * [1/amu] ) --> 0.06465555 [eV]
        units_THz = 98.22906;   % sqrt( [eV/Ang^2] * [1/amu] ) --> 98.22906 [THz=1/ps]
        units_GHz = 98229.06;   % sqrt( [eV/Ang^2] * [1/amu] ) --> 98229.06 [GHz=1/fs]
        
        % mex compiler parameters
        usemex    = false;
        FC        = 'ifort'; 
        FFLAGS    = '-O3 -parallel -fpp -fPIC -lmx -lmex -lmat -nofor_main -bundle -implicitnone -assume realloc_lhs';
        MPATH     = '/Applications/MATLAB_R2017a.app';
        LIBS      = ['-L',am_lib.MPATH,'/bin/maci64 -I',am_lib.MPATH,'/extern/include'];
        EXT       = '.mexmaci64';
        DEBUG     = '-debug'
    end
    

    % general-purpopse functions
    
    methods (Static)
        
        % benchmarking
        
        function t = benchmark_(fnct_,inputs)
            % fnct_ = @(A) outerc_(A,A,A); 
            % x = [1,2,5,10,20,50,100,200,500,1000,2000,5000];
            % for i=1:numel(x); inputs{i}=single(rand(3,x(i))); end
            % loglog(x,benchmark_(fnct_,inputs),'-')
            
            % estimate run times
            ninputs = numel(inputs);
            for i = 1:ninputs
                tic
                    fnct_(inputs{i}); 
                t(i) = toc;
            end
            
            % estimate number of runs which fit in 10 seconds
            nruns = max([round(10./t);ones(size(t))]);
            
            % take avergaes
            for i = 1:ninputs
                tic
                    for j = 1:nruns(i)
                        fnct_(inputs{i});
                    end
                t(i) = toc/nruns(i);
            end
        end


        % numerical precision

        function [C] = mod_(A,tol)

            % set default numerical tolernece
            if nargin < 2; tol = am_lib.tiny; end
            
            C = mod(A+tol,1)-tol;
            
        end

        function [C] = rnd_(A,tol)
            if nargin < 2; tol = am_lib.tiny; end
            C = round(A,-log10(tol));
        end

        function [C] = rndstr_(n)
            set = char(['a':'z','0':'9','_','!','@','#','$','%','^']); nsets = numel(set);
            C = set(ceil(nsets*rand(1,n)));
        end

        function [L] = eq_(a,b,tol)
            if nargin < 3; tol = am_lib.tiny; end
            L = abs(a-b)<tol;
        end
        
        function [L] = lt_(a,b,tol)
            if nargin < 3; tol = am_lib.tiny; end
            L = a<b-tol;
        end
        
        function [L] = gt_(a,b,tol)
            if nargin < 3; tol = am_lib.tiny; end
            L = a>b+tol;
        end
        
        function [x] = wdv_(x,tol)
            % get well-defined value for integers, sqrt(integer), cube-root(integer)
            
            import am_lib.*

            if nargin<2; tol = am_lib.eps; end

            if isreal(x)
                
                for i = 1:numel(x)
                    go=true;
                    % reduce to integer, sqrt, cube-root
                    for j = 1:3; if go
                        if eq_(mod_(x(i).^j),0,tol)
                            x(i) = sign(x(i)) .* round(abs(x(i).^j)).^(1/j); 
                            go=false; break;
                        end
                    end; end
                    % well defined values for cosd(30) = sqrt(3)/2, cosd(60) = 1/2
                    for wdv = [1/2,sqrt(3)/2]; if go
                        if eq_(abs(x(i)),wdv,tol)
                            x(i) = wdv * sign(x(i)); 
                            go=false; break; 
                        end
                    end;end
                end
                
                for i = 1:numel(x)
                    for wdv = [1/2,sqrt(3)/2]
                    if eq_(abs(x(i)),wdv,tol); x(i) = wdv * sign(x(i)); go=false; break; end
                    end
                end
            else
                x = wdv_(real(x)) + wdv_(imag(x))*1i;
            end
        end
        

        % vectorization
        
        function [A] = aug_(A,n)
            % squeeze in 1's in vector A at positions n
            ex_ = any([1:(numel(n)+numel(A))]==n(:),1);
            A(~ex_) = A; 
            A( ex_) = 1;
        end
        
        function [A] = ext_(A,n)
            % remove positions n of vector A
            ex_ = any(1:numel(A)==n(:),1);
            A = A(~ex_);
        end
        
        function [C] = matmul_(A,B,applysqueeze)
            % matrix multiple the first two dimensions and entry-wise the rest
            %
            %       A [(m,n),(a,1,c,d,1,f,...)]   
            %     * B [(n,p),(1,b,c,1,e,f,...)] 
            %   ` = C [(m,p),(a,b,c,d,e,f,...)]
            %
            
            if nargin == 2; applysqueeze = true; end

            aug_ = @(x,i,y) [x(1:(i-1)),y,x(i:end)]; 
            C = sum( reshape(A,aug_(size(A),3,1)) .* reshape(B,aug_(size(B),1,1)) , 2);
            
            if applysqueeze; C = squeeze(C); else
                n = size(C); n(2) = [];
                C = reshape(C,n);
            end
        end

        function [C] = matmulp_(A,B,applysqueeze)
            % high dimensional matrix multiplication: 
            %
            %       A [(m,n),(a,b)]
            %     * B [(n,p),(b,c)]
            %     = C [(m,p),(a,c)]
            %
            % Explicit:
            %     A = rand(4,3,3,1); B = rand(3,2,1,5);
            %     [a{1:4}]=size(A); [b{1:4}]=size(B); C=zeros(a{1},b{2},a{3},b{4});
            %     for b4 = 1:b{4}; for b2 = 1:b{2}
            %     for a4 = 1:a{4}; for a3 = 1:a{3}; for a2 = 1:a{2}; for a1 = 1:a{1}
            %         C(a1,b2,a3,b4) = C(a1,b2,a3,b4) + A(a1,a2,a3,a4)*B(a2,b2,a4,b4);
            %     end; end; end; end
            %     end; end
            %     C - matmulp_(A,B)
            %

            import am_lib.sum_
            
            C = permute(sum_( permute(A,[1,2,6,3,4,5]) ...
                           .* permute(B,[6,1,2,5,3,4]), [2,5]),[1,3,4,6,2,5]);
            
            if nargin == 2; applysqueeze = true; end
            if applysqueeze; C = squeeze(C); end
        end

        function [C] = flatten_(A)
            C = A(:);
        end

        function [C] = kron_(A,B)
            % define high dimensional kronecker self-product: for i=[1:size(A,3)]; C(:,:,i)=kron(A(:,:,i),B(:,:,i)); end
            C = reshape( permute(A,[4 1 5 2 3]) ...
                      .* permute(B,[1 4 2 5 3]), size(A,1)*size(B,1),[],size(A,3));
        end
        
        function [C] = kronpow_(A,n)
            
            import am_lib.kron_
            
            C = A;
            for i = 1:(n-1)
                C = kron_(C,A);
            end
        end

        function [C] = blkdiag_(A,B)
            n = size(A); m = size(B);
            C = zeros(n(1)+m(1),n(2)+m(2),n(3));
            for i = 1:n(3)
                C(:,:,i) = blkdiag(A(:,:,i),B(:,:,i));
            end            
        end
        
        function [C] = osum_(A,B,i)
            % define outer sum of two vector arrays
            % 
            % i = 1: add first dimension
            %       A [ m,a,b,1,d,...]
            %     + B [ p,a,b,c,1,...]
            %   ` = C [(m,p),a,b,c,d,...]
            %
            % i = 2: add second dimension
            %       A [m,n]
            %     + B [m,p]
            %     = C [m,(n,p)]
            %
            
            n=size(A); m=size(B);
            
            if     i == 2
                C = reshape(A,n(1),n(2),1) + reshape(B,m(1),1,m(2));
                C = reshape(C,n(1),n(2)*m(2));
            elseif i == 1
                C = reshape(A,[1,n(1),n(2:end)]) + reshape(B,[m(1),1,m(2:end)]);
                C = reshape(C,[n(1)*m(1),max([n(2:end);m(2:end)])]);
            end
        end

        function [D] = outerc_(A,B,C)
            % outer product of each column of A against each column of B (and C)
            %
            % Explicit (double):
            %     A=rand(3,3);B=rand(3,3); 
            %     D=zeros(size(A,1),size(B,1),size(B,2));
            %     for i = 1:size(A,2)
            %     for m = 1:size(A,1); for n = 1:size(B,1)
            %         D(m,n,i) = D(m,n,i) + A(m,i)*B(n,i); 
            %     end; end
            %     end
            %     outerc_(A,B) - D
            %
            % Explicit (triple):
            %     A=rand(4,5);B=rand(2,5);C=rand(6,5); 
            %     D=zeros(size(A,1),size(B,1),size(C,1),size(B,2));
            %     for i = 1:size(A,2)
            %     for m = 1:size(A,1); for n = 1:size(B,1); for o = 1:size(C,1)
            %         D(m,n,o,i) = D(m,n,o,i) + A(m,i)*B(n,i)*C(o,i); 
            %     end; end; end
            %     end
            %     outerc_(A,B,C) - D
            %
            % Benchmark:
            %     x = [1,2,5,10,20,50,100,200,500,1000,2E4,5E4,1E5,2E5,5E5]; k = max(x)./x;
            %     t1 = zeros(1,numel(x));t2 = zeros(1,numel(x));
            %     for ix = 1:numel(x)
            %         % benchmark with single precision
            %         A=single(rand(3,x(ix)));B=single(rand(3,x(ix)));C=single(rand(3,x(ix))); 
            %         for ik = 1:numel(k)
            %             % explicit
            %             tic;
            %                 D=zeros(size(A,1),size(B,1),size(C,1),size(B,2));
            %                 for i = 1:size(C,2)
            %                 for m = 1:size(A,1); for n = 1:size(B,1); for o = 1:size(C,1)
            %                     D(m,n,o,i) = D(m,n,o,i) + A(m,i)*B(n,i)*C(o,i); 
            %                 end; end; end
            %                 end
            %             t1(ix) = t1(ix) + toc/k(ik);
            %             % vectorized
            %             tic;
            %                 D2=outerc_(A,B,C);
            %             t2(ix) = t2(ix) + toc/k(ik);
            %         end
            %     end
            %     loglog(x,t1,'s-',x,t2,'o:'); legend('explicit','vectorized')
            
            import am_lib.*
            
            if     nargin == 2
                % outer product of two vectors
                % A [ n, 1, m ] * B [ 1, n, m ] = D [ n, n, m ]; [n,m] = size(A)
                D = permute(A,[1,3,2]) ...
                 .* permute(B,[3,1,2]);
            elseif nargin == 3
                % outer product of three vectors
                D = permute(A,[1,3,4,2]) ...
                 .* permute(B,[3,1,4,2]) ...
                 .* permute(C,[3,4,1,2]);
            end
        end
        
        function [C] = tdp_inner_(A,u,w,v)
            % inner triple dot product
            %
            % Explicit:
            %     rng(1); A = rand(3,4,5,6); u = rand(3,6); w = rand(4,6); v = rand(5,6);
            %     % apply triple dot product
            %     n = size(A); B = zeros(n(4),1);
            %     for m = 1:n(4)
            %     for i = 1:n(1); for j = 1:n(2); for k = 1:n(3)
            %         B(m) = B(m) + A(i,j,k,m)*u(i,m)*w(j,m)*v(k,m);
            %     end; end; end
            %     end
            %     B - tdp_(A,u,w,v)
            
            import am_lib.*
            
            if     isempty(v)
                % inner double dot prouct
                C = squeeze(sum_( A .* permute(outerc_(u,w),[1,2,4,3])   , [1,2]  ));
            elseif isempty(w)
                % inner double dot prouct
                C = squeeze(sum_( A .* permute(outerc_(u,v),[1,4,2,3])   , [1,3]  ));
            elseif isempty(u)
                % inner double dot prouct
                C = squeeze(sum_( A .* permute(outerc_(v,w),[4,1,2,3])   , [2,3]  ));
            else
                % inner triple dot product
                C = squeeze(sum_( A .* outerc_(u,w,v) , [1,2,3]));
            end
        end
        
        function [C] = tdp_outer_(A,u,w,v)
            % outer triple dot product (omitting one vector will produce a double dot product)
            
            import am_lib.*
            
            % get sizes
            n = size(A); n(5) = size(u,2);
            
            if     isempty(v)
                % Explicit:
                %     rng(1); A = rand(3,3,5,10); u = rand(3,6); w = rand(3,6); v = rand(5,6);
                %     % apply triple dot product
                %     n = size(A); n(5) = size(u,2); B = zeros(n(3),n(4),n(5));
                %     for o = 1:n(4); for p = 1:n(5)
                %     for i = 1:n(1); for j = 1:n(2); for k = 1:n(3)
                %         B(k,o,p) = B(k,o,p) + A(i,j,k,o)*u(i,p)*w(j,p);
                %     end; end; end
                %     end; end
                %     B - tdp_outer_(A,u,w,[])
                C = squeeze(sum_(reshape(outerc_(u,w),n(1),n(2),1,1,n(5)) .* A,[1,2]));
            elseif isempty(w)
                % Explicit:
                %     rng(1); A = rand(3,7,5,10); u = rand(3,6); w = rand(7,6); v = rand(5,6);
                %     % apply triple dot product
                %     n = size(A); n(5) = size(u,2); B = zeros(n(2),n(4),n(5));
                %     for o = 1:n(4); for p = 1:n(5)
                %     for i = 1:n(1); for j = 1:n(2); for k = 1:n(3)
                %         % B(k,o,p) = B(k,o,p) + A(i,j,k,o)*u(i,p)*w(j,p)*v(k,p);
                %         B(j,o,p) = B(j,o,p) + A(i,j,k,o)*u(i,p)*v(k,p);
                %     end; end; end
                %     end; end
                %     B - tdp_outer_(A,u,[],v)
                C = squeeze(sum_(reshape(outerc_(u,v),n(1),1,n(3),1,n(5)) .* A,[1,3]));
            elseif isempty(u)
                % Explicit:
                %     rng(1); A = rand(7,3,5,10); u = rand(7,6); w = rand(3,6); v = rand(5,6);
                %     % apply triple dot product
                %     n = size(A); n(5) = size(u,2); B = zeros(n(1),n(4),n(5));
                %     for o = 1:n(4); for p = 1:n(5)
                %     for i = 1:n(1); for j = 1:n(2); for k = 1:n(3)
                %         B(i,o,p) = B(i,o,p) + A(i,j,k,o)*w(j,p)*v(k,p);
                %     end; end; end
                %     end; end
                %     B - tdp_outer_(A,[],w,v)
                C = squeeze(sum_(reshape(outerc_(w,v),1,n(2),n(3),1,n(5)) .* A,[2,3]));
            else
                % Explicit: 
                %     rng(1); A = rand(3,3,5,10); u = rand(3,6); w = rand(3,6); v = rand(5,6);
                %     % apply triple dot product
                %     n = size(A); n(5) = size(u,2); B = zeros(n(4),n(5));
                %     for o = 1:n(4); for p = 1:n(5)
                %     for i = 1:n(1); for j = 1:n(2); for k = 1:n(3)
                %         B(o,p) = B(o,p) + A(i,j,k,o)*u(i,p)*w(j,p)*v(k,p);
                %     end; end; end
                %     end; end
                %     B - tdp_outer_(A,u,w,v)
                C = reshape(A,prod(n(1:3)),n(4)).' * reshape(outerc_(u,w,v),prod(n(1:3)),n(5));
            end
        end
        
        function [C] = normc_(A)
            % get length of each column vector
            C = sqrt(sum(abs(A).^2,1));
        end

        function [C] = operm_(A,I)
            % outer permutation operation
            % explicit: 
            %     for i = 1:size(I,2)
            %     for j = 1:size(A,2)
            %         D(:,i,j) = A(I(:,i),j);
            %     end
            %     end
            %     C-D
            
            import am_lib.*
            
            C = reshape(accessc_(repmat(A,1,size(I,2)),repelem(I,1,size(A,2))),size(I,1),size(I,2),size(A,2));
        end
        
        function [C] = accessc_(A,I)
            % permute each column of A according to the indicie matrix I
            % for example: A=randi(10,5,5); [B,I]=sort(A); B-accessc_(A,I)
            % Explicit: 
            % for i = 1:size(A,2)
            %   C(:,i) = C(I(:,i),i);
            % end
            C = A(bsxfun(@plus,I,[0:size(A,2)-1]*size(A,1)));
        end
        
        function [C] = findrow_(A)
            % define function to find the first nonzero value in each row of matrix A
            % returns 0 for rows containing all zeros
            C = (sum(cumsum(A~=0,2)==0,2)+1);
            C = C .* ~all(A==0,2);
        end
       
        function [A] = sum_(A,n,varargin)
            % sum over dimensions n: sum_(A,[2,3])
            for i = 1:numel(n)
                A = sum(A,n(i),varargin{:});
            end
        end

        function B   = changem_(A,newval,oldval)
            B = A;
            [valid,id] = max(bsxfun(@eq,A(:),oldval(:).'),[],2); %//'
            B(valid) = newval(id(valid));
        end


        % matching
        
        function [C] = maxabs_(A)
            C = max(abs(A(:)));
        end
        
        function [C] = minmax_(A)
            C = [min(A(:)),max(A(:))];
        end
        
        function [C] = sortc_(A, tol)
            % column vector-based rank, sort, unique with numeric precision
            import am_lib.rnd_
            
            % set default numerical tolernece
            if nargin<2 || isempty(tol); tol = am_lib.tiny; end
            
            [C] = sortrows(rnd_(A,tol).').'; 
        end

        function [C] = rankc_(A, tol)
            % column vector-based rank, sort, unique with numeric precision
            import am_lib.rnd_
             
            % set default numerical tolernece
            if nargin<2 || isempty(tol); tol = am_lib.tiny; end
            
            C = sortrowsc(rnd_(A,tol).',[1:size(A,1)]).'; 
        end                 

        function [I] = match_(A,B)
            % find the indicies I which permute column entries of A to match B
            % for example: A=[[1;4;3],[1;3;4],[3;1;4]]; B=[3;1;4]; X=match_(A,B); accessc_(A,X)-B
            import am_lib.*
            [~,iA]=sort(A); [~,iB]=sort(B); iB(iB)=[1:numel(iB)]; I = accessc_(iA,iB);
        end
        
        function [c] = member_(A,B,tol)
            % get indicies of column vectors A(:,i,j) in matrix B(:,:)
            % 0 means not A(:,i) is not in B(:,:)
            
            % set default numerical tolernece
            if nargin<3 || isempty(tol); tol = am_lib.tiny; end
            
            [~,m,n] = size(A); c = zeros(m,n);
            for i = 1:m
            for j = 1:n
                c(i,j) = member_engine_(A(:,i,j),B,tol);
            end
            end

            function c = member_engine_(A,B,tol)
                c = 1; r = 1; [d1,d2] = size(B);
                while r <= d1
                    if abs(B(r,c)-A(r))<tol
                        r = r+1;
                    else
                        c = c+1;
                        r = 1;
                    end
                    if c>d2; c=0; return; end
                end
            end
        end
        
        function [C,IA,IC] = uniquec_(A,tol)
            % get unique values with numeric precision
            import am_lib.rnd_
            if nargin<2 || isempty(tol)
                [~,IA,IC] = unique(rnd_(A).'    , 'rows', 'stable'); C=A(:,IA);
            else
                [~,IA,IC] = unique(rnd_(A,tol).', 'rows', 'stable'); C=A(:,IA);
            end
        end
        
        function [C] = uniquemask_(A,tol)
            % returns a matrix with column vectors marking the first
            % occurance of a value in each column
            import am_lib.rnd_
           
            % set default numerical tolernece
            if nargin<2 || isempty(tol); tol = am_lib.tiny; end
            
            C = false(size(A));
            for i = 1:size(A,2)
                [~,b]=unique(rnd_(A(:,i),tol)); C(b,i) = true;
            end
        end

        function c_id = reindex_using_occurances(c_id)
            % relabel based on repetitions 
            import am_lib.*
            % count repetitions/occurances
            occurances = sum(c_id==c_id.',2);
            % sort, lift degeneracies if present using original labeling
            fwd = rankc_([occurances(:).';c_id(:).']);
            % relabel and return indicies to original order
            c_id(fwd) = cumsum([0;diff(c_id(fwd))~=0])+1; 
        end
        
        function [X] = strmatch_(A,B)
            X = zeros(size(A)); 
            for k = 1:numel(B)
                X(strncmp(A,B{k},numel(B{k}))) = k;
            end
        end
        
        function [X] = strmatchi_(A,B)
            % case insensitive match
            X = zeros(size(A)); 
            for k = 1:numel(B)
                X(strncmpi(A,B{k},numel(B{k}))) = k;
            end
        end

        function [A,i2p,p2i] = get_connectivity(PM)

            import am_lib.*

            % exclude rows with all zeros
            PM=PM(~all(PM==0,2),:);
            
            % build connectivity matrix
            for i1 = 1:size(PM,1)
                X=unique(PM(i1,:)); 
                for i2 = 1:numel(X)
                    A(X(i2),X) = 1; 
                end
            end
            
            % merge overlapping rows
            A = merge_(A);

            i2p = round(findrow_(A)).'; p2i = round(([1:size(A,1)]*A));
            
            function [A] = merge_(A)
                % convert to logical
                A = ~am_lib.eq_(A,0);
                
                % loop until number of rows in A stop changing
                m=1; m_last=0;
                while m ~= m_last
                    m_last = m; [m,n] = size(A); i = 1; j = 1; 
                    while (i <= m) && (j <= n)
                        % Find value and index of largest element in the remainder of column j.
                        [p,k] = max(abs(A(i:m,j))); k = k+i-1;
                        if (p == 0)
                           j = j + 1;
                        else
                            % Swap i-th and k-th rows.
                            A([i k],j:n) = A([k i],j:n);
                            % see which rows overlap with the i-th row
                            ex_ = any(and(A(i,j:n), A(:,j:n)),2);
                            % merge overlaps 
                            A(i,j:n) = any(A(ex_,j:n),1);
                            % zero all other rows
                            ex_(i)=false; A(ex_,j:n)=false;
                            i = i + 1;
                            j = j + 1;
                        end
                    end
                    
                    % exclude rows with all zeros
                    A=A(any(A,2),:);
                end
            end 
        end
        
        % file parsing
        
        function t   = extract_token_(str,token,numeric)
            import am_lib.*
            pos = find(strncmp(str,token,numel(token)),1);
            if ~isempty(pos)
                t = strtrim(regexprep(str{pos(1)},[token '|'''],''));
                if ~isempty(t)
                    if nargin>2 && numeric
                        t = sscanf(t,'%f');
                    end
                end
            else
                t = '';
            end
        end

        function [str,nlines] = load_file_(fname)
            import am_lib.count_lines_
            % get number of lines
            nlines = count_lines_(fname); str = cell(nlines,1);
            % reads a file rowise into a cellstr
            fid = fopen(fname,'r'); if fid==-1; error('Failed to open %s. Check that file exists.\n',fname); end
                for j = 1:nlines
                    str{j} = strtrim(fgetl(fid));
                end
            fclose(fid);
        end
        
        
        % special mathematical functions 
        
        function [y] = lorentz_(x)
            % define gaussian function 
            y = 1./(pi*(x.^2+1)); 
        end

        function [y] = gauss_(x)
            y = exp(-abs(x).^2)./sqrt(pi);
        end

        function [y] = delta_(x)
            % define tiny, kronecker delta, and heavside
            y = logical(abs(x)<am_lib.tiny); 
        end

        function [y] = heaviside_(x)
            y = logical(x>0);
        end
                
        function [y] = pvoigt_(x,f)
            % amplitude normalized
            import am_lib.*
            y = (1-f) .* gauss_(x./sqrt(2*log(2))) .* sqrt(pi) + f .* lorentz_(x) .* pi;
        end
        
        function [y] = sinc_(x)
            i=find(x==0); x(i)= 1; y = sin(pi*x)./(pi*x); y(i) = 1;
        end
        
        function [w] = tukeyw_(n,r)
            if nargin == 1; r = 0.5; end
            t = linspace(0,1,n)'; per = r/2; tl = floor(per*(n-1))+1; th = n-tl+1;
            w = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];
        end
        
        function [w] = gaussw_(n,r)
            % [w] = gaussw_(n,r)
            N = n-1; n = (0:N)'-N/2; w = exp(-(1/2)*(r*n/(N/2)).^2);
        end
        
        
        % geometric functions
        
        function R      = rotzd_(t)
            c = cosd(t); s = sind(t);
            R = [c -s 0;
                 s  c 0;
                 0  0 1];
        end
        
        function R      = rotyd_(t)
            c = cosd(t); s = sind(t);
            R = [ c 0 s; 
                  0 1 0; 
                 -s 0 c];
        end
        
        function R      = rotxd_(t)
            c = cosd(t); s = sind(t);
            R = [1 0  0; 
                 0 c -s;
                 0 s  c];
        end
        
        function [A]    = R_axis_(R)
            % [A] = R_axis_(R);
            
            import am_lib.*
            
            nRs = size(R,3);
            
            if nRs == 1
                % convert (im)proper rotation to proper rotation
                R = R*sign(det(R));

                % define basic parameters and functions
                tol = 1E-8; normalize_ = @(v) v/norm(v); 

                % check for identity
                if abs(trace(R)-3)< tol; A=[0;0;1]; return; end

                % get rotation axis
                A = null(R-eye(3));

                % get random point on plane perpendicular to the rotation axis
                v1 = rand(3,1); v1 = normalize_(v1 - dot(v1,A)*A);

                % rotate point on the perpendicular plane
                v2 = R*v1;

                % get cross product (add tiny cross component to deal with 180 deg rotation)
                c = normalize_(cross(v1,v2+cross(v1,A)*tol));

                % adjust sign
                A = sign(dot(c,A))*A;
            else
                
                A = zeros(3,nRs);
                for i = 1:nRs
                    A(:,i) = R_axis_(R(:,:,i));
                end
            end
        end

        function [A]    = R_angle_(R)
            % define conversion of proper rotations to axis & angle representation
            A = acos((trace(R)-1)/2); 
        end
        
        function R      = rot_align_(A,B)
            %
            % Rotation matrix R which aligns unit vector A to unit vector B: B = rotmat(A,B) * A
            %
            % based on Rodrigues' Rotation Formula
            % "A Mathematical Introduction to Robotic Manipulation", Richard M. Murray, Zexiang Li, S. Shankar Sastry, pp. 26-28
            %
            % % quick matlab implementation:
            % A = rand(3,1); A=A./norm(A);
            % B = rand(3,1); B=B./norm(A);
            % v = cross(A,B);
            % s = norm(v);
            % c = dot(A,B);
            % vx(1:3,1) = [0.0,v(3),-v(2)];
            % vx(1:3,2) = [-v(3),0.0,v(1)];
            % vx(1:3,3) = [v(2),-v(1),0.0];
            % R = eye_1(3) + vx + (vx*vx)*(1.0-c)/(s.^2)
            % R*A-B
            %
            A=A./norm(A); B=B./norm(B);
            v = cross(A,B); s = norm(v); c = dot(A,B);
            vx(1:3,1) = [0.0,v(3),-v(2)];
            vx(1:3,2) = [-v(3),0.0,v(1)];
            vx(1:3,3) = [v(2),-v(1),0.0];
            %
            R = eye(3) + vx + vx*vx*(1.0 - c)/(s.^2);
            %
        end

        function v_rot  = rot_vec_(v,k,theta)
            [m,n] = size(v);theta=theta/180*pi;
            if (m ~= 3 && n ~= 3)
                error('input vector is/are not three dimensional'), end
            if (size(v) ~= size(k)) 
                error('rotation vector v and axis k have different dimensions'),end

            k = k/sqrt(k(1)^2 + k(2)^2 + k(3)^2); % normalize rotation axis
            No = numel(v)/3; % number of vectors in array
            v_rot = v; % initialize rotated vector array
            if ( n == 3 )
                crosskv = v(1,:); % initialize cross product k and v with right dim.
                for i = 1:No
                    crosskv(1) = k(2)*v(i,3) - k(3)*v(i,2);
                    crosskv(2) = k(3)*v(i,1) - k(1)*v(i,3); 
                    crosskv(3) = k(1)*v(i,2) - k(2)*v(i,1);
                    v_rot(i,:) = cos(theta)*v(i,:) + (crosskv)*sin(theta)...
                                    + k*(dot(k,v(i,:)))*(1 - cos(theta));
                end
            else % if m == 3 && n ~= 3
                crosskv = v(:,1); % initialize cross product k and v with right dim.
                for i = 1:No
                    crosskv(1) = k(2)*v(3,i) - k(3)*v(2,i);
                    crosskv(2) = k(3)*v(1,i) - k(1)*v(3,i); 
                    crosskv(3) = k(1)*v(2,i) - k(2)*v(1,i);
                    v_rot(:,i) = cos(theta)*v(:,i) + (crosskv)*sin(theta)...
                                    + k*(dot(k,v(:,i)))*(1 - cos(theta));
                end
            end
        end
        
        function [Wtes] = get_wigner(j,R,flag)
            import am_lib.*
            % defaults
            if nargin < 3
                % this will probably cause problems later... if within one hamiltnoian, both
                % tesseral and spherical are used. should stick to spherical some how. but
                % presently, complex values do not work for getting tight binding matrix elements.
                if eq_(mod_(j),0)
                    % j = integer
                    flag='tesseral';
                else
                    % j = half-integer
                    flag='spherical';
                end
            end

            % define tiny, kronecker delta, and heavside
            d_ = @(x,y) logical(x==y); t_ = @(x,y) logical(x>y);
            
            % matrix indices
            [m,mp]=meshgrid([j:-1:-j]);
            
            % define angular momentum operators (Jp raising, Jm lowering, ...)
            Jm = d_(m,mp+1).*sqrt((j+m).*(j-m+1)); Jp = d_(m,mp-1).*sqrt((j-m).*(j+m+1));
            Jx = (Jp+Jm)/2; Jy = (Jp-Jm)/2i; Jz = d_(m,mp).*m; J = cat(3,Jx,Jy,Jz);

            if     contains(flag,'tesseral')
                % define basis change: spherical (complex) to tesseral harmonics (real basis)
                % C. Görller-Walrand and K. Binnemans, in (Elsevier, 1996), pp. 145, eqs. 14-16.
                T = d_(0,m) .* d_(mp,m) + ...
                    t_(m,0) .* sqrt(-1/2) .* ( d_(m,-mp) - (-1).^m.*d_(m, mp) ) + ... % sine   terms
                    t_(0,m) .* sqrt( 1/2) .* ( d_(m, mp) + (-1).^m.*d_(m,-mp) );      % cosine terms

                % if j is a half-integer, tesseral harmonics make no sense, force spherical
                if mod(j*2,1)==1
                    T=1; warning('Tesseral harmonics only make sense for j = integer. Switching to spherical harmonics.'); 
                end
            elseif contains(flag,'spherical')
                % keep it in spherical harmonics, complex basis
                T = 1;
            end

            % batch convert to wigner
            nRs = size(R,3); Wtes = zeros(2*j+1,2*j+1,nRs);
            for i = [1:nRs]; Wtes(:,:,i) = get_wigner_engine(J,T,j,R(:,:,i)); end
            
            
            % engine
            function [Wtes] = get_wigner_engine(J,T,j,R)
                % define wigner function for spherical (complex) and tesseral (real) harmonics
                % Note: for l = 1, Wtes_(R) = R

                % get proper rotation
                d = sign(det(R)); dR = R*d;

                % get rotation axis and angle (the circle shift is required to recover SO(3) rotations)
                an = am_lib.R_angle_(dR); ax = circshift(am_lib.R_axis_(dR),1);

                % define spin-vector dot products [Eq. 1.37 Blundell]
                dotV_ = @(S,V) S(:,:,1)*V(1) + S(:,:,2)*V(2) + S(:,:,3)*V(3);

                if size(R,1)>9 % faster for symbolic and square matrices with dimensions > 9
                    Wsph = expm( -sqrt(-1) * dotV_(J,ax) * an) * d.^j;
                    Wtes = T' * Wsph * T;
                else           % faster for numerical square matrices with dimensions < 9
                    [V,D] = eig( -sqrt(-1) * dotV_(J,ax) * an); 
                    Wsph = V*diag(exp(diag(D)))/V * d.^j;
                    Wtes = T' * Wsph * T ;
                end
            end
        end
        
        function [x,y,z]= sph2cartd_(phi,chi,r)
            z = r .* sind(chi);
            rcoselev = r .* cosd(chi);
            x = rcoselev .* cosd(phi);
            y = rcoselev .* sind(phi);
        end
        
        function [x,y]  = pol2cartd_(th,r)
            x = r.*cosd(th);
            y = r.*sind(th);
        end
        
        function [th,r] = cart2pold_(x,y)
            th = atan2d(y,x);
            r = hypot(x,y);
        end
        
        function [phi,chi,r] = cart2sphd_(x,y,z)
            hypotxy = hypot(x,y); r = hypot(hypotxy,z);
            chi = atan2d(z,hypotxy); phi = atan2d(y,x);
        end
        
        function [C]    = round_(C,tol)
            % this should be essentially the opposite of x./norm(x) for an integer array x
            
            import am_lib.*
            
            if nargin<2; tol=am_lib.tiny; end
            
            s = size(C);
            
            switch ndims_(C)
                case 1
                    % try to convert to all real values if any is imaginary
                    if any(~eq_(imag(C),0)); C = C./1i; end

                    % if imaginary part is small
                    if all( eq_(imag(C),0,tol) )
                        % convert to real
                        C = real(C); 
                        % round values by diving by smallest noninteger
                        for i = 1:100
                            if ~any(~eq_(mod_(C(:),tol),0,tol)); break; end 
                            small=min(C(~eq_(C(:),0,tol))); C=C./small*sign(small); 
                        end
                        C = round(C);
                    end
                    
                case 2
                    for i = 1:s(2)
                        X=C(:,i); C(:,i)=reshape(round_(X(:)),s(1),1);
                    end
                    
                case 3
                    for i = 1:s(3)
                        X=C(:,:,i); C(:,:,i)=reshape(round_(X(:)),s(1),s(2));
                    end
            end
        end
        
        
        % matrix related
        
        function x   = ndims_(A)
            % number of dimensions after squeezing
            x = max(sum(size(A)~=1),1);
        end
        
        function [A] = frref_(A)
            %frref_   Fast reduced row echelon form.
            %   R = frref_(A) produces the reduced row echelon form of A.
            % 
            %   Description: 
            %   For full matrices, the algorithm is based on the vectorization of MATLAB's
            %   RREF function. A typical speed-up range is about 2-4 times of 
            %   the MATLAB's RREF function. However, the actual speed-up depends on the 
            %   size of A. The speed-up is quite considerable if the number of columns in
            %   A is considerably larger than the number of its rows or when A is not dense.
            %
            %   For sparse matrices, the algorithm ignores the tol value and uses sparse
            %   QR to compute the frref_ form, improving the speed by a few orders of 
            %   magnitude.
            %
            %   Authors: Armin Ataei-Esfahani (2008)
            %            Ashish Myles (2012)
            %
            %   Revisions:
            %   25-Sep-2008   Created Function
            %   21-Nov-2012   Added faster algorithm for sparse matrices

            % set size and default tolerences
            [m,n] = size(A); tol=max(m,n)*eps(class(A))*norm(A,'inf');

            % Loop over the entire matrix.
            i = 1; j = 1;
            % t1 = clock;
            while (i <= m) && (j <= n)
               % Find value and index of largest element in the remainder of column j.
               [p,k] = max(abs(A(i:m,j))); k = k+i-1;
               if (p <= tol)
                  % The column is negligible, zero it out.
                  A(i:m,j) = 0; %(faster for sparse) %zeros(m-i+1,1);
                  j = j + 1;
               else
                  % Swap i-th and k-th rows.
                  A([i k],j:n) = A([k i],j:n);
                  % Divide the pivot row by the pivot element.
                  Ai = A(i,j:n)/A(i,j);    
                  % Subtract multiples of the pivot row from all the other rows.
                  A(:,j:n) = A(:,j:n) - A(:,j)*Ai;
                  A(i,j:n) = Ai;
                  i = i + 1;
                  j = j + 1;
               end
            end
        end

        function [F] = fsrref_(A)
            % fast sparse rref based on QR decomposition
            [m,n] = size(A);

            % Non-pivoted Q-less QR decomposition computed by Matlab actually
            % produces the right structure (similar to rref) to identify independent
            % columns.
            R = qr(A);
            % i_dep = pivot columns = dependent variables
            %       = left-most non-zero column (if any) in each row
            % indep_rows (binary vector) = non-zero rows of R
            [indep_rows, i_dep] = max(R ~= 0, [], 2);
            indep_rows = full(indep_rows); % probably more efficient
            i_dep = i_dep(indep_rows);
            i_indep = setdiff(1:n, i_dep);

            % solve R(indep_rows, i_dep) x = R(indep_rows, i_indep)
            %   to eliminate all the i_dep columns
            %   (i.e. we want F(indep_rows, i_dep) = Identity)
            F = sparse([],[],[], m, n);
            F(indep_rows, i_indep) = R(indep_rows, i_dep) \ R(indep_rows, i_indep);
            F(indep_rows, i_dep) = speye(length(i_dep));
        end
        
        function [A] = force_hermiticity_(A)
            A = (A'+A)/2;
        end
        
        function [A] = force_symmetricity_(A)
            A = (A.'+A)/2;
        end
        
        function [a] = trace_(A)
            [m]=size(A,3);[n]=size(A,4);
            a = zeros(m,n);
            for i = 1:m; for j = 1:n
                a(i,j) = trace(A(:,:,i,j));
            end; end
               
        end
        
        function [U] = orth_(U,E)
            % orthonormalize eigencolumns of A within each degenerate
            % subspace using eigenvalues E 
            
            import am_lib.*
            
            if nargin == 2
                % [~,~,a]=unique(rnd_(E));
                a = cumsum([0;diff(sort(rnd_(E(:))))~=0])+1;
                for j = 1:max(a)
                    U(:,a==j)=orth(U(:,a==j));
                end
            else 
                U = orth(U);
            end
        end
        
        function [c] = diagonalicity_(A)
            c = max(max(abs(diag(diag(A))-A)));
        end
        
        function [c] = hermiticity_(A)
            c = max(max(abs(A'-A)));
        end
        
        function [c] = unitaricity_(A)
            c = max(max(abs( A'*A - eye(size(A)) )));
        end
        
        function [c] = symmetricity_(A)
            c = max(max(abs(A-A.')));
        end
        
        function [U,E] = eig_(D)
            % diagonalize and get nice eigenvectors for Hermitian or symmetric matrices D
            % confirm input matrix is hermitian
            
            import am_lib.*
            
            % Assure that D is Hermitian or real-symmetric.
%             assert(hermiticity_(D)<am_lib.eps,'D is not Hermitian.');
%             force_hermiticity_(D)


            % get eigenvectors
            [U,E]=eig(D,'vector');
            
            % maximize vectors 1-norm within degenerate subspaces
            [~,~,unq_]=unique(rnd_(E));
            for j = 1:max(unq_)
                unq_list = find(unq_==j).';
                if numel(unq_list)>1
                for a = unq_list
                    ex_ = unq_list(unq_list~=a);
                    [m,p]=max(abs(U(:,a)));
                    if m > am_lib.eps
                        U(:,ex_) = U(:,ex_) - U(:,a) .* U(p,ex_)/U(p,a);
                    else
                        U(:,a) = 0;
                    end
                end
                end
            end

            % nullify numerically negligible components
            ex_ = abs(U(:))<am_lib.eps; U(ex_) = 0;

            % rotate each column vector in complex space to make the first value be real
            th = angle(U); [~,I]=max(abs(th)); U = U .* exp(-1i*accessc_(th,I));

            % nullify numerically negligible components
            ex_ = abs(imag(U(:)))<am_lib.eps; U(ex_) = real(U(ex_));
            ex_ = abs(real(U(:)))<am_lib.eps; U(ex_) = imag(U(ex_));

            % normalize each colum vector
            U = U./normc_(U);

            % confirm properties of new eigenvectors
%             assert(unitaricity_(U)<am_lib.eps, 'U is not unitary.');
            assert(any(~abs(diag(U'*D*U)-E)<am_lib.eps), 'E is mismatched.');
        end
        
        
        % fitting functions
        
        function [x,r] = lsqnonlin_(cost_,x0,isfixed,varargin)
            % X = lsqnonlin(FUN,X0,LB,UB,OPTIONS)
            % z=[1:500]/250; y=z.^2.3+5.5;
            % isfixed = logical([0 0]); x0=[0 2];
            % fnct_ = @(x) x(1)+z.^x(2);
            % cost_ = @(x) abs(fnct_(x)-y);
            % opts = optimoptions('lsqnonlin','Display','none','MaxIter',7);
            % [x,r]=lsqnonlin_(cost_,x0,isfixed,[],[],opts)
            % plot(z,y,'.',z,fnct_(x),'-')
            
            % fixed and loose indicies
            f = find( isfixed); xf = x0(f);
            l = find(~isfixed); xl = x0(l);

            % Estimate only the non-fixed ones
            [xl,r] = lsqnonlin(@localcost_,xl,varargin{:});

            % Re-create array combining fixed and estimated coefficients
            x([f,l]) = [xf,xl];

            function y = localcost_(x)
               b([f,l]) = [xf,x]; y = cost_(b);
            end
        end

        function [x,r] = ga_(cost_,x0,isfixed,varargin)

            % fixed and loose indicies
            f = find( isfixed); xf = x0(f);
            l = find(~isfixed); xl = x0(l);

            % Estimate only the non-fixed ones
            [xl,r] = ga(@localcost_,sum(isfixed==0),varargin{:});

            % Re-create array combining fixed and estimated coefficients
            x([f,l]) = [xf,xl];

            function y = localcost_(x)
               b([f,l]) = [xf,x]; y = cost_(b);
            end
        end
        
        function [x,r] = sa_(cost_,x0,isfixed,varargin)

            % fixed and loose indicies
            f = find( isfixed); xf = x0(f);
            l = find(~isfixed); xl = x0(l);

            % Estimate only the non-fixed ones
            [xl,r] = simulannealbnd(@localcost_,xl,varargin{:});

            % Re-create array combining fixed and estimated coefficients
            x([f,l]) = [xf,xl];

            function y = localcost_(x)
               b([f,l]) = [xf,x]; y = cost_(b);
            end
        end

        function [c,v] = fit_peak_(x,y,profile)
            import am_lib.*
            
            if nargin < 3; profile='pvoigt'; end
            
            switch profile
                case 'sinc+pvoigt'
                    % define rescaling
                    fscale_= @(c) [log(c(1)),c(2:3),log(c(4)),c(5),log(c(6)),c(7),log(c(8))];
                    rscale_= @(c) [exp(c(1)),c(2:3),exp(c(4)),c(5),exp(c(6)),c(7),exp(c(8))];
                    % define gaussian, sinc, peak, and objective functions
                    label = {'A1','C1','S1','A2','C2','W2','F2','B'};
                    func_ = @(c) c(1).*  sinc_((x-c(2))./c(3)).^2 + ...
                                 c(4).*pvoigt_((x-c(5))./c(6),c(7)) + ...
                                 c(8);
                    cost_ = @(c) sum(abs(log(func_(rscale_(c))) - log(y(:)) ));
                    % lower and upper bounds; starting condition
                    isfixed = [0 0 0 0 0 0 0 0];
                    lb = [0.8*min(y) min(x) 0.02 0.8*min(y) min(x) 1E-4 0 0.8*min(y)]; lb = fscale_(lb); 
                    ub = [1.2*max(y) max(x) 1.00 1.2*max(y) max(x) 1E-2 1 1.2*max(y)]; ub = fscale_(ub);
                    x0 = mean([lb;ub]);
                case 'sinc'
                    % define rescaling
                    fscale_= @(c) [log(c(1)),c(2:3),log(c(4))];
                    rscale_= @(c) [exp(c(1)),c(2:3),exp(c(4))];
                    % define gaussian, sinc, peak, and objective functions
                    label = {'Amp','Center','Width','Background'};
                    func_ = @(c) c(1).*  sinc_((x-c(2))./c(3)).^2 + ...
                                 c(4);
                    cost_ = @(c) sum(abs(log(func_(rscale_(c))) - log(y(:)) ));
                    % lower and upper bounds; starting condition
                    isfixed = [0 0 0 0]; % 0.02 to 2 0.2580
                    lb = [0.8*min(y) min(x) 0.02 0.8*min(y)]; lb = fscale_(lb); 
                    ub = [1.2*max(y) max(x) 1.50 1.2*max(y)]; ub = fscale_(ub);
                    x0 = mean([lb;ub]);
                case 'pvoigt'
                    % define rescaling
                    fscale_= @(c) [log(c(1)),c(2:4),log(c(5))];
                    rscale_= @(c) [exp(c(1)),c(2:4),exp(c(5))];
                    % define gaussian, sinc, peak, and objective functions
                    label = {'Amp','Center','Width','PVoigt','Background'};
                    func_ = @(c) c(1).*pvoigt_((x-c(2))./c(3),c(4)) + ...
                                 c(5);
                    cost_ = @(c) sum(abs(log(func_(rscale_(c))) - log(y(:)) ));
                    % lower and upper bounds; starting condition
                    isfixed = [0 0 0 0 0];
                    lb = [0.5*min(y) min(x) 0.001 0 0.5*min(y)]; lb = fscale_(lb); 
                    ub = [1.5*max(y) max(x) 0.800 1 1.5*max(y)]; ub = fscale_(ub);
                    x0 = mean([lb;ub]);
            end
            
            % optimization options
            opts_hybrid_ = optimoptions(...
                @fmincon,'Display','none',...
                         'MaxFunctionEvaluations',1E2,...
                         'MaxIterations',1E4,...
                         'StepTolerance',1E-18,...
                         'FunctionTolerance',1E-18,...
                         'Algorithm','active-set');
            opts_ = optimoptions(@ga,'PopulationSize',1000, ...
                                     'InitialPopulationMatrix',x0(~isfixed), ...
                                     'MutationFcn',{@mutationadaptfeasible}, ...
                                     'Generations',5, ...
                                     'FunctionTolerance',1E-5, ...
                                     'Display','off',...
                                     'HybridFcn',{@fmincon,opts_hybrid_});
            % plot GA iterations?
            if false; opts_.PlotFcns = {@plot_func_}; end

            % perform optimization
            [c,~] = ga_(cost_,x0,isfixed,[],[],[],[],lb(~isfixed),ub(~isfixed),[],opts_); c = rscale_(c); v = func_(c);
            
            function state = plot_func_(~,state,flag,~)
                switch flag
                    % Plot initialization
                    case 'init'
                        clf; figure(1); set(gcf,'color','w');
                    case {'iter','done'}
                        % find best population
                        [~,j]=min(state.Score);
                        % reconstruct full vector
                        f = find( isfixed); xf = x0(f);
                        l = find(~isfixed); xl = state.Population(j,:);
                        c([f,l]) = [xf,xl];
                        % plot it
                        semilogy(x,func_(rscale_(c)),x,y); axis tight;                        
                        % plot bound info
                        axes('position',[0.4 0.6 0.45 0.25]);
                        a = 1:sum(~isfixed); b = (c(~isfixed)-lb(~isfixed))./(ub(~isfixed)-lb(~isfixed)); plot(a,b,'.-','markersize',20); 
                        h=text(a,b+0.05,strread(sprintf('%0.3g\n',c(~isfixed)),'%s'),'HorizontalAlignment','left'); set(h,'rotation',90);
                        ylim([0 1]); set(gca,'YTick',[0 1],'YTickLabel',{'LB','UB'}); set(gca,'XTickLabel',label(~isfixed)); 
                end
            end
        end
        
        
        % combinatorial 
        
        function [x]    = nchoosek_(n,k)
            % much faster than matlab's nchoosek
            kk=min(k,n-k);
            if kk<2
               if kk<1
                  if k==n
                     x=1:n;
                  else
                     x=[];
                  end
               else
                  if k==1
                     x=(1:n)';
                  else
                     x=1:n;
                     x=reshape(x(ones(n-1,1),:),n,n-1);
                  end
               end   
            else
               n1=n+1;
               m=prod(n1-kk:n)/prod(1:kk);
               x=zeros(m,k);
               f=n1-k;
               x(1:f,k)=(k:n)';
               for a=k-1:-1:1
                  d=f;
                  h=f;
                  x(1:f,a)=a;
                  for b=a+1:a+n-k
                     d=d*(n1+a-b-k)/(n1-b);
                     e=f+1;
                     f=e+d-1;
                     x(e:f,a)=b;
                     x(e:f,a+1:k)=x(h-d+1:h,a+1:k);
                  end
               end
            end
            x=x.';
        end
        
        function [x]    = nchoosrk_(n,k)
            % with duplications allowed
            
            import am_lib.nchoosek_
            
            x=nchoosek_(n+k-1,k).';
            x=x-repmat(0:k-1,size(x,1),1);
            x=x.';
        end
        
        function PN     = perm_norep_(N,K)
            % Subfunction: permutations without replacement.
            % Uses the algorithm in combs_no_rep as a basis, then permutes each row.
            % pn = @(N,K) prod(1:N)/(prod(1:(N-K)));  Number of rows.

            if N==K
                PN = perms_loop(N);  % Call helper function.
                return
            elseif K==1
                PN = (1:N);  % Easy case.
                return
            end

            if K>N  % Since there is no replacement, this cannot happen.
                error(['When no repetitions are allowed, '...
                       'K must be less than or equal to N'])
            end

            M = double(N);  % Single will give us trouble on indexing.
            WV = 1:K;  % Working vector.
            lim = K;   % Sets the limit for working index.
            inc = 1;   % Controls which element of WV is being worked on.
            BC = prod(M-K+1:M);  % Pre-allocation of return arg.
            BC1 = BC / ( prod(1:K)); % Number of comb blocks.
            PN = zeros(round(BC),K,class(N));
            L = prod(1:K) ;  % To get the size of the blocks.
            cnt = 1+L;
            P = perms_loop(K);  % Only need to use this once.
            PN(1:(1+L-1),:) = WV(P);  % The first row.

            for ii = 2:(BC1 - 1)
                if logical((inc+lim)-N)  % The logical is nec. for class single(?)
                    stp = inc;  % This is where the for loop below stops.
                    flg = 0;  % Used for resetting inc.
                else
                    stp = 1;
                    flg = 1;
                end

                for jj = 1:stp
                    WV(K  + jj - inc) = lim + jj;  % Faster than a vector assignment!
                end

                PN(cnt:(cnt+L-1),:) = WV(P);  % Assign block.
                cnt = cnt + L;  % Increment base index.    
                inc = inc*flg + 1;  % Increment the counter.
                lim = WV(K - inc + 1 );  % lim for next run.
            end

            V = (N-K+1):N;  % Final vector.
            PN(cnt:(cnt+L-1),:) = V(P);  % Fill final block.
            % The sorting below is NOT necessary.  If you prefer this nice
            % order, the next two lines can be un-commented.
            % [id,id] = sort(PN(:,1));  %#ok  This is not necessary!
            % PN = PN(id,:);  % Return values.
            
            % rotate for consistency with other functions
            PN = PN.';

            function P = perms_loop(N)
                % Helper function to perms_no_rep.  This is basically the same as the
                % MATLAB function perms.  It has been un-recursed for a runtime of around  
                % half the recursive version found in perms.m  For example:
                %
                %      tic,Tp = perms(1:9);toc
                %      %Elapsed time is 0.222111 seconds.  Allow Ctrl+T+C+R on block
                %      tic,Tc = combinator(9,9,'p');toc  
                %      %Elapsed time is 0.143219 seconds.
                %      isequal(Tc,Tp)  % Yes

                Q = double(N); % Single will give us trouble on indexing.
                P = 1;  % Initializer.
                G = cumprod(1:(Q-1));  % Holds the sizes of P.
                CN = class(N);

                for n = 2:Q
                    q = P;
                    m = G(n-1);
                    P = zeros(n*m,n,CN);
                    P(1:m, 1) = n;
                    P(1:m, 2:n) = q;
                    a = m + 1;

                    for kk = n-1:-1:1
                        t = q;
                        t(t == kk) = n;
                        b = a + m - 1;
                        P(a:b, 1) = kk;
                        P(a:b, 2:n) = t;
                        a = b + 1;
                    end 
                end
            end
        end
        
        function [M, I] = permn_(V, N, K)
            % PERMN - permutations with repetition
            %   Using two input variables V and N, M = PERMN(V,N) returns all
            %   permutations of N elements taken from the vector V, with repetitions.
            %   V can be any type of array (numbers, cells etc.) and M will be of the
            %   same type as V.  If V is empty or N is 0, M will be empty.  M has the
            %   size numel(V).^N-by-N. 
            %
            %   When only a subset of these permutations is needed, you can call PERMN
            %   with 3 input variables: M = PERMN(V,N,K) returns only the K-ths
            %   permutations.  The output is the same as M = PERMN(V,N) ; M = M(K,:),
            %   but it avoids memory issues that may occur when there are too many
            %   combinations.  This is particulary useful when you only need a few
            %   permutations at a given time. If V or K is empty, or N is zero, M will
            %   be empty. M has the size numel(K)-by-N. 
            %
            %   [M, I] = PERMN(...) also returns an index matrix I so that M = V(I).
            %
            %   Examples:
            %     M = permn([1 2 3],2) % returns the 9-by-2 matrix:
            %              1     1
            %              1     2
            %              1     3
            %              2     1
            %              2     2
            %              2     3
            %              3     1
            %              3     2
            %              3     3
            
            narginchk(2,3) ;

            if fix(N) ~= N || N < 0 || numel(N) ~= 1 
                error('permn:negativeN','Second argument should be a positive integer') ;
            end
            nV = numel(V) ;

            if nargin==2 % PERMN(V,N) - return all permutations

                if nV==0 || N == 0
                    M = zeros(nV,N) ;
                    I = zeros(nV,N) ;

                elseif N == 1
                    % return column vectors
                    M = V(:) ;
                    I = (1:nV).' ;
                else
                    % this is faster than the math trick used for the call with three
                    % arguments.
                    [Y{N:-1:1}] = ndgrid(1:nV) ;
                    I = reshape(cat(N+1,Y{:}),[],N) ;
                    M = V(I) ;
                end
            else % PERMN(V,N,K) - return a subset of all permutations
                nK = numel(K) ;
                if nV == 0 || N == 0 || nK == 0
                    M = zeros(numel(K), N) ;
                    I = zeros(numel(K), N) ;
                elseif nK < 1 || any(K<1) || any(K ~= fix(K))
                    error('permn:InvalidIndex','Third argument should contain positive integers.') ;
                else

                    V = reshape(V,1,[]) ; % v1.1 make input a row vector
                    nV = numel(V) ;
                    Npos = nV^N ;
                    if any(K > Npos)
                        warning('permn:IndexOverflow', ...
                            'Values of K exceeding the total number of combinations are saturated.')
                        K = min(K, Npos) ;
                    end

                    % The engine is based on version 3.2 with the correction
                    % suggested by Roger Stafford. This approach uses a single matrix
                    % multiplication.
                    B = nV.^(1-N:0) ;
                    I = ((K(:)-.5) * B) ; % matrix multiplication
                    I = rem(floor(I),nV) + 1 ;
                    M = V(I) ;
                end
            end
        end

        function [A,c]  = perm_heap_(A,n,c)
            % [A,c] = heap_perm_(A,n,c)
            % generates on permutation at a time. check it with:
            % n=4; A = zeros(factorial(n),n);
            % % initialize
            % [A(1,:),c] = heap_perm_([1:n],n);
            % % generate permutations
            % for i = 2:factorial(n)
            %     [A(i,:),c] = heap_perm_(A(i-1,:),n,c);
            % end
            % sortrows(perms([1:n]))-sortrows(A)

            if nargin<3
                % c is used to continue generating permutations
                c = ones(1,n);
            end
            i=1;
            while i <= n
                if  c(i) < i
                    if mod(i,2)==1
                        A([1,i])=A([i,1]);
                    else
                        A([i,c(i)])=A([c(i),i]);
                    end
                    if nargout > 1
                        % everytime a new permutation is generated, reset i to 1
                        c(i) = c(i) + 1; i = 1;
                    end
                    % return A and c for continuation
                    return
                else
                    c(i) = 1;
                    i = i + 1;
                end
            end
        end
        
        function [a]    = perm_narayana_(a)
            % [a]   = narayana_perm_(a)
            % given a permutation, generates the next permutation in
            % lexicographical order. start generations with a = [1:n]
            % for example:      
            %     clear;clc
            %     n = 5;            
            %     i=0;   a=zeros(factorial(n),n);
            %     i=i+1; a(i,1:n)=1:n; tic
            %     t = narayana_perm_(a(1,:));
            %     while ~isempty(t)
            %         i=i+1; a(i,:) = t;
            %         t = narayana_perm_(a(i,:));
            %     end
            %

            % get number of elements
            n = numel(a);
            % one element
            if n==1; return; end
            % find largest k such that a(k) < a(k+1)
            k=n-1;
            while a(k)>=a(k+1)
                k=k-1;
                if k==0
                    % last permutation
                    a = []; return;
                end
            end
            % find largest l > k such that a(k) < a(l)
            l = n;
            while a(l)<a(k)
                l=l-1;
            end
            % swap k and l positions
            a([k,l])=a([l,k]);
            % reverse k+1 to n positions
            a([(k+1):n])=a([n:-1:(k+1)]);
        end
        
        function [x]    = ind2bas_(index,base)
            % [dhms] = ind2bas_(seconds,[7,24,60,60])
            % use this to convert from seconds to weeks,days,hours minutes,seconds 
            n = numel(base)+1; x = zeros(1,n-1); b = cumprod([base,1],'reverse');
            for i = 2:n
                x(i-1)  = floor(index./b(i));
                index = index - x(i-1).*b(i);
            end
        end

        function [a]    = group_perm_(a,label,j)
            % [a] = group_perm_(a,label)
            % call like this: n=9; a = 1:n; label = [1,1,1,2,2,2,2,2,3];
            
            import am_lib.*
            
            %  are in monotonically increasing order?
            if any(diff(label)<0) % yes
                dosort = true; [~,fwd] = sort(label); rev(fwd)=[1:n]; a=a(fwd); 
            else % no
                dosort = false;
            end

                % sorted
                ulabels = unique(label);
                % nlabels = numel(ulabels);
                mdigits = factorial(sum(label(:)==ulabels,1)); 
                % ncounts = prod(mdigits);

                S = find([1,diff(label)]);
                E = S + sum(ulabels==label.',1) - 1;

                % initialize
                % fprintf('%3i',a); fprintf('\n');
                % for j = 2:ncounts
                    % update
                    k = find(ind2bas_(j-1,mdigits)-ind2bas_(j-2,mdigits)>0);
                    % permute
                    t = narayana_perm_(a(S(k):E(k)));
                    % restart
                    if isempty(t); t = S(k):E(k); end
                    % update
                    a(S(k):E(k)) = t;
                    % print
                    % fprintf('%3i',a); fprintf('\n');
                % end
                
            if dosort
                a = a(rev);
            end
        end
        
        
        % general plotting
        
        function [th] = assign_color_(V)
            % assigns a number between [0,1] based on how close vectors V
            % in are to identity vectors in that same basis. It essentially
            % reduces the dimensionality of the data

            % set number of points
            n = size(V,1);

            % define points around circle on complex plane
            p = exp(2*pi*1i*[1:n]/n);

            % get complex points
            c = (p*abs(V)).';

            % get linear mapping between [0,1]
            th = mod( (atan2d(imag(c),real(c))/180+1)/2 + 3/8, 1);
        end
        
        function [h] = plotc_(x,y,c)
            x = x(:).'; y=y(:).'; c=c(:).'; z=zeros(size(x));
            % col = x;  % This is the color, vary with x in this case.
            h = surface([x;x],[y;y],[z;z],[c;c],'facecol','no','edgecol','interp','linew',1);
        end
        
        function [h] = plot3_(A,varargin)
           h = plot3(A(1,:),A(2,:),A(3,:),varargin{:});
        end
        
        function [h] = scatter3_(A,varargin)
           h = scatter3(A(1,:),A(2,:),A(3,:),varargin{:});
        end
        
        function [h] = plothull_(A,varargin)
            
            import am_lib.*
            
            % initialize figure
            set(gcf,'color','w'); hold on;
            
            % get points in convex hull
            DT = delaunayTriangulation(A(1,:).',A(2,:).',A(3,:).'); 
            CH = convexHull(DT); A=A(:,unique(CH));
            
            % get and plot faces
            DT = delaunayTriangulation(A(1,:).',A(2,:).',A(3,:).');
            CH = convexHull(DT); TR = triangulation(CH,A.');
            h = trisurf(TR,'FaceColor','black','EdgeColor','none','FaceAlpha',0.01); 

            % get and plot edges
            FE = featureEdges(TR,pi/100).'; 
            for i = 1:size(FE,2); plot3_(A(:,FE(:,i)),'-','color','k','linewidth',2); end

            hold off; daspect([1 1 1])
        end
        
        function [h] = plotv3_(A,varargin)
           A = repelem(A,1,2); A(:,1:2:end)=0;
           h = plot3(A(1,:),A(2,:),A(3,:),varargin{:});
        end

        function [h] = spyc_(A)
            [x,y] = find(A);
            h = scatter(y,x,200,A(A~=0),'.');
            set(gca,'YDir','rev'); box on;
            ylim([0 size(A,1)+1]); xlim([0 size(A,2)+1]); 
            set(gca,'XTick',[]); set(gca,'YTick',[]);
            daspect([1 1 1])
        end       

        function [cmap] = color_(n,palette)
            % color_(n,palette)
            
            if nargin < 2; palette='spectral'; end
            %   
            brighten_ = @(x,alpha) x*alpha + (1-alpha);
            dim_ = @(x,beta) x*beta;
            % predefine some colormaps
            set1 = brighten_([[ 55, 126, 184]*.85;[228, 26, 28];[ 77, 175, 74];[ 255, 127, 0];[ 152, 78, 163]],0.8);
            set2 = brighten_([[228, 26, 28];[ 55, 126, 184]; [ 77, 175, 74];[ 255, 127, 0];[ 255, 237, 111]*.85;[ 166, 86, 40];[ 247, 129, 191];[ 153, 153, 153];[ 152, 78, 163]],0.9);
            set3 = dim_([[141, 211, 199];[ 255, 237, 111];[ 190, 186, 218];[ 251, 128, 114];[ 128, 177, 211];[ 253, 180, 98];[ 179, 222, 105];[ 188, 128, 189];[ 217, 217, 217];[ 204, 235, 197];[ 252, 205, 229];[ 255, 255, 179]],0.93);
            % interpolating function
            map_ = @(n,cmap) interp1([0:(size(cmap,1)-1)]./(size(cmap,1)-1),cmap,linspace(0,1,n));
            % switch based on N
            switch n
                case {1}
                    cmap = [55, 126, 184]./255;
                case {2, 3, 4, 5 }
                    cmap = set1(1:n,:)./255;
                case {6 , 7, 8, 9}
                    cmap = set2(1:n,:)./255;
                case {10, 11, 12}
                    cmap = set3(1:n,:)./255;
                otherwise 
                    switch palette
                        case 'virdis'
                            cmap =[...
                            0.267004,0.004874,0.329415; 0.268510,0.009605,0.335427; 0.269944,0.014625,0.341379; 0.271305,0.019942,0.347269; 0.272594,0.025563,0.353093; ...
                            0.273809,0.031497,0.358853; 0.274952,0.037752,0.364543; 0.276022,0.044167,0.370164; 0.277018,0.050344,0.375715; 0.277941,0.056324,0.381191; ...
                            0.278791,0.062145,0.386592; 0.279566,0.067836,0.391917; 0.280267,0.073417,0.397163; 0.280894,0.078907,0.402329; 0.281446,0.084320,0.407414; ...
                            0.281924,0.089666,0.412415; 0.282327,0.094955,0.417331; 0.282656,0.100196,0.422160; 0.282910,0.105393,0.426902; 0.283091,0.110553,0.431554; ...
                            0.283197,0.115680,0.436115; 0.283229,0.120777,0.440584; 0.283187,0.125848,0.444960; 0.283072,0.130895,0.449241; 0.282884,0.135920,0.453427; ...
                            0.282623,0.140926,0.457517; 0.282290,0.145912,0.461510; 0.281887,0.150881,0.465405; 0.281412,0.155834,0.469201; 0.280868,0.160771,0.472899; ...
                            0.280255,0.165693,0.476498; 0.279574,0.170599,0.479997; 0.278826,0.175490,0.483397; 0.278012,0.180367,0.486697; 0.277134,0.185228,0.489898; ...
                            0.276194,0.190074,0.493001; 0.275191,0.194905,0.496005; 0.274128,0.199721,0.498911; 0.273006,0.204520,0.501721; 0.271828,0.209303,0.504434; ...
                            0.270595,0.214069,0.507052; 0.269308,0.218818,0.509577; 0.267968,0.223549,0.512008; 0.266580,0.228262,0.514349; 0.265145,0.232956,0.516599; ...
                            0.263663,0.237631,0.518762; 0.262138,0.242286,0.520837; 0.260571,0.246922,0.522828; 0.258965,0.251537,0.524736; 0.257322,0.256130,0.526563; ...
                            0.255645,0.260703,0.528312; 0.253935,0.265254,0.529983; 0.252194,0.269783,0.531579; 0.250425,0.274290,0.533103; 0.248629,0.278775,0.534556; ...
                            0.246811,0.283237,0.535941; 0.244972,0.287675,0.537260; 0.243113,0.292092,0.538516; 0.241237,0.296485,0.539709; 0.239346,0.300855,0.540844; ...
                            0.237441,0.305202,0.541921; 0.235526,0.309527,0.542944; 0.233603,0.313828,0.543914; 0.231674,0.318106,0.544834; 0.229739,0.322361,0.545706; ...
                            0.227802,0.326594,0.546532; 0.225863,0.330805,0.547314; 0.223925,0.334994,0.548053; 0.221989,0.339161,0.548752; 0.220057,0.343307,0.549413; ...
                            0.218130,0.347432,0.550038; 0.216210,0.351535,0.550627; 0.214298,0.355619,0.551184; 0.212395,0.359683,0.551710; 0.210503,0.363727,0.552206; ...
                            0.208623,0.367752,0.552675; 0.206756,0.371758,0.553117; 0.204903,0.375746,0.553533; 0.203063,0.379716,0.553925; 0.201239,0.383670,0.554294; ...
                            0.199430,0.387607,0.554642; 0.197636,0.391528,0.554969; 0.195860,0.395433,0.555276; 0.194100,0.399323,0.555565; 0.192357,0.403199,0.555836; ...
                            0.190631,0.407061,0.556089; 0.188923,0.410910,0.556326; 0.187231,0.414746,0.556547; 0.185556,0.418570,0.556753; 0.183898,0.422383,0.556944; ...
                            0.182256,0.426184,0.557120; 0.180629,0.429975,0.557282; 0.179019,0.433756,0.557430; 0.177423,0.437527,0.557565; 0.175841,0.441290,0.557685; ...
                            0.174274,0.445044,0.557792; 0.172719,0.448791,0.557885; 0.171176,0.452530,0.557965; 0.169646,0.456262,0.558030; 0.168126,0.459988,0.558082; ...
                            0.166617,0.463708,0.558119; 0.165117,0.467423,0.558141; 0.163625,0.471133,0.558148; 0.162142,0.474838,0.558140; 0.160665,0.478540,0.558115; ...
                            0.159194,0.482237,0.558073; 0.157729,0.485932,0.558013; 0.156270,0.489624,0.557936; 0.154815,0.493313,0.557840; 0.153364,0.497000,0.557724; ...
                            0.151918,0.500685,0.557587; 0.150476,0.504369,0.557430; 0.149039,0.508051,0.557250; 0.147607,0.511733,0.557049; 0.146180,0.515413,0.556823; ...
                            0.144759,0.519093,0.556572; 0.143343,0.522773,0.556295; 0.141935,0.526453,0.555991; 0.140536,0.530132,0.555659; 0.139147,0.533812,0.555298; ...
                            0.137770,0.537492,0.554906; 0.136408,0.541173,0.554483; 0.135066,0.544853,0.554029; 0.133743,0.548535,0.553541; 0.132444,0.552216,0.553018; ...
                            0.131172,0.555899,0.552459; 0.129933,0.559582,0.551864; 0.128729,0.563265,0.551229; 0.127568,0.566949,0.550556; 0.126453,0.570633,0.549841; ...
                            0.125394,0.574318,0.549086; 0.124395,0.578002,0.548287; 0.123463,0.581687,0.547445; 0.122606,0.585371,0.546557; 0.121831,0.589055,0.545623; ...
                            0.121148,0.592739,0.544641; 0.120565,0.596422,0.543611; 0.120092,0.600104,0.542530; 0.119738,0.603785,0.541400; 0.119512,0.607464,0.540218; ...
                            0.119423,0.611141,0.538982; 0.119483,0.614817,0.537692; 0.119699,0.618490,0.536347; 0.120081,0.622161,0.534946; 0.120638,0.625828,0.533488; ...
                            0.121380,0.629492,0.531973; 0.122312,0.633153,0.530398; 0.123444,0.636809,0.528763; 0.124780,0.640461,0.527068; 0.126326,0.644107,0.525311; ...
                            0.128087,0.647749,0.523491; 0.130067,0.651384,0.521608; 0.132268,0.655014,0.519661; 0.134692,0.658636,0.517649; 0.137339,0.662252,0.515571; ...
                            0.140210,0.665859,0.513427; 0.143303,0.669459,0.511215; 0.146616,0.673050,0.508936; 0.150148,0.676631,0.506589; 0.153894,0.680203,0.504172; ...
                            0.157851,0.683765,0.501686; 0.162016,0.687316,0.499129; 0.166383,0.690856,0.496502; 0.170948,0.694384,0.493803; 0.175707,0.697900,0.491033; ...
                            0.180653,0.701402,0.488189; 0.185783,0.704891,0.485273; 0.191090,0.708366,0.482284; 0.196571,0.711827,0.479221; 0.202219,0.715272,0.476084; ...
                            0.208030,0.718701,0.472873; 0.214000,0.722114,0.469588; 0.220124,0.725509,0.466226; 0.226397,0.728888,0.462789; 0.232815,0.732247,0.459277; ...
                            0.239374,0.735588,0.455688; 0.246070,0.738910,0.452024; 0.252899,0.742211,0.448284; 0.259857,0.745492,0.444467; 0.266941,0.748751,0.440573; ...
                            0.274149,0.751988,0.436601; 0.281477,0.755203,0.432552; 0.288921,0.758394,0.428426; 0.296479,0.761561,0.424223; 0.304148,0.764704,0.419943; ...
                            0.311925,0.767822,0.415586; 0.319809,0.770914,0.411152; 0.327796,0.773980,0.406640; 0.335885,0.777018,0.402049; 0.344074,0.780029,0.397381; ...
                            0.352360,0.783011,0.392636; 0.360741,0.785964,0.387814; 0.369214,0.788888,0.382914; 0.377779,0.791781,0.377939; 0.386433,0.794644,0.372886; ...
                            0.395174,0.797475,0.367757; 0.404001,0.800275,0.362552; 0.412913,0.803041,0.357269; 0.421908,0.805774,0.351910; 0.430983,0.808473,0.346476; ...
                            0.440137,0.811138,0.340967; 0.449368,0.813768,0.335384; 0.458674,0.816363,0.329727; 0.468053,0.818921,0.323998; 0.477504,0.821444,0.318195; ...
                            0.487026,0.823929,0.312321; 0.496615,0.826376,0.306377; 0.506271,0.828786,0.300362; 0.515992,0.831158,0.294279; 0.525776,0.833491,0.288127; ...
                            0.535621,0.835785,0.281908; 0.545524,0.838039,0.275626; 0.555484,0.840254,0.269281; 0.565498,0.842430,0.262877; 0.575563,0.844566,0.256415; ...
                            0.585678,0.846661,0.249897; 0.595839,0.848717,0.243329; 0.606045,0.850733,0.236712; 0.616293,0.852709,0.230052; 0.626579,0.854645,0.223353; ...
                            0.636902,0.856542,0.216620; 0.647257,0.858400,0.209861; 0.657642,0.860219,0.203082; 0.668054,0.861999,0.196293; 0.678489,0.863742,0.189503; ...
                            0.688944,0.865448,0.182725; 0.699415,0.867117,0.175971; 0.709898,0.868751,0.169257; 0.720391,0.870350,0.162603; 0.730889,0.871916,0.156029; ...
                            0.741388,0.873449,0.149561; 0.751884,0.874951,0.143228; 0.762373,0.876424,0.137064; 0.772852,0.877868,0.131109; 0.783315,0.879285,0.125405; ...
                            0.793760,0.880678,0.120005; 0.804182,0.882046,0.114965; 0.814576,0.883393,0.110347; 0.824940,0.884720,0.106217; 0.835270,0.886029,0.102646; ...
                            0.845561,0.887322,0.099702; 0.855810,0.888601,0.097452; 0.866013,0.889868,0.095953; 0.876168,0.891125,0.095250; 0.886271,0.892374,0.095374; ...
                            0.896320,0.893616,0.096335; 0.906311,0.894855,0.098125; 0.916242,0.896091,0.100717; 0.926106,0.897330,0.104071; 0.935904,0.898570,0.108131; ...
                            0.945636,0.899815,0.112838; 0.955300,0.901065,0.118128; 0.964894,0.902323,0.123941; 0.974417,0.903590,0.130215; 0.983868,0.904867,0.136897; ...
                            0.993248,0.906157,0.143936];
                        case 'magma'
                            cmap = [ ...
                            0.001462,0.000466,0.013866; 0.002258,0.001295,0.018331; 0.003279,0.002305,0.023708; 0.004512,0.003490,0.029965; 0.005950,0.004843,0.037130; ...
                            0.007588,0.006356,0.044973; 0.009426,0.008022,0.052844; 0.011465,0.009828,0.060750; 0.013708,0.011771,0.068667; 0.016156,0.013840,0.076603; ...
                            0.018815,0.016026,0.084584; 0.021692,0.018320,0.092610; 0.024792,0.020715,0.100676; 0.028123,0.023201,0.108787; 0.031696,0.025765,0.116965; ...
                            0.035520,0.028397,0.125209; 0.039608,0.031090,0.133515; 0.043830,0.033830,0.141886; 0.048062,0.036607,0.150327; 0.052320,0.039407,0.158841; ...
                            0.056615,0.042160,0.167446; 0.060949,0.044794,0.176129; 0.065330,0.047318,0.184892; 0.069764,0.049726,0.193735; 0.074257,0.052017,0.202660; ...
                            0.078815,0.054184,0.211667; 0.083446,0.056225,0.220755; 0.088155,0.058133,0.229922; 0.092949,0.059904,0.239164; 0.097833,0.061531,0.248477; ...
                            0.102815,0.063010,0.257854; 0.107899,0.064335,0.267289; 0.113094,0.065492,0.276784; 0.118405,0.066479,0.286321; 0.123833,0.067295,0.295879; ...
                            0.129380,0.067935,0.305443; 0.135053,0.068391,0.315000; 0.140858,0.068654,0.324538; 0.146785,0.068738,0.334011; 0.152839,0.068637,0.343404; ...
                            0.159018,0.068354,0.352688; 0.165308,0.067911,0.361816; 0.171713,0.067305,0.370771; 0.178212,0.066576,0.379497; 0.184801,0.065732,0.387973; ...
                            0.191460,0.064818,0.396152; 0.198177,0.063862,0.404009; 0.204935,0.062907,0.411514; 0.211718,0.061992,0.418647; 0.218512,0.061158,0.425392; ...
                            0.225302,0.060445,0.431742; 0.232077,0.059889,0.437695; 0.238826,0.059517,0.443256; 0.245543,0.059352,0.448436; 0.252220,0.059415,0.453248; ...
                            0.258857,0.059706,0.457710; 0.265447,0.060237,0.461840; 0.271994,0.060994,0.465660; 0.278493,0.061978,0.469190; 0.284951,0.063168,0.472451; ...
                            0.291366,0.064553,0.475462; 0.297740,0.066117,0.478243; 0.304081,0.067835,0.480812; 0.310382,0.069702,0.483186; 0.316654,0.071690,0.485380; ...
                            0.322899,0.073782,0.487408; 0.329114,0.075972,0.489287; 0.335308,0.078236,0.491024; 0.341482,0.080564,0.492631; 0.347636,0.082946,0.494121; ...
                            0.353773,0.085373,0.495501; 0.359898,0.087831,0.496778; 0.366012,0.090314,0.497960; 0.372116,0.092816,0.499053; 0.378211,0.095332,0.500067; ...
                            0.384299,0.097855,0.501002; 0.390384,0.100379,0.501864; 0.396467,0.102902,0.502658; 0.402548,0.105420,0.503386; 0.408629,0.107930,0.504052; ...
                            0.414709,0.110431,0.504662; 0.420791,0.112920,0.505215; 0.426877,0.115395,0.505714; 0.432967,0.117855,0.506160; 0.439062,0.120298,0.506555; ...
                            0.445163,0.122724,0.506901; 0.451271,0.125132,0.507198; 0.457386,0.127522,0.507448; 0.463508,0.129893,0.507652; 0.469640,0.132245,0.507809; ...
                            0.475780,0.134577,0.507921; 0.481929,0.136891,0.507989; 0.488088,0.139186,0.508011; 0.494258,0.141462,0.507988; 0.500438,0.143719,0.507920; ...
                            0.506629,0.145958,0.507806; 0.512831,0.148179,0.507648; 0.519045,0.150383,0.507443; 0.525270,0.152569,0.507192; 0.531507,0.154739,0.506895; ...
                            0.537755,0.156894,0.506551; 0.544015,0.159033,0.506159; 0.550287,0.161158,0.505719; 0.556571,0.163269,0.505230; 0.562866,0.165368,0.504692; ...
                            0.569172,0.167454,0.504105; 0.575490,0.169530,0.503466; 0.581819,0.171596,0.502777; 0.588158,0.173652,0.502035; 0.594508,0.175701,0.501241; ...
                            0.600868,0.177743,0.500394; 0.607238,0.179779,0.499492; 0.613617,0.181811,0.498536; 0.620005,0.183840,0.497524; 0.626401,0.185867,0.496456; ...
                            0.632805,0.187893,0.495332; 0.639216,0.189921,0.494150; 0.645633,0.191952,0.492910; 0.652056,0.193986,0.491611; 0.658483,0.196027,0.490253; ...
                            0.664915,0.198075,0.488836; 0.671349,0.200133,0.487358; 0.677786,0.202203,0.485819; 0.684224,0.204286,0.484219; 0.690661,0.206384,0.482558; ...
                            0.697098,0.208501,0.480835; 0.703532,0.210638,0.479049; 0.709962,0.212797,0.477201; 0.716387,0.214982,0.475290; 0.722805,0.217194,0.473316; ...
                            0.729216,0.219437,0.471279; 0.735616,0.221713,0.469180; 0.742004,0.224025,0.467018; 0.748378,0.226377,0.464794; 0.754737,0.228772,0.462509; ...
                            0.761077,0.231214,0.460162; 0.767398,0.233705,0.457755; 0.773695,0.236249,0.455289; 0.779968,0.238851,0.452765; 0.786212,0.241514,0.450184; ...
                            0.792427,0.244242,0.447543; 0.798608,0.247040,0.444848; 0.804752,0.249911,0.442102; 0.810855,0.252861,0.439305; 0.816914,0.255895,0.436461; ...
                            0.822926,0.259016,0.433573; 0.828886,0.262229,0.430644; 0.834791,0.265540,0.427671; 0.840636,0.268953,0.424666; 0.846416,0.272473,0.421631; ...
                            0.852126,0.276106,0.418573; 0.857763,0.279857,0.415496; 0.863320,0.283729,0.412403; 0.868793,0.287728,0.409303; 0.874176,0.291859,0.406205; ...
                            0.879464,0.296125,0.403118; 0.884651,0.300530,0.400047; 0.889731,0.305079,0.397002; 0.894700,0.309773,0.393995; 0.899552,0.314616,0.391037; ...
                            0.904281,0.319610,0.388137; 0.908884,0.324755,0.385308; 0.913354,0.330052,0.382563; 0.917689,0.335500,0.379915; 0.921884,0.341098,0.377376; ...
                            0.925937,0.346844,0.374959; 0.929845,0.352734,0.372677; 0.933606,0.358764,0.370541; 0.937221,0.364929,0.368567; 0.940687,0.371224,0.366762; ...
                            0.944006,0.377643,0.365136; 0.947180,0.384178,0.363701; 0.950210,0.390820,0.362468; 0.953099,0.397563,0.361438; 0.955849,0.404400,0.360619; ...
                            0.958464,0.411324,0.360014; 0.960949,0.418323,0.359630; 0.963310,0.425390,0.359469; 0.965549,0.432519,0.359529; 0.967671,0.439703,0.359810; ...
                            0.969680,0.446936,0.360311; 0.971582,0.454210,0.361030; 0.973381,0.461520,0.361965; 0.975082,0.468861,0.363111; 0.976690,0.476226,0.364466; ...
                            0.978210,0.483612,0.366025; 0.979645,0.491014,0.367783; 0.981000,0.498428,0.369734; 0.982279,0.505851,0.371874; 0.983485,0.513280,0.374198; ...
                            0.984622,0.520713,0.376698; 0.985693,0.528148,0.379371; 0.986700,0.535582,0.382210; 0.987646,0.543015,0.385210; 0.988533,0.550446,0.388365; ...
                            0.989363,0.557873,0.391671; 0.990138,0.565296,0.395122; 0.990871,0.572706,0.398714; 0.991558,0.580107,0.402441; 0.992196,0.587502,0.406299; ...
                            0.992785,0.594891,0.410283; 0.993326,0.602275,0.414390; 0.993834,0.609644,0.418613; 0.994309,0.616999,0.422950; 0.994738,0.624350,0.427397; ...
                            0.995122,0.631696,0.431951; 0.995480,0.639027,0.436607; 0.995810,0.646344,0.441361; 0.996096,0.653659,0.446213; 0.996341,0.660969,0.451160; ...
                            0.996580,0.668256,0.456192; 0.996775,0.675541,0.461314; 0.996925,0.682828,0.466526; 0.997077,0.690088,0.471811; 0.997186,0.697349,0.477182; ...
                            0.997254,0.704611,0.482635; 0.997325,0.711848,0.488154; 0.997351,0.719089,0.493755; 0.997351,0.726324,0.499428; 0.997341,0.733545,0.505167; ...
                            0.997285,0.740772,0.510983; 0.997228,0.747981,0.516859; 0.997138,0.755190,0.522806; 0.997019,0.762398,0.528821; 0.996898,0.769591,0.534892; ...
                            0.996727,0.776795,0.541039; 0.996571,0.783977,0.547233; 0.996369,0.791167,0.553499; 0.996162,0.798348,0.559820; 0.995932,0.805527,0.566202; ...
                            0.995680,0.812706,0.572645; 0.995424,0.819875,0.579140; 0.995131,0.827052,0.585701; 0.994851,0.834213,0.592307; 0.994524,0.841387,0.598983; ...
                            0.994222,0.848540,0.605696; 0.993866,0.855711,0.612482; 0.993545,0.862859,0.619299; 0.993170,0.870024,0.626189; 0.992831,0.877168,0.633109; ...
                            0.992440,0.884330,0.640099; 0.992089,0.891470,0.647116; 0.991688,0.898627,0.654202; 0.991332,0.905763,0.661309; 0.990930,0.912915,0.668481; ...
                            0.990570,0.920049,0.675675; 0.990175,0.927196,0.682926; 0.989815,0.934329,0.690198; 0.989434,0.941470,0.697519; 0.989077,0.948604,0.704863; ...
                            0.988717,0.955742,0.712242; 0.988367,0.962878,0.719649; 0.988033,0.970012,0.727077; 0.987691,0.977154,0.734536; 0.987387,0.984288,0.742002; ...
                            0.987053,0.991438,0.749504];
                        case 'red2blue' % sns.color_palette("RdBu_r", 7)
                            cmap = [ ... 
                            0.16339870177063293, 0.44498270983789490, 0.6975009791991290;
                            0.42068437209316328, 0.67643216077019186, 0.8186851319144753;
                            0.76147636946509856, 0.86851211856393251, 0.9245674785445717;
                            0.96908881383783674, 0.96647443490869855, 0.9649365649503820;
                            0.98246828247519102, 0.80069205340217142, 0.7061130509657018;
                            0.89457901435739851, 0.50380624217145586, 0.3997693394913390;
                            0.72848905885920801, 0.15501730406985564, 0.1973856272650700];
                        otherwise % spectral
                            cmap = [ ...
                            0.81414841553744144, 0.21968473703143937, 0.30480585554066825;
                            0.93302576331531295, 0.39131103777417953, 0.27197233193060932;
                            0.98177624099394856, 0.60738179087638855, 0.34579008992980509;
                            0.99469434864380779, 0.80922723167082844, 0.48696657138712268;
                            0.99823144954793597, 0.94517493598601399, 0.65705499929540301;
                            0.95578623869839840, 0.98231449547935934, 0.68004615517223588;
                            0.82029989537070780, 0.92756632496328917, 0.61268745099796973;
                            0.59100347582031698, 0.83552480795804196, 0.64429067864137535;
                            0.36001538412243711, 0.71618609919267540, 0.66551328406614418;
                            0.21299500558890549, 0.51141871132102668, 0.73079586379668293];
                    end
                    % interpolate and return
                    cmap = map_(n,cmap);
            end
        end
        
        
        % correlation and polynomial fitting
        
        function       plotcorr_(x,y)
            % plot correlation x vs y
            
            import am_lib.*
            
            % plot correlation for dft vs bvk forces on atoms
            plot(x(:),y(:),'.','color',[1 1 1]*0.70); daspect([1 1 1]); box on;
            title(sprintf('R^2 = %f',pcorr_(x(:),y(:)).^2)); 
            
            % linear regression
            % mv = [-1 1]*max(abs(axis));
            mv = minmax_([x(:);y(:)]);
            line( mv, mv, 'linewidth',1.5); 
            line( mv, pinterp_(mv,x(:),y(:),1) ,'linestyle','--','linewidth',1.5,'color','r'); 
            axis([mv,mv]);
        end
        
        function [R] = pcorr_(A,B)
            % pearson's correlation coefficient
            xcor_ = @(A,B) numel(A)*A(:).'*B(:)-sum(A(:))*sum(B(:));
            R = xcor_(A,B)/(sqrt(xcor_(A,A))*sqrt(xcor_(B,B)));
        end
        
        function [C] = pfit_(x,y,n)
            % least squares fit polynomial of degree n to data
            C = (x(:).^[0:n])\y(:);
        end
        
        function [y] = peval_(x,C)
            n = numel(C)-1;
            % evaluate least squares fit polynomial of degree n to data
            y = (x(:).^[0:n])*C(:);
        end
        
        function [y] = pinterp_(x,x0,y0,n)
            % interpolate x0 vs y0 data at points x using a polynomial of degree n
            
            import am_lib.*
            
            y = peval_(x(:),pfit_(x0(:),y0(:),n));
        end

        
        % fft related
        
        function [r,gr] = fft_(k,gk)
            % interpolate gk on equidistant points
            kp= linspace(min(k),max(k),numel(k));
            g = interp1(k,gk,kp); k=kp; clear ke;
            % apply fft
            N = numel(k);
            r = fftshift(([0:(N-1)]'-floor(N/2))/(k(end)-k(1))); 
            gr= fft(g); gr = gr(1:floor(N/2)); r = r(1:floor(N/2));
        end
        
        function          plot_power_spectrum_(t,y,leg)
            import am_lib.*
            % count number of scans
            nys=size(y,2);
            % apply fft
            for i = 1:nys; [f(:,i),yf(:,i)] = fft_(t(:,i).',y(:,i).'); end
            % plot results
            figure(1); set(gcf,'color','w'); clf
            a=arrayfun(@(i){t(:,i),y(:,i)},[1:nys],'UniformOutput',false);a=[a{:}];
            axes('position',[0.1 0.1+0.5 0.85 0.35]); plot(a{:}); axis tight;
            xlabel('t'); ylabel('g(t)');
            a=arrayfun(@(i){1./f(:,i),abs(yf(:,i)).^2},[1:nys],'UniformOutput',false);a=[a{:}];
            axes('position',[0.1 0.1 0.85 0.35]); loglog(a{:}); axis tight;
            xlabel('1/f'); ylabel('abs( F[g(t)](f) )^2');
            if nargin<3
            legend(arrayfun(@(i){num2str(i)},[1:nys]),'location','southeast');
            else
            legend(leg,'location','southeast');
            end
        end
        
        
        % interpolation

        function [fq] = fftinterp_(f,q,n,algo)
            % fourier interpolate f(k) at points q; f must be periodic over [0,1)
            %
            % generate f using like this:
            %
            % mpgrid_ = @(N) [0:(N-1)]/N;
            % [xk,yk,zk]=meshgrid(mpgrid_(n(1)),mpgrid_(n(2)),mpgrid_(n(3)));
            % k = [xk(:),yk(:),zk(:)];
            % f = cos(2*pi*xk)+cos(2*pi*yk)+cos(2*pi*zk);
            %

            import am_lib.*
            
            % input formats:
            %
            %       f [nbands, kx*ky*kz ]   or   f [ kx,  ky,  kz ]   or   f [ nbands, kx, ky, kz ]
            %
            % the first format requires n to be input as well
            m = size(f);
            switch numel(m)
                case 2
                    if nargin ~= 3
                        error('fbz dimensions n are required');
                    else
                        n = [n,m(1)];
                    end
                case 3
                    n = [m,1];
                    f = reshape(f,1,[]);
                case 4
                    n = [m(2:4),m(1)];
                    f = reshape(f,m(1),[]);
            end
            
            % set default algorithm to cos interpolation
            if ~exist('algo','var'); algo = 'cos'; end

            % select an algorithm
            switch algo
                case 'fft'
                    % define function to get kernel
                    ifft_kernel_ = @(k,r,n) exp(+2i*pi*matmul_(reshape(k,1,3,[],1),reshape(r,3,1,1,[])))./sqrt(prod(n(1:3)));
                    % generate direct mesh r
                    m_ = @(i) [0:(n(i)-1)]-floor(n(i)/2); [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); r = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    % get fft expansion coefficients: c [ nrs, n ]
                    c = zeros(prod(n(1:3)),n(4)); f = permute(reshape(f,[n(4),n(1:3)]),[2,3,4,1]);
                    for j = 1:n(4); c(:,j) = reshape(fftshift(fftn( f(:,:,:,j) )),[],1) ./ sqrt(prod(n(1:3))); end
                    % interpolate f on q
                    fq = real( ifft_kernel_(q,r,n) * c ).';
                case 'dft'
                    % define function to get kernel
                     fft_kernel_ = @(k,r,n) exp(-2i*pi*matmul_(reshape(r,1,3,[],1),reshape(k,3,1,1,[])))./sqrt(prod(n(1:3)));
                    ifft_kernel_ = @(k,r,n) exp(+2i*pi*matmul_(reshape(k,1,3,[],1),reshape(r,3,1,1,[])))./sqrt(prod(n(1:3)));
                    % generate direct and reciprocal meshes r and k
                    m_ = @(i) [0:(n(i)-1)]-floor(n(i)/2); [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); r = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    m_ = @(i) [0:(n(i)-1)]./n(i);         [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); k = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    % get dft expansion coefficients: c [ nrs, n ]
                    c = fft_kernel_(k,r,n) * f.';
                    % interpolate f on q
                    fq = real( ifft_kernel_(q,r,n) * c ).';
                case 'cos'
                    % cosine transform kernel
                     cos_kernel_ = @(k,r,n) cos( 2*pi*matmul_(reshape(r,1,3,[],1),reshape(k,3,1,1,[])))./sqrt(prod(n(1:3)));
                    icos_kernel_ = @(k,r,n) cos( 2*pi*matmul_(reshape(k,1,3,[],1),reshape(r,3,1,1,[])))./sqrt(prod(n(1:3)));
                    % generate direct and reciprocal meshes r and k
                    m_ = @(i) [0:(n(i)-1)]-floor(n(i)/2); [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); r = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    m_ = @(i) [0:(n(i)-1)]./n(i);         [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); k = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    % get dft expansion coefficients: c [ nrs, n ]
                    c = cos_kernel_(k,r,n) * f.';
                    % interpolate f on q
                    fq = ( icos_kernel_(q,r,n) * c ).';
                case 'star'
                    % testing star function: DOES NOT WORK.
                    % DOES NOT WORK. DOES NOT WORK. DOES NOT WORK. 
                    % DOES NOT WORK. DOES NOT WORK. DOES NOT WORK. 
                    % DOES NOT WORK. DOES NOT WORK. DOES NOT WORK. 
                    % DOES NOT WORK. DOES NOT WORK. DOES NOT WORK. 
                    x=1;
                    % define function to get kernel
                     fft_kernel_ = @(k,r,n) exp(-2i*pi*matmul_(reshape(r,1,3,[],1),reshape(k,3,1,1,[])))./sqrt(prod(x*n(1:3)));
                    ifft_kernel_ = @(k,r,n) exp(+2i*pi*matmul_(reshape(k,1,3,[],1),reshape(r,3,1,1,[])))./sqrt(prod(x*n(1:3)));
                    % generate direct and reciprocal meshes r and k
                    m_ = @(i) [0:(x*n(i)-1)]-floor(x*n(i)/2); [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); r = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    m_ = @(i) [0:(n(i)-1)]./n(i);             [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); k = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    % conver to [cart]  
                    r = uc.bas*r;     r = uc2ws( r, diag(n(1:3))*uc.bas     );
                    k = fbz.recbas*k; k = uc2ws( k, diag(n(1:3))*fbz.recbas );
                    q = fbz.recbas*q; q = uc2ws( q, diag(n(1:3))*fbz.recbas );
                    % get symmetrically equivalent r-points (stars)
                    r4sym = uc2ws( uc.bas*r, diag(n(1:3))*uc.bas     );
                    [~,~,~,R] = get_symmetries(pc); R=matmul_(pc.bas,matmul_(R,inv(pc.bas))); nRs=size(R,3);
                    PM = member_(uc2ws(matmul_(R,r4sym),diag(n(1:3))*uc.bas),r4sym); A = get_connectivity(PM);
                    i2f = round(findrow_(A)).'; f2i = [1:numel(i2f)]*A; nstars=numel(i2f); 
            %         % get star function kernel
                    K=fft_kernel_(k,r,n); Ki=ifft_kernel_(q,r,n); Z = [f2i==[1:max(f2i)].']; K = Z*K; Ki=Ki/Z;


                    K = sum(exp(+2i*pi*matmul_(reshape(                  k, 1, 3, 1,[]    ), ....
                                               reshape(matmul_(R,r(:,i2f)), 3, 1,[], 1,nRs))) ,3)./nRs./sqrt(prod(x*n(1:3))) .* sum(f2i==[1:numel(i2f)].',2);

                    Ki= sum(exp(-2i*pi*matmul_(reshape(                  q, 1, 3,[], 1    ), ...
                                               reshape(matmul_(R,r(:,i2f)), 3, 1, 1,[],nRs))),3)./nRs./sqrt(prod(x*n(1:3)));% .* sqrt(sum(f2i==[1:numel(i2f)].',2)).';
            % r = matmul_(R,r(:,i2f(3)));
            %         % get factorization matrix (to accumlate rows of matrix)
            %         K=fft_kernel_(k,r,n); Ki=ifft_kernel_(q,r,n); %Z = [f2i==[1:max(f2i)].']; K = Z*K; Ki=Ki/Z;
                    % get dft expansion coefficients: c [ nrs, n ]
                    c = K * f.';
                    % interpolate f on q
                    fq = ( Ki * c ).';
                    
                    [~,a,b]=unique(rnd_(c),'rows','stable'); a=a.'; b=b.';
            end
            
        end

        function y    = linspacen_(d1, d2, n)
            n  = double(n); d1 = squeeze(d1); d2 = squeeze(d2);
            NDim = ndims(d1);
            if NDim==2 && any(size(d1)==1)
                NDim = NDim-1;
                if all(size(d1)==1)
                    NDim = 0;
                end
            end
            pp      = (0:n-2)./(floor(n)-1);
            Sum1    = tensor_product(d1, ones(1,n-1));
            Sum2    = tensor_product((d2-d1), pp);
            y = cat(NDim+1, Sum1  + Sum2, shiftdim(d2, size(d1, 1)==1 ));

            function Z = tensor_product(X,Y)
                sX=size(X);sY=size(Y); ndim1=ndims(X);ndim2=ndims(Y); indperm=[ndim2+1:ndim1+ndim2,1:ndim2];
                Z=squeeze(repmat(X,[ones(1,ndims(X)),sY]).*permute(repmat(Y,[ones(1,ndims(Y)),sX]),indperm));
            end
        end
        
        function n = gcd_(n)
            % operates on vector n = [n1,n2,n3 ...]
            x=1; p=n;
            while(size(n,2))>=2
                p= n(:,size(n,2)-1:size(n,2));
                n=n(1,1:size(n,2)-2);
                x=1;
                while(x~=0)
                    x= max(p)-min(p);
                    p = [x,min(p)];
                end    
                n = [n,max(p)];
                p = [];
            end
        end
    end

    % unix functions and scripts
    
    methods (Static)
        
        function [str] = verbatim_()
            %VERBATIM  Get the text that appears in the next comment block.
            %  Returns the text of the first comment block following the call.  The
            %  block comment delimiters, %{ and %}, must appear on lines by themselves
            %  (optionally preceded by white space).
            %
            %  If you want a final end-of-line character then leave a blank line before
            %  the %}.
            %
            %  If both comment delimiters are preceded by the same white space (same
            %  combination of spaces and tabs) then that white space will be deleted
            %  (if possible) from the beginning of each line of the commented text.
            %  This is so the whole block can be indented.
            %
            %  Example,
            %
            %      str = verbatim;
            %          %{
            %          This is the text
            %          that will be returned by verbatim.
            %          %}
            %
            %  VERBATIM can only be used in an m-file.

            % Get the function call stack.
            [dbs,thisWorkspace] = dbstack('-completenames');
            assert(length(dbs) > 1,'VERBATIM must be called from an M-file.')
            dbs = dbs(thisWorkspace + 1);
            lines = repmat({''},1,100);

            % Open the file.
            fid = fopen(dbs.file);

            try
                % Skip lines up to the current line in the calling function.
                textscan(fid,'%*s',0,'HeaderLines',dbs.line,'Delimiter','');

                % Read lines until one begins with '%{'.
                line = '';
                while ~isequal(line,-1) && isempty(regexp(line,'^\s*%{','once'))
                    line = fgetl(fid);
                end

                % Read and save lines until one begins with '%}'.  Leave first cell
                % empty.
                k = 1;
                while true
                    k = k + 1;
                    lines{k} = fgetl(fid);
                    if isequal(lines{k},-1) || ...
                            ~isempty(regexp(lines{k},'^\s*%}','once'))
                        break
                    end
                end

                % Close the file.
                fclose(fid);

            catch err
                % Close the file and rethrow the error.
                fclose(fid);
                rethrow(err)
            end

            % If white space preceeding '%{' and '%}' is the same then delete it from
            % the beginning of each line of the commented text so you can have indented
            % code.
            white_space1 = regexp(line,'^\s*','match','once');
            white_space2 = regexp(lines{k},'^\s*','match','once');
            if strcmp(white_space1,white_space2)
                lines = regexprep(lines,['^',white_space1],'');
            end

            % Construct the output string.
            str = [sprintf('%s\n',lines{2:k-2}),lines{k-1}];
        end
        
        
        % matlab-integrated
        
        function [n] = count_lines_(fname)
            if ispc
                [~,a] = system(sprintf('type %s | find /c /v ""',fname));
            elseif or(ismac,isunix)
                [~,a] = system(sprintf('wc -l %s',fname));
            end
            n = sscanf(a,'%i');
        end
        
        function hostName = get_host_()
            if ispc
                hostName = getenv('COMPUTERNAME');
            else
                hostName = getenv('HOSTNAME');
            end
        end
        
        function userName = get_user_()
            if ispc
                userName = getenv('username');
            else
                userName = getenv('USER');
            end
        end
        
        function rel_path = relativepath( tgt_path, act_path )

                if nargin<2; act_path = pwd; end
            
                % Predefine return string:
                rel_path = '';

                % Make sure strings end by a filesep character:
                if  isempty(act_path)   ||   ~isequal(act_path(end),filesep)
                   act_path = [act_path filesep];
                end
                if  isempty(tgt_path)   ||   ~isequal(tgt_path(end),filesep)
                   tgt_path = [tgt_path filesep];
                end

                % Convert to all lowercase:
                [act_path] = fileparts( lower(act_path) );
                [tgt_path] = fileparts( lower(tgt_path) );

                % Create a cell-array containing the directory levels:
                act_path_cell = pathparts(act_path);
                tgt_path_cell = pathparts(tgt_path);

                % If volumes are different, return absolute path:
                if  isempty(act_path_cell)   ||   isempty(tgt_path_cell)
                   return  % rel_path = ''
                else
                   if  ~isequal( act_path_cell{1} , tgt_path_cell{1} )
                      rel_path = tgt_path;
                      return
                   end
                end

                % Remove level by level, as long as both are equal:
                while  ~isempty(act_path_cell)   &&   ~isempty(tgt_path_cell)
                   if  isequal( act_path_cell{1}, tgt_path_cell{1} )
                      act_path_cell(1) = [];
                      tgt_path_cell(1) = [];
                   else
                      break
                   end
                end

                % As much levels down ('../') as levels are remaining in "act_path":
                for  i = 1 : length(act_path_cell)
                   rel_path = ['..' filesep rel_path];
                end

                % Relative directory levels to target directory:
                for  i = 1 : length(tgt_path_cell)
                   rel_path = [rel_path tgt_path_cell{i} filesep];
                end

                % Start with '.' or '..' :
                if  isempty(rel_path)
                   rel_path = ['.' filesep];
                elseif  ~isequal(rel_path(1),'.')
                   rel_path = ['.' filesep rel_path];
                end
                
                % clean up
                rel_path = strrep(rel_path,'//','/');

            function  path_cell = pathparts(path_str)
                path_str = [filesep path_str filesep];
                path_cell = {};
                sep_pos = findstr( path_str, filesep );
                for j = 1 : length(sep_pos)-1
                   path_cell{j} = path_str( sep_pos(j)+1 : sep_pos(j+1)-1 );
                end
            end
        end
        
    end
end








