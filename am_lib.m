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
        MPATH     = '/Applications/MATLAB_R2016b.app';
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

        function [C] = mod_(A)
            C = mod(A+am_lib.tiny,1)-am_lib.tiny;
        end

        function [C] = rnd_(A)
            C = round(A,-log10(am_lib.tiny));
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
            %     * B [ p,a,b,c,1,...]
            %   ` = C [(m,p),a,b,c,d,...]
            %
            % i = 2: add second dimension
            %       A [m,n]
            %     * B [m,p]
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

        % matching
        
        function [C] = maxabs_(A)
            C = max(abs(A(:)));
        end
        
        function [C] = minmax_(A)
            C = [min(A(:)),max(A(:))];
        end
        
        function [C] = sortc_(A)
            % column vector-based rank, sort, unique with numeric precision
            import am_lib.rnd_
            [C] = sortrows(rnd_(A).').'; 
        end

        function [C] = rankc_(A)
            % column vector-based rank, sort, unique with numeric precision
            import am_lib.rnd_
            C = sortrowsc(rnd_(A).',[1:size(A,1)]).'; 
        end                 

        function [I] = match_(A,B)
            % find the indicies I which permute column entries of A to match B
            % for example: A=[[1;4;3],[1;3;4],[3;1;4]]; B=[3;1;4]; X=match_(A,B); accessc_(A,X)-B
            import am_lib.*
            [~,iA]=sort(A); [~,iB]=sort(B); iB(iB)=[1:numel(iB)]; I = accessc_(iA,iB);
        end
        
        function [c] = member_(A,B)
            % get indicies of column vectors A(:,i,j) in matrix B(:,:)
            % 0 means not A(:,i) is not in B(:,:)

            [~,m,n] = size(A); c = zeros(m,n);
            for i = 1:m
            for j = 1:n
                c(i,j) = member_engine_(A(:,i,j),B,am_lib.tiny);
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
        
        function [C,inds] = uniquecol_(A)
            % get unique values with numeric precision
            import am_lib.rnd_
            [~,inds] = unique(rnd_(A).','rows','stable'); C=A(:,inds);
        end
        
        function [C] = uniquemask_(A)
            % returns a matrix with column vectors marking the first
            % occurance of a value in each column
            import am_lib.rnd_
            
            C = false(size(A));
            for i = 1:size(A,2)
                [~,b]=unique(rnd_(A(:,i))); C(b,i) = true;
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
        
        % special mathematical functions 
        
        function [C] = lorentz_(A)
            % define gaussian function 
            C = 1./(pi*(A.^2+1)); 
            
        end

        function [C] = gauss_(A)
            C = exp(-abs(A).^2)./sqrt(pi);
        end

        function [C] = delta_(A)
            % define tiny, kronecker delta, and heavside
            C = logical(abs(A)<am_lib.tiny); 
        end

        function [C] = heaviside_(A)
            C = logical(A>0);
        end

        
        % geometric functions
        
        function [A] = R_axis_(R)
            % define basic parameters and functions
            tiny = 1E-8; normalize_ = @(v) v/norm(v); 

            % check for identity
            if abs(trace(R)-3)< tiny; A=[0;0;1]; return; end

            % get rotation axis
            A = null(R-eye(3));

            % get random point on plane perpendicular to the rotation axis
            v1 = rand(3,1); v1 = normalize_(v1 - dot(v1,A)*A);

            % rotate point on the perpendicular plane
            v2 = R*v1;

            % get cross product (add tiny cross component to deal with 180 deg rotation)
            c = normalize_(cross(v1,v2+cross(v1,A)*tiny));

            % adjust sign
            A = sign(dot(c,A))*A;
        end

        function [A] = R_angle_(R)
            % define conversion of proper rotations to axis & angle representation
            A = acos((trace(R)-1)/2); 
        end

        function [Wtes] = get_wigner(j,R)

            import am_lib.get_wigner_engine

            % define tiny, kronecker delta, and heavside
            d_ = @(x,y) logical(x==y); t_ = @(x,y) logical(x>y);
            
            % matrix indices
            [m,mp]=meshgrid([j:-1:-j]);

            % define angular momentum operators (Jp raising, Jm lowering, ...)
            Jm = d_(m,mp+1).*sqrt((j+m).*(j-m+1)); Jp = d_(m,mp-1).*sqrt((j-m).*(j+m+1));
            Jx = (Jp+Jm)/2; Jy = (Jp-Jm)/2i; Jz = d_(m,mp).*m; J = cat(3,Jx,Jy,Jz);

            % define basis change: spherical to tesseral harmonics (real basis)
            T = d_(0,m) .* d_(mp,m) + ...
                t_(m,0) .* - sqrt(-1/2) .* ( (-1).^m.*d_(m,-mp) - d_(m, mp) ) + ...
                t_(0,m) .* - sqrt( 1/2) .* ( (-1).^m.*d_(m, mp) + d_(m,-mp) );

            % batch convert to wigner
            nRs = size(R,3); Wtes = zeros(2*j+1,2*j+1,nRs);
            for i = [1:nRs]; Wtes(:,:,i) = get_wigner_engine(J,T,j,R(:,:,i)); end
        end
        
        function [Wtes] = get_wigner_engine(J,T,j,R)
            % define wigner function for spherical (complex) and tesseral (real) harmonics
            % Note: for l = 1, Wtes_(R) = R

            import am_lib.*

            % get proper rotation
            d = det(R); dR = R*d;

            % get rotation axis and angle
            an = R_angle_(dR); ax = circshift(R_axis_(dR),1);

            % define spin-vector dot products [Eq. 1.37 Blundell]
            dotV_ = @(S,V) S(:,:,1)*V(1) + S(:,:,2)*V(2) + S(:,:,3)*V(3);

            if size(R,1)>9 % faster for symbolic and square matrices with dimensions > 9
                Wsph = expm( -sqrt(-1)*dotV_(J,ax)*an ) * d^j;
                Wtes = T' * Wsph * T;
            else % faster for numerical square matrices with dimensions < 9
                [V,D] = eig( -sqrt(-1) * dotV_(J,ax) * an); 
                Wsph = V*diag(exp(diag(D)))/V * d^j;
                Wtes = T' * Wsph * T ;
            end
        end
        
        
        % matrix related
        
        function [A] = merge_(A)

            [m,n] = size(A); tol=max(m,n)*eps(class(A))*norm(A,'inf');

            i = 1; j = 1; go=true;
            while go
                go = false;
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
                        % see which rows overlap with the i-th row
                        ex_ = (A(i,:)*A~=0);
                        % merge overlaps and zero all other rows
                        A(ex_,:)=0; A(i,ex_)=1;
                        i = i + 1;
                        j = j + 1;
                    end
                end
                % is one loop enough?! added go loop here to try and help.
                % maybe this is not right.
                AA = A*A';
                if any(any( diag(diag(AA))-AA >tol ,2),1)
                    go = true;
                end
            end
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
        
        
        % combinatorial 
        
        function [x] = nchoosek_(n,k)
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
        
        function [x] = nchoosrk_(n,k)
            % with duplications allowed
            
            import am_lib.nchoosek_
            
            x=nchoosek_(n+k-1,k).';
            x=x-repmat(0:k-1,size(x,1),1);
            x=x.';
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
        
        
        % general plotting
        
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
            C = (x(:).^[0:n])\y;
        end
        
        function [y] = peval_(x,C,n)
            % evaluate least squares fit polynomial of degree n to data
            y = (x(:).^[0:n])*C;
        end
        
        function [y] = pinterp_(x,x0,y0,n)
            % interpolate x0 vs y0 data at points x using a polynomial of degree n
            
            import am_lib.*
            
            y = peval_(x(:),pfit_(x0(:),y0(:),n),n);
        end
        
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

    end

    % unix functions and scripts
    
    methods (Static)
        
        % stand-alone
        
        function generate_scripts()
            
            import am_lib.*
            
            % initialize counter
            i=0; fprintf('Generating scripts:\n')
            
            % scripts functions
            fname='outcar2fp.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
            %{
            #!/bin/bash
            # Parse outcar to extract forces and displacements at each timestep.
            # Antonio Mei Nov/2014
            # Antonio Mei Jan/2017
            usage_ () {
                echo "Creates infile.force_position based on supplied outcar files."
                echo ""
                echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <natoms> [-c <compress_name>]"
                echo ""
                echo "Example: $0 -f -t -o \"\$(find . -name "OUTCAR*" | grep 4.00-300)\" -n 250"
                echo ""
                echo "-h : prints this message"
                echo "-n : [REQUIRED] number of atoms in the simulation cell"
                echo "-o : [REQUIRED] list of outcar files to parse"
                echo "-t : trims the last md run (useful for removing runs which have not completed)"
                echo "-f : overwrites existing infile.force_position"
                echo "-c : compresses infile.force_position to a tar.gz file"
                echo ""
                echo "infile.force_position file contents:"
                echo "   x position   y position   z position     x force      y force      z force"
                exit 1
            }
            
            main_ () {
                # trim the last md run which may not have completed
                trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
                # get position and forces
                get_  () { cat $2 | grep -h -A $(($1+1)) POSITION  ; }
                # cut header lines
                cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/POSITION/d' ; }
                # compress produced infile.force_position
                compress_ () { tar -zcvf infile.force_position.tar.gz infile.force_position ; }
                #
                if ${ISFORCE}; then
                    if [ -f "./infile.force_position" ]; then
                        rm ./infile.force_position
                        printf " ... ./infile.force_position overwritten\n"
                    fi
                fi
                #
                if ${ISTRIM}; then
                    printf " ... trim:\n"
                    for F in "${FLIST}"; do
                        printf " ...     %-100s\n" "${F}"
                        trim_ ${F} | get_ ${NATOMS} | cut_ >> infile.force_position
                    done
                else
                    printf " ... batch parsing without trim\n"
                    get_ ${NATOMS} "${FLIST}" | cut_ >> infile.force_position
                fi
                #
                printf " ... infile.force_position created\n"
                #
                if ${ISCOMPRESS}; then
                    printf " ... infile.force_position.tar.gz compressed\n"
                    compress_ 
                fi
            }
            
            ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
            if (($# == 0)); then usage_; exit 1; fi
            while getopts "n:o:htfc" o; do
                case "${o}" in
                    o)  FLIST=${OPTARG} ;;
                    n)  NATOMS=${OPTARG} ;;
                    c)  ISCOMPRESS=true ;;
                    t)  ISTRIM=true ;;
                    f)  ISFORCE=true ;;
                    h)  usage_; exit 0 ;;
                    *)  usage_; exit 1 ;;
                esac
            done
            main_
            %}
            
            fname='outcar2en.sh'; fid=fopen([fname],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); fprintf(' ... %s (succeeded)\n',fname);
            %{
            #!/bin/bash
            # Preprocess outcar to remove the last run, which may not have finished
            # Antonio Mei Nov/2014
            # Antonio Mei Jan/2017
            usage_ () {
                echo "Creates infile.electron_energies based on supplied outcar files."
                echo ""
                echo "Usage: $0 [-h] [-t] [-f] -o <outcar_list> -n <nbands> [-c <compress_name>]"
                echo ""
                echo "Example: $0 -f -t -o \"\$(find . -name "OUTCAR*" | grep 4.00-300)\" -n 751"
                echo ""
                echo "-h : prints this message"
                echo "-n : [REQUIRED] number of bands"
                echo "-o : [REQUIRED] list of outcar files to parse"
                echo "-t : trims the last md run (useful for removing runs which have not completed)"
                echo "-f : overwrites existing infile.electron_energies"
                echo "-c : compresses infile.electron_energies to a tar.gz file"
                echo ""
                echo "infile.electron_energies file contents:"
                echo "   n index    En energy    fn occupation"
                exit 1
            }
            main_ () {
                # trim the last md run which may not have completed
                trim_ () { tac $1 | awk '!found && /POSITION/{found=1;next}1' | tac ; }
                # get energies
                get_  () { cat $2 | grep -h -A ${1} occupation  ; }
                # cut header lines
                cut_  () { cat $1 | sed '/^--$/d' | sed '/--/d' | sed '/occupation/d' ; }
                # compress produced infile.electron_energies
                compress_ () { tar -zcvf infile.electron_energies.tar.gz infile.electron_energies ; }
                #
                if ${ISFORCE}; then
                    if [ -f "./infile.electron_energies" ]; then
                        rm ./infile.electron_energies
                        printf " ... ./infile.electron_energies overwritten\n"
                    fi
                fi
                # 
                if ${ISTRIM}; then
                    printf " ... trim:\n"
                    for F in "${FLIST}"; do
                        printf " ...     %-100s\n" "${F}"
                        trim_ ${F} | get_ ${NBANDS} | cut_ >> infile.electron_energies
                    done
                else
                    printf " ... batch parsing without trim\n"
                    get_ ${NBANDS} "${FLIST}" | cut_ >> infile.electron_energies
                fi
                #
                awk '{ print $2 }' infile.electron_energies > infile.electron_energies.tmp && mv infile.electron_energies.tmp infile.electron_energies
                #
                printf " ... infile.electron_energies created\n"
                #
                if ${ISCOMPRESS}; then
                    printf " ... infile.electron_energies.tar.gz compressed\n"
                    compress_ 
                fi
            }
            ISCOMPRESS=false; ISTRIM=false; ISFORCE=false;
            if (($# == 0)); then usage_; exit 1; fi
            while getopts "n:o:htfc" o; do
                case "${o}" in
                    o)  FLIST=${OPTARG} ;;
                    n)  NBANDS=${OPTARG} ;;
                    c)  ISCOMPRESS=true ;;
                    t)  ISTRIM=true ;;
                    f)  ISFORCE=true ;;
                    h)  usage_; exit 0 ;;
                    *)  usage_; exit 1 ;;
                esac
            done
            main_
            %}
            
                    
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
        
    end    
    

    % mex functions
    
    methods (Static)
        
        function compile_mex()
            % compile mex functions
            %
            % to compile the library manually
            %
            %       fort -O3 -parallel -fpp -c am_mex_lib.f90
            %
            % to compile individual interfaces manually
            %
            %       ifort -O3 -parallel -fpp -fPIC -L/Applications/MATLAB_R2016b.app/bin/maci64
            %       -I/Applications/MATLAB_R2016b.app/extern/include -lmx -lmex -lmat -nofor_main 
            %       -bundle ./uc2ws_mex.f90 -o uc2ws_mex.mexmaci64    
            %
            
            import am_lib.*
            
            % mex library
            flib = 'am_mex_lib'; fid=fopen([flib,'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
module am_mex_lib

    implicit none

    public

    !  $$$$$$\   $$$$$$\  $$\   $$\  $$$$$$\ $$$$$$$$\  $$$$$$\  $$\   $$\ $$$$$$$$\  $$$$$$\
    ! $$  __$$\ $$  __$$\ $$$\  $$ |$$  __$$\\__$$  __|$$  __$$\ $$$\  $$ |\__$$  __|$$  __$$\
    ! $$ /  \__|$$ /  $$ |$$$$\ $$ |$$ /  \__|  $$ |   $$ /  $$ |$$$$\ $$ |   $$ |   $$ /  \__|
    ! $$ |      $$ |  $$ |$$ $$\$$ |\$$$$$$\    $$ |   $$$$$$$$ |$$ $$\$$ |   $$ |   \$$$$$$\
    ! $$ |      $$ |  $$ |$$ \$$$$ | \____$$\   $$ |   $$  __$$ |$$ \$$$$ |   $$ |    \____$$\
    ! $$ |  $$\ $$ |  $$ |$$ |\$$$ |$$\   $$ |  $$ |   $$ |  $$ |$$ |\$$$ |   $$ |   $$\   $$ |
    ! \$$$$$$  | $$$$$$  |$$ | \$$ |\$$$$$$  |  $$ |   $$ |  $$ |$$ | \$$ |   $$ |   \$$$$$$  |
    !  \______/  \______/ \__|  \__| \______/   \__|   \__|  \__|\__|  \__|   \__|    \______/

    integer    , parameter :: sp      = kind(0.0E0) !> single precision
    integer    , parameter :: dp      = kind(0.0D0) !> double precision
    real(dp)   , parameter :: tiny    =  1.0D-5
    real(dp)   , parameter :: halfpi  =  1.570796326794897_dp
    real(dp)   , parameter :: sqrtpi  =  1.772453850905516_dp
    real(dp)   , parameter :: pi      =  3.141592653589793_dp
    real(dp)   , parameter :: twopi   =  6.283185307179586_dp
    real(dp)   , parameter :: fourpi  = 12.566370614359172_dp
    complex(dp), parameter :: cmplx_i = cmplx(0,1,dp)
    complex(dp), parameter :: itwopi  = cmplx_i*twopi

    ! $$$$$$\ $$\   $$\ $$$$$$$$\ $$$$$$$$\ $$$$$$$\  $$$$$$$$\  $$$$$$\   $$$$$$\  $$$$$$$$\  $$$$$$\
    ! \_$$  _|$$$\  $$ |\__$$  __|$$  _____|$$  __$$\ $$  _____|$$  __$$\ $$  __$$\ $$  _____|$$  __$$\
    !   $$ |  $$$$\ $$ |   $$ |   $$ |      $$ |  $$ |$$ |      $$ /  $$ |$$ /  \__|$$ |      $$ /  \__|
    !   $$ |  $$ $$\$$ |   $$ |   $$$$$\    $$$$$$$  |$$$$$\    $$$$$$$$ |$$ |      $$$$$\    \$$$$$$\
    !   $$ |  $$ \$$$$ |   $$ |   $$  __|   $$  __$$< $$  __|   $$  __$$ |$$ |      $$  __|    \____$$\
    !   $$ |  $$ |\$$$ |   $$ |   $$ |      $$ |  $$ |$$ |      $$ |  $$ |$$ |  $$\ $$ |      $$\   $$ |
    ! $$$$$$\ $$ | \$$ |   $$ |   $$$$$$$$\ $$ |  $$ |$$ |      $$ |  $$ |\$$$$$$  |$$$$$$$$\ \$$$$$$  |
    ! \______|\__|  \__|   \__|   \________|\__|  \__|\__|      \__|  \__| \______/ \________| \______/

    interface foursort
        module procedure d_foursort, i_foursort
    end interface ! foursort

    interface fourrank
        module procedure d_fourrank, i_fourrank
    end interface ! fourrank

    contains

    !  $$$$$$\ $$$$$$$$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\   $$$$$$\
    ! $$  __$$\\__$$  __|$$  __$$\\__$$  __|\_$$  _|$$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$  __$$\
    ! $$ /  \__|  $$ |   $$ /  $$ |  $$ |     $$ |  $$ /  \__|  $$ |     $$ |  $$ /  \__|$$ /  \__|
    ! \$$$$$$\    $$ |   $$$$$$$$ |  $$ |     $$ |  \$$$$$$\    $$ |     $$ |  $$ |      \$$$$$$\
    !  \____$$\   $$ |   $$  __$$ |  $$ |     $$ |   \____$$\   $$ |     $$ |  $$ |       \____$$\
    ! $$\   $$ |  $$ |   $$ |  $$ |  $$ |     $$ |  $$\   $$ |  $$ |     $$ |  $$ |  $$\ $$\   $$ |
    ! \$$$$$$  |  $$ |   $$ |  $$ |  $$ |   $$$$$$\ \$$$$$$  |  $$ |   $$$$$$\ \$$$$$$  |\$$$$$$  |
    !  \______/   \__|   \__|  \__|  \__|   \______| \______/   \__|   \______| \______/  \______/

    pure function factorial(n) result(y)
        !
        implicit none
        !
        integer, intent(in) :: n
        real(dp) :: y
        !
        if (n.ge.1) then
            y = product([1:n])
        else
            y = 1
        endif
    end function  factorial

    pure function nchoosek(n,k) result (res)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, intent(in) :: k
        integer :: res
        !
        res=factorial(n)/(factorial(k)*factorial(n-k))
        !
    end function  nchoosek

    function      perms(n) result(PT)
        !
        implicit none
        !
        integer, intent(in) :: n
        integer, allocatable :: P(:,:)
        integer, allocatable :: PT(:,:)
        integer :: m
        !
        m = product([1:n])
        !
        allocate(P(m,n))
        !
        call permutate([1:n],P)
        !
        allocate(PT,source=transpose(P))
        !
        contains
        recursive subroutine permutate(E, P)
            !
            implicit none
            !
            integer, intent(in)  :: E(:) ! array of objects 
            integer, intent(out) :: P(:,:) ! permutations of E 
            integer :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
            N = size(E); Nfac = size(P,1); 
            do i = 1, N
              if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S) 
              forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/) 
            enddo 
        end subroutine permutate 
    end function  perms

    ! $$\   $$\ $$\   $$\ $$\      $$\ $$$$$$$\  $$$$$$$$\ $$$$$$$\        $$$$$$$$\ $$\   $$\ $$$$$$$$\  $$$$$$\  $$$$$$$\ $$\     $$\
    ! $$$\  $$ |$$ |  $$ |$$$\    $$$ |$$  __$$\ $$  _____|$$  __$$\       \__$$  __|$$ |  $$ |$$  _____|$$  __$$\ $$  __$$\\$$\   $$  |
    ! $$$$\ $$ |$$ |  $$ |$$$$\  $$$$ |$$ |  $$ |$$ |      $$ |  $$ |         $$ |   $$ |  $$ |$$ |      $$ /  $$ |$$ |  $$ |\$$\ $$  /
    ! $$ $$\$$ |$$ |  $$ |$$\$$\$$ $$ |$$$$$$$\ |$$$$$\    $$$$$$$  |         $$ |   $$$$$$$$ |$$$$$\    $$ |  $$ |$$$$$$$  | \$$$$  /
    ! $$ \$$$$ |$$ |  $$ |$$ \$$$  $$ |$$  __$$\ $$  __|   $$  __$$<          $$ |   $$  __$$ |$$  __|   $$ |  $$ |$$  __$$<   \$$  /
    ! $$ |\$$$ |$$ |  $$ |$$ |\$  /$$ |$$ |  $$ |$$ |      $$ |  $$ |         $$ |   $$ |  $$ |$$ |      $$ |  $$ |$$ |  $$ |   $$ |
    ! $$ | \$$ |\$$$$$$  |$$ | \_/ $$ |$$$$$$$  |$$$$$$$$\ $$ |  $$ |         $$ |   $$ |  $$ |$$$$$$$$\  $$$$$$  |$$ |  $$ |   $$ |
    ! \__|  \__| \______/ \__|     \__|\_______/ \________|\__|  \__|         \__|   \__|  \__|\________| \______/ \__|  \__|   \__|

    pure function primes(nprimes)
        !
        ! naive approach to generating the first n prime numbers
        !
        implicit none
        !
        integer, intent(in)  :: nprimes
        integer, allocatable :: primes(:) ! array that will hold the primes
        integer :: at, found, i
        logical :: is_prime
        !
        allocate (primes(nprimes))
        !
        primes(1) = 2
        at = 2
        found = 1
        do
            is_prime = .true. ! assume prime
            do i = 1, found
                if (modulo(at,primes(i)).eq.0) then ! if divisible by any other element
                    is_prime = .false.               ! in the array, then not prime.
                    at = at + 1
                    continue
                end if
            end do
            found = found + 1
            primes(found) = at
            at = at + 1
            if (found == nprimes) then ! stop when all primes are found
                exit
            endif
        end do
        !
    end function  primes

    !  $$$$$$\   $$$$$$\  $$$$$$$\ $$$$$$$$\ $$$$$$\ $$\   $$\  $$$$$$\
    ! $$  __$$\ $$  __$$\ $$  __$$\\__$$  __|\_$$  _|$$$\  $$ |$$  __$$\
    ! $$ /  \__|$$ /  $$ |$$ |  $$ |  $$ |     $$ |  $$$$\ $$ |$$ /  \__|
    ! \$$$$$$\  $$ |  $$ |$$$$$$$  |  $$ |     $$ |  $$ $$\$$ |$$ |$$$$\
    !  \____$$\ $$ |  $$ |$$  __$$<   $$ |     $$ |  $$ \$$$$ |$$ |\_$$ |
    ! $$\   $$ |$$ |  $$ |$$ |  $$ |  $$ |     $$ |  $$ |\$$$ |$$ |  $$ |
    ! \$$$$$$  | $$$$$$  |$$ |  $$ |  $$ |   $$$$$$\ $$ | \$$ |\$$$$$$  |
    !  \______/  \______/ \__|  \__|  \__|   \______|\__|  \__| \______/


    pure function d_foursort(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x(4)
        real(dp) :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  d_foursort

    pure function i_foursort(x) result(y)
        !
        implicit none
        !
        integer, intent(in) :: x(4)
        integer :: y(4)
        !
        y = x
        if (y(1).gt.y(2)) then; y([1,2]) = y([2,1]); endif
        if (y(3).gt.y(4)) then; y([3,4]) = y([4,3]); endif
        if (y(1).gt.y(3)) then; y([1,3]) = y([3,1]); endif
        if (y(2).gt.y(4)) then; y([2,4]) = y([4,2]); endif
        if (y(2).gt.y(3)) then; y([2,3]) = y([3,2]); endif
        !
    end function  i_foursort

    pure function d_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        real(dp),intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !    
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else 
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  d_fourrank

    pure function i_fourrank(n) result(ind)
        ! borrowed from olle, who took it from stack overflow
        implicit none
        !
        integer,intent(in) :: n(4)
        integer :: ind(4)
        integer :: low1,high1,low2,high2,highest,lowest,middle1,middle2
        !    
        if ( n(1) <= n(2) ) then
            low1 = 1
            high1 = 2
        else 
            low1 = 2
            high1 = 1
        endif

        if ( n(3) <= n(4) ) then
            low2 = 3
            high2 = 4
        else
            low2 = 4
            high2 = 3
        endif

        if ( n(low1) <= n(low2) ) then
            lowest = low1
            middle1 = low2
        else
            lowest = low2
            middle1 = low1
        endif

        if ( n(high1) >= n(high2) ) then
            highest = high1
            middle2 = high2
        else
            highest = high2
            middle2 = high1
        endif

        if ( n(middle1) < n(middle2) ) then
            ind=(/lowest,middle1,middle2,highest/)
        else
            ind=(/lowest,middle2,middle1,highest/)
        endif
        !
    end function  i_fourrank

    pure function d_rank(v) result(ind)
      !
      ! based on dlasrt2
      !
      implicit none
      !
      real(dp), intent(in) :: v(:)
      real(dp), allocatable :: d(:)
      integer, allocatable :: ind(:)
      integer :: select, endd, i, j, start, stkpnt, tmpind, stack(2,32), n
      real(dp) :: d1, d2, d3, dmnmx, tmp
      !
      select = 20
      !
      n = size(v)
      !
      allocate(d,source=v)
      !
      allocate(ind,source=[1:n])
      !
      if (n.le.1) return
      !
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
      !
      do
         start = stack(1, stkpnt)
         endd = stack(2, stkpnt)
         stkpnt = stkpnt - 1
         if(endd-start.gt.0) then
            do i = start + 1, endd
               do j = i, start + 1, -1
                  if(d(j).lt.d(j-1)) then
                     dmnmx = d(j)
                     d(j) = d(j-1)
                     d(j-1) = dmnmx
                     tmpind = ind(j)
                     ind(j) = ind(j-1)
                     ind(j-1) = tmpind
                  else
                     exit
                  end if
               enddo
            enddo
         else if(endd-start.gt.select) then
            d1 = d(start)
            d2 = d(endd)
            i = (start+endd) / 2
            d3 = d(i)
            if(d1.lt.d2) then
               if(d3.lt.d1) then
                  dmnmx = d1
               else if(d3.lt.d2) then
                  dmnmx = d3
               else
                  dmnmx = d2
               end if
            else
               if(d3.lt.d2) then
                  dmnmx = d2
               else if(d3.lt.d1) then
                  dmnmx = d3
               else
                  dmnmx = d1
               end if
            end if
            i = start - 1
            j = endd + 1
            do
               do
                  j = j - 1
                  if (.not.(d(j).gt.dmnmx)) exit
               enddo
               do
                  i = i + 1
                  if (.not.(d(i).lt.dmnmx)) exit
               enddo
               if(i.gt.j) then
                  exit   
               else 
                  tmp = d(i)
                  d(i) = d(j)
                  d(j) = tmp
                  tmpind = ind(j)
                  ind(j) = ind(i)
                  ind(i) = tmpind
               end if
            enddo
            if(j-start.gt.endd-j-1) then
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
            else
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
            end if
         end if
         if (.not.(stkpnt.gt.0)) exit
      enddo
    end function  d_rank

    pure function i_rank(v) result(ind)
      !
      ! based on dlasrt2
      !
      implicit none
      !
      integer, intent(in) :: v(:)
      integer, allocatable :: d(:)
      integer, allocatable :: ind(:)
      integer :: select, endd, i, j, start, stkpnt, tmpind, stack(2,32), n
      integer :: d1, d2, d3, dmnmx, tmp
      !
      select = 20
      !
      n = size(v)
      !
      allocate(d,source=v)
      !
      allocate(ind,source=[1:n])
      !
      if (n.le.1) return
      !
      stkpnt = 1
      stack(1,1) = 1
      stack(2,1) = n
      !
      do
         start = stack(1, stkpnt)
         endd = stack(2, stkpnt)
         stkpnt = stkpnt - 1
         if(endd-start.gt.0) then
            do i = start + 1, endd
               do j = i, start + 1, -1
                  if(d(j).lt.d(j-1)) then
                     dmnmx = d(j)
                     d(j) = d(j-1)
                     d(j-1) = dmnmx
                     tmpind = ind(j)
                     ind(j) = ind(j-1)
                     ind(j-1) = tmpind
                  else
                     exit
                  end if
               enddo
            enddo
         else if(endd-start.gt.select) then
            d1 = d(start)
            d2 = d(endd)
            i = (start+endd) / 2
            d3 = d(i)
            if(d1.lt.d2) then
               if(d3.lt.d1) then
                  dmnmx = d1
               else if(d3.lt.d2) then
                  dmnmx = d3
               else
                  dmnmx = d2
               end if
            else
               if(d3.lt.d2) then
                  dmnmx = d2
               else if(d3.lt.d1) then
                  dmnmx = d3
               else
                  dmnmx = d1
               end if
            end if
            i = start - 1
            j = endd + 1
            do
               do
                  j = j - 1
                  if (.not.(d(j).gt.dmnmx)) exit
               enddo
               do
                  i = i + 1
                  if (.not.(d(i).lt.dmnmx)) exit
               enddo
               if(i.gt.j) then
                  exit   
               else 
                  tmp = d(i)
                  d(i) = d(j)
                  d(j) = tmp
                  tmpind = ind(j)
                  ind(j) = ind(i)
                  ind(i) = tmpind
               end if
            enddo
            if(j-start.gt.endd-j-1) then
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
            else
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = j + 1
               stack(2, stkpnt) = endd
               stkpnt = stkpnt + 1
               stack(1, stkpnt) = start
               stack(2, stkpnt) = j
            end if
         end if
         if (.not.(stkpnt.gt.0)) exit
      enddo
    end function  i_rank

    ! $$\      $$\  $$$$$$\ $$$$$$$$\ $$\   $$\       $$$$$$$$\ $$\   $$\ $$\   $$\  $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\  $$$$$$\
    ! $$$\    $$$ |$$  __$$\\__$$  __|$$ |  $$ |      $$  _____|$$ |  $$ |$$$\  $$ |$$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$$\  $$ |$$  __$$\
    ! $$$$\  $$$$ |$$ /  $$ |  $$ |   $$ |  $$ |      $$ |      $$ |  $$ |$$$$\ $$ |$$ /  \__|  $$ |     $$ |  $$ /  $$ |$$$$\ $$ |$$ /  \__|
    ! $$\$$\$$ $$ |$$$$$$$$ |  $$ |   $$$$$$$$ |      $$$$$\    $$ |  $$ |$$ $$\$$ |$$ |        $$ |     $$ |  $$ |  $$ |$$ $$\$$ |\$$$$$$\
    ! $$ \$$$  $$ |$$  __$$ |  $$ |   $$  __$$ |      $$  __|   $$ |  $$ |$$ \$$$$ |$$ |        $$ |     $$ |  $$ |  $$ |$$ \$$$$ | \____$$\
    ! $$ |\$  /$$ |$$ |  $$ |  $$ |   $$ |  $$ |      $$ |      $$ |  $$ |$$ |\$$$ |$$ |  $$\   $$ |     $$ |  $$ |  $$ |$$ |\$$$ |$$\   $$ |
    ! $$ | \_/ $$ |$$ |  $$ |  $$ |   $$ |  $$ |      $$ |      \$$$$$$  |$$ | \$$ |\$$$$$$  |  $$ |   $$$$$$\  $$$$$$  |$$ | \$$ |\$$$$$$  |
    ! \__|     \__|\__|  \__|  \__|   \__|  \__|      \__|       \______/ \__|  \__| \______/   \__|   \______| \______/ \__|  \__| \______/

    ! special functions

    pure function legendre(l,m,x) result(y)
        !
        ! Associated legrende polynomial
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England?; New York, 1992), p 246.
        !
        implicit none
        !
        integer , intent(in) :: l,m
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer  :: ll
        real(dp) :: pll(size(x))
        real(dp) :: pmm(size(x))
        real(dp) :: pmmp1(size(x))
        real(dp) :: somx2(size(x))
        !
        allocate(y(size(x)))
        !
        if ((m.lt.0).or.(m.gt.l).or.(all(abs(x).gt.1.0_dp))) then
            y = 0
        endif
        !
        pmm=1.0
        if (m > 0) then
            somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
            pmm=product(arth(1.0_dp,2.0_dp,m))*somx2**m
            if (mod(m,2) == 1) pmm=-pmm
        end if
        if (l == m) then
            y=pmm
        else
            pmmp1=x*(2*m+1)*pmm
            if (l == m+1) then
                y=pmmp1
            else
                do ll=m+2,l
                    pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
                    pmm=pmmp1
                    pmmp1=pll
                end do
                y=pll
            end if
        end if
        contains
        pure function arth(first,increment,n)
            !
            implicit none
            !
            real(dp), intent(in) :: first
            real(dp), intent(in) :: increment
            integer , intent(in) :: n
            real(dp) :: arth(n)
            integer  :: k,k2
            real(dp) :: temp
            !
            if (n > 0) arth(1)=first
            if (n <= 16) then
                do k=2,n
                    arth(k)=arth(k-1)+increment
                enddo
            else
                do k=2,8
                    arth(k)=arth(k-1)+increment
                end do
                temp=increment*8
                k=8
                do
                    if (k >= n) exit
                    k2=k+k
                    arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
                    temp=temp+temp
                    k=k2
                enddo
            end if
        end function arth
    end function  legendre

    pure function laguerre(k,p,x) result(y)
        ! associated Laguerre polynomial L_k^p(x)(k,p,x)
        ! Using the expression from Samuel Shaw Ming Wong "Computational Methods in Physics and Engineering", Eq 4-86 p 139
        ! Note, there is a typographical error on http://mathworld.wolfram.com/AssociatedLaguerrePolynomial.html Eq 10
        ! Also see Linus Pauling "Introuction to Quantum Mechancs"
        implicit none 
        !
        integer , intent(in) :: k, p
        real(dp), intent(in) :: x(:)
        real(dp), allocatable :: y(:)
        integer :: j
        !
        allocate(y,mold=x)
        y=0.0_dp
        !
        do j = 0, k
            y = y + nchoosek(k+p,k-j) * (-x)**j / factorial(j)
        enddo
        !
    end function  laguerre

    pure function Ylm(l,m,theta,phi)
        !
        ! computes spherical harmonics. Theta and phi in radians. size(theta) must equal size(phi); together they define the radial coordinates.
        !
        ! W. H. Press, B. P. Flannery, S. A. Teukolsky, and W. T. Vetterling, Numerical Recipes in Fortran 77: The Art
        ! of Scientific Computing, 2 edition (Cambridge University Press, Cambridge England?; New York, 1992), p 246.
        ! 
        implicit none
        !
        integer , intent(in) :: l, m
        real(dp), intent(in) :: theta(:), phi(:)
        complex(dp), allocatable :: Ylm(:)
        integer :: n
        ! 
        n = size(theta)
        !
        allocate(Ylm(n))
        !
        Ylm = sqrt( real(2*l+1,dp)/fourpi * factorial(l-m)/real(factorial(l+m),dp) ) * legendre(l,m,cos(theta)) * exp(cmplx_i*m*phi)
        !
    end function  Ylm

    pure function heavi(m)
        !
        ! A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 15
        !
        implicit none
        !
        integer, intent(in) :: m
        real(dp) :: heavi
        !
        if (m.ge.0) then
            heavi = 1.0_dp
        else
            heavi = 0.0_dp
        endif
        !
    end function  heavi

    ! heaviside theta functions

    pure function fermi_dirac(x) result(y)
        !
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Fermi-Dirac smearing
        if (x.lt.-maxarg) then
            y = 0.0_dp
        elseif (x.gt.maxarg) then
            y = 1.0_dp
        else
            y = 1.0_dp/(1.0_dp + exp(-x))
        endif
    end function  fermi_dirac
    
    pure function methfessel_paxton(x,n) result(y)
        ! theta function : PRB 40, 3616 (1989).
        implicit none
        !
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y
        
        real(dp) :: a, hp, arg, hd
        integer :: i, ni
        real(dp), parameter :: maxarg = 200.0_dp
        
        ! Methfessel-Paxton
        y = gauss_freq(x * sqrt(2.0_dp) )
        if (n.eq.0) return
        hd = 0.0_dp
        arg = min(maxarg, x**2)
        hp = exp(- arg)
        ni = 0
        a = 1.0_dp / sqrt (pi)
        do i=1,n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a = -a/(real(i,dp)*4.0_dp)
            y = y-a*hd
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
        enddo
        contains
        pure function gauss_freq(x)
            !     gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
            implicit none
            !
            real(dp),intent(in) :: x
            real(dp)            :: gauss_freq
            !
            gauss_freq = 0.5_dp * erfc(-x*0.7071067811865475_dp)
            !
        end function  gauss_freq
    end function  methfessel_paxton

    pure function marzari_vanderbilt(x) result(y)
        ! theta function: PRL 82, 3296 (1999)
        ! 1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
        implicit none
        !
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg, xp
        real(dp) :: maxarg
        maxarg = 200.0_dp
        ! Cold smearing
         xp = x - 1.0_dp / sqrt (2.0_dp)
         arg = min (maxarg, xp**2)
         y = 0.5d0 * erf (xp) + 1.0_dp / sqrt (2.0_dp * pi) * exp (- arg) + 0.5d0
    end function  marzari_vanderbilt

    pure function erf(x)
      !---------------------------------------------------------------------
      !
      !     Error function - computed from the rational approximations of
      !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
      !
      !     for abs(x) le 0.47 erf is calculated directly
      !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
      !
      implicit none  
      real(dp), intent(in) :: x
      real(dp) :: x2, p1(4), q1(4)
      real(dp) :: erf
      !
      p1 = [2.426679552305318E2_dp, 2.197926161829415E1_dp, 6.996383488619136_dp,  -3.560984370181538E-2_dp]
      q1 = [2.150588758698612E2_dp, 9.116490540451490E1_dp, 1.508279763040779E1_dp, 1.000000000000000_dp]
      !
      if (abs(x).gt.6.0_dp) then  
         !  erf(6)=1-10^(-17) cannot be distinguished from 1
         erf = sign (1.0_dp, x)  
      else  
         if (abs(x).le.0.47_dp) then  
            x2 = x**2  
            erf = x*(p1(1)+x2*(p1(2)+x2*(p1(3)+x2*p1(4))))/(q1(1)+x2*(q1(2)+x2*(q1(3)+x2*q1(4))))
         else  
            erf = 1.0_dp - erfc(x)  
         endif
      endif
    end function  erf

    pure function erfc(x)
      !     erfc(x) = 1-erf(x)  - See comments in erf
      implicit none  
      real(dp),intent(in) :: x
      real(dp)            :: erfc
      real(dp) :: ax, xm2, p2(8), q2(8), p3(5), q3(5), pim1
      !
      p2 = [ 3.004592610201616E2_dp,  4.519189537118719E2_dp,  3.393208167343437E2_dp,  1.529892850469404E2_dp,  4.316222722205674E1_dp,  7.211758250883094_dp,    5.641955174789740E-1_dp,-1.368648573827167E-7_dp]
      q2 = [ 3.004592609569833E2_dp,  7.909509253278980E2_dp,  9.313540948506096E2_dp,  6.389802644656312E2_dp,  2.775854447439876E2_dp,  7.700015293522947E1_dp,  1.278272731962942E1_dp,  1.000000000000000_dp]
      p3 = [-2.996107077035422E-3_dp,-4.947309106232507E-2_dp,  -2.269565935396869E-1_dp,-2.786613086096478E-1_dp,  -2.231924597341847E-2_dp]
      q3 = [ 1.062092305284679E-2_dp, 1.913089261078298E-1_dp,  1.051675107067932_dp,    1.987332018171353_dp,     1.000000000000000_dp]
      pim1 = 0.56418958354775629_dp  ! sqrt(1/pi)
      ax = abs(x)  
      if (ax > 26.0_dp) then  
         !  erfc(26.0)=10^(-296); erfc(9.0)=10^(-37);
         erfc = 0.0_dp  
      elseif (ax.gt.4.0_dp) then
         xm2=(1.0_dp/ax)**2
         erfc=(1.0_dp/ax)*exp(-2.0_dp)*(pim1+xm2*(p3(1)+xm2*(p3(2)+xm2*(p3(3)+xm2*(p3(4)+xm2*p3(5)))))/(q3(1)+xm2*(q3(2)+xm2*(q3(3)+xm2*(q3(4)+xm2*q3(5))))))
      elseif(ax.gt.0.47_dp)then
         erfc=exp(-x**2)*(p2(1)+ax*(p2(2)+ax*(p2(3)+ax*(p2(4)+ax*(p2(5)+ax*(p2(6)+ax*(p2(7)+ax*p2(8))))))))/(q2(1)+ax*(q2(2)+ax*(q2(3)+ax*(q2(4)+ax*(q2(5)+ax*(q2(6)+ax*(q2(7)+ax*q2(8))))))))
      else
         erfc=1.0_dp-erf(ax)
      endif
      ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
      if (x < 0.0_dp) erfc = 2.0_dp - erfc 
    end function  erfc

    ! delta functions and theta function derivatives

    pure function fermi_dirac_dydx(x) result(y)
        ! derivative of Fermi-Dirac function: 0.5/(1.0+cosh(x))
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        if (abs(x).le.36.0_dp) then
            y = 1.0_dp/(2.0_dp+exp(-x)+exp(+x))
        else
            y = 0.0_dp
        endif
    end function  fermi_dirac_dydx

    pure function marzari_vanderbilt_dydx(x) result(y)
        ! 1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        real(dp) :: arg
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp,(x-1.0_dp/sqrt(2.0_dp))**2)
        y = sqrtpm1*exp(-arg)*(2.0_dp-sqrt(2.0_dp)*x)
    end function  marzari_vanderbilt_dydx

    pure function methfessel_paxton_dydx(x,n) result(y)
        ! derivative of the corresponding Methfessel-Paxton wgauss
        implicit none
        real(dp), intent(in) :: x
        integer , intent(in) :: n
        real(dp) :: y
        real(dp) :: a, arg, hp, hd
        integer :: i, ni
        real(dp) :: sqrtpm1
        !
        sqrtpm1 = 0.564189583547756_dp ! 1/sqrt(pi)
        arg = min(200.0_dp, x**2)
        y = exp(-arg)*sqrtpm1
        if (n.eq.0) return
        hd = 0.0_dp
        hp = exp(-arg)
        ni = 0
        a  = sqrtpm1
        do i = 1, n
            hd = 2.0_dp*x*hp-2.0_dp*real(ni,dp)*hd
            ni = ni+1
            a  = -a/(real(i,dp)*4.0_dp)
            hp = 2.0_dp*x*hd-2.0_dp*real(ni,dp)*hp
            ni = ni+1
            y  = y+a*hp
        enddo
    end function  methfessel_paxton_dydx

    pure function gauss(x) result(y)
        ! derivative of Fermi-Dirac function: exp( - (x-x_o)^2 / (2*sigma^2) ) / sigma * sqrt( 2 * pi )
        ! degauss = (sqrt(2)*sigma) => exp(-(x-xo)**2/degauss**2)/(degauss*sqrt(pi))
        ! set degauss = 1 and center xo = 0
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        if (abs(x).le.7.0_dp) then
            y = exp(-x**2.0_dp)/sqrtpi
            ! in order to avoid problems for large values of x in the e
        else
            y = 0.0_dp
        endif
    end function  gauss

    pure function lorentz(x) result(y)
        ! x = ( E(j,k)-Ep(i) )/degauss 
        ! derivative of Fermi-Dirac function: 1/( ((x-xo)/degauss)**2 + 1) * 1/(pi * degauss), 
        ! notice that the 1/degauss factor is taken outside of the equation 
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: y
        !
        y = 1.0_dp/( pi * ( x**2 + 1 ) )
        !
    end function  lorentz

    ! $$$$$$\ $$\   $$\ $$$$$$$$\ $$$$$$$$\  $$$$$$\  $$$$$$$\   $$$$$$\ $$$$$$$$\ $$$$$$\  $$$$$$\  $$\   $$\
    !  \_$$  _|$$$\  $$ |\__$$  __|$$  _____|$$  __$$\ $$  __$$\ $$  __$$\\__$$  __|\_$$  _|$$  __$$\ $$$\  $$ |
    !    $$ |  $$$$\ $$ |   $$ |   $$ |      $$ /  \__|$$ |  $$ |$$ /  $$ |  $$ |     $$ |  $$ /  $$ |$$$$\ $$ |
    !    $$ |  $$ $$\$$ |   $$ |   $$$$$\    $$ |$$$$\ $$$$$$$  |$$$$$$$$ |  $$ |     $$ |  $$ |  $$ |$$ $$\$$ |
    !    $$ |  $$ \$$$$ |   $$ |   $$  __|   $$ |\_$$ |$$  __$$< $$  __$$ |  $$ |     $$ |  $$ |  $$ |$$ \$$$$ |
    !    $$ |  $$ |\$$$ |   $$ |   $$ |      $$ |  $$ |$$ |  $$ |$$ |  $$ |  $$ |     $$ |  $$ |  $$ |$$ |\$$$ |
    !  $$$$$$\ $$ | \$$ |   $$ |   $$$$$$$$\ \$$$$$$  |$$ |  $$ |$$ |  $$ |  $$ |   $$$$$$\  $$$$$$  |$$ | \$$ |
    !  \______|\__|  \__|   \__|   \________| \______/ \__|  \__|\__|  \__|  \__|   \______| \______/ \__|  \__|

    ! gaussian method

    function       get_dos_quick(Ep,E,kptw,degauss,flags) result(D)
        ! flags = fermi, mp, mv, gauss, lorentz
        implicit none
        !
        real(dp), intent(in) :: Ep(:)       ! probing energies
        real(dp), intent(in) :: E(:,:)      ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: kptw(:)     ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: degauss
        character(*), intent(in) :: flags
        real(dp), allocatable :: D(:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: i,j,k
        real(dp) :: xp
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i,j,k) SHARED(D,E,Ep,nEs,kptw)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            D(i) = 0.0_dp
            do j = 1, nbands
            do k = 1, nkpts
                xp = ( E(j,k)-Ep(i) )/degauss
                ! get contribution
                if     (index(flags,'fermi').ne.0) then
                    D(i) = D(i) + kptw(k) * fermi_dirac_dydx(x=xp)
                elseif (index(flags,'mp').ne.0) then
                    D(i) = D(i) + kptw(k) * methfessel_paxton_dydx(x=xp,n=1)
                elseif (index(flags,'mv').ne.0) then
                    D(i) = D(i) + kptw(k) * marzari_vanderbilt_dydx(x=xp)
                elseif (index(flags,'gauss').ne.0) then
                    D(i) = D(i) + kptw(k) * gauss(x=xp)
                elseif (index(flags,'lorentz').ne.0) then
                    D(i) = D(i) + kptw(k) * lorentz(x=xp)
                else
                    stop 'ERROR [get_dos_quick]: flag not recognized'
                endif
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        D = D / degauss
    end function   get_dos_quick

    function       get_pdos_quick(Ep,E,kptw,degauss,W,flags) result(pD)
        ! flags = fermi, mp, mv, gauss, lorentz
        implicit none
        !
        real(dp), intent(in) :: Ep(:)    ! probing energies
        real(dp), intent(in) :: E(:,:)   ! E(nbands,nkpts) band energies
        real(dp), intent(in) :: W(:,:,:) ! W(nprojections,nbands,nkpts) band weights
        real(dp), intent(in) :: kptw(:)  ! kptw(nkpts) normalized kpoint weights
        real(dp), intent(in) :: degauss
        character(*), intent(in) :: flags
        real(dp), allocatable :: pD(:,:)
        integer  :: nEs
        integer  :: nbands
        integer  :: nkpts
        integer  :: nprojections
        integer  :: i,j,k,l
        real(dp) :: xp,yp
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get n of kpoints
        nkpts = size(kptw)
        ! numbe of projections
        nprojections = size(W,1)
        ! allocate space for dos
        allocate(pD(nprojections,nEs))
        ! initialize
        pD = 0.0_dp
        !$OMP PARALLEL PRIVATE(i,j,k,l) SHARED(pD,E,Ep,nEs,kptw)
        !$OMP DO
        do i = 1, nEs
            do j = 1, nbands
            do k = 1, nkpts
                xp = ( E(j,k)-Ep(i) )/degauss
                ! get contribution from this band
                if     (index(flags,'fermi').ne.0) then
                    yp = fermi_dirac_dydx(x=xp)
                elseif (index(flags,'mp').ne.0) then
                    yp = methfessel_paxton_dydx(x=xp,n=1)
                elseif (index(flags,'mv').ne.0) then
                    yp = marzari_vanderbilt_dydx(x=xp)
                elseif (index(flags,'gauss').ne.0) then
                    yp = gauss(x=xp)
                elseif (index(flags,'lorentz').ne.0) then
                    yp = lorentz(x=xp)
                else
                    stop 'ERROR [get_dos_quick]: flag not recognized'
                endif
                ! multiply by kpoint weight
                yp = yp * kptw(k)
                ! multiply by band weights
                do l = 1, nprojections
                    pD(l,i) = pD(l,i) + yp * W(l,j,k)
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        pD = pD / degauss
    end function   get_pdos_quick

    ! DOS properties

    function       get_Ef(E,tet,tetw,nelecs) result(Ef)
        ! adapted from qe
        implicit none
        !
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:)
        real(dp), intent(in) :: nelecs
        real(dp) :: Ef
        integer  :: i
        integer  :: ntets
        integer  :: nbands
        integer  :: nkpts
        real(dp) :: Emin    ! absolute lowest band energy
        real(dp) :: Emax    ! absolute highest band energy
        real(dp) :: elw     ! lower bound for EF
        real(dp) :: eup     ! upper bound for EF
        real(dp) :: idos_up ! number of states for Ef=eup
        real(dp) :: idos_lw ! number of states for Ef=elw
        real(dp) :: idos_ep ! number of states for Ef=(elw+eup)/2
        real(dp) :: try_Ef  ! new Ef to try
        real(dp) :: crit    ! stopping criterion
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        nkpts  = size(E,2)
        ntets  = sizE(tet,2)
        ! get lower and upper bracketing bounds
        Emin = minval(E)
        Emax = maxval(E)
        elw = Emin
        eup = Emax
        ! initailize bracketing 
        idos_up = get_tet_idos_engine(Ep=eup,E=E,tet=tet,tetw=tetw)
        idos_lw = get_tet_idos_engine(Ep=elw,E=E,tet=tet,tetw=tetw)
        crit = 1.0d+10
        ! search for Ef by gradualling reducing bracketted region
        search : do i = 1, 1000 ! hard-coded 1000 maximum iterations
            try_Ef = (eup+elw) / 2.0_dp
            idos_ep = get_tet_idos_engine(Ep=try_Ef,E=E,tet=tet,tetw=tetw)
            if (abs(idos_ep-nelecs).lt.crit) then
                crit = abs(idos_ep-nelecs)
                Ef = try_Ef
            endif
            ! converged
            if (abs(idos_ep-nelecs).lt.tiny) then
                exit search
            elseif((idos_ep-nelecs).lt.-tiny) then
                elw = try_Ef
            else
                eup = try_Ef
            endif
        enddo search
        ! check convergence
        if (abs(idos_ep-nelecs).gt.tiny) stop 'ERROR [XX] failed to converge on Ef.'
        ! check that Ef < highest band
        if (Ef.gt.Emax) stop 'ERROR [XX]: Ef is above the highest band energy'
        !
    end function   get_Ef

    ! tetrahedra methods

    function       get_dos_tet(Ep,E,tet,tetw) result(D)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: D(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(D(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(D,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! get dos
            D(i) = get_tet_dos_engine(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_dos_tet

    function       get_idos_tet(Ep,E,tet,tetw) result(iD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), allocatable :: iD(:)
        integer :: nEs
        integer :: i
        ! get number of probing energies
        nEs = size(Ep)
        ! allocate space for dos
        allocate(iD(nEs))
        !$OMP PARALLEL PRIVATE(i) SHARED(iD,nEs,E,Ep,tet,tetw)
        !$OMP DO
        do i = 1, nEs
            ! get dos
            iD(i) = get_tet_idos_engine(Ep=Ep(i),E=E,tet=tet,tetw=tetw)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_idos_tet

    function       get_pdos_tet(Ep,E,tet,tetw,weight) result(pD)
        !
        implicit none
        !
        real(dp), intent(in) :: Ep(:) ! probing energies
        real(dp), intent(in) :: E(:,:) ! E(nbands,nkpts) band energies
        integer , intent(in) :: tet(:,:) ! tet(:,ntets) tetrahedra conenctivity
        real(dp), intent(in) :: tetw(:) ! tetw(ntets) weight of each tetrahedron
        real(dp), intent(in) :: weight(:,:,:) ! weights(nprojections,nbands,nkpts), weights (per band per kpoint) for projected dos: can be absolute square of TB or BvK eigenvectors
        real(dp), allocatable :: pD(:,:)
        integer  :: nEs
        integer  :: nbands
        integer  :: ntets
        integer  :: nprojections
        real(dp) :: wc(4) ! corner tetrahedron weights
        integer  :: i,j,k,m
        ! get number of probing energies
        nEs = size(Ep)
        ! get n of band
        nbands = size(E,1)
        ! get number of tetrahedra
        ntets = size(tet,2)
        ! get number of projections
        nprojections = size(weight,1)
        ! allocate space 
        allocate(pD(nprojections,nEs))
        ! loop over energies, bands, tetrahedra
        !$OMP PARALLEL PRIVATE(i,j,k,m,wc) SHARED(pD,E,Ep,nEs,tetw,weight)
        !$OMP DO
        do i = 1, nEs
            ! initialize D(i)
            pD(:,i) = 0.0_dp
            do j = 1, nbands
            do k = 1, ntets
                ! get tetrahedron corner weights
                wc = get_delta_wc(Ep=Ep(i),Ec=E(j,tet(:,k)))
                ! loop over projections, increment pDOS with contributions from band j in tetrahedron k
                do m = 1, nprojections
                    pD(m,i) = pD(m,i) + tetw(k) * sum(wc * weight(m,j,tet(:,k)) )
                enddo
            enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end function   get_pdos_tet

    pure function  get_theta_wc(Ep,Ec,ntets) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        integer , intent(in) :: ntets
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: etot
        real(dp) :: c1,c2,c3,c4
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        integer  :: inds(4)
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.le.e2) then
            ! Eq B6 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
            ! Eq C2 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(Ep-e1)**2/(e2-e1)/(e3-e1)/(e4-e1)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B2 from Blochl's PhysRevB.49.16223
            wc(k1)=c4*(4.0_dp-(Ep-e1)*(1.0_dp/(e2-e1)+1.0_dp/(e3-e1)+1.0_dp/(e4-e1)))+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B3 from Blochl's PhysRevB.49.16223
            wc(k2)=c4*(Ep-e1)/(e2-e1)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B4 from Blochl's PhysRevB.49.16223
            wc(k3)=c4*(Ep-e1)/(e3-e1)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B5 from Blochl's PhysRevB.49.16223
            wc(k4)=c4*(Ep-e1)/(e4-e1)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e3) then
            ! Eq B11 from Blochl's PhysRevB.49.16223
            c1=0.25_dp/ntets*(Ep-e1)**2/(e4-e1)/(e3-e1)
            ! Eq B12 from Blochl's PhysRevB.49.16223
            c2=0.25_dp/ntets*(Ep-e1)*(Ep-e2)*(e3-Ep)/(e4-e1)/(e3-e2)/(e3-e1)
            ! Eq B13 from Blochl's PhysRevB.49.16223
            c3=0.25_dp/ntets*(Ep-e2)**2*(e4-Ep)/(e4-e2)/(e3-e2)/(e4-e1)
            ! Eq C3 from Blochl's PhysRevB.49.16223
            dosEp=1.0_dp/ntets/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)*(Ep-e2)**2/(e3-e2)/(e4-e2))
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B7 from Blochl's PhysRevB.49.16223
            wc(k1)=c1+(c1+c2)*(e3-Ep)/(e3-e1)+(c1+c2+c3)*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B8 from Blochl's PhysRevB.49.16223
            wc(k2)=c1+c2+c3+(c2+c3)*(e3-Ep)/(e3-e2)+c3*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B9 from Blochl's PhysRevB.49.16223
            wc(k3)=(c1+c2)*(Ep-e1)/(e3-e1)+(c2+c3)*(Ep-e2)/(e3-e2)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B10 from Blochl's PhysRevB.49.16223
            wc(k4)=(c1+c2+c3)*(Ep-e1)/(e4-e1)+c3*(Ep-e2)/(e4-e2)+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.le.e4) then
            ! Eq B18 from Blochl's PhysRevB.49.16223
            c4=0.25_dp/ntets*(e4-Ep)**3/(e4-e1)/(e4-e2)/(e4-e3)
            ! Eq C4 from Blochl's PhysRevB.49.16223
            dosEp=3.0_dp/ntets*(e4-Ep)**2/(e4-e1)/(e4-e2)/(e4-e3)
            ! shortcut for Blochl corrections (Eq.22 PhysRevB.49.16223)
            etot=e1+e2+e3+e4
            ! Eq B14 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e1)+dosEp*(etot-4.0_dp*e1)*0.025_dp
            ! Eq B15 from Blochl's PhysRevB.49.16223
            wc(k2)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e2)+dosEp*(etot-4.0_dp*e2)*0.025_dp
            ! Eq B16 from Blochl's PhysRevB.49.16223
            wc(k3)=0.25_dp/ntets-c4*(e4-Ep)/(e4-e3)+dosEp*(etot-4.0_dp*e3)*0.025_dp
            ! Eq B17 from Blochl's PhysRevB.49.16223
            wc(k4)=0.25_dp/ntets-c4*(4.0_dp-(e4-Ep)*(1.0_dp/(e4-e1)+1.0_dp/(e4-e2)+1.0_dp/(e4-e3)))+dosEp*(etot-4.0_dp*e4)*0.025_dp
        elseif (Ep.ge.e4) then
            ! Eq B19 from Blochl's PhysRevB.49.16223
            wc(k1)=0.25_dp/ntets
            wc(k2)=0.25_dp/ntets
            wc(k3)=0.25_dp/ntets
            wc(k4)=0.25_dp/ntets
        endif
        ! 
        wc = wc * ntets
        !
    end function   get_theta_wc

    pure function  get_delta_wc(Ep,Ec) result(wc)
        ! get linear tetrahedron corner weights for delta function with Blochl corrections
        implicit none
        !
        real(dp), intent(in) :: Ep    ! probing energy
        real(dp), intent(in) :: Ec(4) ! corner energies
        real(dp) :: wc(4)
        real(dp) :: dosEp
        real(dp) :: e1,e2,e3,e4
        integer  :: k1,k2,k3,k4
        real(dp) :: o13,f12,f13,f14,f21,f23,f31,f32,f34,f24,f41,f42
        integer  :: inds(4)
        ! initialize shortcut
        o13  = 1.0_dp/3.0_dp
        ! rank corner weights in increasing order
        inds = fourrank(Ec)
        ! sort weights
        e1 = Ec(inds(1))
        e2 = Ec(inds(2))
        e3 = Ec(inds(3))
        e4 = Ec(inds(4))
        ! k1-k4 are the irreducible k-points corresponding to e1-e4
        k1 = inds(1)
        k2 = inds(2)
        k3 = inds(3)
        k4 = inds(4)
        ! calculate weights wc
        if     (Ep.le.e1) then
            ! Eq B1 from Blochl's PhysRevB.49.16223
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
        elseif (Ep.lt.e2) then
            f12 = (Ep-e2)/(e1-e2)
            f21 = 1.0_dp - f12
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp - f14
            dosEp  = 3.0_dp * f21 * f31 * f41 / (Ep-e1)
            wc(k1) = o13 * (f12 + f13 + f14)
            wc(k2) = o13 * f21
            wc(k3) = o13 * f31
            wc(k4) = o13 * f41
            wc = wc * dosEp
        elseif (Ep.lt.e3) then
            f13 = (Ep-e3)/(e1-e3)
            f31 = 1.0_dp - f13
            f14 = (Ep-e4)/(e1-e4)
            f41 = 1.0_dp-f14
            f23 = (Ep-e3)/(e2-e3)
            f32 = 1.0_dp - f23
            f24 = (Ep-e4)/(e2-e4)
            f42 = 1.0_dp - f24
            dosEp  = 3.0_dp * (f23*f31 + f32*f24)
            wc(k1) = f14 * o13 + f13*f31*f23 / dosEp
            wc(k2) = f23 * o13 + f24*f24*f32 / dosEp
            wc(k3) = f32 * o13 + f31*f31*f23 / dosEp
            wc(k4) = f41 * o13 + f42*f24*f32 / dosEp
            dosEp  = dosEp / (e4-e1)
            wc = wc * dosEp
        elseif (Ep.lt.e4) then
            f14 = (Ep-e4)/(e1-e4)
            f24 = (Ep-e4)/(e2-e4)
            f34 = (Ep-e4)/(e3-e4)
            dosEp  = 3.0_dp * f14 * f24 * f34 / (e4-Ep)
            wc(k1) = f14 * o13
            wc(k2) = f24 * o13
            wc(k3) = f34 * o13
            wc(k4) = (3.0_dp - f14 - f24 - f34 ) * o13
            wc = wc * dosEp
        elseif (Ep.gt.e4) then
            wc(k1)=0.0_dp
            wc(k2)=0.0_dp
            wc(k3)=0.0_dp
            wc(k4)=0.0_dp
       endif
       !
    end function   get_delta_wc

    pure function  get_tet_idos_engine(Ep,E,tet,tetw) result(idos)
        ! return integrated DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: idos
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer  :: i,j
        integer  :: nbands
        integer  :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        idos = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the vertices of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to 
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq A1 from Blochl's PhysRevB.49.16223
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq A2 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(Ep-e1)**3/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq A3 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)/(e3-e1)/(e4-e1)*((e2-e1)**2+3.0_dp*(e2-e1)*(Ep-e2)+3.0_dp*(Ep-e2)**2.0_dp-(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**3.0_dp)
                elseif (Ep.le.e4) then
                    ! Eq A4 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j)*(1.0_dp-(e4-Ep)**3.0_dp/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq A5 from Blochl's PhysRevB.49.16223
                    idos = idos + tetw(j) 
                endif
            enddo
        enddo
    end function   get_tet_idos_engine

    pure function  get_tet_dos_engine(Ep,E,tet,tetw) result(dosEp)
        ! returns DOS at Ep
        implicit none
        !
        real(dp), intent(in) :: Ep
        real(dp), intent(in) :: E(:,:)
        integer , intent(in) :: tet(:,:)
        real(dp), intent(in) :: tetw(:) ! weight of irreducible tetrahedron
        real(dp) :: dosEp
        real(dp) :: Ec(4)
        real(dp) :: e1,e2,e3,e4
        integer :: i,j
        integer :: nbands
        integer :: ntets
        ! get bands, kpts, tetrahedra
        nbands = size(E,1)
        ntets  = sizE(tet,2)
        !
        dosEp = 0.0_dp
        do j = 1, ntets
            do i = 1, nbands
                ! get energies at the verticies of the j-th tethedron in ascending order
                Ec(1:4) = foursort(E(i,tet(:,j)))
                ! map to 
                e1 = Ec(1)
                e2 = Ec(2)
                e3 = Ec(3)
                e4 = Ec(4)
                ! calculate sum over k of the integrated charge
                if     (Ep.le.e1) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp*1.0_dp
                    ! do nothing
                elseif (Ep.le.e2) then
                    ! Eq C2 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*3.0_dp*(Ep-e1)**2.0_dp/(e2-e1)/(e3-e1)/(e4-e1)
                elseif (Ep.le.e3) then
                    ! Eq C3 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)/(e3-e1)/(e4-e1)*(3.0_dp*(e2-e1)+6.0_dp*(Ep-e2)-3.0_dp*(e3-e1+e4-e2)/(e3-e2)/(e4-e2)*(Ep-e2)**2.0_dp)
                elseif (Ep.le.e4) then
                    ! Eq C4 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp + tetw(j)*(3.0_dp*(e4-Ep)**2.0_dp/(e4-e1)/(e4-e2)/(e4-e3))
                elseif (Ep.ge.e4) then
                    ! Eq C1 from Blochl's PhysRevB.49.16223
                    dosEp = dosEp*1.0_dp
                    ! do nothing
                endif
            enddo
        enddo
        !
    end function   get_tet_dos_engine

    ! $$\      $$\ $$$$$$\  $$$$$$\   $$$$$$\  $$$$$$$$\ $$\       $$\        $$$$$$\  $$\   $$\ $$$$$$$$\  $$$$$$\  $$\   $$\  $$$$$$\
    ! $$$\    $$$ |\_$$  _|$$  __$$\ $$  __$$\ $$  _____|$$ |      $$ |      $$  __$$\ $$$\  $$ |$$  _____|$$  __$$\ $$ |  $$ |$$  __$$\
    ! $$$$\  $$$$ |  $$ |  $$ /  \__|$$ /  \__|$$ |      $$ |      $$ |      $$ /  $$ |$$$$\ $$ |$$ |      $$ /  $$ |$$ |  $$ |$$ /  \__|
    ! $$\$$\$$ $$ |  $$ |  \$$$$$$\  $$ |      $$$$$\    $$ |      $$ |      $$$$$$$$ |$$ $$\$$ |$$$$$\    $$ |  $$ |$$ |  $$ |\$$$$$$\
    ! $$ \$$$  $$ |  $$ |   \____$$\ $$ |      $$  __|   $$ |      $$ |      $$  __$$ |$$ \$$$$ |$$  __|   $$ |  $$ |$$ |  $$ | \____$$\
    ! $$ |\$  /$$ |  $$ |  $$\   $$ |$$ |  $$\ $$ |      $$ |      $$ |      $$ |  $$ |$$ |\$$$ |$$ |      $$ |  $$ |$$ |  $$ |$$\   $$ |
    ! $$ | \_/ $$ |$$$$$$\ \$$$$$$  |\$$$$$$  |$$$$$$$$\ $$$$$$$$\ $$$$$$$$\ $$ |  $$ |$$ | \$$ |$$$$$$$$\  $$$$$$  |\$$$$$$  |\$$$$$$  |
    ! \__|     \__|\______| \______/  \______/ \________|\________|\________|\__|  \__|\__|  \__|\________| \______/  \______/  \______/

    pure subroutine uc2ws(K,M,tiny)

        implicit none

        integer :: j
        real(dp), allocatable, intent(inout) :: K(:,:)
        real(dp), intent(in) :: M(3,3)
        real(dp), intent(in) :: tiny
        real(dp) :: G(3,26) , G2(26)

        ! generate mesh
        G(1,:) = real([-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1])
        G(2,:) = real([-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1])
        G(3,:) = real([-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1])
        G = matmul(M,G)
        G2 = sum(G**2,1)

        ! perform calc
        ! do j = 1, size(K,2)
        do concurrent (j = 1:size(K,2))
            call uc2ws_engine(K(:,j),G,G2,tiny)
        enddo
    end subroutine  uc2ws

    pure subroutine uc2ws_engine(K,G,G2,tiny)

        implicit none

        real(dp), intent(inout) :: K(3)
        real(dp), intent(in) :: G(3,26), G2(26), tiny
        real(dp)  :: P
        logical :: go
        integer :: i

        go = .true.
        do while (go)
            go = .false.
            do i = 1, 26
                P = 2*(K(1)*G(1,i) + K(2)*G(2,i) + K(3)*G(3,i)) - G2(i)
                if (P .gt. 1D-8) then
                    K = K - G(:,i)
                    go = .true.
                endif
            enddo
        enddo
    end subroutine  uc2ws_engine

end module
%}
            
            % mex interfaces
            i=0;
            i=i+1; f{i} = 'uc2ws_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : uc2ws

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: K(:,:)
    real*8 :: tiny, G(3,26) , G2(26) , M(3,3)

    ! [K] = uc2ws(K,M,tiny)
    if(nrhs .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Three inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','One output produced.')
    if(mxGetM(prhs(1)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 1 requires three-dimensional column vectors.')
    if(mxGetM(prhs(2)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 2 requires a three-dimensional transformation matrix.')
    if(mxGetN(prhs(2)) .ne. 3) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 2 requires a three-dimensional transformation matrix.')
    if(mxGetM(prhs(3)) .ne. 1 .or. &
     & mxGetN(prhs(3)) .ne. 1) call mexErrMsgIdAndTxt ('MATLAB:uc2ws_mex','Input 3 requires a scalar.')

    ! inputs 
    ! K
    allocate( K(3,mxGetN(prhs(1))) )
    call mxCopyPtrToReal8(mxGetPr(prhs(1)),K,mxGetM(prhs(1))*mxGetN(prhs(1)))
    ! M(3,3)
    call mxCopyPtrToReal8(mxGetPr(prhs(2)),M,mxGetM(prhs(2))*mxGetN(prhs(2)))
    ! tiny
    call mxCopyPtrToReal8(mxGetPr(prhs(3)),tiny,mxGetM(prhs(3)))

    ! execute code
    call uc2ws(K,M,tiny)

    ! output results
    plhs(1) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(K,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_dos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_dos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: D(:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Four inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_dos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')

    ! inputs 
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))

    ! execute code
    D = get_dos_tet(Ep,E,nint(tet),tetw)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_idos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_idos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: D(:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Four inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_idos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')

    ! inputs 
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))

    ! execute code
    D = get_idos_tet(Ep,E,nint(tet),tetw)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(1)),mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            i=i+1; f{i} = 'get_pdos_tet_mex'; fid=fopen([f{i},'.f90'],'w'); fprintf(fid,'%s',verbatim_()); fclose(fid); 
%{
#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    use am_mex_lib , only : get_pdos_tet

    implicit none

    mwPointer plhs(*), prhs(*)
    integer nlhs, nrhs
    mwPointer mxGetPr
    mwPointer mxCreateDoubleMatrix
    mwPointer mxGetM, mxGetN

    real*8, allocatable :: Ep(:)
    real*8, allocatable :: E(:,:)
    real*8, allocatable :: tetw(:)
    real*8, allocatable :: tet(:,:)
    real*8, allocatable :: weight(:,:,:)
    real*8, allocatable :: D(:,:)

    integer :: i

    ! [D] = get_dos_tet(Ep,E,tet,tetw)
    if(nrhs .ne. 5) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Five inputs required.')
    if(nlhs .gt. 1) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','One output produced.')
    if(mxGetM(prhs(3)) .ne. 4) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Connectivity list requires four numbers in each column.')
    if(mxGetN(prhs(3)) .ne. mxGetN(prhs(4))) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Mismatched number of tetrahedra in weights and connectivity list.')
    if(mxGetN(prhs(2)) .ne. mxGetN(prhs(5))) call mexErrMsgIdAndTxt ('MATLAB:get_pdos_tet_mex','Mismatched number of kpoints in projections weights and connectivity list.')

    ! inputs 
    i=0
    ! 1) Ep
    i=i+1
    allocate(   Ep( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),Ep, mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 2) E(nbands,nkpts)
    i=i+1
    allocate(   E( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),E,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 3) tet(1:4,ntets)
    i=i+1
    allocate( tet( mxGetM(prhs(i)), mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tet,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 4) tetw(ntets)
    i=i+1
    allocate( tetw( mxGetM(prhs(i))*mxGetN(prhs(i)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),tetw,mxGetM(prhs(i))*mxGetN(prhs(i)))
    ! 5) weight(nprojections,nbands,nkpts)
    i=i+1 
    allocate( weight( mxGetM(prhs(i)) , mxGetM(prhs(2)) , mxGetN(prhs(2)) ) )
    call mxCopyPtrToReal8(mxGetPr(prhs(i)),weight,mxGetM(prhs(i))*mxGetM(prhs(2))*mxGetN(prhs(2)))

    ! execute code
    D = get_pdos_tet(Ep,E,nint(tet),tetw,weight)

    ! outputs
    i=0
    ! output results
    i=i+1
    plhs(i) = mxCreateDoubleMatrix(mxGetM(prhs(5)),mxGetM(prhs(1))*mxGetN(prhs(1)),0)
    call mxCopyReal8ToPtr(D,mxGetPr(plhs(1)),mxGetM(prhs(5))*mxGetM(prhs(1))*mxGetN(prhs(1)))
end subroutine mexFunction
%}
            
            % compile everything
            fprintf('Building mex library:\n')
                if system([am_lib.FC,' ',am_lib.FFLAGS,' ',am_lib.DEBUG,' ',am_lib.LIBS,' -c ',flib,'.f90'])
                    error(' ... %s (failed)\n',flib)
                else
                    fprintf(' ... %s (succeeded)\n',flib);
                end
            
            fprintf('Building mex interfaces:\n')
            for i = 1:numel(f)
                if system([am_lib.FC,' ',am_lib.FFLAGS,' ',am_lib.DEBUG,' ',am_lib.LIBS,' ',flib,'.o ',f{i},'.f90 -o ',f{i},am_lib.EXT])
                    warning(' ... %s (failed)\n',f{i})
                else
                    fprintf(' ... %s (succeeded)\n',f{i});
                end
            end
        
        end
        
        function [K] = uc2ws(K,M)
            % uc2ws uses M real (reciprocal) lattice vectors to reduces K(1:3,:) vectors 
            % in cartesian (reciprocal) coordinates to the definiging Wigner-Seitz cell.

            if  and( am_lib.usemex , which('uc2ws_mex') )
                % use mex function if available
                K = uc2ws_mex(K,M,am_lib.eps);
            else
                % generate mesh
                G = M * [-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1;
                         -1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1;
                         -1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1];
                G2=sum(G.^2,1);
                % call engine
                for j = 1:size(K,2); K(:,j) = uc2ws_engine(K(:,j),G,G2,am_lib.eps); end
            end

            function [K] = uc2ws_engine(K,G,G2,tiny)
                go = true;
                while go; go=false;
                    for i = 1:26
                        P = 2*K.'*G(:,i)-G2(i);
                        if P > + tiny
                            K = K - G(:,i); go=true;
                        end
                    end
                end
            end
        end

        function [D] = get_dos_tet(Ep,E,tet,tetw)
            %
            if  and( am_lib.usemex , which('get_dos_tet_mex') )
                % use mex function if available
                D = get_dos_tet_mex(Ep,E,tet,tetw);
            else
                nEps = numel(Ep); D = zeros(nEps,1);
                for m = 1:nEps
                    D(m) = get_dos_tet_engine(Ep(m),E,tet,tetw);
                end
                D = reshape(D,size(Ep));
            end

            function [dosEp] = get_dos_tet_engine(Ep,E,tet,tetw)
                nbands = size(E,1);
                ntets = size(tet,2);
                dosEp = 0;
                for j = 1:ntets
                    for i = 1:nbands
                        % get energies at the verticies of the j-th tethedron in ascending order
                        Ec(1:4) = sort(E(i,tet(:,j)));
                        % calculate sum over k of the integrated charge
                        if     Ep<=Ec(1)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        elseif Ep<=Ec(2)
                            % Eq C2 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)*3.0*(Ep-Ec(1)).^2/(Ec(2)-Ec(1))/(Ec(3)-Ec(1))/(Ec(4)-Ec(1));
                        elseif Ep<=Ec(3)
                            % Eq C3 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)/(Ec(3)-Ec(1))/(Ec(4)-Ec(1))*(3.0*(Ec(2)-Ec(1))+6.0*(Ep-Ec(2))-3.0*(Ec(3)-Ec(1)+Ec(4)-Ec(2))/(Ec(3)-Ec(2))/(Ec(4)-Ec(2))*(Ep-Ec(2)).^2);
                        elseif Ep<=Ec(4)
                            % Eq C4 from Blochl's PhysRevB.49.16223
                            dosEp = dosEp + tetw(j)*(3.0*(Ec(4)-Ep).^2/(Ec(4)-Ec(1))/(Ec(4)-Ec(2))/(Ec(4)-Ec(3)));
                        elseif Ep>=Ec(4)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        end
                    end
                end
            end
        end

        function [D] = get_idos_tet(Ep,E,tet,tetw)
            %
            if  and( am_lib.usemex , which('get_idos_tet_mex') )
                % use mex function if available
                D = get_idos_tet_mex(Ep,E,tet,tetw);
            else
                nEps = numel(Ep); D = zeros(nEps,1);
                for m = 1:nEps
                    D(m) = get_dos_tet_engine(Ep(m),E,tet,tetw);
                end
                D = reshape(D,size(Ep));
            end

            function [idos] = get_dos_tet_engine(Ep,E,tet,tetw)
                nbands = size(E,1);
                ntets = size(tet,2);
                idos = 0;
                for j = 1:ntets
                    for i = 1:nbands
                        % get energies at the verticies of the j-th tethedron in ascending order
                        Ec(1:4) = sort(E(i,tet(:,j)));
                        % calculate sum over k of the integrated charge
                        if     Ep<=Ec(1)
                            % Eq C1 from Blochl's PhysRevB.49.16223
                            % do nothing
                        elseif Ep<=Ec(2)
                            % Eq A2 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)*(Ep-Ec(1)).^3/(Ec(2)-Ec(1))/(Ec(3)-Ec(1))/(Ec(4)-Ec(1));
                        elseif Ep<=Ec(3)
                            % Eq A3 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)/(Ec(3)-Ec(1))/(Ec(4)-Ec(1))*((Ec(2)-Ec(1)).^2+3.0*(Ec(2)-Ec(1))*(Ep-Ec(2))+3.0*(Ep-Ec(2)).^2.0-(Ec(3)-Ec(1)+Ec(4)-Ec(2))/(Ec(3)-Ec(2))/(Ec(4)-Ec(2))*(Ep-Ec(2)).^3.0);
                        elseif Ep<=Ec(4)
                            % Eq A4 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j)*(1.0-(Ec(4)-Ep).^3.0/(Ec(4)-Ec(1))/(Ec(4)-Ec(2))/(Ec(4)-Ec(3)));
                        elseif Ep>=Ec(4)
                            % Eq A5 from Blochl's PhysRevB.49.16223
                            idos = idos + tetw(j);
                        end
                    end
                end
            end
        end

        function [D] = get_pdos_tet(Ep,E,tet,tetw)
            %
            if  and( am_lib.usemex , which('get_pdos_tet_mex') )
                % use mex function if available
                D = get_pdos_tet_mex(Ep,E,tet,tetw);
            else
                error('matlab version not yet implemented. mex is required for pdos');
            end
        end
        
    end
    
    % color map
    
    methods (Static)
        
        % aesthetic

        function [cmap] =  get_colormap(palette,n)
            % color palettes
            switch palette
            case 'spectral' % sns.color_palette("Spectral", 10))
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
            case 'RdGy' % sns.color_palette("RdGy", 10)
            cmap = [ ... 
             0.66920416904430757, 0.08489042841920666, 0.16401384522517523;
             0.81153403660830326, 0.32110727271612954, 0.27581700390460445;
             0.92226067360709696, 0.56747406545807333, 0.44867361117811766;
             0.97970011655022116, 0.78408305785235233, 0.68489044554093303;
             0.99646289909587182, 0.93633218372569371, 0.90096117468441239;
             0.94517493598601399, 0.94517493598601399, 0.94517493598601399;
             0.82583622722064742, 0.82583622722064742, 0.82583622722064742;
             0.67058825492858887, 0.67058825492858887, 0.67058825492858887;
             0.48481355812035354, 0.48481355812035354, 0.48481355812035354;
             0.28235295195789900, 0.28235295195789900, 0.28235295195789900];
            case 'OrRd' % sns.color_palette("OrRd", 10)
            cmap = [ ... 
             0.99717031927669753, 0.92618224200080423, 0.82362169448067157;
             0.99434063855339494, 0.87504806588677797, 0.71132643012439500;
             0.99215686321258545, 0.81522492450826312, 0.60281432586557726;
             0.99215686321258545, 0.74140716650906735, 0.52604385754641370;
             0.98965013518052947, 0.61802385975332819, 0.40985776314548417;
             0.96984237011741192, 0.49634757789911010, 0.32496733069419859;
             0.92950404111076801, 0.37896194370353925, 0.26911189196740881;
             0.85863899343154015, 0.22246828552554634, 0.14805075201918097;
             0.76452135198256554, 0.08341407314235088, 0.05387158794145958;
             0.64518263491929750, 0.00000000000000000, 0.00000000000000000];
            case 'GnBu' % sns.color_palette("GnBu", 10)
            cmap = [ ... 
             0.90354479621438422, 0.96276816830915568, 0.88175317890503824;
             0.84367551873711977, 0.93903883485233086, 0.82059209066278793;
             0.77674741815118231, 0.91252595677095305, 0.76221454704509062;
             0.67044984663234042, 0.87118801229140341, 0.71497118192560527;
             0.54602076925483400, 0.82405229806900027, 0.74740485934650192;
             0.41868512595401092, 0.76462900287964763, 0.78985007019603959;
             0.29457901909070855, 0.68936564711963433, 0.82066898416070377;
             0.19123414667213665, 0.57420994464088881, 0.75866206744137932;
             0.09219531267881393, 0.47040370585871677, 0.70579009827445538;
             0.03137255087494850, 0.36416763593168822, 0.62755865349489104];
            case 'YlGn' % sns.color_palette("YlGn", 10)
            cmap = [ ... 
             0.97736255421357998, 0.99151095783009247, 0.77353326993830063;
             0.91649366126340981, 0.96738177818410542, 0.68725876527674057;
             0.82256056420943313, 0.92890427182702462, 0.62565169474657845;
             0.69264131013084862, 0.87280277574763576, 0.56364477802725399;
             0.54557478696692230, 0.80901192987666415, 0.50422146378778943;
             0.39277201715637655, 0.73826991179410151, 0.43489427531466762;
             0.24521339342874640, 0.65799309996997613, 0.35630912044469049;
             0.15663207261001361, 0.54283739749123061, 0.27953865212552687;
             0.06082276345468031, 0.45650136143553488, 0.23653979979309381;
             0.00000000000000000, 0.36962707428371205, 0.20039984916939454];
            case 'YlGnBu' % In [11]: sns.color_palette("YlGnBu", 10)
            cmap = [ ... 
             0.94906574698055490, 0.98019223493688246, 0.73779317210702333;
             0.86337563290315522, 0.94648212545058308, 0.69933104444952576;
             0.73388697750428145, 0.89564014462863695, 0.71040370814940512;
             0.52129181202720187, 0.81296425566953767, 0.73107268038917994;
             0.34262207781567294, 0.74626683557734774, 0.75589390151640945;
             0.20396771430969238, 0.66137641387827251, 0.76296810332466569;
             0.11534025493790122, 0.55215688698432031, 0.74519032660652607;
             0.13010381214758929, 0.40156863822656519, 0.67432527892729821;
             0.13988466630963720, 0.27690888619890397, 0.61514804293127623;
             0.11343329991487897, 0.17880815288015439, 0.51487891335113378];
            case 'GnBu_d' % sns.color_palette("GnBu_d")
            cmap = [ ...
             0.21697808798621680, 0.32733564601225013, 0.36941176807179171;
             0.23442778952760630, 0.45820839330261826, 0.54352941859002213;
             0.25140587751382310, 0.58554403931486831, 0.71294118666181383;
             0.32480841754308700, 0.68493145540648814, 0.78994746862673293;
             0.45066770474895150, 0.75099834881576832, 0.77038576275694604;
             0.58002308326608990, 0.81890043370863974, 0.75028067616855398];
            case 'cubehelix_greenblue' % sns.cubehelix_palette(8, start=.5, rot=-.75)
            cmap = [ ...
             0.84232988177938480, 0.87374044279641840, 0.75249540307310370;
             0.68251876357072430, 0.81069190728320800, 0.63524701801182060;
             0.51093657786460060, 0.73671986965753190, 0.57304087941263320;
             0.37208664465749840, 0.63786334195260290, 0.55503689058379240;
             0.28846276635927040, 0.51638144597481420, 0.54342177164221150;
             0.24670155725826660, 0.37340824813905654, 0.49725690696587516;
             0.22179626547231540, 0.23841378594571613, 0.39797674055755683;
             0.17250549177124480, 0.11951843162770594, 0.24320155229883056];
            case 'cubehelix_purple' % sns.cubehelix_palette(8)
            cmap = [ ... 
             0.93126922233253720, 0.82019217960821180, 0.79714809746635920;
             0.88228981687371890, 0.69582086670574200, 0.70654571194854310;
             0.81353802547006760, 0.57050551823578220, 0.63928085946815500;
             0.71958007083491190, 0.45537982893127477, 0.58610629958109260;
             0.60469068026344690, 0.35739308184976665, 0.53374078536924060;
             0.46496993672552040, 0.26868986121314253, 0.46365277636406470;
             0.32101947432593470, 0.19303051265196464, 0.37078816777247920;
             0.17508656489522050, 0.11840023306916837, 0.24215989137836502];
            case 'red2blue' % sns.color_palette("RdBu_r", 7)
            cmap = [ ... 
             0.16339870177063293, 0.44498270983789490, 0.6975009791991290;
             0.42068437209316328, 0.67643216077019186, 0.8186851319144753;
             0.76147636946509856, 0.86851211856393251, 0.9245674785445717;
             0.96908881383783674, 0.96647443490869855, 0.9649365649503820;
             0.98246828247519102, 0.80069205340217142, 0.7061130509657018;
             0.89457901435739851, 0.50380624217145586, 0.3997693394913390;
             0.72848905885920801, 0.15501730406985564, 0.1973856272650700];
            case 'colorbrewer' % http://colorbrewer2.org/#type=diverging&scheme=Spectral&n=11
            cmap = [ ...
             0.22656250000000000, 0.00390625000000000, 0.2578125000000000;
             0.83203125000000000, 0.24218750000000000, 0.3085937500000000;
             0.95312500000000000, 0.42578125000000000, 0.2617187500000000;
             0.98828125000000000, 0.67968750000000000, 0.3789062500000000;
             0.99218750000000000, 0.87500000000000000, 0.5429687500000000;
             0.99609375000000000, 0.99609375000000000, 0.7460937500000000;
             0.89843750000000000, 0.95703125000000000, 0.5937500000000000;
             0.66796875000000000, 0.86328125000000000, 0.6406250000000000;
             0.39843750000000000, 0.75781250000000000, 0.6445312500000000;
             0.19531250000000000, 0.53125000000000000, 0.7382812500000000;
             0.36718750000000000, 0.30859375000000000, 0.6328125000000000];
            case 'magma'
            cmap = [ ...
             1.46159096e-03,   4.66127766e-04,   1.38655200e-02;
             2.25764007e-03,   1.29495431e-03,   1.83311461e-02;
             3.27943222e-03,   2.30452991e-03,   2.37083291e-02;
             4.51230222e-03,   3.49037666e-03,   2.99647059e-02;
             5.94976987e-03,   4.84285000e-03,   3.71296695e-02;
             7.58798550e-03,   6.35613622e-03,   4.49730774e-02;
             9.42604390e-03,   8.02185006e-03,   5.28443561e-02;
             1.14654337e-02,   9.82831486e-03,   6.07496380e-02;
             1.37075706e-02,   1.17705913e-02,   6.86665843e-02;
             1.61557566e-02,   1.38404966e-02,   7.66026660e-02;
             1.88153670e-02,   1.60262753e-02,   8.45844897e-02;
             2.16919340e-02,   1.83201254e-02,   9.26101050e-02;
             2.47917814e-02,   2.07147875e-02,   1.00675555e-01;
             2.81228154e-02,   2.32009284e-02,   1.08786954e-01;
             3.16955304e-02,   2.57651161e-02,   1.16964722e-01;
             3.55204468e-02,   2.83974570e-02,   1.25209396e-01;
             3.96084872e-02,   3.10895652e-02,   1.33515085e-01;
             4.38295350e-02,   3.38299885e-02,   1.41886249e-01;
             4.80616391e-02,   3.66066101e-02,   1.50326989e-01;
             5.23204388e-02,   3.94066020e-02,   1.58841025e-01;
             5.66148978e-02,   4.21598925e-02,   1.67445592e-01;
             6.09493930e-02,   4.47944924e-02,   1.76128834e-01;
             6.53301801e-02,   4.73177796e-02,   1.84891506e-01;
             6.97637296e-02,   4.97264666e-02,   1.93735088e-01;
             7.42565152e-02,   5.20167766e-02,   2.02660374e-01;
             7.88150034e-02,   5.41844801e-02,   2.11667355e-01;
             8.34456313e-02,   5.62249365e-02,   2.20755099e-01;
             8.81547730e-02,   5.81331465e-02,   2.29921611e-01;
             9.29486914e-02,   5.99038167e-02,   2.39163669e-01;
             9.78334770e-02,   6.15314414e-02,   2.48476662e-01;
             1.02814972e-01,   6.30104053e-02,   2.57854400e-01;
             1.07898679e-01,   6.43351102e-02,   2.67288933e-01;
             1.13094451e-01,   6.54920358e-02,   2.76783978e-01;
             1.18405035e-01,   6.64791593e-02,   2.86320656e-01;
             1.23832651e-01,   6.72946449e-02,   2.95879431e-01;
             1.29380192e-01,   6.79349264e-02,   3.05442931e-01;
             1.35053322e-01,   6.83912798e-02,   3.14999890e-01;
             1.40857952e-01,   6.86540710e-02,   3.24537640e-01;
             1.46785234e-01,   6.87382323e-02,   3.34011109e-01;
             1.52839217e-01,   6.86368599e-02,   3.43404450e-01;
             1.59017511e-01,   6.83540225e-02,   3.52688028e-01;
             1.65308131e-01,   6.79108689e-02,   3.61816426e-01;
             1.71713033e-01,   6.73053260e-02,   3.70770827e-01;
             1.78211730e-01,   6.65758073e-02,   3.79497161e-01;
             1.84800877e-01,   6.57324381e-02,   3.87972507e-01;
             1.91459745e-01,   6.48183312e-02,   3.96151969e-01;
             1.98176877e-01,   6.38624166e-02,   4.04008953e-01;
             2.04934882e-01,   6.29066192e-02,   4.11514273e-01;
             2.11718061e-01,   6.19917876e-02,   4.18646741e-01;
             2.18511590e-01,   6.11584918e-02,   4.25391816e-01;
             2.25302032e-01,   6.04451843e-02,   4.31741767e-01;
             2.32076515e-01,   5.98886855e-02,   4.37694665e-01;
             2.38825991e-01,   5.95170384e-02,   4.43255999e-01;
             2.45543175e-01,   5.93524384e-02,   4.48435938e-01;
             2.52220252e-01,   5.94147119e-02,   4.53247729e-01;
             2.58857304e-01,   5.97055998e-02,   4.57709924e-01;
             2.65446744e-01,   6.02368754e-02,   4.61840297e-01;
             2.71994089e-01,   6.09935552e-02,   4.65660375e-01;
             2.78493300e-01,   6.19778136e-02,   4.69190328e-01;
             2.84951097e-01,   6.31676261e-02,   4.72450879e-01;
             2.91365817e-01,   6.45534486e-02,   4.75462193e-01;
             2.97740413e-01,   6.61170432e-02,   4.78243482e-01;
             3.04080941e-01,   6.78353452e-02,   4.80811572e-01;
             3.10382027e-01,   6.97024767e-02,   4.83186340e-01;
             3.16654235e-01,   7.16895272e-02,   4.85380429e-01;
             3.22899126e-01,   7.37819504e-02,   4.87408399e-01;
             3.29114038e-01,   7.59715081e-02,   4.89286796e-01;
             3.35307503e-01,   7.82361045e-02,   4.91024144e-01;
             3.41481725e-01,   8.05635079e-02,   4.92631321e-01;
             3.47635742e-01,   8.29463512e-02,   4.94120923e-01;
             3.53773161e-01,   8.53726329e-02,   4.95501096e-01;
             3.59897941e-01,   8.78311772e-02,   4.96778331e-01;
             3.66011928e-01,   9.03143031e-02,   4.97959963e-01;
             3.72116205e-01,   9.28159917e-02,   4.99053326e-01;
             3.78210547e-01,   9.53322947e-02,   5.00066568e-01;
             3.84299445e-01,   9.78549106e-02,   5.01001964e-01;
             3.90384361e-01,   1.00379466e-01,   5.01864236e-01;
             3.96466670e-01,   1.02902194e-01,   5.02657590e-01;
             4.02547663e-01,   1.05419865e-01,   5.03385761e-01;
             4.08628505e-01,   1.07929771e-01,   5.04052118e-01;
             4.14708664e-01,   1.10431177e-01,   5.04661843e-01;
             4.20791157e-01,   1.12920210e-01,   5.05214935e-01;
             4.26876965e-01,   1.15395258e-01,   5.05713602e-01;
             4.32967001e-01,   1.17854987e-01,   5.06159754e-01;
             4.39062114e-01,   1.20298314e-01,   5.06555026e-01;
             4.45163096e-01,   1.22724371e-01,   5.06900806e-01;
             4.51270678e-01,   1.25132484e-01,   5.07198258e-01;
             4.57385535e-01,   1.27522145e-01,   5.07448336e-01;
             4.63508291e-01,   1.29892998e-01,   5.07651812e-01;
             4.69639514e-01,   1.32244819e-01,   5.07809282e-01;
             4.75779723e-01,   1.34577500e-01,   5.07921193e-01;
             4.81928997e-01,   1.36891390e-01,   5.07988509e-01;
             4.88088169e-01,   1.39186217e-01,   5.08010737e-01;
             4.94257673e-01,   1.41462106e-01,   5.07987836e-01;
             5.00437834e-01,   1.43719323e-01,   5.07919772e-01;
             5.06628929e-01,   1.45958202e-01,   5.07806420e-01;
             5.12831195e-01,   1.48179144e-01,   5.07647570e-01;
             5.19044825e-01,   1.50382611e-01,   5.07442938e-01;
             5.25269968e-01,   1.52569121e-01,   5.07192172e-01;
             5.31506735e-01,   1.54739247e-01,   5.06894860e-01;
             5.37755194e-01,   1.56893613e-01,   5.06550538e-01;
             5.44015371e-01,   1.59032895e-01,   5.06158696e-01;
             5.50287252e-01,   1.61157816e-01,   5.05718782e-01;
             5.56570783e-01,   1.63269149e-01,   5.05230210e-01;
             5.62865867e-01,   1.65367714e-01,   5.04692365e-01;
             5.69172368e-01,   1.67454379e-01,   5.04104606e-01;
             5.75490107e-01,   1.69530062e-01,   5.03466273e-01;
             5.81818864e-01,   1.71595728e-01,   5.02776690e-01;
             5.88158375e-01,   1.73652392e-01,   5.02035167e-01;
             5.94508337e-01,   1.75701122e-01,   5.01241011e-01;
             6.00868399e-01,   1.77743036e-01,   5.00393522e-01;
             6.07238169e-01,   1.79779309e-01,   4.99491999e-01;
             6.13617209e-01,   1.81811170e-01,   4.98535746e-01;
             6.20005032e-01,   1.83839907e-01,   4.97524075e-01;
             6.26401108e-01,   1.85866869e-01,   4.96456304e-01;
             6.32804854e-01,   1.87893468e-01,   4.95331769e-01;
             6.39215638e-01,   1.89921182e-01,   4.94149821e-01;
             6.45632778e-01,   1.91951556e-01,   4.92909832e-01;
             6.52055535e-01,   1.93986210e-01,   4.91611196e-01;
             6.58483116e-01,   1.96026835e-01,   4.90253338e-01;
             6.64914668e-01,   1.98075202e-01,   4.88835712e-01;
             6.71349279e-01,   2.00133166e-01,   4.87357807e-01;
             6.77785975e-01,   2.02202663e-01,   4.85819154e-01;
             6.84223712e-01,   2.04285721e-01,   4.84219325e-01;
             6.90661380e-01,   2.06384461e-01,   4.82557941e-01;
             6.97097796e-01,   2.08501100e-01,   4.80834678e-01;
             7.03531700e-01,   2.10637956e-01,   4.79049270e-01;
             7.09961888e-01,   2.12797337e-01,   4.77201121e-01;
             7.16387038e-01,   2.14981693e-01,   4.75289780e-01;
             7.22805451e-01,   2.17193831e-01,   4.73315708e-01;
             7.29215521e-01,   2.19436516e-01,   4.71278924e-01;
             7.35615545e-01,   2.21712634e-01,   4.69179541e-01;
             7.42003713e-01,   2.24025196e-01,   4.67017774e-01;
             7.48378107e-01,   2.26377345e-01,   4.64793954e-01;
             7.54736692e-01,   2.28772352e-01,   4.62508534e-01;
             7.61077312e-01,   2.31213625e-01,   4.60162106e-01;
             7.67397681e-01,   2.33704708e-01,   4.57755411e-01;
             7.73695380e-01,   2.36249283e-01,   4.55289354e-01;
             7.79967847e-01,   2.38851170e-01,   4.52765022e-01;
             7.86212372e-01,   2.41514325e-01,   4.50183695e-01;
             7.92426972e-01,   2.44242250e-01,   4.47543155e-01;
             7.98607760e-01,   2.47039798e-01,   4.44848441e-01;
             8.04751511e-01,   2.49911350e-01,   4.42101615e-01;
             8.10854841e-01,   2.52861399e-01,   4.39304963e-01;
             8.16914186e-01,   2.55894550e-01,   4.36461074e-01;
             8.22925797e-01,   2.59015505e-01,   4.33572874e-01;
             8.28885740e-01,   2.62229049e-01,   4.30643647e-01;
             8.34790818e-01,   2.65539703e-01,   4.27671352e-01;
             8.40635680e-01,   2.68952874e-01,   4.24665620e-01;
             8.46415804e-01,   2.72473491e-01,   4.21631064e-01;
             8.52126490e-01,   2.76106469e-01,   4.18572767e-01;
             8.57762870e-01,   2.79856666e-01,   4.15496319e-01;
             8.63320397e-01,   2.83729003e-01,   4.12402889e-01;
             8.68793368e-01,   2.87728205e-01,   4.09303002e-01;
             8.74176342e-01,   2.91858679e-01,   4.06205397e-01;
             8.79463944e-01,   2.96124596e-01,   4.03118034e-01;
             8.84650824e-01,   3.00530090e-01,   4.00047060e-01;
             8.89731418e-01,   3.05078817e-01,   3.97001559e-01;
             8.94700194e-01,   3.09773445e-01,   3.93994634e-01;
             8.99551884e-01,   3.14616425e-01,   3.91036674e-01;
             9.04281297e-01,   3.19609981e-01,   3.88136889e-01;
             9.08883524e-01,   3.24755126e-01,   3.85308008e-01;
             9.13354091e-01,   3.30051947e-01,   3.82563414e-01;
             9.17688852e-01,   3.35500068e-01,   3.79915138e-01;
             9.21884187e-01,   3.41098112e-01,   3.77375977e-01;
             9.25937102e-01,   3.46843685e-01,   3.74959077e-01;
             9.29845090e-01,   3.52733817e-01,   3.72676513e-01;
             9.33606454e-01,   3.58764377e-01,   3.70540883e-01;
             9.37220874e-01,   3.64929312e-01,   3.68566525e-01;
             9.40687443e-01,   3.71224168e-01,   3.66761699e-01;
             9.44006448e-01,   3.77642889e-01,   3.65136328e-01;
             9.47179528e-01,   3.84177874e-01,   3.63701130e-01;
             9.50210150e-01,   3.90819546e-01,   3.62467694e-01;
             9.53099077e-01,   3.97562894e-01,   3.61438431e-01;
             9.55849237e-01,   4.04400213e-01,   3.60619076e-01;
             9.58464079e-01,   4.11323666e-01,   3.60014232e-01;
             9.60949221e-01,   4.18323245e-01,   3.59629789e-01;
             9.63310281e-01,   4.25389724e-01,   3.59469020e-01;
             9.65549351e-01,   4.32518707e-01,   3.59529151e-01;
             9.67671128e-01,   4.39702976e-01,   3.59810172e-01;
             9.69680441e-01,   4.46935635e-01,   3.60311120e-01;
             9.71582181e-01,   4.54210170e-01,   3.61030156e-01;
             9.73381238e-01,   4.61520484e-01,   3.61964652e-01;
             9.75082439e-01,   4.68860936e-01,   3.63111292e-01;
             9.76690494e-01,   4.76226350e-01,   3.64466162e-01;
             9.78209957e-01,   4.83612031e-01,   3.66024854e-01;
             9.79645181e-01,   4.91013764e-01,   3.67782559e-01;
             9.81000291e-01,   4.98427800e-01,   3.69734157e-01;
             9.82279159e-01,   5.05850848e-01,   3.71874301e-01;
             9.83485387e-01,   5.13280054e-01,   3.74197501e-01;
             9.84622298e-01,   5.20712972e-01,   3.76698186e-01;
             9.85692925e-01,   5.28147545e-01,   3.79370774e-01;
             9.86700017e-01,   5.35582070e-01,   3.82209724e-01;
             9.87646038e-01,   5.43015173e-01,   3.85209578e-01;
             9.88533173e-01,   5.50445778e-01,   3.88365009e-01;
             9.89363341e-01,   5.57873075e-01,   3.91670846e-01;
             9.90138201e-01,   5.65296495e-01,   3.95122099e-01;
             9.90871208e-01,   5.72706259e-01,   3.98713971e-01;
             9.91558165e-01,   5.80106828e-01,   4.02441058e-01;
             9.92195728e-01,   5.87501706e-01,   4.06298792e-01;
             9.92784669e-01,   5.94891088e-01,   4.10282976e-01;
             9.93325561e-01,   6.02275297e-01,   4.14389658e-01;
             9.93834412e-01,   6.09643540e-01,   4.18613221e-01;
             9.94308514e-01,   6.16998953e-01,   4.22949672e-01;
             9.94737698e-01,   6.24349657e-01,   4.27396771e-01;
             9.95121854e-01,   6.31696376e-01,   4.31951492e-01;
             9.95480469e-01,   6.39026596e-01,   4.36607159e-01;
             9.95809924e-01,   6.46343897e-01,   4.41360951e-01;
             9.96095703e-01,   6.53658756e-01,   4.46213021e-01;
             9.96341406e-01,   6.60969379e-01,   4.51160201e-01;
             9.96579803e-01,   6.68255621e-01,   4.56191814e-01;
             9.96774784e-01,   6.75541484e-01,   4.61314158e-01;
             9.96925427e-01,   6.82827953e-01,   4.66525689e-01;
             9.97077185e-01,   6.90087897e-01,   4.71811461e-01;
             9.97186253e-01,   6.97348991e-01,   4.77181727e-01;
             9.97253982e-01,   7.04610791e-01,   4.82634651e-01;
             9.97325180e-01,   7.11847714e-01,   4.88154375e-01;
             9.97350983e-01,   7.19089119e-01,   4.93754665e-01;
             9.97350583e-01,   7.26324415e-01,   4.99427972e-01;
             9.97341259e-01,   7.33544671e-01,   5.05166839e-01;
             9.97284689e-01,   7.40771893e-01,   5.10983331e-01;
             9.97228367e-01,   7.47980563e-01,   5.16859378e-01;
             9.97138480e-01,   7.55189852e-01,   5.22805996e-01;
             9.97019342e-01,   7.62397883e-01,   5.28820775e-01;
             9.96898254e-01,   7.69590975e-01,   5.34892341e-01;
             9.96726862e-01,   7.76794860e-01,   5.41038571e-01;
             9.96570645e-01,   7.83976508e-01,   5.47232992e-01;
             9.96369065e-01,   7.91167346e-01,   5.53498939e-01;
             9.96162309e-01,   7.98347709e-01,   5.59819643e-01;
             9.95932448e-01,   8.05527126e-01,   5.66201824e-01;
             9.95680107e-01,   8.12705773e-01,   5.72644795e-01;
             9.95423973e-01,   8.19875302e-01,   5.79140130e-01;
             9.95131288e-01,   8.27051773e-01,   5.85701463e-01;
             9.94851089e-01,   8.34212826e-01,   5.92307093e-01;
             9.94523666e-01,   8.41386618e-01,   5.98982818e-01;
             9.94221900e-01,   8.48540474e-01,   6.05695903e-01;
             9.93865767e-01,   8.55711038e-01,   6.12481798e-01;
             9.93545285e-01,   8.62858846e-01,   6.19299300e-01;
             9.93169558e-01,   8.70024467e-01,   6.26189463e-01;
             9.92830963e-01,   8.77168404e-01,   6.33109148e-01;
             9.92439881e-01,   8.84329694e-01,   6.40099465e-01;
             9.92089454e-01,   8.91469549e-01,   6.47116021e-01;
             9.91687744e-01,   8.98627050e-01,   6.54201544e-01;
             9.91331929e-01,   9.05762748e-01,   6.61308839e-01;
             9.90929685e-01,   9.12915010e-01,   6.68481201e-01;
             9.90569914e-01,   9.20048699e-01,   6.75674592e-01;
             9.90174637e-01,   9.27195612e-01,   6.82925602e-01;
             9.89814839e-01,   9.34328540e-01,   6.90198194e-01;
             9.89433736e-01,   9.41470354e-01,   6.97518628e-01;
             9.89077438e-01,   9.48604077e-01,   7.04862519e-01;
             9.88717064e-01,   9.55741520e-01,   7.12242232e-01;
             9.88367028e-01,   9.62878026e-01,   7.19648627e-01;
             9.88032885e-01,   9.70012413e-01,   7.27076773e-01;
             9.87690702e-01,   9.77154231e-01,   7.34536205e-01;
             9.87386827e-01,   9.84287561e-01,   7.42001547e-01;
             9.87052509e-01,   9.91437853e-01,   7.49504188e-01];
            end

            % interpolating function
            map_ = @(n,cmap) interp1([0:(size(cmap,1)-1)]./(size(cmap,1)-1),cmap,linspace(0,1,n));

            % interpolate and return
            cmap = map_(n,cmap);
        end
        
    end
end








