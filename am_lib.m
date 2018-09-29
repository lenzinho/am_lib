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
        
        % constants
        r_0       = 2.81794032E-6; % [nm]       classical electron radius
        N_A       = 6.022141E23;   % [mol]      Avogadro's number
    end
    
    % unit conversion
    
    methods (Static)
        
        function [varargout] = kxkz2angle(kx,kz,hv,flag)
            % [alpha_i,alpha_f] = kxkz2angle(kx,kz,hv,'in/exit')
            % [   w   ,  th2  ] = kxkz2angle(kx,kz,hv,'w2th')
            %
            % test 
            % kx = rand(1)*10;
            % kz = rand(1)*10;
            % [w,th2] = kxkz2angle(kx,kz,hv,'w2th')
            % [kxp,kzp]     = angle2kxkz(w,th2,hv)
            % kx-kxp 
            % kz-kzp
            
            import am_lib.*
            
            if isempty(hv)
                Cu = am_atom.define('Cu'); 
                hv = Cu.get_emission_line('kalpha1'); 
            end
            lambda = get_photon_wavelength(hv);
            
            th2_ = @(kx,kz,lambda) 2.*atan2d(...
                 real(sqrt(    (kx.^2 + kz.^2).*lambda.^2 )),...
                 real(sqrt(4 - (kx.^2 + kz.^2).*lambda.^2 )));
            w_   = @(kx,kz,lambda)  atan2d(...
                 real( +(kz.*kx.^2.*lambda + kz.^3.*lambda + kx.*sqrt(kx.^2+kz.^2).*sqrt(4-kx.^2.*lambda.^2-kz.^2.*lambda.^2))./(kx.^2+kz.^2) ), ...
                 real( -(kx.*kz.^2.*lambda + kx.^3.*lambda - kz.*sqrt(kx.^2+kz.^2).*sqrt(4-kx.^2.*lambda.^2-kz.^2.*lambda.^2))./(kx.^2+kz.^2) ));
             
            switch flag
                case 'in/exit'
                    % recover th2 and w values to determine which points are accessible
                    alpha_i =   w_(kx,kz,lambda); 
                    alpha_f = th2_(kx,kz,lambda)-alpha_i;
                    varargout{1} = alpha_i;
                    varargout{2} = alpha_f;
                case 'w2th'
                    varargout{1} =   w_(kx,kz,lambda);
                    varargout{2} = th2_(kx,kz,lambda);
                case 'test'
                    kx = rand(1)*10; kz = rand(1)*10; hv = 1;
                    [w,th2] = kxkz2angle(kx,kz,hv,'w2th');
                    [kxp,kzp] = angle2kxkz(w,th2,hv);
                    fprintf('%f \t %f\n',kx-kxp, kz-kzp)
            end
        end
        
        function [kx,kz]     = angle2kxkz(w,th2,hv,flag)
            % [kx,kz] = angle2kxkz(w,th2,hv)
            
            import am_lib.*
            
            if isempty(hv)
                Cu = am_atom.define('Cu'); 
                hv = Cu.get_emission_line('kalpha1'); 
            end
            lambda = get_photon_wavelength(hv);
            
            switch flag
                case 'in/exit'
                    % in this case
                    % w = incident angle
                    % 2th = exit angle
                    th2 = th2 + w;
                    [kx,kz]  = angle2kxkz(w,th2,hv,'w2th');
                case 'w2th'        
                    % convert to reciprocal coordinates
                    kx_ = @(w,th2,lambda)  2/(lambda).*sind(th2/2).*sind(th2/2-w);
                    kz_ = @(w,th2,lambda)  2/(lambda).*sind(th2/2).*cosd(th2/2-w);
                    kx  = kx_(w,th2,lambda); kz = kz_(w,th2,lambda);
            end
        end

        function [kz]        = get_kz(th,hv)
            
            import am_lib.*
            
            kz = sind(th)/get_photon_wavelength(hv);
            
        end

        function [th]        = get_th(kz,hv)
            
            import am_lib.*
            
            th = asind(get_photon_wavelength(hv)*kz);
            
        end
          
        function [hv]        = get_photon_energy(lambda)
            %
            % E = get_photon_energy(lambda)
            % 
            % E      [eV]           photon energy
            % lambda [nm]           photon wavelength
            %
            % Conversion factor = 
            %       Plank's constant * speed of light / nm / eV
            %
            
            hv = 1239.842 ./ lambda;
            
        end
        
        function [lambda]    = get_photon_wavelength(hv)
            %
            % E = get_photon_energy(lambda)
            % 
            % E      [eV]           photon energy
            % lambda [nm]          photon wavelength
            %
            % Conversion factor = 
            %       Plank's constant * speed of light / nm / eV
            %
            
            lambda = 1239.842 ./ hv;
            
        end
        
        function [x]         = rebase_(index,base)
            % [dhms] = rebase_(seconds,[7,24,60,60])
            % use this to convert from seconds to weeks,days,hours minutes,seconds 
            n = numel(base)+1; x = zeros(1,n-1); b = cumprod([base,1],'reverse');
            for i = 2:n
                x(i-1)  = floor(index./b(i));
                index = index - x(i-1).*b(i);
            end
        end
        
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

        function     test()

            z=[];
            syms z; n = 50;

            % [-1,+1]
            F = {@(n) am_lib.chebyshevTr_(n,'edge'), ...
                 @(n) am_lib.chebyshevUr_(n,'edge'), ...
                 @(n) am_lib.legendrer_(n), ...
                 @(n) am_lib.clenshawcurtisr_(n)};

            f_ = @(x) x.^2; w_ = @(x) 1;
            exact = double(int(f_(z).*w_(z),z,-1,+1));
            dfdx_ = matlabFunction(diff(f_(z)));

            label = {'Chebyshev T','Chebyshev U','Legendre','Clenshaw-Curtis'};
            for i = 1:numel(F)
                [x,D,w]= F{i}(n); D=D(:,:,1);
                L(1) = exact-sum(f_(x).*w);
                L(2) = max(abs(dfdx_(x)-D*f_(x)));
                criteria = L<am_lib.tiny;
                test_(all(criteria), sprintf('%s',label{i}), ...
                                     sprintf('failed (%g)',L(~criteria)));
            end

            % [0,2*pi]
            F = {@(n) am_lib.fourierr_(n)};

            f_ = @(x) sin(12*pi*x); w_ = @(x) 1;
            exact = double(int(f_(z).*w_(z),z,-1,+1));
            dfdx_ = matlabFunction(diff(f_(z)));

            label = {'Fourier'};
            for i = 1:numel(F)
                [x,D,w]= F{i}(n); D=D(:,:,1);
                L(1) = exact-sum(f_(x).*w);
                L(2) = max(abs(dfdx_(x)-D*f_(x)));
                criteria = L<am_lib.tiny;
                test_(all(criteria), sprintf('%s',label{i}), ...
                                     sprintf('failed (%g)',L(~criteria)));
            end

            % [0,+Inf]
            F = {@(n) am_lib.laguerrer_(n,'edge')};

            f_ = @(x) exp(-3*x.^2); w_ = @(x) exp(-x);
            exact = double(int(f_(z).*w_(z),z,0,Inf));
            dfdx_ = matlabFunction(diff(f_(z)));

            label = {'Laguerre'};
            for i = 1:numel(F)
                [x,D,w]= F{i}(n); D=D(:,:,1);
                L(1) = exact-sum(f_(x).*w);
                L(2) = max(abs(dfdx_(x)-D*f_(x)));
                criteria = L<am_lib.tiny;
                test_(all(criteria), sprintf('%s',label{i}), ...
                                     sprintf('failed (%g)',L(~criteria)));
            end

            % [-Inf,+Inf]
            F = {@(n) am_lib.hermiter_(n)};

            f_ = @(x) exp(-3*x.^2); w_ = @(x) exp(-x.^2);
            exact = double(int(f_(z).*w_(z),z,-Inf,Inf));
            dfdx_ = matlabFunction(diff(f_(z)));

            label = {'Hermite'};
            for i = 1:numel(F)
                [x,D,w]= F{i}(n); D=D(:,:,1);
                L(1) = exact-sum(f_(x).*w);
                L(2) = max(abs(dfdx_(x)-D*f_(x)));
                criteria = L<am_lib.tiny;
                test_(all(criteria), sprintf('%s',label{i}), ...
                                     sprintf('failed (%g)',L(~criteria)));
            end

            % finite difference (with pbc)
            n = 10*n; x = [0:n-1]'/n; f_ = @(x) sin(2*pi*x);
            label = {'Central Difference [-1, 0,+1]','Central Difference [-1,+1]   ','Forward Difference [ 0,+1,+2]','Backward Difference [-2,-1, 0]'};
            q = {[-1,0,1],[-1,1],[0,1,2],[-2,-1,0]}; 
            for i = 1:numel(q)
                c = am_lib.get_differentiation_weights(q{i},1);
                f    = f_(x);
                dfdx = squeeze(am_lib.circshift_(f(:),-q{i})) *c(:) ./(x(2)-x(1));
                dfdx_= matlabFunction(diff(f_(z))); % exact
                
                L(1) = max(abs(dfdx_(x)-dfdx(:)));
                criteria = L<am_lib.tiny;
                test_(all(criteria), sprintf('%s',label{i}), ...
                                     sprintf('failed (%g)',L(~criteria)));
            end

            plot(x,dfdx_(x),'.',x,dfdx,'-')
            
            function test_(logical,test_name,fail_msg)
                if logical
                    fprintf('      %s: pass\n',test_name);
                else
                    fprintf('      %s: %s\n',test_name,fail_msg);
                end 
            end
        end

        
        % numerical precision

        function [L] = iseven_(A,tol)
            % set default numerical tolernece
            if nargin < 2; tol = am_lib.tiny; end
            % evaluate
            L = abs(mod(A+tol,2)-tol)<tol;
        end
        
        function [L] = isodd_(A,tol)
            % set default numerical tolernece
            if nargin < 2; tol = am_lib.tiny; end
            % evaluate
            L = abs(mod(A+tol,2)-tol-1)<tol;
        end
        
        function [C] = mod_(A,tol)
            % set default numerical tolernece
            if nargin < 2; tol = am_lib.tiny; end
            % evaluate
            C = mod(A+tol,1)-tol;
        end

        function [C] = rnd_(A,tol)
            if nargin < 2; tol = am_lib.tiny; end
            C = round(A,-log10(tol));
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
                    % well defined values for 1/3, 2/3, 1/4, 3/4, 1/6, 5/6, 1/8, 3/8, 5/8, 7/8
                    for wdv = [1/3, 2/3, 1/4, 3/4, 1/6, 5/6, 1/8, 3/8, 5/8, 7/8]; if go
                        if eq_(abs(x(i)),wdv,tol)
                            x(i) = wdv * sign(x(i)); 
                            go=false; break; 
                        end
                    end;end
                end
            else
                x = wdv_(real(x)) + wdv_(imag(x))*1i;
            end
        end

        function [L] = isdiag_(A)
            L = am_lib.all_(am_lib.eq_(diag(diag(A)),A));
        end
        
        function [N] = null_(A,tol)
            if nargin < 2; tol = am_lib.eps; end
            [V,D]=eig(A,'vector'); N = V(:,am_lib.eq_(D,0,tol));
        end
        
        function [C,IA,IC] = uniquetol_(X,tol)
            if nargin < 2; tol = am_lib.eps; end
            [~,IA,IC] = unique(am_lib.rnd_(X,tol),'stable'); C = X(:,IA);
        end

        
        % data-type conversion
        
        function [B] = str2num_(A)
            if iscell(A)
                nAs = numel(A);
                B = zeros(1,nAs);
                for i = 1:nAs
                    B(i) = str2num(A{i});
                end
            else isstring(A)
                B = str2num(A);
            end
        end
        
        
        % symbolic
        
        function M = subs_(M,var,val)
            % converts list of symbolic variables to double
            %     var = sort(symvar([ip.v_fc{:}]));
            %     val = [ip.fc{:}];
            if iscell(M)
                for j = 1:numel(M)
                    M{j} = am_lib.subs_(M{j},var,val);
                end
            else
                for i = 1:numel(val)
                    M = subs(M,var(i),val(i));
                end
                M = double(M);
            end
        end
        
        
        % vectorization
        
        function [C] = field2array_(A,field)
            n = numel(A); C = zeros(1,n);
            for i = 1:n
                C(i) = A(i).(field);
            end
        end

        function [C] = field2cell_(A,field)
            n = numel(A); C = cell(1,n);
            for i = 1:n
                C{i} = A(i).(field){:};
            end
        end

        function [C] = flatten_(A)
            C = A(:);
        end
        
        function X   = outer_(A,B,op) % Define outer operations (generalizes outer product to any operation).
            X = op(repmat(A,1,length(B)), repmat(B,length(A),1)); 
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
        
        function [C] = circshift_(A,B)
            for i = 1:size(B,2); C{i} = circshift(A,B(:,i)); end
            C = cat(ndims(C{1})+1,C{:});
        end
        
        function [c] = if_(L,a,b)
            if L
                c = a;
            else
                c = b;
            end
        end
        
        function [A] = aug_(A,n)
            % squeeze in 1's in vector A at positions n
            ex_ = any([1:(numel(n)+numel(A))]==n(:),1);
            A(~ex_) = A; A( ex_) = 1;
        end
        
        function [A] = ext_(A,n)
            % remove positions n of vector A
            ex_ = any(1:numel(A)==n(:),1);
            A = A(~ex_);
        end
        
        function [L] = all_(x)
            L = all(x(:));
        end
        
        function [L] = any_(x)
            L = any(x(:));
        end

        function [d] = det_(A)
            n = size(A,3); 
            d = zeros(1,n);
            for i = 1:n
                d(i) = det(A(:,:,i));
            end
        end
        
        function [d] = diag_(A)
            % get diagonal elements:  diag(i,:,:,...) = A(i,i,:,:,...)
            n = size(A); m = min(n(1),n(2));
            d = zeros(m,prod(n(3:end)));
            for i = 1:m; d(i,:) = A(i,i,:); end
            d = reshape(d,[m,n(3:end)]);
        end
        
        function [T] = trace_(A)
            % sum over diagonal elements: tr(1,:,:,...) = sum_i A(i,i,:,:,...) 
            T = sum(am_lib.diag_(A),1);
        end
        
        function [C] = mtimes_(varargin)
            switch nargin
                case 1
                    C = varargin{1};
                otherwise
                    C = varargin{1};                    
                    for i = 2:nargin; C = C * varargin{i}; end
            end
        end
        
        function [A] = wrap_(A,dim,n)
            % wraps matrix
            % e.g. am_lib.wrap_(magic(3),[1,2],[2,2])
            nds = numel(dim);
            if nds>1
                if nds ~= numel(n); error('incorrect input'); end
                for i = 1:nds; A = am_lib.wrap_(A,dim(i),n(i)); end
            else
                if n==0; return; end
                ind = [1:ndims(A)]; fwd = circshift(ind,dim-1); rev(fwd) = ind;
                A = permute(A,fwd); s = size(A); A = reshape(A,s(1),[]);
                    if n > 0
                        A = cat(1,A,A(1:n,:));
                    else
                        A = cat(1,A(end+n+1:end,:),A);
                    end
                A = reshape(A,[s(1)+abs(n),s(2:end)]); A = permute(A,rev);
            end
        end
        
        function [A] = reflect_(A,dim,n)
            % reflects matrix
            % e.g. am_lib.reflect_(reshape([1:16],4,4).',2,-3)
            nds = numel(dim);
            if nds>1
                if nds ~= numel(n); error('incorrect input'); end
                for i = 1:nds; A = am_lib.wrap_(A,dim(i),n(i)); end
            else
                if n==0; return; end
                ind = [1:ndims(A)]; fwd = circshift(ind,dim-1); rev(fwd) = ind;
                A = permute(A,fwd); s = size(A); A = reshape(A,s(1),[]);
                    if n > 0
                        A = cat(1,A,A(end-1:-1:end-abs(n),:));
                    else
                        A = cat(1,A(abs(n)+1:-1:2,:),A);
                    end
                A = reshape(A,[s(1)+abs(n),s(2:end)]); A = permute(A,rev);
            end
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

        function [C] = kronsum_(A,B)
            % % only works for 2d matrices A and B
            % [m,n]=size(A); [p,q]=size(B);
            % C = repelem(A,p,q)+repmat(B,m,n);
            
            % works for all cases:
            % define high dimensional kronecker self-sum: for i=[1:size(A,3)]; C(:,:,i)=kron(A(:,:,i),B(:,:,i)); end
            C = reshape( permute(A,[4 1 5 2 3]) ...
                       + permute(B,[1 4 2 5 3]), size(A,1)*size(B,1),[],size(A,3));
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
        
        function [C] = ssum_(A)
           % sums over the third dimension of a symmetry tensor because sum(A,3) does not currently exist in matlab
           [l,m,n] = size(A); C = sym(zeros(l,m)); for i = 1:n; C = C + A(:,:,i); end
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

        function T   = transpose_(d,p)
            % T = transpose_(d,p)
            %       p = permutation 
            %       d = dimensions
            % For 3-dimensional transpose super operator:
            %       transpose_([3,3],[2,1])
            %
            % test with:
            %     clear;clc;
            %     p = [1,3,4,2];
            %     d = [4,5,6,2];
            %     T = get_transpose_superop_(p,d);
            %     X = rand(d);
            %     Y = T * X(:);
            %     reshape(Y,d(p)) - permute(X,p)

            % get total dimensions
            pd = prod(d); 
            % create indicies
            A = reshape([1:pd],d); B = A; B = permute(B,p);
            % create permutation super operator
            T = sparse( A(:), B(:), ones(1,pd), pd, pd );
        end

        function Dq  = hist_(varargin) % hist_(x,D,xq), hist_(x,y,D,xq,yq)
            switch nargin
                case 3 % hist_(x,D,xq)
                    x = varargin{1}(:); 
                    D = varargin{2}(:); 
                    xq= varargin{3}(:); [m,n] = size(varargin{3});
                    x = x.'; xq = xq.';
                case 5 % hist_(x,y,D,xq,yq)
                    x = varargin{1}(:); 
                    y = varargin{2}(:); 
                    D = varargin{3}(:); 
                    xq= varargin{4}(:); 
                    yq= varargin{5}(:); [m,n] = size(varargin{4});
                    x = [x,y].'; xq = [xq,yq].';
                otherwise
                    error('invalid input');
            end
            % find index of closest xq point for each x
            i = knnsearch(xq.',x.'); 
            % integrate
            Dq = accumarray(i(:),D,[],@sum).';
            % normalize
            nvalues = sum(i(:)==[1:max(i(:))],1); Dq = Dq./nvalues; 
            % ????
            Dq(numel(Dq)+1:m*n) = NaN; 
            % reshape
            Dq = reshape(Dq,[m,n]);
            % laplacian interpolation on nan values with Dirichlet b.c.
            Dq = regionfill(Dq,isnan(Dq));
        end
        
        
        % functions operating on matrix column
        
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
            
            % deals with real first and then imaginary
            C =    sortrowsc(rnd_(real(A     ),tol).',[1:size(A,1)]).'; 
            C = C( sortrowsc(rnd_(imag(A(:,C)),tol).',[1:size(A,1)]) ); 
        end                 

        function [C] = findc_(A)
            % define function to find the first nonzero value in each row of matrix A
            % returns 0 for rows containing all zeros
            C = sum(cumsum(A~=0,1)==0,1)+1;
            C = C .* ~all(A==0,1);
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
        
        function [C] = normc_(A)
            % get length of each column vector
            C = sqrt(sum(abs(A).^2,1));
        end
        
        function [C] = norm_(A)
            % normalize columns
            C = A./am_lib.normc_(A);
        end
        
        function [C] = intercalate_(A,B)
            % intercalates columns of A and B into C
            C = B(:,[1;1]*(1:size(B,2)));
            C(:,1:2:end) = A;
        end
        

        % matching
        
        function [i,j] = max2_(A)
            [~,ind] = max(A(:)); [i,j]=ind2sub(size(A),ind);
        end
        
        function [i,j] = min2_(A)
            [~,ind] = max(A(:)); 
            [i,j]=ind2sub(size(A),ind);
        end
        
        function [C] = maxabs_(A)
            C = max(abs(A(:)));
        end
        
        function [varargout] = minmax_(A)
            switch nargout
                case 1
                    varargout{1} = [min(A(:)),max(A(:))];
                case 2
                    varargout{1} = min(A(:));
                    varargout{2} = max(A(:));
                case 4
                    [varargout{1:2}] = min(A(:));
                    [varargout{3:4}] = max(A(:));
                otherwise
                    error('invalid number of outputs');
            end
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
        
        function [C] = setdiffi_(A,B)
            % return at most n values
            % integer set diff
            C = false(1,max(numel(A),numel(B)));
            try
            C(A(:)) = true; C(B(:)) = false; C = C(A(:));
            catch
                asdf
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
            % reshape
            [n,m]=size(c_id); if n < m; c_id = c_id.'; end
                % count repetitions/occurances
                occurances = sum(c_id==c_id.',2);
                % sort, lift degeneracies if present using original labeling
                fwd = rankc_([occurances(:).';c_id(:).']);
                % relabel and return indicies to original order
                c_id(fwd) = cumsum([0;diff(c_id(fwd))~=0])+1; 
            c_id = reshape(c_id,[n,m]);
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

        function [A] = merge_(A)
            
            import am_lib.*
            
            % convert to logical
            A = ~eq_(A,0);

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
        end
        
        function [x] = sort4_(x) % ultra-fast sorting of four numbers 
            % ultra-fast sorting of four numbers
            if x(1)>y(2); x([1,2]) = x([2,1]); end
            if x(3)>y(4); x([3,4]) = x([4,3]); end
            if x(1)>y(3); x([1,3]) = x([3,1]); end
            if x(2)>y(4); x([2,4]) = x([4,2]); end
            if x(2)>y(3); x([2,3]) = x([3,2]); end
        end
        
        function [a] = rank4_(x) % ultra-fast ranking of four numbers 
            % ultra-fast ranking of four numbers
            a = zeros(1,8);
            if x(1) <= x(2);  a(1) = 1; a(2) = 2; 
            else;             a(1) = 2; a(2) = 1; end
            if x(3) <= x(4);  a(3) = 3; a(4) = 4; 
            else;             a(3) = 4; a(4) = 3; end
            if x(a(1)) <= x(a(3)); a(6) = a(1); a(7) = a(3); 
            else;                  a(6) = a(3); a(7) = a(1); end
            if x(a(2)) >= x(a(4)); a(5) = a(2); a(8) = a(4); 
            else;                  a(5) = a(4); a(8) = a(2); end
            if x(a(7)) < x(a(8));  a=[a(6),a(7),a(8),a(5)]; 
            else;                  a=[a(6),a(8),a(7),a(5)]; end
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
        
        function t   = extract_token_below(str,token,numeric)
            import am_lib.*
            l = find(contains(str,token),1);
            if ~isempty(l)
                b = strsplit(str{l},' ');
                w = find(strcmp(b,token))-1;
                b = strsplit(str{l+1},' ');
                t = b(w);
                if nargin>2 && numeric
                    t = sscanf(t,'%f');
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
        
        function [Z] = zeros_(varargin)
            n = numel([varargin{:}]);
            if n == 1
                Z = zeros(varargin{1},1);
            else
                Z = zeros(varargin{:});
            end
        end
        
        function [y] = rand_(varargin) % uniformly-distributed random number 
            % uniformly distributed random number between 0 and 1
            % Garcia p 347
            seed = round(clock*flipud(cumprod([1 60 60 24 31 12].')));
            y = am_lib.zeros_(varargin{:}); n = prod([varargin{:}]);
            a = 7^5; c = 0; M = 2^31-1;
            % generate the seed and throw the first five values away
            y(1) = mod(seed,M); for i = 1:5; y(1) = mod(a*y(1)+c,M); end
            if n > 2; for i = 2:n
                y(i) = mod(a*y(i-1)+c,M);
            end; end
            y = y./M;
        end
        
        function [y] = rande_(varargin) % exponentially distributed random number with decay 1 
            % exponentially distributed random number with decay 
            % Garcia p 349
            y = am_lib.rand_([varargin{:}]); y = -log(1-y);
        end
        
        function [y] = randn_(varargin) % normally-distributed random number with std 1, mean 0 
            % Box-Muller normally distributed number with std 1, mean 0
            % Garcia p 349
            % generate an even number of uniformly distributed random variables between 0 and 1
            n = prod([varargin{:}]); y = am_lib.rand_(n+mod(n,2));
            % apply box muller
            a = sqrt(-2*log(y(1:2:end))); b = 2*pi*y(2:2:end);
            y(1:2:end) = a.*sin(b); y(2:2:end) = a.*cos(b);
            % reshape output variable
            m = numel([varargin{:}]); 
            if m==1; y = reshape(y(1:n),[varargin{1},1]);
            else;    y = reshape(y(1:n),[varargin{:}]); end
        end
        
        function [p] = arg_(cmplx) % argument (phase) of complex number 
            p = atan2(imag(cmplx),real(cmplx));
        end 

        function [y] = sinc_(x) % sinc function 
            y = sin(x)./x; y(x==0) = 1;
        end
  
        function [y] = expsum_(x,N) % sum_{n=0}^{N} exp(i n x)
            if N < 1E8 % finite
                y = (1 - exp(1i*N*x))./(1-exp(1i*x)); y(x==0)=1;
            else
                y = 1./(1-exp(1i*x));
            end
        end
        
        function [y] = heaviside_(x) % heaviside function
            % A. V. Podolskiy and P. Vogl, Phys. Rev. B 69, 233101 (2004). Eq 15
            y = logical(x>=0);
        end
        
        function [y] = fermi_dirac_(x) % fermi-dirac funcrtion 
            y = 1.0./(1.0 + exp(-x)); 
        end
        
        function [y] = methfessel_paxton_(x,n) % theta function PRB 40, 3616 (1989).
            maxarg = 200;
            % Methfessel-Paxton
            y = gauss_freq(x .* sqrt(2.0));
            if n==0; return; end
            hd = 0; arg = min(maxarg, x.^2); hp = exp(-arg);
            ni = 0; a = 1.0./sqrt(pi); 
            for i=1:n
                hd = 2.0.*x.*hp-2.0.*ni.*hd; 
                ni = ni+1; a = -a./(i.*4.0); y = y-a.*hd;
                hp = 2.0.*x.*hd-2.0.*ni.*hp; ni = ni+1;
            end
            function g = gauss_freq(x)
                %     gauss_freq(x) = (1+erf(x./sqrt(2)))./2 = erfc(-x./sqrt(2))./2
                g = 0.5 .* am_lib.erfc_(-x.*0.7071067811865475);
                %
            end
        end
        
        function [y] = marzari_vanderbilt_(x) % Marzari-Vanderbilt theta-function PRL 82, 3296 (1999)
            % theta function: PRL 82, 3296 (1999)
            % 1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
            %, xp
            maxarg = 200;
            % Cold smearing
             xp = x - 1.0./sqrt(2.0); arg = min(maxarg, xp.^2);
             y = 0.5.*am_lib.erf_(xp) + 1.0./sqrt(2.0.*pi).*exp(-arg) + 0.5d0;
        end
        
        function [y] = erf_(x) % error function 
            y = arrayfun(@erf_engine_,x);
            function y = erf_engine_(x)
                p1 = [2.426679552305318E2, 2.197926161829415E1, 6.996383488619136,  -3.560984370181538E-2];
                q1 = [2.150588758698612E2, 9.116490540451490E1, 1.508279763040779E1, 1.000000000000000];
                %
                if abs(x)>6.0
                    y = sign(x);
                else
                    if abs(x)<0.47
                        x2 = x.^2;
                        y = x.*(p1(1)+x2.*(p1(2)+x2.*(p1(3)+x2.*p1(4))))/(q1(1)+x2.*(q1(2)+x2.*(q1(3)+x2.*q1(4))));
                    else
                        y = 1.0 - am_lib.erfc_(x);
                    end
                end
            end
        end

        function [y] = erfc_(x) % complementary error function 
            
            y = arrayfun(@erfc_engine_,x);
            
            function y = erfc_engine_(x)
                p2 = [ 3.004592610201616E2,  4.519189537118719E2,  3.393208167343437E2,  1.529892850469404E2,  4.316222722205674E1,  7.211758250883094,    5.641955174789740E-1,-1.368648573827167E-7];
                q2 = [ 3.004592609569833E2,  7.909509253278980E2,  9.313540948506096E2,  6.389802644656312E2,  2.775854447439876E2,  7.700015293522947E1,  1.278272731962942E1,  1.000000000000000];
                p3 = [-2.996107077035422E-3,-4.947309106232507E-2, -2.269565935396869E-1,-2.786613086096478E-1,  -2.231924597341847E-2];
                q3 = [ 1.062092305284679E-2, 1.913089261078298E-1,  1.051675107067932,    1.987332018171353,     1.000000000000000];
                pim1 = 0.56418958354775629; % sqrt(1./pi)
                ax = abs(x);
                if ax > 26.0
                    y = 0.0;
                elseif ax >4.0
                    x2=x.^2; xm2=(1.0/ax).^2;
                    y=(1.0/ax)*exp(-x2)*(pim1+xm2*(p3(1)+xm2*(p3(2)+xm2*(p3(3)+xm2*(p3(4)+xm2*p3(5)))))/(q3(1)+xm2*(q3(2)+xm2*(q3(3)+xm2*(q3(4)+xm2*q3(5))))));
                elseif ax>0.47
                    y = exp(-x.^2).*(p2(1)+ax.*(p2(2)+ax.*(p2(3)+ax.*(p2(4)+ax.*(p2(5)+ax.*(p2(6)+ax.*(p2(7)+ax.*p2(8))))))))./(q2(1)+ax.*(q2(2)+ax.*(q2(3)+ax.*(q2(4)+ax.*(q2(5)+ax.*(q2(6)+ax.*(q2(7)+ax.*q2(8))))))));
                else
                    y = 1.0-am_lib.erf_(ax);
                end
                if x < 0.0;  y = 2.0 - y; end
            end
        end
        
        function [y] = lorentz_(x) % lorentz peak function
            y = 1./(pi.*(x.^2+1)); 
        end

        function [y] = gauss_(x) % gaussian normal distribution peak function
            y = exp(-abs(x).^2)./sqrt(pi);
        end

        function [y] = delta_(x) % kronecker delta
            y = logical(abs(x)<am_lib.tiny); 
        end

        function [y] = fermi_dirac_dydx_(x) % derivative of fermi dirac function
            % derivative of Fermi-Dirac function: 0.5./(1.0+cosh(x))
            y = zeros(size(x)); ex_ = abs(x)<36.0; y(ex_) = 1.0./(2.0+exp(-x(ex_))+exp(+x(ex_)));
        end
        
        function [y] = marzari_vanderbilt_dydx_(x) % derivative of Marzari-Vanderbilt theta-function 
            % 1./sqrt(pi).*exp(-(x-1./sqrt(2)).^2).*(2-sqrt(2).*x)
            %
            sqrtpm1 = 0.564189583547756; % 1./sqrt(pi)
            arg = min(200,(x-1.0./sqrt(2.0)).^2);
            y = sqrtpm1.*exp(-arg).*(2.0-sqrt(2.0).*x);
        end
        
        function [y] = methfessel_paxton_dydx_(x,n) % derivative of Methfessel-Paxton 
            sqrtpm1 = 0.564189583547756; % 1./sqrt(pi)
            arg = min(200, x.^2);
            y = exp(-arg).*sqrtpm1;
            if n==0 return; end
            hd = 0; hp = exp(-arg); ni = 0; a = sqrtpm1;
            for i = 1:n
                hd = 2.0.*x.*hp-2.0.*ni.*hd;
                ni = ni+1; a = -a./(i.*4.0);
                hp = 2.0.*x.*hd-2.0.*ni.*hp;
                ni = ni+1; y = y+a.*hp;
            end
        end
        
        function [y] = pvoigt_(x,f) % pseudovoigt peak function
            y = sqrt(pi) .* (1-f) .* am_lib.gauss_(x./sqrt(2*log(2))) + f .* am_lib.lorentz_(x) .* pi;
        end
        
        function [w] = tukeyw_(n,r) % tukey window
            if nargin == 1; r = 0.5; end
            if r <= 0 ; w = ones(1,n); elseif r >= 1; w = am_lib.hannw_(n); else
                t = linspace(0,1,n)'; per = r/2; tl = floor(per*(n-1))+1; th = n-tl+1;
                w = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];
            end
        end
        
        function [w] = hannw_(n)
            w = 0.5-0.5*cos(2*pi*[0:n-1].'/(n-1));
        end % hanning window
        
        function [w] = gaussw_(n,r) % gaussian window
            N = n-1; n = (0:N)'-N/2; w = exp(-(1/2)*(r*n/(N/2)).^2);
        end
        
        function [y] = rect_(x,dx) % rectangular function 

          y = mod(floor(x/dx),2);

        end

        function y = saw_(x,dx,a,b) % saw function 

        %y = sawfct(x,dx,a,b)
        % saw tooth function on R with period dx onto [a,b]

        x = x/dx-floor(x/dx);
        y = a + (b-a)*x;

        end
        
        function y = step_(x,nmax) % step function 

            %y = stepfct(x,nmax)
            % integer step function of x with period 1 such that [0,1] --> [1,nmax]

            x = x-floor(x);
            y = floor(nmax*x)+1;

            y = int16(y);
            y = y+int16(y==0);

        end

        
        % integration

        function [X,W] = simplex_quad_(N,vert)
            % example: integrate (x+y+z).^2 over tetrahedron
            % V=[1/sqrt(3) 0 0; -sqrt(3)/6,1/2,0;-sqrt(3)/6,-1/2,0;0 0 sqrt(6)/3];
            % n = 3; V=eye(n+1,n); 
            % [X,W]=simplex_quad_(n,V); 
            % Q=W'*(X(:,1)+X(:,2)+X(:,3)).^2
            % 
            % % plot stuff
            % figure(1);set(gcf,'color','w'); hold on;
            % % plot points
            % scatter3(X(:,1),X(:,2),X(:,3),'o');
            % % plot convex hull
            % scatter3(V(:,1),V(:,2),V(:,3),'.k');
            % CH=convhull(V(:,1),V(:,2),V(:,3));
            % h=trisurf(CH,V(:,1),V(:,2),V(:,3));
            % h.FaceColor='k'; h.FaceAlpha=0.02;
            % daspect([1 1 1]);
            % hold off;

            [m,n]=size(vert); Nn=N^n;
            if n==1 % The 1-D simplex is only an interval
                [q,w]=rquad(N,0); len=diff(vert);
                X=vert(1)+len*q;  W=abs(len)*w;
            else % Find quadrature rules for higher dimensional domains     
                for k=1:n 
                    [q{k},w{k}]=rquad(N,n-k);
                end
                [Q{1:n}]=ndgrid(q{:}); q=reshape(cat(n,Q{:}),Nn,n);
                [W{1:n}]=ndgrid(w{:}); w=reshape(cat(n,W{:}),Nn,n);
                map=eye(m); map(2:m,1)=-1; c=map*vert;
                W=abs(det(c(2:m,:)))*prod(w,2);
                qp=cumprod(q,2); e=ones(Nn,1);
                X=[e [(1-q(:,1:n-1)),e].*[e,qp(:,1:n-2),qp(:,n)]]*c;
            end
            function [x,w]=rquad(N,k)
                k1=k+1; k2=k+2; v=1:N;  nnk=2*v+k;
                A=[k/k2 repmat(k^2,1,N)./(nnk.*(nnk+2))];
                v=2:N; nnk=nnk(v);
                B1=4*k1/(k2*k2*(k+3)); nk=v+k; nnk2=nnk.*nnk;
                B=4*(v.*nk).^2./(nnk2.*nnk2-nnk2);
                ab=[A' [(2^k1)/k1; B1; B']]; s=sqrt(ab(2:N,2));
                [V,Y]=eig(diag(ab(1:N,1),0)+diag(s,-1)+diag(s,1));
                [Y,I]=sort(diag(Y));    
                x=(Y+1)/2; w=(1/2)^(k1)*ab(1,2)*V(1,I)'.^2;
            end
        end
        
        
        % functions that should of existed
        
        function n   = gcd_(n)
            if any(n==0)
                n = am_lib.gcd_(n(n~=0));
            else
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
         
        function [C] = rndstr_(n)
            set = char(['a':'z','0':'9','_','!','@','#','$','%','^']); nsets = numel(set);
            C = set(ceil(nsets*rand(1,n)));
        end

        function x   = ndims_(A)
            % number of dimensions after squeezing
            x = max(sum(size(A)~=1),1);
        end
        
        
        % geometric functions
        
        function L      = isorthogonal_(v1,v2,tol)
            if nargin<3;tol=am_lib.tiny;end
            L = abs(dot(v1,v2)./norm(v1)./norm(v2))<tol;
        end
        
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
        
        function R      = rotz_(t)
            c = cos(t); s = sin(t);
            R = [c -s 0;
                 s  c 0;
                 0  0 1];
        end
        
        function R      = roty_(t)
            c = cos(t); s = sin(t);
            R = [ c 0 s; 
                  0 1 0; 
                 -s 0 c];
        end
        
        function R      = rotx_(t)
            c = cos(t); s = sin(t);
            R = [1 0  0; 
                 0 c -s;
                 0 s  c];
        end
        
        function [A]    = R_axis_(R)
            % [A] = R_axis_(R);
            
            import am_lib.*
            
            nRs = size(R,3);
            
            % convert U(2) double group to O(3) representation
            if size(R,1)==2
                R = SU2_to_SO3(R); R = wdv_(R);
            end
            
            % get rotation axis
            if nRs == 1
                % convert (im)proper rotation to proper rotation
                R = R*sign(det(R));

                % define basic parameters and functions
                tol = 1E-8; normalize_ = @(v) v/norm(v); 

                % check for identity
                if abs(trace(R)-3)< tol; A=[0;0;1]; return; end

                % get rotation axis
                A = null_(R-eye(3));

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
            import am_lib.*
            
            % convert U(2) double group to O(3) representation
            if size(R,1)==2
                R = SU2_to_SO3(R); R = wdv_(R);
            end
            
            n = size(R,3);
            if n==1
                A = acos((trace(R)-1)/2); 
            else
                for i = 1:n
                    A(i) = am_lib.R_angle_(R(:,:,i));
                end
            end
        end
        
        function R      = R_align_(A,B)
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

        function v_rot  = R_vec_(v,k,theta)
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

        function J      = J_(j)
            % Angular momentum operators in 2j+1 dimensions. For j = 1/2, J are the generators of SU(2). For j = 1, J are the generators of SO(3).            
            % The number of matrices needed to generate a Lie group is the same as the dimension of the group; SU(n) groups have dimension n^2 - 1; SO(n) groups have dimension n(n-1)/2
            
            % define kronecker delta
            d_ = @(x,y) logical(x==y); 
            
            % matrix indices
            [m,mp]=meshgrid([j:-1:-j]);
            
            % define angular momentum operators (Jp raising, Jm lowering, ...)
            Jm = d_(m,mp+1).*sqrt((j+m).*(j-m+1)); Jp = d_(m,mp-1).*sqrt((j-m).*(j+m+1));
            Jx = (Jp+Jm)/2; Jy = (Jp-Jm)/2i; Jz = d_(m,mp).*m; J = cat(3,Jx,Jy,Jz); % diagonal in Jz basis
        end
        
        function [W]    = get_wigner(j,R,flag)
            %
            % The D in Wigner D-matrices stands for Darstellung, which means "representation" in German.
            %
            % The group SO(3) is the set of all three dimensional, real orthogonal matrices with unit
            % determinant. [Requiring that the determinant equal 1 and not -1 excludes inversion.]
            % The group O is obtained from the direct product of SO(3) and the the group Ci = {E,I}. 
            % Inversions matrix representatives in j are formed using the parity operator, (-1)^j, see:
            %
            %    - P. Jacobs, Group Theory with Applications in Chemical Physics
            %       (Cambridge University Press, 2005), p 208, eq 11.8.2. 
            %
            % The group SO(3) represents the set of all possible rotations of a three dimensional real
            % vector; it is essentially the collection of all proper rotations. The orthogonal group O(3) 
            % includes improper rotations and it is given by the direct product of SO(3) and the inversion
            % group i.
            %
            %   - R. M. Martin, Electronic Structure: Basic Theory and Practical Methods, 1 edition
            %        (Cambridge University Press, Cambridge, UK?; New York, 2008), p 573.
            %   - Romero Nichols "Density Functional Study of Fullerene-Based Solids: Crystal Structure,
            %        Doping, and Electron- Phonon Interaction", Ph.D. Thesis UIUC, p 96.
            %   - C. Cohen-Tannoudji, B. Diu, and F. Laloe, Quantum Mechanics, 1 edition (Wiley-VCH, New
            %        York; Paris, 1992), p 666.
            %   - J. J. Sakurai, Modern Quantum Mechanics, Revised edition (Addison Wesley, Reading, Mass,
            %        1993). p 207
            %
            % SU(2) is in the Cartan gauge, namely the inversion operator is I =  -1i*eye(2). In the Pauli gauge,
            % the inversion operator is I = -eye(2). Excellent discussion about gauges in:   
            % 
            %   - O. Chalaev, arXiv cond-mat.mtrl-sci, (2012).
            %
                
            import am_lib.*

            % get angular momentum operators [Jx,Jy,Jz]
            J = J_(j);

            % batch convert all rotation matrices to wigner functions
            nRs = size(R,3); W = zeros(2*j+1,2*j+1,nRs);
            for i = [1:nRs]; W(:,:,i) = get_wigner_engine(J,j,R(:,:,i)); end

            % convert to real harmonics
            if     contains(flag,'tesseral') || contains(flag,'real')
                % define basis change: spherical (complex) to tesseral harmonics (real basis)
                % C. Görller-Walrand and K. Binnemans, in (Elsevier, 1996), pp. 145, eqs. 14-16.
                [m,mp]=meshgrid([j:-1:-j]); d_ = @(x,y) logical(x==y); t_ = @(x,y) logical(x>y);
                T = d_(0,m) .* d_(mp,m) + ...
                    t_(m,0) .* sqrt(-1/2) .* ( d_(m,-mp) - (-1).^m .* d_(m, mp) ) + ... % sine   terms
                    t_(0,m) .* sqrt( 1/2) .* ( d_(m, mp) + (-1).^m .* d_(m,-mp) );      % cosine terms

                % reorder to recover SO(3) rotations
                O = circshift([1:(2*j+1)],floor((2*j+1)/2)); T = T(:,O);

                % do the actual conversion
                W = matmul_(matmul_(T',W),T);

                % make sure the teseral harmonics are real
                if any_(~eq_(imag(W),0))
                    if ~eq_(mod_(j),0)
                        error('Real harmonics only make sense for j = integer.'); 
                    else
                        error('Real harmonics have non-negligible imaginary components. Something is wrong');
                    end
                end
                
                W = real(W);
            end
            
            % engine
            function [W_sph] = get_wigner_engine(J,j,R)
                % Note: for l = 1, Wtes_(R) = R

                % remove inversion compontent to get pure rotation
                d = sign(det(R)); dR = R*d;

                % get rotation axis and angle
                an = am_lib.R_angle_(dR); ax = am_lib.R_axis_(dR);
                
                % define spin-vector dot products [Eq. 1.37 Blundell]
                dotV_ = @(S,V) S(:,:,1)*V(1) + S(:,:,2)*V(2) + S(:,:,3)*V(3);

                % compute spherical harmonics
                %   explicit matrix exponent version:
                %     [V,D] = eig( -1i * dotV_(J,ax) * an); 
                %     W_sph = V*diag(exp(diag(D)))/V;
                W_sph = expm(-1i * dotV_(J,ax) * an);
                
                % reinstate inversion to recover improper rotation if required
                W_sph = (d).^(j) * W_sph;
            end
        end
        
        function [R]    = SU2_to_SO3(W)
            
            import am_lib.det_
            
            R_ = @(x,y) [ ...
                 real(x^2-y^2),     imag(x^2+y^2),      -2*real(x*y);
                -imag(x^2-y^2),     real(x^2+y^2),       2*imag(x*y);
             2*real(x*conj(y)), 2*imag(x*conj(y)), abs(x)^2-abs(y)^2];
             n = size(W,3); R = zeros(3,3,n); 
             % convert to proper rotation
             d = permute(det_(W),[1,3,2]); W = d.^(1/2) .* W;
             % convert to SO(3)
             for i = 1:n; R(:,:,i) = R_(W(1,1,i),W(1,2,i)); end    
             % restore inversion
             R = d.^(1) .* R;
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

        function [x,y]  = pol2cart_(th,r)
            x = r.*cos(th);
            y = r.*sin(th);
        end
        
        function [th,r] = cart2pold_(x,y)
            th = atan2d(y,x); r = hypot(x,y);
        end
        
        function [th,r] = cart2pol_(x,y)
            th = atan2(y,x); r = hypot(x,y);
        end
        
        function [phi,chi,r] = cart2sphd_(x,y,z)
            hypotxy = hypot(x,y); r = hypot(hypotxy,z);
            chi = atan2d(z,hypotxy); phi = atan2d(y,x);
        end
        
        
        % matrix related
        
        function [A,A_inv] = circulant_(v)
            % get the circulant matrix corresponding to vector v and its inverse
            switch 1 % algo
                case 1; A = toeplitz(v(:),circshift(flipud(v(:)),1).');
                case 2; n = numel(v); F = dftmtx(n); A = F*diag(fft(v))*F'./n;
            end
            % get inverse of A if requested
            if nargout==2
                switch 0 % algo
                    case 0; A_inv = am_lib.am_lib.circulant_(ifft(1./fft(v))); % complete fft-based solution
                    case 1; n = numel(v); F = dftmtx(n); A_inv = F'*1./(fft(v)*n)*F; % using FFT
                    case 2; n = numel(v); F = dftmtx(n); A_inv = F*diag(1./diag(F'*A*F))*F'; % explicitly (diagonalize, flip, revert to original basis)
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
        
        function [C] = conv_(A,B,flag)
            % A = rand(6); B = rand(5);
            % conv_(A,B,'full')  - conv2(A,B,'full')
            % conv_(A,B,'valid') - conv2(A,B,'valid')
            %
            % test:
            %
            % clear;clc;
            % % fft convolution demo
            % n = 50; m = 60;
            % A=zeros(n);B=zeros(m);
            % A(20:30,25:26) = 1;
            % B(23:26,23:26) = 1;
            % C = convn(A,B,'same'); size(A)
            % % A = padarray(A,(m-size(A))/2); size(A)
            % % B = padarray(B,(m-size(B))/2); size(B)
            % % D = fftn(  fftn(A).*fftn(B) );
            % % F = fftshift((ifftn(  fftn(A).*fftn(B) )));
            % F = am_lib.conv_(A,B,'same');
            % norm(C(:)-F(:))
            % subplot(2,2,1);imagesc(A); daspect([1 1 1]);
            % subplot(2,2,2);imagesc(B); daspect([1 1 1]);
            % subplot(2,2,3);imagesc(C); daspect([1 1 1]);
            % subplot(2,2,4);imagesc(F); daspect([1 1 1]);
            % 

            switch flag
                case 'full'
                    [m1]=size(A); [m2]=size(B); m = m1+m2-1;
                    C = ifftn(fftn(A,m).*fftn(B,m));
                case 'same'
                    [m1] = size(A); [m2] = size(B); 
                    m = m1 + m2 - 1; a = ceil((m2-1)./2);
                    % pad and fft
                    C = ifftn(fftn(A,m).* fftn(B,m)); 
                    % unpad
                    for i = 1:numel(a); x_{i} = a(i)+1:m1(i)+a(i); end
                    C = C(x_{:});
                case 'valid'
                    [m1]=size(A); [m2]=size(B); m = m1+m2-1;
                    C = ifftn(fftn(A,m).*fftn(B,m));
                    C = C(m2(1):m1(1),m2(2):m1(2));
                otherwise
                    error('unknown flag');
            end
        end
        
        
        % matrix properties
        
        function [L] = isvector_(A)
            L = any(numel(A)==size(A));
        end
        
        function [L] = ismatrix_(A)
            [n,m]=size(A); L = any(numel(A)==n*m);
        end
        
        function [L] = isdiagdom_(A)
            L = all((2*abs(diag(A))) >= sum(abs(A),2));
        end
        
        function [L] = ishermitian_(A,tol)
            if nargin==1; tol=am_lib.tiny; end
            L = max(max(abs(A'-A))) < tol;
        end
        
        function [L] = isunitary_(A,tol)
            if nargin==1; tol=am_lib.tiny; end
            L = max(max(abs( A'*A - eye(size(A)) ))) < tol;
        end
        
        function [L] = issymmetric_(A,tol)
            if nargin==1; tol=am_lib.tiny; end
            L = max(max(abs(A-A.'))) < tol;
        end
        
        function [L] = isdiagonal_(A,tol)
            if nargin==1; tol=am_lib.tiny; end
            L = max(max(abs(diag(diag(A))-A))) < tol;
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
        
        
        % A x = b solvers
        
        function x = tridiag_(A,b)
            %  Solve the  n x n  tridiagonal system for y:
            %
            %  [ a(1)  c(1)                                  ] [  y(1)  ]   [  f(1)  ]
            %  [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
            %  [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
            %  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
            %  [                    ...    ...    ...        ] [        ]   [        ]
            %  [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
            %  [                                 b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]
            %
            %  f must be a vector (row or column) of length n
            %  a, b, c must be vectors of length n (note that b(1) and c(n) are not used)s
            n = length(b);
            % get tridiagonal vectores
            a = diag(A); b = diag(A,-1); c = diag(A,+1);
            % initialize loop
            v = zeros(n,1); y = v; w = a(1); y(1) = b(1)/w;
            % solve tridiagonal system of equations
            for i=2:n
                v(i-1) = c(i-1)/w;
                w = a(i) - b(i)*v(i-1);
                y(i) = ( b(i) - b(i)*y(i-1) )/w;
            end
            for j=n-1:-1:1
               y(j) = y(j) - v(j)*y(j+1);
            end
        end
        
        function [x,converged] = ddsolver_(A,x,b,m,algo)
            % break condition m = [tol, max iterations]
            if ~am_lib.isdiagdom_(A); error('ddsolver_ only works on diagonally dominant matrices'); end
            %
            if isempty(x); x = zeros(size(A,2),1); end
            if isempty(b); b = zeros(size(A,2),1); end
            % set break condition m = [tol, max iterations]
            if numel(m) == 1
                if m < 1
                    [tol,maxiter] = deal(m,Inf);
                else  
                    [tol,maxiter] = deal(0,m);
                end
            else
                [tol,maxiter] = deal(m(1),m(2));
            end
            % A = D + L + U; Notice the different defifinition from Multigrid Tutorial
            LDU_ = @(A,n) deal( (tril(A)-spdiags(diag(A),0,n,n)), ...
                                         spdiags(diag(A),0,n,n) , ...
                                (triu(A)-spdiags(diag(A),0,n,n))); 
            % initialize variables
            i = 0; n = size(A,1); 
            % define iterative A x = b solver of the form: x(j+1) = R*x(j) + B;
            switch algo
                case {'J','jacobi'}
                    [L,D,U] = LDU_(A,n); R = -D\(L+U); B = D\b;
                case {'Jw','weighted-jacobi'}
                    [L,D,U] = LDU_(A,n); w = 2/3;
                    R = -D\(L+U); R = (1-w)*speye(n) + w*R; B = (D/w)\b;
                case {'SOR','successive-over-relaxation'}
                    [L,D,U] = LDU_(A,n); w = 2/3;
                    if am_lib.issymmetric_(A) 
                        % Symmetric Successive Over-Relaxation
                        % Templates p 12
                        B1 = (D+w*U)\(-w*L+(1-w)*D);
                        B2 = (D+w*L)\(-w*U+(1-w)*D);
                        R = B1*B2; B = w*(2-w) * ((D+w*U)\D) * ((D+w*L)\b);
                    else
                        % Successive Over-Relaxation
                        % Templates p 11
                        B = (D+w*L)\(w*b); R = (D+w*L)\((1-w)*D-w*U);
                    end
                case {'GS','gauss-seidel'}
                    [L,D,U] = LDU_(A,n); R = -(D+L)\U; B = (D+L)\b; 
            end
            % iterate
            while true
                xp = x; x = R*xp + B;
                if     norm(x-xp)           < tol; converged = true;  return;
                elseif i > maxiter;                converged = false; return;
                else   i = i+1; end
            end
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

        function [x,r] = ga_(cost_,x0,isfixed,islinked,varargin)

            nvars = numel(x0);
            if isempty(isfixed);  isfixed  = zeros(1,nvars); end
            if isempty(islinked); islinked = zeros(1,nvars); end
            
            % check input dimensions
            if numel(x0)~=numel(isfixed); error('x0 and isfixed size mismatch'); end
            if numel(x0)~=numel(islinked); error('x0 and islinked size mismatch'); end
            if (nargin-4)>5; if numel(x0)~=numel(varargin{5}); error('x0 and lb size mismatch'); end; end
            if (nargin-4)>6; if numel(x0)~=numel(varargin{6}); error('x0 and ub size mismatch'); end; end
            
            % unlinked variables linked to 0
            islinked(islinked==0) = max([0,islinked(islinked~=0)])+[1:sum(islinked==0)];
            
            % find which variables are linked
            [~,i2u,u2i]=unique(islinked,'stable');
            
            % if one linked variable is fixed, make all fixed
            for i = 1:max(u2i); if any(isfixed(u2i==i)); isfixed(u2i==i)=1; end; end
            
            % 1) use only 1 variable per linked variables
            x0 = x0(i2u); isfixed = isfixed(i2u); 
            if (nargin-4)>4; varargin{5} = varargin{5}(i2u); end % lb
            if (nargin-4)>5; varargin{6} = varargin{6}(i2u); end % ub
            if (nargin-4)>7; varargin{8}.InitialPopulationMatrix = varargin{8}.InitialPopulationMatrix(i2u); end % x0
            
            % 2) use only floating variables
            if (nargin-4)>4; varargin{5} = varargin{5}(~isfixed); end % lb
            if (nargin-4)>5; varargin{6} = varargin{6}(~isfixed); end % ub
            if (nargin-4)>7; varargin{8}.InitialPopulationMatrix = varargin{8}.InitialPopulationMatrix(~isfixed); end % x0
            
                % fixed and loose indicies
                f = find( isfixed); xf = x0(f);
                l = find(~isfixed); xl = x0(l);

                % Estimate only the non-fixed ones
                [xl,r] = ga(@localcost_,sum(isfixed(:)==0),varargin{:});

                % Re-create array combining fixed and estimated coefficients
                x([f,l]) = [xf,xl];

            % restore all varaibles (removing links)
            x = x(u2i);

            function y = localcost_(x)
                b([f,l]) = [xf,x]; b = b(u2i); y = cost_(b);
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
        
        function [c,f,l] = fit_(x,y,flag)
            import am_lib.*
            if contains(flag,'lorentz')
                % define function
                func_ = @(c,x) c(1).*lorentz_((x-c(2)).*c(3)) + c(5);
                FWHM_ = @(c) 2/c(3); 
                % name parameter
                label = {'Amp','Center','Width''Background'};
                % estimate parameters
                [~,j] = max(conv(y,ones(1,20)/20,'same'));
                x0    = [ max(y), x(j), 1E-1, 0.5, 0]; 
                % define rescaling
                fscale_= @(c) [log(c(1)),c(2:4),log(c(5))];
                rscale_= @(c) [exp(c(1)),c(2:4),exp(c(5))];
            elseif contains(flag,'gauss')
                % define function
                func_ = @(c,x) c(1).*gauss_((x-c(2)).*c(3)) + c(5);
                FWHM_ = @(c) 2*sqrt(log(2))*c(3); 
                % name parameter
                label = {'Amp','Center','Width''Background'};
                % estimate parameters
                [~,j] = max(conv(y,ones(1,20)/20,'same'));
                x0    = [ max(y), x(j), 1E-1, 0.5, 0]; 
                % define rescaling
                fscale_= @(c) [log(c(1)),c(2:4),log(c(5))];
                rscale_= @(c) [exp(c(1)),c(2:4),exp(c(5))];
            elseif contains(flag,'pvoigt')
                % define function
                func_ = @(c,x) c(1).*pvoigt_((x-c(2)).*c(3),c(4)) + c(5);
                FWHM_ = @(c) 2/c(3); 
                % name parameter
                label = {'Amp','Center','Width','Lorentzian Fraction','Background'};
                % estimate parameters
                [~,j] = max(conv(y,ones(1,20)/20,'same'));
                x0    = [ max(y), x(j), 1E-1, 0.5, 0]; 
                % define rescaling
                fscale_= @(c) [log(c(1)),c(2:4),log(c(5))];
                rscale_= @(c) [exp(c(1)),c(2:4),exp(c(5))];
            else
                error('unknown profile');
            end

            % define cost function
%             cost_ = @(c) abs(log(func_(rscale_(c),x)) - log(y(:))).*y(:).^5; % weigh very top heavy
            cost_ = @(c) abs(log(func_(rscale_(c),x)) - log(y(:))).*y(:); % weigh top heavy
%             cost_ = @(c) abs(log(func_(rscale_(c),x)) - log(y(:))); % uniform weight

            % optimization options
            opts_ = optimoptions(@lsqnonlin,'Display','none','MaxIterations',1E10,'StepTolerance',1E-18,'FunctionTolerance',1E-18);

            % optimize
            c = lsqnonlin(cost_,fscale_(x0),[-Inf -Inf -Inf -Inf -Inf],[Inf Inf Inf 1 Inf],opts_); c = rscale_(c); f = func_; l = label;
            
            
            % plot if no output is requested
            if nargout == 0
                figure(1); set(gcf,'color','w'); plot(x,func_(c,x),'-k',x,y,'.-','linewidth',1); 
                title(sprintf('FWHM = %g',FWHM_(c))); % set(gca,'yscale','log');
            end
            
            % only report the FWHM
            if contains(flag,'FWHM')
                c = FWHM_(c); 
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
                    k = find(rebase_(j-1,mdigits)-rebase_(j-2,mdigits)>0);
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
        
        
        % image processing
        
        function [Ir,th]= level_horizon_(Ir,flag)
            % rotate image?
            switch flag
                case 'interactively'
                    am_lib.imagesc_(Ir); title('Select two points to level horizon'); p = ginput(2); 
                    if isempty(p); th=0; else; dp = diff(p.',1,2); th = atan2d(dp(2),dp(1))+90; end
                case 'fft';   th = -fft_(Ir,0.1);
                case 'hough'; th = -hough_(Ir,0.1);
                otherwise; error('unknown flag');
            end
            Ir = imrotate(Ir,th,'bicubic');
            if nargout==0; am_lib.imagesc_(Ir); title(['th = ', num2str(th)]); end

            function th = fft_(image, precision)
                % FFT.
                maxsize = max(size(image));
                T = fftshift(fft2(image, maxsize, maxsize)); % create rectangular transform
                T = log(abs(T)+1);                           % get magnitude in <0..inf)  

                % Combine two FFT quadrants together (another two quadrants are symetric).
                center = ceil((maxsize+1)/2); evenS = mod(maxsize+1, 2);
                T = (rot90(T(center:end, 1+evenS:center), 1) + T(center:end, center:end)); 
                T = T(2:end, 2:end); T(1,:)=0; T(:,1)=0;
                
                % Find the dominant orientation
                angles = floor(90/precision); score = zeros(angles, 1); maxDist = maxsize/2-1;
                for th = 0:angles-1
                    [y,x] = pol2cart(deg2rad(th*precision), 0:maxDist-1); % all [x,y]
                    for i = 1:maxDist
                        score(th+1) = score(th+1) + T(round(y(i)+1), round(x(i)+1));
                    end
                end

                % Return the most dominant direction.
                [~, position] = max(score); th = (position-1)*precision;
            end

            function th = hough_(image, precision)
                % Detect edges.
                BW = edge(image,'prewitt');

                % Perform the Hough transform.
                [H, T, ~] = hough(BW,'Theta',-90:precision:90-precision);  

                % Find the most dominant line direction.
                data=var(H);                      % measure variance at each angle 
                fold=floor(90/precision);         % assume right angles & fold data
                data=data(1:fold) + data(end-fold+1:end);
                [~, column] = max(data);          % the column with the crispiest peaks
                th = -T(column);               % column to degrees 
            end
        end
        
        function [I]    = low_pass_fourier_filter(I,wn,wm)
            % [wn, wm]=deal(21,19);
%             w = am_lib.hannw_(wn).*am_lib.hannw_(wm).';
            w = am_lib.gaussw_(size(I,1),wn).*am_lib.gaussw_(size(I,2),wm).';
            I = real(ifftn(fftn(I) .* fftshift(w))); % low pass filter
        end
        
        
        % general plotting
        
        function        set_plot_defaults_()
            set(groot,'defaultFigureColor','w');
            set(groot,'defaultLineLineWidth',1.5);
            set(groot,'defaultAxesLineWidth',1.5);
            set(groot,'defaultLineMarkerSize',9);
            set(groot,'defaultAxesFontSize',16);
            set(groot,'defaultAxesFontName','Helvetica');
            set(gcf,'color','w');
        end
        
        function variabilityplot_(x,y,varargin) 
            % x = randn(400,1);
            % y1 = nominal(randi(2,400,1),{'little','lots'});
            % y2 = nominal(randi(3,400,1),{'large','medium','small'});
            % y3 = nominal(randi(2,400,1),{'aardvark','potato'},[1,2]);
            % y = [y1,y2,y3];
            % variabilityplot_(x,y)

            % get sizes
            n = size(y,2); numgrps = zeros(1,n);
            for k = 1:n; numgrps(k) = numel(unique(y(:,k))); end
            numgrps = cumprod(numgrps); N = numgrps(n); y = fliplr(y);
            % plot 
            boxplot(x,y,...
                'plotstyle','compact','labelorientation','horizontal',...
                'factorseparator',1:n,varargin{:}); % set(findobj(gca,'Type','text'),'FontSize',16)
            % get handels
            hbxplt = get(gca,'children'); hall = get(hbxplt,'children'); halltype = get(hall,'type'); hsepln = hall(end-n+1:end);
            htxt = hall(strcmpi('text',halltype)); set(htxt,'units','data')
            txtpos = get(htxt,'position'); txtpos = cat(1,txtpos{:}); txtpos(:,2) = flipud(txtpos(:,2));

            % position label text appropriately
            x = reshape(txtpos(:,1),N,n); 
            for k = 2:n
                m = numgrps(k-1);
                for j = 1:N
                    ii = floor((j-1)/m);
                    i1 = 1 + m*ii;
                    i2 = m*(1+ii);
                    x(j,k) = mean(x(i1:i2,1));
                end
            end
            txtpos(:,1) = x(:);
            for k = 1:length(htxt); set(htxt(k),'position',txtpos(k,:)); end
            % draw label lines 
            tlcol = 0.5*[1,1,1]; txtpos = get(htxt,'extent'); txtpos = cat(1,txtpos{:});
            xl = xlim; yl = ylim; y1 = min(yl); y2 = min(txtpos(:,2)); y = linspace(y1,y2,n+1);
            % draw label box's horizontal lines
            for k = 2:(n+1)
                line(xl,[y(k),y(k)],'parent',gca,'clipping','off','color',tlcol);
            end
            % draw label box's vertical lines
            line(xl(1)*[1,1],[y1,y2],'parent',gca,'clipping','off','color',tlcol);
            line(xl(2)*[1,1],[y1,y2],'parent',gca,'clipping','off','color',tlcol);
            % draw vertical seperators
            for j = 1:n
                newy = get(hsepln(j),'YData');
                newy(newy==yl(2)) = y(j+1);
                line(get(hsepln(j),'XData'),newy,'parent',gca,'clipping','off','color',tlcol);
            end
            delete(hsepln(1))
        end
        
        function [h,ax] = overlay_(I,varargin) % overlay_(I,x,y): plot y(x) ontop of image I; overlay_(I,C,V): plot V contours of C ontop of image I 
            % plots x,y on top of a background image I 
            
            if     am_lib.isvector_(varargin{1})
                [x,y]=deal(varargin{:});
                % plot field as background
                h(1) = am_lib.imagesc_(flipud(I.')); hold on; drawnow;
                % overlay statistical function on field
                ax(1) = gca; ax(2) = axes('position',get(gca,'position'));
                % plot function
                h(2) = loglog(x,y,'-w'); 
                % remove white background
                set(gca,'color','none'); 
                % make sure it's the correct size
                ar = size(I); pbaspect([ar(2:-1:1),1]); 
            elseif am_lib.ismatrix_(varargin{1}) && am_lib.isvector_(varargin{2})
                switch numel(varargin)
                    case 2; [C,V]=deal(varargin{:});
                    case 1; [C]=deal(varargin{:}); V=[];
                    otherwise; error('invalid input');
                end
                % over lay the strain on the image
                h(1) = am_lib.imagesc_(flipud(I.')); hold on; drawnow; 
                % overlay statistical function on field
                ax(1) = gca; ax(2) = axes('position',get(gca,'position')); linkaxes(ax); colormap(ax(1),colormap('gray'));
                % plot isostrain
%                 if isempty(V); [~,h(2)]=am_lib.imcontour_(flipud(C.')); else; [~,h(2)]=am_lib.imcontour_(flipud(C.'),V); end
%                 if isempty(V); [~,h(2)]=contourf(flipud(C.')); else; [~,h(2)]=contourf(flipud(C.'),V); end
                am_lib.imagesc(flipud(C.'),'alphadata',rescale(I));
                % remove white background
                set(gca,'color','none'); set(h(2),'linewidth',1); colormap(ax(2),flipud(am_lib.colormap_('red2blue'))); 
            else
                error('invalid input');
            end
            
        end  
        
        function [h] = plot_isosurface_(X,dX,V,Vp,C,w,m,flag) 
            %  X = [(x,y,z),1:n(1),1:n(2),1:n(3)] coordinates
            % dX = translations for periodic boundary conditions (if X is defined from [0,1) then dX is probably eye(3))
            %  V = [(1:m)  ,1:n(1),1:n(2),1:n(3)] volumetric data
            % Vp = [1:q] list of points to sample V at
            %  w = number of times to extend the boundary
            %  m = number of times to refine the mesh
            % flag = 'cspline/lspline, center, extend, color'

            % extend weights?
            if ~isempty(w); flag=[flag,',extend']; end
            
            % add color to plot?
            if isempty(C); C = [1:size(V,1)].'.*ones(1,size(V,2),size(V,3),size(V,4)); end
            
            % loop over bands, removing ones that don't matter. (speeds things up >10000x)
            ex_ = false(size(V,1),1);
            for i = 1:size(V,1); for j = 1:numel(Vp)
                if max(V(i,:))>Vp(j) && min(V(i,:))<Vp(j)
                    ex_(i) = true;
                end
            end; end
            V = V(ex_,:,:,:); C = C(ex_,:,:,:);
    
            % up-sample?
            if     contains(flag,'cspline'); upflag = 'cspline';
            elseif contains(flag,'lspline'); upflag = 'lspline'; 
            else;  m = 1; upflag = 'none'; end 
            switch upflag
                case 'none'
                    % do nothing
                case {'lspline','cspline'}
                    % apply p.b.c to E and C
                    ww = [1,1,1]; 
                    X = cat(2,X,X(:,1:ww(1),:,:)+dX(:,1)); V = cat(2,V,V(:,1:ww(1),:,:)); C = cat(2,C,C(:,1:ww(1),:,:)); 
                    X = cat(3,X,X(:,:,1:ww(2),:)+dX(:,2)); V = cat(3,V,V(:,:,1:ww(2),:)); C = cat(3,C,C(:,:,1:ww(2),:)); 
                    X = cat(4,X,X(:,:,:,1:ww(3))+dX(:,3)); V = cat(4,V,V(:,:,:,1:ww(3))); C = cat(4,C,C(:,:,:,1:ww(3)));
                    % upscale
                    n = [size(V,2),size(V,3),size(V,4)]-1; 
                    Kup = zeros([size(X,1),n*2^(m-1)+1]); Eup = zeros([size(V,1),n*2^(m-1)+1]); Cup = zeros([size(C,1),n*2^(m-1)+1]);
                    for i = 1:size(X,1); Kup(i,:,:,:) = interpn(permute(X(i,:,:,:),[2,3,4,1]),m-1,upflag); end; X = Kup;
                    for i = 1:size(V,1); Eup(i,:,:,:) = interpn(permute(V(i,:,:,:),[2,3,4,1]),m-1,upflag); end; V = Eup;
                    for i = 1:size(C,1); Cup(i,:,:,:) = interpn(permute(C(i,:,:,:),[2,3,4,1]),m-1,upflag); end; C = Cup;
                otherwise
                    error('unknown interpolation method');
            end

            % move gamma to center?
            if contains(flag,'center')
                V = circshift(V,floor(size(V)/2)); 
            end

            % extend on both sides?
            if contains(flag,'extend')
                w = w*m;
                V = cat(2,cat(2,V(:,end-w(1):end,:,:),V),V(:,1:w(1),:,:)); C = cat(2,cat(2,C(:,end-w(1):end,:,:),C),C(:,1:w(1),:,:)); X = cat(2,cat(2,X(:,end-w(1):end,:,:)-dX(:,1),X),X(:,1:w(1),:,:)+dX(:,1)); 
                V = cat(3,cat(3,V(:,:,end-w(2):end,:),V),V(:,:,1:w(2),:)); C = cat(3,cat(3,C(:,:,end-w(2):end,:),C),C(:,:,1:w(2),:)); X = cat(3,cat(3,X(:,:,end-w(2):end,:)-dX(:,2),X),X(:,:,1:w(2),:)+dX(:,2)); 
                V = cat(4,cat(4,V(:,:,:,end-w(3):end),V),V(:,:,:,1:w(3))); C = cat(4,cat(4,C(:,:,:,end-w(3):end),C),C(:,:,:,1:w(3))); X = cat(4,cat(4,X(:,:,:,end-w(3):end)-dX(:,3),X),X(:,:,:,1:w(3))+dX(:,3)); 
            end

            % convert ndgrid to meshgrid [meshgrid([1:2],[3:5],[6:10]) - permute(ndgrid([1:2],[3:5],[6:10]),[2,1,3])]
            X = permute(X,[1,3,2,4]); V = permute(V,[1,3,2,4]); C = permute(C,[1,3,2,4]);

            % plot isosurface
            figure(1); clf; set(gcf,'color','w'); delete(findall(gcf,'Type','light')); hold on; b=0; h=[];
            for j = 1:numel(Vp); for i = 1:size(V,1); if max(V(i,:))>Vp(j) && min(V(i,:))<Vp(j); b=b+1;
                h(b) = patch(isosurface(squeeze(X(1,:,:,:)),squeeze(X(2,:,:,:)),squeeze(X(3,:,:,:)),squeeze(V(i,:,:,:)),Vp(j), squeeze(C(i,:,:,:))                   ));
                set(h(b),'FaceColor','interp','EdgeColor','none','AmbientStrength',0.3,'DiffuseStrength',1,'SpecularStrengt',0.4,'SpecularExponent',30);
            end; end; end
        end
        
        function [h] = plotc_(x,y,c,w) % line plot which changes color (and width)
            if nargin<4
                x = x(:).'; y=y(:).'; c=c(:).'; z=zeros(size(x));
                h = surface([x;x],[y;y],[z;z],[c;c],'facecol','no','edgecol','interp','linew',1);
            else
                x = x(:).'; y=y(:).'; c=c(:).'; w=w(:).'/2; z=zeros(size(x));
                h = surface([x;x],[y-w;y+w],[z;z],[c;c],'edgecol','interp','linew',1);
            end
            % rotate sideways? 
            % view([90 -90])
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

        function [h] = quiverc_(varargin)
            
            if     nargin == 8
                [X,Y,Z,U,V,W,C,clist]=deal(varargin{:});
                h = quiver3(X, Y, Z, U, V, W);
            elseif nargin == 5
                [U,V,W,C,clist]=deal(varargin{:});
                h = quiver3(U, V, W);
            elseif nargin == 6
                [X,Y,U,V,C,clist]=deal(varargin{:});
                % center
                    [T,R] = am_lib.cart2pol_(U,V); I=abs(X(:).'-X(:)); S = min(I(I(:)~=0))./max(R(R(:)~=0));
                    [U,V] = am_lib.pol2cart_(T,R*S);
                    % X = X - U/2; Y = Y - V/2;
                h = quiver(X, Y, U, V, 0.5, 'linewidth', 1.5);
            elseif nargin == 4
                [U,V,C,clist]=deal(varargin{:});
                h = quiver(U, V);
            elseif nargin == 3
                % just update the color
                [h,C,clist]=deal(varargin{:});
            end

            %// Now determine the color to make each arrow using a colormap
            [~, ~, ind] = histcounts(C, size(clist, 1));

            %// Now map this to a colormap to get RGB
            cmap = uint8(ind2rgb(ind(:), clist) * 255);
            cmap(:,:,4) = 255;
            cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

            %// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
            set(h.Head, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

            %// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
            set(h.Tail, ...
                'ColorBinding', 'interpolated', ...
                'ColorData', reshape(cmap(1:2,:,:), [], 4).');
        end
        
        function [h] = spyc_(A)
            [x,y] = find(A);
            h = scatter(y,x,200,A(A~=0),'.');
            set(gca,'YDir','rev'); box on;
            ylim([0 size(A,1)+1]); xlim([0 size(A,2)+1]); 
            set(gca,'XTick',[]); set(gca,'YTick',[]);
            daspect([1 1 1])
        end
        
        function [h] = imagesc_(A,varargin)
            h=imagesc(A,varargin{:}); axis tight; daspect([1 1 1]); axis off; box on; view([0 0 1]); daspect([1 1 1]);
        end

        function [h] = contour_(A,V,threshold,varargin) % contour_(matrix,levels,threshold,varargin)
            % histogram(c(3,:));set(gca,'yscale','log');

            [a] = contourc(A,V); axis tight; daspect([1 1 1]); axis off; box on; view([0 0 1]); daspect([1 1 1]);
            
            c = get_contour_line_properties(a);
            
            for n = 1:size(c,2); if c(3,n)>threshold
                sl_ = c(4,n)+1:c(5,n); h(n) = patch(a(1,sl_),a(2,sl_),c(1,n));
            end; end

            function c = get_contour_line_properties(a)
                % quick pass to get number of contours
                i=1; j=0; while i < size(a,2); j=j+1; i=i+a(2,i)+1; end; ncontours=j;
                % c(level,nnodes,length,start,end)
                i=1; j=0; c = zeros( 5, ncontours ); 
                while i < size(a,2) 
                    j=j+1; c(1,j)=a(1,i); c(2,j)=a(2,i); c(4,j)=i; c(5,j)=i+c(2,j);
                    c(3,j)=sum(am_lib.normc_(diff(a(:,i+[1:c(2,j)]),1,2)));
                    i=c(5,j)+1;
                end
                [~,inds]=sort(-c(3,:)); c=c(:,inds);
            end
        end
        
        function [varargout] = imcontour_(A,varargin)
            [varargout{1:2}]=imcontour(A,varargin{:}); axis tight; daspect([1 1 1]); axis off; 
        end
        
        function [h] = hist2_(x,y)
            [v,e] = hist3([x(:) y(:)],[1,1]*200); [e{:}]=meshgrid(e{1:2});
            h = contourf(e{1},e{2},log(v.'),100,'edgecolor','none');
        end

        function       frame2gif_(F,fname)
            for i = 1:numel(F)
                I = frame2im(F(i)); 
                [imind,cm] = rgb2ind(I,256); 
                if i == 1; imwrite(imind,cm,fname,'gif','DelayTime',0,'Loopcount',inf); 
                else;      imwrite(imind,cm,fname,'gif','DelayTime',0,'WriteMode','append'); 
                end
            end
            % loop back
            for i = numel(F)-1:-1:2
                I = frame2im(F(i));  [imind,cm] = rgb2ind(I,256); 
                imwrite(imind,cm,fname,'gif','DelayTime',0,'WriteMode','append'); 
            end
        end
        
        
        % coloring
        
        function hsl = hsv2hsl_(hsv)
            hsl = hsv;
            hsl(:,3) = 0.5 .* hsv(:,3) .* (2 - hsv(:,2));
            hsl(:,2) = hsv(:,3) .* hsv(:,2) ./ (1 - abs(2.*hsl(:,3)-1));
        end

        function hsv = hsl2hsv_(hsl)
            hsv = hsl; 
            hsv(:,3) = (2.*hsl(:,3) + hsv(:,2).*(1-abs(2.*hsl(:,3)-1)))/2;
            hsv(:,2) = 2.*(hsv(:,3) - hsl(:,3))./hsv(:,3);
        end
        
        function hsl = rgb2hsl_(rgb)
            hsl = am_lib.hsv2hsl_(rgb2hsv(rgb));
        end
        
        function rgb = hsl2rgb_(hsl)
            rgb = hsv2rgb(am_lib.hsl2hsv_(hsl));
        end
        
        function [cmap] = colormap_(flag,n)
            % switch based on N
            switch flag
                case 'discrete'
                    switch n
                        case {1}
                            cmap = [55, 126, 184]./255;
                        case {2, 3, 4, 5}
                            cmap = [[   47, 107, 156];
                                    [  228,  26,  28];
                                    [   77, 175,  74];
                                    [  255, 127,   0];
                                    [  152,  78, 163]]./255;
                            cmap = am_lib.clight_(cmap(1:n,:),0.20);
                        case {6 , 7, 8, 9}
                            cmap = [[ 228,  26,  28];
                                    [  55, 126, 184];
                                    [  77, 175,  74];
                                    [ 255, 127,   0];
                                    [ 217, 202,  94];
                                    [ 166,  86,  40];
                                    [ 247, 129, 191];
                                    [ 153, 153, 153];
                                    [ 152,  78, 163]]./255;
                            cmap = am_lib.clight_(cmap(1:n,:),0.10);
                        case {10, 11, 12}
                            cmap = [[ 141, 211, 199];
                                    [ 255, 237, 111];
                                    [ 190, 186, 218];
                                    [ 251, 128, 114];
                                    [ 128, 177, 211];
                                    [ 253, 180,  98];
                                    [ 179, 222, 105];
                                    [ 188, 128, 189];
                                    [ 217, 217, 217];
                                    [ 204, 235, 197];
                                    [ 252, 205, 229];
                                    [ 255, 255, 179]]./255;
                            cmap = am_lib.clight_(cmap(1:n,:),-0.07);
                        otherwise 
                            error('n is too big');
                    end
                case 'jet'
                    cmap = jet(5);
                case 'hsv' % for plotting phase, based on nist recommendation, puts primary colors at 90 degs.
                  x = hsv(900);
                  cmap(1:150,:)  =x(1:150,:);
                  cmap(151:300,:)=x(151:2:450,:);
                  cmap(301:450,:)=x(451:600,:);
                  cmap(451:600,:)=x(601:2:900,:);
                case 'spectral'
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
                case {'r2b','red2blue'} % sns.cmap_palette("RdBu_r", 7)
                    cmap = [ ... 
                    0.16339870177063293, 0.44498270983789490, 0.6975009791991290;
                    0.42068437209316328, 0.67643216077019186, 0.8186851319144753;
                    0.76147636946509856, 0.86851211856393251, 0.9245674785445717;
                    0.96908881383783674, 0.96647443490869855, 0.9649365649503820;
                    0.98246828247519102, 0.80069205340217142, 0.7061130509657018;
                    0.89457901435739851, 0.50380624217145586, 0.3997693394913390;
                    0.72848905885920801, 0.15501730406985564, 0.1973856272650700];
                case 'thermal' 
                    cmap = [1.555601333154079877e-02 1.382442454646408414e-01 2.018108864558305071e-01
                    1.620183633850513089e-02 1.410507428866217272e-01 2.089765125440807836e-01
                    1.685648942708358952e-02 1.438270143621834252e-01 2.162386804476043589e-01
                    1.752640064782528401e-02 1.465717250667996630e-01 2.235996996833259920e-01
                    1.821871873545745021e-02 1.492834638238061673e-01 2.310618693528040390e-01
                    1.894137836902154426e-02 1.519607349643580241e-01 2.386274839825403005e-01
                    1.969967580211434005e-02 1.546014513385217670e-01 2.463049741539924953e-01
                    2.050331512714091350e-02 1.572037794856419868e-01 2.540971092002056730e-01
                    2.136720981691051999e-02 1.597664500864496573e-01 2.619991459036808412e-01
                    2.230340677460975612e-02 1.622875536944013708e-01 2.700132114112391291e-01
                    2.332520459934065912e-02 1.647650548255009395e-01 2.781413941166567261e-01
                    2.444727823684482437e-02 1.671967824349823439e-01 2.863857270216258466e-01
                    2.568581629165060318e-02 1.695804203313354408e-01 2.947481661145329168e-01
                    2.705867185137392564e-02 1.719134976584987262e-01 3.032305630924797546e-01
                    2.858552764394725262e-02 1.741933796245943022e-01 3.118346316168387755e-01
                    3.028807626872338093e-02 1.764172587161002559e-01 3.205619061457981589e-01
                    3.219021609777513587e-02 1.785821467123839268e-01 3.294136922283489310e-01
                    3.431826321592867240e-02 1.806848679100377109e-01 3.383910069732901649e-01
                    3.670117943457189280e-02 1.827220540829147533e-01 3.474945082255425644e-01
                    3.937081594703384368e-02 1.846901418460877020e-01 3.567244107935257369e-01
                    4.230474116360182640e-02 1.865853732643254770e-01 3.660803878810173217e-01
                    4.544128431656369732e-02 1.884038007524381497e-01 3.755614556934170345e-01
                    4.879889459256544354e-02 1.901412975600472455e-01 3.851658390246698871e-01
                    5.238564999909287728e-02 1.917935754201579024e-01 3.948908155057742619e-01
                    5.620896927989652708e-02 1.933562112707872260e-01 4.047325361348014794e-01
                    6.027560971317643540e-02 1.948246853300575621e-01 4.146858197468099028e-01
                    6.459519137924149557e-02 1.961877478810097886e-01 4.247714594384868758e-01
                    6.917293885235362150e-02 1.974458334676339466e-01 4.349572770994127868e-01
                    7.401397924875777190e-02 1.985943734629224688e-01 4.452322549865626589e-01
                    7.912632580773251711e-02 1.996251445642635014e-01 4.555965555444902448e-01
                    8.452074659970570947e-02 2.005284212595413451e-01 4.660508732253974551e-01
                    9.019392390291772199e-02 2.013079448693051998e-01 4.765478779135125520e-01
                    9.616430834359415702e-02 2.019472529514862447e-01 4.871044455010897223e-01
                    1.024253987466881288e-01 2.024520213616345932e-01 4.976646226688278829e-01
                    1.089944270249512126e-01 2.028088855832538839e-01 5.082270885561682716e-01
                    1.158597351626668159e-01 2.030273547727350913e-01 5.187245256305992314e-01
                    1.230424263443430921e-01 2.030937970306834206e-01 5.291483757425788914e-01
                    1.305276694613609623e-01 2.030217816424806365e-01 5.394201195013325068e-01
                    1.383099125215292158e-01 2.028195613734770086e-01 5.494767764101642360e-01
                    1.463797060150685281e-01 2.025001769277504637e-01 5.592461312000727158e-01
                    1.547186326280874380e-01 2.020850680852618320e-01 5.686422620965899677e-01
                    1.632970499077517346e-01 2.016054641193363861e-01 5.775676861911989146e-01
                    1.720728231749732162e-01 2.011027291278179585e-01 5.859183845130113699e-01
                    1.809917622181062835e-01 2.006270524284321510e-01 5.935918186078640302e-01
                    1.899902199574210471e-01 2.002341583564870575e-01 6.004970669758479263e-01
                    1.989997441940680178e-01 1.999802927301375655e-01 6.065651516258429021e-01
                    2.079529779081165097e-01 1.999164256738651391e-01 6.117570630961892686e-01
                    2.167895160216459782e-01 2.000830340237900185e-01 6.160673924886361785e-01
                    2.254604257150365498e-01 2.005067081654643424e-01 6.195227918364496489e-01
                    2.339306302023402284e-01 2.011991984057693306e-01 6.221760825969698816e-01
                    2.421790683949211487e-01 2.021587162723321729e-01 6.240979346009031259e-01
                    2.501971277034634733e-01 2.033727314506766359e-01 6.253682433520739714e-01
                    2.579861074778473928e-01 2.048213417591974728e-01 6.260687943604991146e-01
                    2.655544225081862275e-01 2.064804694302365684e-01 6.262779744769025880e-01
                    2.729150389729214088e-01 2.083244557603688429e-01 6.260675767301320249e-01
                    2.800833915559314269e-01 2.103279194298604549e-01 6.255013265069014894e-01
                    2.870675074033041674e-01 2.124691386407265292e-01 6.246420009551492125e-01
                    2.938851360627945386e-01 2.147254304183257023e-01 6.235380183242871244e-01
                    3.005577213461362307e-01 2.170756514619715527e-01 6.222251903623999825e-01
                    3.070843773107185259e-01 2.195059802015210115e-01 6.207554984633752992e-01
                    3.134916285549604886e-01 2.219983100391352271e-01 6.191452159520036691e-01
                    3.197798380866583856e-01 2.245420356602937373e-01 6.174343435010863912e-01
                    3.259695527216341926e-01 2.271241237863139140e-01 6.156329324235468858e-01
                    3.320579089479565038e-01 2.297373625924510887e-01 6.137771467682858750e-01
                    3.380660052498126178e-01 2.323715331127139128e-01 6.118660583919267593e-01
                    3.439917420525137604e-01 2.350215076527270575e-01 6.099280657295180763e-01
                    3.498460685618723365e-01 2.376808945869600675e-01 6.079694872801727490e-01
                    3.556384580181179977e-01 2.403444043084654869e-01 6.059953497614897211e-01
                    3.613686295558476425e-01 2.430089153836439420e-01 6.040241385974032262e-01
                    3.670434468015162932e-01 2.456707832167414618e-01 6.020607616318822686e-01
                    3.726708785615932551e-01 2.483266964722255499e-01 6.001059364008434205e-01
                    3.782547954216434749e-01 2.509743622345630976e-01 5.981656274234778969e-01
                    3.837960822688916696e-01 2.536122834643630419e-01 5.962498766549472196e-01
                    3.892988345534371120e-01 2.562387736058999166e-01 5.943618049976240325e-01
                    3.947690991030736729e-01 2.588520700720948198e-01 5.924996351953921714e-01
                    4.002098946086877773e-01 2.614510659770345469e-01 5.906660913059953444e-01
                    4.056240909379103532e-01 2.640348396856425084e-01 5.888633104220148962e-01
                    4.110144223076953041e-01 2.666026246434279878e-01 5.870929096766683841e-01
                    4.163832172801885667e-01 2.691538344975527575e-01 5.853565951388559618e-01
                    4.217313231007376317e-01 2.716882525666767800e-01 5.836583632684495537e-01
                    4.270634345309600177e-01 2.742050306264073312e-01 5.819941042354170868e-01
                    4.323818104897113601e-01 2.767038027770706843e-01 5.803639256565888971e-01
                    4.376886243188534142e-01 2.791842625143060586e-01 5.787676418773346487e-01
                    4.429859696643737577e-01 2.816461519867726748e-01 5.772048027503156042e-01
                    4.482758654902720408e-01 2.840892531259122666e-01 5.756747186760452495e-01
                    4.535602602863915700e-01 2.865133804028727749e-01 5.741764823913096949e-01
                    4.588410355259065487e-01 2.889183750000149931e-01 5.727089879419176022e-01
                    4.641200084236148382e-01 2.913041002126984247e-01 5.712709472315609105e-01
                    4.693989340425177015e-01 2.936704379215548943e-01 5.698609044983501404e-01
                    4.746795067933508583e-01 2.960172859965477521e-01 5.684772490345112450e-01
                    4.799633613697112389e-01 2.983445565121622955e-01 5.671182264323321176e-01
                    4.852520731601697723e-01 3.006521746683436525e-01 5.657819486102849682e-01
                    4.905471581781303825e-01 3.029400783247205298e-01 5.644664028469250638e-01
                    4.958500725502145712e-01 3.052082180664184574e-01 5.631694600262432404e-01
                    5.011622116043696895e-01 3.074565577287760587e-01 5.618888822762294621e-01
                    5.064849085998188727e-01 3.096850753156841218e-01 5.606223301621769961e-01
                    5.118194331420761189e-01 3.118937642524056697e-01 5.593673695773668797e-01
                    5.171669893275794294e-01 3.140826349187781363e-01 5.581214784559522801e-01
                    5.225287136638799845e-01 3.162517164128481606e-01 5.568820534159194535e-01
                    5.279056728126164666e-01 3.184010584984696690e-01 5.556464164236831760e-01
                    5.332988612036093645e-01 3.205307336933521101e-01 5.544118215561065766e-01
                    5.387091985692399332e-01 3.226408394566178117e-01 5.531754619204082291e-01
                    5.441375274486375258e-01 3.247315004373090841e-01 5.519344767774559957e-01
                    5.495846107110635703e-01 3.268028707475315597e-01 5.506859588994191812e-01
                    5.550511291471261766e-01 3.288551362262008837e-01 5.494269621786239677e-01
                    5.605376791749926424e-01 3.308885166617086537e-01 5.481545094908906179e-01
                    5.660447559653561944e-01 3.329032737201366166e-01 5.468656290966650291e-01
                    5.715723572408055730e-01 3.348998704862335973e-01 5.455581199822507887e-01
                    5.771213104289429907e-01 3.368784486432290226e-01 5.442280104688511644e-01
                    5.826918609018115758e-01 3.388393806259411556e-01 5.428722728323092106e-01
                    5.882841610050709713e-01 3.407830824811499126e-01 5.414878910241607279e-01
                    5.938982692324897839e-01 3.427100153754324974e-01 5.400718694374621043e-01
                    5.995341497366565298e-01 3.446206869260000083e-01 5.386212416598933350e-01
                    6.051916721861749782e-01 3.465156523493675422e-01 5.371330791454219655e-01
                    6.108706119725975103e-01 3.483955154267821541e-01 5.356044997333838653e-01
                    6.165706507631435462e-01 3.502609292895154658e-01 5.340326759432469927e-01
                    6.222913773881607602e-01 3.521125970311436149e-01 5.324148429745680922e-01
                    6.280322890453092777e-01 3.539512721578863541e-01 5.307483063447554494e-01
                    6.337927927958686425e-01 3.557777588917046541e-01 5.290304491020394462e-01
                    6.395721854014132512e-01 3.575929224264804973e-01 5.272587829694529438e-01
                    6.453697498455218673e-01 3.593976455700171879e-01 5.254307632896803026e-01
                    6.511846984697030605e-01 3.611928549527803622e-01 5.235439096149940852e-01
                    6.570160963803246545e-01 3.629795579684892415e-01 5.215959730192853971e-01
                    6.628629283562820218e-01 3.647588132682739737e-01 5.195848109675399451e-01
                    6.687241012633005077e-01 3.665317306708469891e-01 5.175083912494201632e-01
                    6.745984464659754432e-01 3.682994710922303794e-01 5.153647953022981731e-01
                    6.804847221871188623e-01 3.700632465224873990e-01 5.131522209386339961e-01
                    6.863816157645579175e-01 3.718243200760849021e-01 5.108689845040222943e-01
                    6.922877457570544291e-01 3.735840061411623281e-01 5.085135225034096429e-01
                    6.982016638532075881e-01 3.753436706510840937e-01 5.060843927435032530e-01
                    7.041218565401857754e-01 3.771047314992945210e-01 5.035802750493091340e-01
                    7.100468845018954589e-01 3.788685911907690995e-01 5.009996451582281463e-01
                    7.159752518731292703e-01 3.806367017978663503e-01 4.983410633235164089e-01
                    7.219050254760910335e-01 3.824107548202065887e-01 4.956041045679995816e-01
                    7.278344332015346252e-01 3.841923879692217270e-01 4.927879483140180095e-01
                    7.337616407032357957e-01 3.859832960516177969e-01 4.898918977085310877e-01
                    7.396847510296571393e-01 3.877852324461281697e-01 4.869153794944616753e-01
                    7.456018038508585022e-01 3.896000107579556393e-01 4.838579440403712462e-01
                    7.515107742785041012e-01 3.914295066368485010e-01 4.807192656291401911e-01
                    7.574095712835625660e-01 3.932756597383701980e-01 4.774991431076147097e-01
                    7.632960357235062387e-01 3.951404758009958162e-01 4.741975009994948143e-01
                    7.691679379984673881e-01 3.970260288041012053e-01 4.708143911830061090e-01
                    7.750229753641669772e-01 3.989344631636088656e-01 4.673499952331199858e-01
                    7.808587689384376418e-01 4.008679959130103665e-01 4.638046275250513051e-01
                    7.866728604480874854e-01 4.028289188075652172e-01 4.601787391913814695e-01
                    7.924627087737056153e-01 4.048196002786729752e-01 4.564729230191116316e-01
                    7.982256863620216247e-01 4.068424871537119070e-01 4.526879193650572009e-01
                    8.039590755885861473e-01 4.089001060440320967e-01 4.488246231577860956e-01
                    8.096600651680347926e-01 4.109950642903736351e-01 4.448840920414381395e-01
                    8.153258420561760866e-01 4.131300072582320126e-01 4.408672005707100494e-01
                    8.209533846301936277e-01 4.153077124896901728e-01 4.367753702901882029e-01
                    8.265394587133882975e-01 4.175310847974926798e-01 4.326106584748535266e-01
                    8.320808476618816174e-01 4.198030498289601065e-01 4.283748912056176694e-01
                    8.375742298072511582e-01 4.221266063894242304e-01 4.240701127101388912e-01
                    8.430161781421201539e-01 4.245048216497959159e-01 4.196985997897492715e-01
                    8.484031611801808870e-01 4.269408247931844591e-01 4.152628771038663902e-01
                    8.537315452228415591e-01 4.294377989032187592e-01 4.107657331305429871e-01
                    8.589975982777037222e-01 4.319989708985695898e-01 4.062102365684474026e-01
                    8.641974958823694930e-01 4.346275993272235572e-01 4.015997528869852951e-01
                    8.693273290890853877e-01 4.373269598520185819e-01 3.969379606684441675e-01
                    8.743831148591909574e-01 4.401003282876529976e-01 3.922288673203608855e-01
                    8.793608090992056647e-01 4.429509610905758565e-01 3.874768236699652202e-01
                    8.842563225402819693e-01 4.458820732584555802e-01 3.826865368883226592e-01
                    8.890655396175203284e-01 4.488968136663485375e-01 3.778630811333164030e-01
                    8.937843693689974112e-01 4.519982404973394430e-01 3.730116014425700621e-01
                    8.984086715918495614e-01 4.551892900487020666e-01 3.681382748140558103e-01
                    9.029343899688814234e-01 4.584727367754222738e-01 3.632494215164466245e-01
                    9.073575586033109097e-01 4.618511663085561048e-01 3.583516614288673741e-01
                    9.116743277656134126e-01 4.653269421603676848e-01 3.534519722899168159e-01
                    9.158809992051784032e-01 4.689021723136149178e-01 3.485576699619067353e-01
                    9.199740625655835613e-01 4.725786767270368505e-01 3.436763818638891022e-01
                    9.239502320744639174e-01 4.763579567194722864e-01 3.388160133952918263e-01
                    9.278064825504417357e-01 4.802411672570275347e-01 3.339847073943023048e-01
                    9.315400836682339314e-01 4.842290931875970483e-01 3.291907969273191181e-01
                    9.351486313601243827e-01 4.883221304365225612e-01 3.244427519772630775e-01
                    9.386300752174523421e-01 4.925202730905588466e-01 3.197491208755643965e-01
                    9.419827407976890665e-01 4.968231071526066356e-01 3.151184675888816233e-01
                    9.452053458453827384e-01 5.012298115494860928e-01 3.105593062092800172e-01
                    9.482970095990234105e-01 5.057391667267315816e-01 3.060800341879720277e-01
                    9.512572545755481057e-01 5.103495708806106146e-01 3.016888659817511531e-01
                    9.540860004899643920e-01 5.150590635745326828e-01 2.973937688352701891e-01
                    9.567835502644689294e-01 5.198653561843505910e-01 2.932024023936338208e-01
                    9.593505683914355098e-01 5.247658683349885056e-01 2.891220637265379811e-01
                    9.617880522174245828e-01 5.297577692490467172e-01 2.851596391521897811e-01
                    9.640972969908582213e-01 5.348380227429672118e-01 2.813215639873033469e-01
                    9.662798557460373639e-01 5.400034344900018768e-01 2.776137910349112392e-01
                    9.683374952663609259e-01 5.452507001277228094e-01 2.740417682745930894e-01
                    9.702721494708056449e-01 5.505764528210163045e-01 2.706104258619907443e-01
                    9.720889605050128113e-01 5.559744240639272750e-01 2.673323057285149629e-01
                    9.737882391976736551e-01 5.614431075234689317e-01 2.642051375197135843e-01
                    9.753722517954447335e-01 5.669791681305653697e-01 2.612321043607351290e-01
                    9.768433321430202154e-01 5.725792962307277856e-01 2.584161030939570725e-01
                    9.782037550763590383e-01 5.782403184217296266e-01 2.557595328984585970e-01
                    9.794556966363962003e-01 5.839592236834146854e-01 2.532643241222839459e-01
                    9.806011989875685897e-01 5.897331838357899869e-01 2.509319710696591987e-01
                    9.816445119805891073e-01 5.955576445021050214e-01 2.487674614053176358e-01
                    9.825957966513861885e-01 6.014235301335163486e-01 2.467838886220111161e-01
                    9.834484089207010671e-01 6.073354334336248384e-01 2.449669167474177733e-01
                    9.842037798764916579e-01 6.132913094304376367e-01 2.433164315861336413e-01
                    9.848630901763529844e-01 6.192893396808167861e-01 2.418320393357376030e-01
                    9.854335855591603854e-01 6.253232095344845032e-01 2.405207443125048083e-01
                    9.859266054527725531e-01 6.313839151833324781e-01 2.393924955683320310e-01
                    9.863290788081210403e-01 6.374805817839866995e-01 2.384285005271009616e-01
                    9.866412575432023102e-01 6.436122274808632193e-01 2.376272705850792089e-01
                    9.868743031654042541e-01 6.497702193582076680e-01 2.369977788720153966e-01
                    9.870387195949218428e-01 6.559468364085242476e-01 2.365460008515470891e-01
                    9.871160198111447182e-01 6.621544902577213287e-01 2.362516434857291903e-01
                    9.871054832314208882e-01 6.683929123881801049e-01 2.361125093170109435e-01
                    9.870452158568817635e-01 6.746361467497961062e-01 2.361551640855951706e-01
                    9.869008071972434903e-01 6.809069196703931848e-01 2.363489608363473771e-01
                    9.866692878063002548e-01 6.872064512263130753e-01 2.366900769944741967e-01
                    9.863929517876753872e-01 6.935072010148546351e-01 2.372025181179537867e-01
                    9.860358081545986808e-01 6.998320059248218650e-01 2.378591308273384497e-01
                    9.855936903889445100e-01 7.061828277187244263e-01 2.386553467137359497e-01
                    9.851166760681809853e-01 7.125285611733751523e-01 2.396137846497558566e-01
                    9.845515404596055786e-01 7.189016348706930293e-01 2.407027465239837682e-01
                    9.839219208542443473e-01 7.252873073536283410e-01 2.419304971970357987e-01
                    9.832411577231264799e-01 7.316776102841332508e-01 2.432981102176038357e-01
                    9.824686360853925882e-01 7.380960010628512258e-01 2.447853973017022899e-01
                    9.816675490572134288e-01 7.445058075421049359e-01 2.464118597749720974e-01
                    9.807811261674898029e-01 7.509394073137781733e-01 2.481527624504580309e-01
                    9.798377198237012697e-01 7.573804700508821597e-01 2.500137146942527089e-01
                    9.788392654586509645e-01 7.638279008720277874e-01 2.519905378906656113e-01
                    9.777574675384567149e-01 7.702969047505190403e-01 2.540716600569991046e-01
                    9.766479350869351483e-01 7.767573771217611833e-01 2.562664070834250185e-01
                    9.754356977304895482e-01 7.832491288133777152e-01 2.585534361426328198e-01
                    9.742107025732685832e-01 7.897246981741067318e-01 2.609472103669700505e-01
                    9.728835679779427315e-01 7.962305102434532600e-01 2.634261983272434549e-01
                    9.715297353896644728e-01 8.027276438790406088e-01 2.659994820979044161e-01
                    9.700861348373609472e-01 8.092479763321948072e-01 2.686528986078252079e-01
                    9.686057099481873989e-01 8.157648637457888263e-01 2.713899293880046582e-01
                    9.670430691608270513e-01 8.223006128282058791e-01 2.742007328066633498e-01
                    9.654377156641301694e-01 8.288358297838601674e-01 2.770858440538391254e-01
                    9.637524161478506768e-01 8.353882936883590959e-01 2.800376262549509332e-01
                    9.620231537308009395e-01 8.419407901977008502e-01 2.830554448577937698e-01
                    9.602104885023859948e-01 8.485116373841136150e-01 2.861325388392048086e-01
                    9.583576301078738924e-01 8.550807238511934916e-01 2.892681513706822360e-01
                    9.564117528707279936e-01 8.616719638797333269e-01 2.924557053792812833e-01
                    9.544348419615947821e-01 8.682572954460299197e-01 2.956947937778980906e-01
                    9.523487086523741985e-01 8.748712522095133393e-01 2.989788164145827376e-01
                    9.502464564071694264e-01 8.814728114137494464e-01 3.023077648025994102e-01
                    9.480141307284709606e-01 8.881110895508603775e-01 3.056751963404329420e-01
                    9.457819798963852387e-01 8.947301765899248194e-01 3.090811208408512645e-01
                    9.434198646152930356e-01 9.013849121338212145e-01 3.125200216900483885e-01
                    9.410286159691669816e-01 9.080328518084608280e-01 3.159906396320988908e-01
                    9.385331275251967975e-01 9.147047487143539213e-01 3.194893313796646206e-01
                    9.359711083432949996e-01 9.213848126149789541e-01 3.230138415459749002e-01
                    9.333375985220991877e-01 9.280748341844049509e-01 3.265615469169362850e-01
                    9.305915654403119630e-01 9.347905093445134650e-01 3.301299810381577715e-01
                    9.278142496620381818e-01 9.414998557883345054e-01 3.337168327083647745e-01
                    9.248899165263079203e-01 9.482470197415310276e-01 3.373199677733378365e-01
                    9.219411401959984875e-01 9.549849147932994997e-01 3.409370492608514991e-01
                    9.188613878011584468e-01 9.617532494963948464e-01 3.445662997232426528e-01
                    9.156931782520092433e-01 9.685354904254356301e-01 3.482056900726151483e-01
                    9.124490701578419349e-01 9.753266872784461805e-01 3.518533597970244786e-01
                    9.090418416674036495e-01 9.821574063216705897e-01 3.555078064299531104e-01];
               case 'deep'
                    cmap = [9.928371765383620096e-01 9.943734553013935384e-01 8.001361955494933342e-01
                    9.849374457410008388e-01 9.913545172197536504e-01 7.953271573982337861e-01
                    9.770418482420034634e-01 9.883418276673759939e-01 7.905717472557142189e-01
                    9.691488355955474310e-01 9.853357520972024775e-01 7.858697059005903540e-01
                    9.612382061888059548e-01 9.823439661710483550e-01 7.812181455049906909e-01
                    9.533233708513803029e-01 9.793607550577946297e-01 7.766197776697654209e-01
                    9.454069776374501854e-01 9.763847481120443428e-01 7.720753845210462929e-01
                    9.374799177499437697e-01 9.734191074102575003e-01 7.675843461643131471e-01
                    9.295353395208965086e-01 9.704660061861888343e-01 7.631467470398153319e-01
                    9.215850348376466439e-01 9.675205198521025229e-01 7.587646801057830181e-01
                    9.136274030668488644e-01 9.645828701859077148e-01 7.544386825530340346e-01
                    9.056514559384564178e-01 9.616567079910524063e-01 7.501688017241657791e-01
                    8.976569922484586295e-01 9.587415918264282633e-01 7.459562027416458685e-01
                    8.896512126622908578e-01 9.558344072358454513e-01 7.418023698547455691e-01
                    8.816325974742758032e-01 9.529352564923027069e-01 7.377082218108865774e-01
                    8.735986658557798323e-01 9.500445573198598170e-01 7.336747582572585857e-01
                    8.655355018678880796e-01 9.471666805105979359e-01 7.297031220835451526e-01
                    8.574559441871676402e-01 9.442965454208677167e-01 7.257948454644865821e-01
                    8.493586413080321806e-01 9.414341108550433601e-01 7.219511807036801398e-01
                    8.412422934868647451e-01 9.385792958890862847e-01 7.181734629573161000e-01
                    8.331056602292961077e-01 9.357319779845385543e-01 7.144631077114116380e-01
                    8.249438571064254822e-01 9.328932140048489252e-01 7.108218561594963347e-01
                    8.167536485550084269e-01 9.300634438432798801e-01 7.072516022427790539e-01
                    8.085407925860420564e-01 9.272401607530683654e-01 7.037535913222745521e-01
                    8.003043988510905038e-01 9.244230580301131539e-01 7.003294929842703853e-01
                    7.920436855303579771e-01 9.216117772405966191e-01 6.969810377273069069e-01
                    7.837579923068553889e-01 9.188059058340631857e-01 6.937100117646464170e-01
                    7.754467945117471395e-01 9.160049747102424478e-01 6.905182510191696377e-01
                    7.671097184662585278e-01 9.132084557727712104e-01 6.874076341947011892e-01
                    7.587465580267132026e-01 9.104157595099175992e-01 6.843800748081978469e-01
                    7.503572923166539343e-01 9.076262326497597233e-01 6.814375120698907828e-01
                    7.419421046029871514e-01 9.048391559449654453e-01 6.785819005039074314e-01
                    7.335014022414247936e-01 9.020537421501462205e-01 6.758151982106467281e-01
                    7.250358375800118882e-01 8.992691342626382145e-01 6.731393536848668813e-01
                    7.165463296678448168e-01 8.964844041051395207e-01 6.705562911206142118e-01
                    7.080340865696405084e-01 8.936985513356989763e-01 6.680678941563225059e-01
                    6.995006280355556827e-01 8.909105029766873907e-01 6.656759880411459163e-01
                    6.909478082204841831e-01 8.881191135591984809e-01 6.633823202371271766e-01
                    6.823778380886152961e-01 8.853231659823359578e-01 6.611885395113676900e-01
                    6.737933070789577927e-01 8.825213731875580780e-01 6.590961736178440056e-01
                    6.651972035471812594e-01 8.797123807461469935e-01 6.571066057195337207e-01
                    6.565929334411800822e-01 8.768947704523611941e-01 6.552210497573525139e-01
                    6.479843366143701600e-01 8.740670650055498703e-01 6.534405250319036407e-01
                    6.393757001353049807e-01 8.712277338509245572e-01 6.517658303256151919e-01
                    6.307693113306602761e-01 8.683757593089920235e-01 6.501986071208971651e-01
                    6.221681809999368706e-01 8.655099786991408140e-01 6.487403091120448329e-01
                    6.135828162115329887e-01 8.626276255310663110e-01 6.473888753222819537e-01
                    6.050194699551968425e-01 8.597270135144634562e-01 6.461439080647282118e-01
                    5.964848373853082197e-01 8.568064537576968176e-01 6.450046593443751197e-01
                    5.879840953182007279e-01 8.538646600753627691e-01 6.439710905713145195e-01
                    5.795218127572832056e-01 8.509005270271372545e-01 6.430435569229606685e-01
                    5.711115213696185133e-01 8.479112775740021979e-01 6.422171169895777298e-01
                    5.627614837100215484e-01 8.448953489705807174e-01 6.414893854468299850e-01
                    5.544783665092802849e-01 8.418515860497588488e-01 6.408587965640016870e-01
                    5.462695843545278818e-01 8.387787842708503971e-01 6.403232223573597226e-01
                    5.381481232862911357e-01 8.356748924828519831e-01 6.398763669277579558e-01
                    5.301228248165350543e-01 8.325387580839158641e-01 6.395142187190840932e-01
                    5.221997330825246530e-01 8.293697946738645133e-01 6.392346216632595057e-01
                    5.143908532894702068e-01 8.261665755050031645e-01 6.390304309659169402e-01
                    5.067050118087177424e-01 8.229283397233977393e-01 6.388963162980869637e-01
                    4.991495064662723746e-01 8.196546689563362076e-01 6.388277739100475250e-01
                    4.917330703220189059e-01 8.163450533427777378e-01 6.388184799685107107e-01
                    4.844633050708381794e-01 8.129992646664957467e-01 6.388624455765857801e-01
                    4.773469739538908629e-01 8.096172922533247940e-01 6.389538165444101914e-01
                    4.703904182695624603e-01 8.061992786917829834e-01 6.390864379965265352e-01
                    4.635986539031755060e-01 8.027456073395251579e-01 6.392548734619343254e-01
                    4.569766618393767410e-01 7.992567478306369377e-01 6.394529319719780558e-01
                    4.505287643111087204e-01 7.957333261733031682e-01 6.396742566585853496e-01
                    4.442569334258392177e-01 7.921762343010834151e-01 6.399148201157260907e-01
                    4.381634653667463852e-01 7.885863788123392837e-01 6.401694292235835526e-01
                    4.322513883314121341e-01 7.849647158245569578e-01 6.404304190942151642e-01
                    4.265195639210723755e-01 7.813124732411278472e-01 6.406960731925235297e-01
                    4.209679602144861810e-01 7.776308288475485275e-01 6.409622276506293792e-01
                    4.155968407812616339e-01 7.739210344346504344e-01 6.412227317683711902e-01
                    4.104045242070183952e-01 7.701844253858858291e-01 6.414742347181462412e-01
                    4.053883482522956383e-01 7.664223146281371468e-01 6.417151578859574546e-01
                    4.005459328058830759e-01 7.626360316845265386e-01 6.419424749255564500e-01
                    3.958750262046074608e-01 7.588270315806920907e-01 6.421508285001665817e-01
                    3.913718273240804346e-01 7.549966666219802836e-01 6.423387337484185444e-01
                    3.870322615152646528e-01 7.511462095764531721e-01 6.425055472082883412e-01
                    3.828523536620307977e-01 7.472769743775978801e-01 6.426493131526664904e-01
                    3.788278449989629926e-01 7.433902546660990929e-01 6.427683061328856029e-01
                    3.749542811367644335e-01 7.394876286687900313e-01 6.428573615340694714e-01
                    3.712266123106295335e-01 7.355700686179795778e-01 6.429189618774799886e-01
                    3.676399929257251897e-01 7.316387718079003788e-01 6.429520487195032885e-01
                    3.641894512092864744e-01 7.276949011320993366e-01 6.429557270674552960e-01
                    3.608699207518579755e-01 7.237395816489646805e-01 6.429292709767442382e-01
                    3.576762479926827720e-01 7.197739138067411613e-01 6.428719883285997083e-01
                    3.546027248366774298e-01 7.157992458648064771e-01 6.427810759039616073e-01
                    3.516446005650431528e-01 7.118162697085735902e-01 6.426589560156208414e-01
                    3.487967412412643076e-01 7.078259414050183107e-01 6.425054728664090220e-01
                    3.460540467004642462e-01 7.038291754672099110e-01 6.423205839967465192e-01
                    3.434114694865876838e-01 6.998268449671132263e-01 6.421043512357870187e-01
                    3.408640318201912600e-01 6.958197818437616977e-01 6.418569321473704958e-01
                    3.384068406811324148e-01 6.918087773726725453e-01 6.415785719753461791e-01
                    3.360351011011302735e-01 6.877945827659479594e-01 6.412695960850323118e-01
                    3.337441277691141628e-01 6.837779098758932639e-01 6.409304028912300444e-01
                    3.315291322901112725e-01 6.797595093162658308e-01 6.405610222271946874e-01
                    3.293855731203558790e-01 6.757400463409192204e-01 6.401618521317685717e-01
                    3.273095072934468774e-01 6.717199943156246800e-01 6.397341828681051279e-01
                    3.252967835603527980e-01 6.676999151191671533e-01 6.392786349130804568e-01
                    3.233433949344308722e-01 6.636803359579184214e-01 6.387958727683986648e-01
                    3.214454825944703664e-01 6.596617502568347113e-01 6.382866005148117861e-01
                    3.195993388733257556e-01 6.556446185422359907e-01 6.377515576836865208e-01
                    3.178014094235615539e-01 6.516293693095515094e-01 6.371915154212616228e-01
                    3.160482946464280851e-01 6.476163998706647718e-01 6.366072729215248582e-01
                    3.143367504654775990e-01 6.436060771766987099e-01 6.359996541044424800e-01
                    3.126636885203824545e-01 6.395987386132337971e-01 6.353695045172913503e-01
                    3.110261758511674857e-01 6.355946927658273626e-01 6.347176884379418516e-01
                    3.094214341373541788e-01 6.315942201545601264e-01 6.340450861601970578e-01
                    3.078468385510386152e-01 6.275975739369872297e-01 6.333525914425652825e-01
                    3.062999162775806861e-01 6.236049805794583456e-01 6.326411091031564071e-01
                    3.047783447524610168e-01 6.196166404972232034e-01 6.319115527447410896e-01
                    3.032799496578738041e-01 6.156327286641396501e-01 6.311648425953061414e-01
                    3.018027027180503197e-01 6.116533951930802626e-01 6.304019034507418739e-01
                    3.003447193279823457e-01 6.076787658883653354e-01 6.296236627075371128e-01
                    2.989042560461285802e-01 6.037089427717098333e-01 6.288310484745529561e-01
                    2.974797079780254205e-01 5.997440045832753697e-01 6.280249877540684533e-01
                    2.960696060742640245e-01 5.957840072594756675e-01 6.272064046833601969e-01
                    2.946726143633343065e-01 5.918289843891730850e-01 6.263762188290394883e-01
                    2.932875271368584058e-01 5.878789476499071132e-01 6.255353435272825724e-01
                    2.919132661023912667e-01 5.839338872257091584e-01 6.246846842638799080e-01
                    2.905488775167375803e-01 5.799937722079637759e-01 6.238251370887518688e-01
                    2.891935293105973304e-01 5.760585509806988025e-01 6.229575870601862242e-01
                    2.878465082138687570e-01 5.721281515915097593e-01 6.220829067145717817e-01
                    2.865071747347951447e-01 5.682024945496734203e-01 6.212019079466474247e-01
                    2.851749861947763254e-01 5.642814857106783766e-01 6.203153697409601319e-01
                    2.838497050469647731e-01 5.603649539888247988e-01 6.194242689076802089e-01
                    2.825310677598035780e-01 5.564527485481920444e-01 6.185294158631611250e-01
                    2.812189151446758406e-01 5.525446999518639490e-01 6.176316008625786225e-01
                    2.799131896111074491e-01 5.486406204645675189e-01 6.167315925457064196e-01
                    2.786139324558262742e-01 5.447403043452749838e-01 6.158301364812402978e-01
                    2.773212811865227168e-01 5.408435281299215358e-01 6.149279537054621603e-01
                    2.760354668808263079e-01 5.369500509042623992e-01 6.140257392505246159e-01
                    2.747568115804877587e-01 5.330596145668072827e-01 6.131241606570347891e-01
                    2.734857257205313141e-01 5.291719440816851083e-01 6.122238564648385672e-01
                    2.722227055926091377e-01 5.252867477212617153e-01 6.113254346750160995e-01
                    2.709683308416185876e-01 5.214037172983150281e-01 6.104294711750373192e-01
                    2.697232619941725140e-01 5.175225283875982685e-01 6.095365081178055755e-01
                    2.684882380171457750e-01 5.136428405367017280e-01 6.086470522439324515e-01
                    2.672640739041831637e-01 5.097642974662250914e-01 6.077615731349965689e-01
                    2.660516582874621339e-01 5.058865272594796902e-01 6.068805013837572648e-01
                    2.648519510715739433e-01 5.020091425421505660e-01 6.060042266652778675e-01
                    2.636659810857353015e-01 4.981317406526836744e-01 6.051330956906521008e-01
                    2.624948437498288434e-01 4.942539038045714039e-01 6.042674100224807443e-01
                    2.613396987489819967e-01 4.903751992421883643e-01 6.034074237284018372e-01
                    2.602017811880152354e-01 4.864951750164123179e-01 6.025533553474511361e-01
                    2.590825816132030779e-01 4.826133005480903182e-01 6.017055804809932074e-01
                    2.579832374637801018e-01 4.787291639776031782e-01 6.008639770499764055e-01
                    2.569051367797046126e-01 4.748422758923547815e-01 6.000285749995881712e-01
                    2.558497193309375306e-01 4.709521329962904068e-01 5.991993407256713811e-01
                    2.548184740794410263e-01 4.670582183399953347e-01 5.983761730113198452e-01
                    2.538129365659009817e-01 4.631600015854966945e-01 5.975588985257048735e-01
                    2.528346862021874641e-01 4.592569393176054171e-01 5.967472668238706923e-01
                    2.518853434473967146e-01 4.553484754161292725e-01 5.959409447786497838e-01
                    2.509665668417216944e-01 4.514340415062205181e-01 5.951395103673859932e-01
                    2.500800498681828299e-01 4.475130575075850214e-01 5.943424457267542094e-01
                    2.492275176071283571e-01 4.435849323073600692e-01 5.935491293785424283e-01
                    2.484107231427581941e-01 4.396490645862106139e-01 5.927588275177814170e-01
                    2.476314436738927260e-01 4.357048438328084417e-01 5.919706842419948378e-01
                    2.468914762732849488e-01 4.317516515883728090e-01 5.911837105866416531e-01
                    2.461926332303804310e-01 4.277888629705645096e-01 5.903967722170088139e-01
                    2.455367369014649914e-01 4.238158485348740845e-01 5.896085756110598375e-01
                    2.449256139784778408e-01 4.198319765418777605e-01 5.888176525512366366e-01
                    2.443611399380560822e-01 4.158365980702596332e-01 5.880223891276964432e-01
                    2.438461792188084676e-01 4.118287207942545325e-01 5.872218507211212080e-01
                    2.433817108226277726e-01 4.078080077755985022e-01 5.864131553289283483e-01
                    2.429695407602900925e-01 4.037738393736669540e-01 5.855939512690606641e-01
                    2.426114394907277760e-01 3.997256137994158465e-01 5.847615988448080504e-01
                    2.423091294682429564e-01 3.956627542073260506e-01 5.839131390484204598e-01
                    2.420642707838099872e-01 3.915847171721096864e-01 5.830452587571144374e-01
                    2.418784445905285962e-01 3.874910027793410650e-01 5.821542521747482546e-01
                    2.417531339648231747e-01 3.833811665857123074e-01 5.812359783131098023e-01
                    2.416897018162703636e-01 3.792548337300427064e-01 5.802858143751018494e-01
                    2.416907565718227624e-01 3.751112253980636857e-01 5.792995847596160708e-01
                    2.417595622535183009e-01 3.709493751547713325e-01 5.782729391945388153e-01
                    2.418940089237943680e-01 3.667702659833258494e-01 5.771972660307567171e-01
                    2.420946965383475313e-01 3.625739944720036134e-01 5.760653956197256953e-01
                    2.423619198524908369e-01 3.583608468185064400e-01 5.748693688323588402e-01
                    2.426978313653707087e-01 3.541305556383861908e-01 5.736016010478480753e-01
                    2.431056167697961956e-01 3.498826161637011989e-01 5.722540340252517677e-01
                    2.435790408918375727e-01 3.456199709880173332e-01 5.708127541670769967e-01
                    2.441164681518231960e-01 3.413441041408621368e-01 5.692658625959313712e-01
                    2.447220255685461088e-01 3.370546900675063240e-01 5.676029832475817383e-01
                    2.453897790312707938e-01 3.327551263503047418e-01 5.658082210414479007e-01
                    2.461121481740833894e-01 3.284495884603169102e-01 5.638648470620498676e-01
                    2.468912126448456479e-01 3.241391993124722593e-01 5.617586120246663706e-01
                    2.477157002362878613e-01 3.198299095156284522e-01 5.594704462476630669e-01
                    2.485761162021524751e-01 3.155272777788998839e-01 5.569822129976372826e-01
                    2.494682375867288693e-01 3.112353605277175528e-01 5.542766318947622839e-01
                    2.503729202057394798e-01 3.069632816204772574e-01 5.513352372309007210e-01
                    2.512817651002536290e-01 3.027167370222472731e-01 5.481422664464320471e-01
                    2.521753740527091225e-01 2.985048904754492582e-01 5.446839245346258851e-01
                    2.530378388279626023e-01 2.943355684042912035e-01 5.409499814033098541e-01
                    2.538514672799476179e-01 2.902167617845944347e-01 5.369345375383627328e-01
                    2.545978636988186494e-01 2.861560770732184955e-01 5.326370186681820273e-01
                    2.552597334794323158e-01 2.821600341185547811e-01 5.280624026249696179e-01
                    2.558254018044066047e-01 2.782328077909240194e-01 5.232189983935346955e-01
                    2.562737222001675308e-01 2.743799915974786674e-01 5.181261586082565040e-01
                    2.566051651711689918e-01 2.706008890121900934e-01 5.127963599343589030e-01
                    2.568071232679735028e-01 2.668971447705413280e-01 5.072541841677337127e-01
                    2.568729082739400482e-01 2.632680806498599591e-01 5.015251413221865073e-01
                    2.568088786232387566e-01 2.597098435540172168e-01 4.956253550777797723e-01
                    2.566151212870083631e-01 2.562195453869526296e-01 4.895770754997314511e-01
                    2.562879598478859378e-01 2.527944349049762174e-01 4.834087259028412853e-01
                    2.558322794870629413e-01 2.494300496133808887e-01 4.771399728195321877e-01
                    2.552572079534596305e-01 2.461214417402471932e-01 4.707834286138410929e-01
                    2.545675711557803811e-01 2.428643843006248471e-01 4.643557571030281772e-01
                    2.537686502971923108e-01 2.396546885092911139e-01 4.578718126285494239e-01
                    2.528659964243425984e-01 2.364882912834416206e-01 4.513446141213755536e-01
                    2.518652755598765891e-01 2.333613180085921113e-01 4.447853910070804773e-01
                    2.507721436866522935e-01 2.302701237479483076e-01 4.382036814126420432e-01
                    2.495921493408267966e-01 2.272113165747541852e-01 4.316074651476269897e-01
                    2.483306606950595463e-01 2.241817667402334346e-01 4.250033167249620547e-01
                    2.469928136315570621e-01 2.211786051187170643e-01 4.183965667275956202e-01
                    2.455817950644241798e-01 2.181991605590304917e-01 4.117955160707330586e-01
                    2.440963016438912891e-01 2.152407079180826410e-01 4.052187485307188752e-01
                    2.425474613143962510e-01 2.123011804696942062e-01 3.986531809377649171e-01
                    2.409394174709658110e-01 2.093786088479052676e-01 3.921004515047561423e-01
                    2.392760085101002243e-01 2.064712127546999287e-01 3.855616338648150121e-01
                    2.375559510487991743e-01 2.035768990151702318e-01 3.790513023869773179e-01
                    2.357791582066215419e-01 2.006936504843463420e-01 3.725810406048237766e-01
                    2.339566833930070699e-01 1.978208439176197264e-01 3.661282774694271658e-01
                    2.320913707908846546e-01 1.949572807641995476e-01 3.596924727656225507e-01
                    2.301734316844865069e-01 1.921000566693645550e-01 3.533131710447778295e-01
                    2.282155151615900546e-01 1.892493734782149939e-01 3.469585312854190362e-01
                    2.262221013532286218e-01 1.864046213456004575e-01 3.406203233993785884e-01
                    2.241821302892643142e-01 1.835626593748518609e-01 3.343431537784405938e-01
                    2.221083620693462546e-01 1.807243380946543521e-01 3.280900612914751102e-01
                    2.200027865362182422e-01 1.778889096627078448e-01 3.218593457818205161e-01
                    2.178554139765243036e-01 1.750533261852178779e-01 3.156932341101472139e-01
                    2.156819223383368844e-01 1.722195412332822306e-01 3.095395593647118360e-01
                    2.134730781003147948e-01 1.693846918206770857e-01 3.034377234324792671e-01
                    2.112350828916275125e-01 1.665490510713014960e-01 2.973690772834624574e-01
                    2.089702734855372057e-01 1.637121940546225618e-01 2.913286801723345421e-01
                    2.066728560428846562e-01 1.608719311731062473e-01 2.853437530392856081e-01
                    2.043551063834350701e-01 1.580301135357929931e-01 2.793689192405632293e-01
                    2.020028382233076125e-01 1.551826890051333785e-01 2.734655011163185101e-01
                    1.996333580091984583e-01 1.523327614839000976e-01 2.675659276139747411e-01
                    1.972319620276792862e-01 1.494761585027865047e-01 2.617344119314569673e-01
                    1.948137533226730334e-01 1.466155558553730864e-01 2.559107021731873988e-01
                    1.923665760685599468e-01 1.437473901381048913e-01 2.501489823140411461e-01
                    1.899031789417767457e-01 1.408738040042367690e-01 2.443975812138057813e-01
                    1.874125186315797886e-01 1.379915476831030108e-01 2.387062880680187460e-01
                    1.849071752981109873e-01 1.351026683367795023e-01 2.330231779340583564e-01
                    1.823751828761793481e-01 1.322038111876145949e-01 2.274020989868827392e-01
                    1.798308834661569988e-01 1.292972741720705698e-01 2.217828373102059270e-01
                    1.772595702951434704e-01 1.263792859047284667e-01 2.162309361702130506e-01
                    1.746769766714075522e-01 1.234522307754042925e-01 2.106797904214151584e-01
                    1.720703275074147443e-01 1.205129174047056828e-01 2.051860710532931176e-01
                    1.694509066570421274e-01 1.175626362090423926e-01 1.997023294326863430e-01
                    1.668117590877864764e-01 1.145994139402238821e-01 1.942594655847655616e-01
                    1.641571774416621943e-01 1.116231764992316466e-01 1.888399155026514453e-01
                    1.614878059346168959e-01 1.086331782717616379e-01 1.834416540200605461e-01
                    1.587995102818766935e-01 1.056281331759973130e-01 1.780821476813713722e-01
                    1.561019746507273376e-01 1.026082525875711138e-01 1.727215696232307918e-01]; 
                case 'curl'
                    cmap = [8.225559928700268419e-02 1.149244079727295142e-01 2.647901677800857390e-01
                    8.312616532498406929e-02 1.190383729463048712e-01 2.668628892216621806e-01
                    8.400180885962132971e-02 1.231074880892656653e-01 2.689526699064171411e-01
                    8.487294239495335457e-02 1.271387529060027943e-01 2.710541708402016137e-01
                    8.574385298640457842e-02 1.311333174761502018e-01 2.731691209373900975e-01
                    8.661249189260347703e-02 1.350944971238551839e-01 2.752961432065319514e-01
                    8.747533041314431435e-02 1.390258052165279645e-01 2.774332852121961235e-01
                    8.833858505105957049e-02 1.429270910011002649e-01 2.795831537842536352e-01
                    8.919012906146844832e-02 1.468043594975814992e-01 2.817400195447572475e-01
                    9.004099984169053328e-02 1.506555099870153513e-01 2.839086654207542693e-01
                    9.088231952195491292e-02 1.544850037627045203e-01 2.860850125083750362e-01
                    9.171714479257989105e-02 1.582931942356169963e-01 2.882702798507874586e-01
                    9.254607948203208423e-02 1.620811665705463311e-01 2.904645825457518593e-01
                    9.336173420340779239e-02 1.658523274558137695e-01 2.926648040593683997e-01
                    9.417284157369981701e-02 1.696050767223750977e-01 2.948744162588470830e-01
                    9.496899572502048859e-02 1.733434725855039771e-01 2.970892156823352059e-01
                    9.575619444438937533e-02 1.770666290075273708e-01 2.993115284087380368e-01
                    9.653316478613682694e-02 1.807757562460599599e-01 3.015407883893915231e-01
                    9.729328810023782359e-02 1.844734436312414905e-01 3.037745103921251633e-01
                    9.804493118338306057e-02 1.881580898967480098e-01 3.060157402408537064e-01
                    9.877832247043097369e-02 1.918329889556127654e-01 3.082609202522414993e-01
                    9.949803783218733044e-02 1.954975058870322413e-01 3.105116763778897337e-01
                    1.002054858430543316e-01 1.991518634625495388e-01 3.127684262427183892e-01
                    1.008900241588315538e-01 2.027992883367119026e-01 3.150275553737127421e-01
                    1.015620277820400430e-01 2.064376400065636719e-01 3.172925125709148420e-01
                    1.022150454689789434e-01 2.100689913462399083e-01 3.195611195376787395e-01
                    1.028468845169810686e-01 2.136942776676366007e-01 3.218326666201724029e-01
                    1.034637486817424901e-01 2.173124203806320875e-01 3.241090419310820314e-01
                    1.040556247342821483e-01 2.209261472703716311e-01 3.263871052200472134e-01
                    1.046272664909374817e-01 2.245346804307608024e-01 3.286682704446855507e-01
                    1.051814556007013568e-01 2.281377376539149293e-01 3.309532535866355762e-01
                    1.057055066316897329e-01 2.317384538919968207e-01 3.332383062494437276e-01
                    1.062093640276061124e-01 2.353348932541681759e-01 3.355262248019345583e-01
                    1.066927192498384747e-01 2.389274243968986799e-01 3.378167764830257158e-01
                    1.071428862568302165e-01 2.425189920877232619e-01 3.401063918909076889e-01
                    1.075713140889580088e-01 2.461073946949010050e-01 3.423980926413837667e-01
                    1.079763885418836000e-01 2.496932209057829977e-01 3.446912766824791197e-01
                    1.083460004297124302e-01 2.532791299120676354e-01 3.469826761312421182e-01
                    1.086913329090623825e-01 2.568630402753477870e-01 3.492750549670107785e-01
                    1.090116287119651528e-01 2.604453159192556266e-01 3.515680214913912693e-01
                    1.092932808921200649e-01 2.640287615628073015e-01 3.538580573547516761e-01
                    1.095478763802813504e-01 2.676112773022403801e-01 3.561478615212890775e-01
                    1.097748910210358808e-01 2.711931376089152246e-01 3.584370856671010852e-01
                    1.099635342633853707e-01 2.747764822051754763e-01 3.607229892872333421e-01
                    1.101198294511886444e-01 2.783603021914852760e-01 3.630068099573608986e-01
                    1.102460152182262176e-01 2.819443225282238785e-01 3.652888319764409086e-01
                    1.103361572207214036e-01 2.855297163446686715e-01 3.675674754730824945e-01
                    1.103867349279119559e-01 2.891171759115234718e-01 3.698417464033080804e-01
                    1.104047114339987423e-01 2.927055780381122019e-01 3.721129484539063559e-01
                    1.103893448894039397e-01 2.962951553198432397e-01 3.743806401368951486e-01
                    1.103290785873510815e-01 2.998879083652261635e-01 3.766420819050353419e-01
                    1.102317288825286901e-01 3.034825829921489193e-01 3.788986922133656954e-01
                    1.100986315611906241e-01 3.070790351009552999e-01 3.811504545068737926e-01
                    1.099286556573810802e-01 3.106775192369479188e-01 3.833968205472055302e-01
                    1.097092559883022234e-01 3.142800415550652815e-01 3.856349205402823110e-01
                    1.094518269158459012e-01 3.178848431435130073e-01 3.878667847148202785e-01
                    1.091556968975302966e-01 3.214920812428978536e-01 3.900919340676701208e-01
                    1.088202096165073185e-01 3.251019032868823211e-01 3.923098843692563453e-01
                    1.084337294201929702e-01 3.287160421828903556e-01 3.945179917556369542e-01
                    1.080052802278752000e-01 3.323331728550899533e-01 3.967176863054378000e-01
                    1.075355258379605550e-01 3.359532282299735328e-01 3.989087127976621017e-01
                    1.070239159802940931e-01 3.395763147211236510e-01 4.010905648833596460e-01
                    1.064666158592602885e-01 3.432029795356145718e-01 4.032620996037361571e-01
                    1.058582376509708545e-01 3.468339423292918222e-01 4.054218838433404359e-01
                    1.052065591856250482e-01 3.504681464569319171e-01 4.075709632165070984e-01
                    1.045112201334900126e-01 3.541056562885774861e-01 4.097088092956981398e-01
                    1.037719209811681642e-01 3.577465259466254266e-01 4.118348877541359032e-01
                    1.029884320307208612e-01 3.613907991323849767e-01 4.139486582680417803e-01
                    1.021513134274243950e-01 3.650396452261166491e-01 4.160478671051655586e-01
                    1.012689167529539358e-01 3.686920150348446112e-01 4.181335233490640069e-01
                    1.003423119365809690e-01 3.723477936835040136e-01 4.202052548968338574e-01
                    9.937169919107982641e-02 3.760069784468507703e-01 4.222625006231348066e-01
                    9.835741884413440328e-02 3.796695548525758634e-01 4.243046936775203282e-01
                    9.729997123759509536e-02 3.833354963232004087e-01 4.263312615677605777e-01
                    9.620003933222870396e-02 3.870047637840662302e-01 4.283416262899890081e-01
                    9.505172587236093706e-02 3.906780327832828914e-01 4.303339883746660766e-01
                    9.386233813144639893e-02 3.943545459798681874e-01 4.323088724355508838e-01
                    9.263427266200571775e-02 3.980341179549674036e-01 4.342658626828024837e-01
                    9.136927887187012987e-02 4.017166535612498035e-01 4.362043674792987491e-01
                    9.006945702453716951e-02 4.054020428715323643e-01 4.381237915805220595e-01
                    8.873730500447504776e-02 4.090901606132017476e-01 4.400235366709063789e-01
                    8.737577058920215078e-02 4.127808655736769361e-01 4.419030019956885491e-01
                    8.598830961259482097e-02 4.164739999785817548e-01 4.437615850971945997e-01
                    8.457895032121595658e-02 4.201693888446985103e-01 4.455986826648858368e-01
                    8.315236408743081897e-02 4.238668393102046350e-01 4.474136915088433031e-01
                    8.171394242970148047e-02 4.275661399451892164e-01 4.492060096666721791e-01
                    8.026987997778653461e-02 4.312670600459782566e-01 4.509750376540915817e-01
                    7.882726258105521300e-02 4.349693489173894201e-01 4.527201798696426915e-01
                    7.739415915996109008e-02 4.386727351476746306e-01 4.544408461640829233e-01
                    7.597971511267428979e-02 4.423769258816079852e-01 4.561364535850265800e-01
                    7.459424408415230023e-02 4.460816060979165276e-01 4.578064283072901808e-01
                    7.324931366840245484e-02 4.497864378980456768e-01 4.594502077591564038e-01
                    7.195781915786156335e-02 4.534910598140922677e-01 4.610672429543411499e-01
                    7.073403782879061907e-02 4.571950861446216208e-01 4.626570010388615928e-01
                    6.959365457471067273e-02 4.608981063279830592e-01 4.642189680611767399e-01
                    6.855374817251533304e-02 4.645996843636931994e-01 4.657526519729214276e-01
                    6.763272639280798471e-02 4.682993582933911436e-01 4.672575858662252890e-01
                    6.685019795786503738e-02 4.719966397538147285e-01 4.687333314519904204e-01
                    6.622677049492370349e-02 4.756910136151853430e-01 4.701794827815656830e-01
                    6.577673858968982601e-02 4.793824370720338179e-01 4.715941633832017588e-01
                    6.552566064603559948e-02 4.830700650743400826e-01 4.729777438223161101e-01
                    6.549820807888259711e-02 4.867530881918184504e-01 4.743305399824718216e-01
                    6.571558023976076246e-02 4.904308747270616498e-01 4.756523204648347991e-01
                    6.619782080433814220e-02 4.941027658684106760e-01 4.769429100625798834e-01
                    6.696293828126623215e-02 4.977680989992983585e-01 4.782021118153178540e-01
                    6.801610562745827315e-02 5.014269919387287500e-01 4.794267081520149909e-01
                    6.938287907938911481e-02 5.050778255979795350e-01 4.806197953737798012e-01
                    7.107372065756567547e-02 5.087198306495862576e-01 4.817814896929000779e-01
                    7.309574065371191032e-02 5.123522156200341904e-01 4.829119931483520367e-01
                    7.544312997722565917e-02 5.159750173542567708e-01 4.840076729265158639e-01
                    7.812721320542403980e-02 5.195865216476734938e-01 4.850725749743695636e-01
                    8.114572394888117102e-02 5.231858403590214923e-01 4.861073630469811557e-01
                    8.448832972524200624e-02 5.267726392289812098e-01 4.871097239159773440e-01
                    8.815110303261089464e-02 5.303457163263245455e-01 4.880817131293228583e-01
                    9.212706115366336990e-02 5.339038791645212001e-01 4.890257976801760109e-01
                    9.640117372078826907e-02 5.374467824904106683e-01 4.899392661237967350e-01
                    1.009651127375016111e-01 5.409730397620680087e-01 4.908260102866368046e-01
                    1.058058868411170528e-01 5.444817357603629615e-01 4.916873206428047927e-01
                    1.109091719419828259e-01 5.479722539284174188e-01 4.925220286879308795e-01
                    1.162631754482068569e-01 5.514432364311973034e-01 4.933355557142476422e-01
                    1.218535592637135234e-01 5.548942058020517321e-01 4.941256571865651481e-01
                    1.276672732366210816e-01 5.583239554235875923e-01 4.948978636889357907e-01
                    1.336912106530599997e-01 5.617318570854707982e-01 4.956519172255148820e-01
                    1.399119531735926458e-01 5.651170059000043544e-01 4.963918396953513335e-01
                    1.463171276735238113e-01 5.684787369209851615e-01 4.971190242696297834e-01
                    1.528937099908951325e-01 5.718163393279495077e-01 4.978371217540351057e-01
                    1.596300036906863895e-01 5.751292412901443107e-01 4.985481743532478860e-01
                    1.665133003434038084e-01 5.784169089899582339e-01 4.992561455145328453e-01
                    1.735328783097749850e-01 5.816789154432452369e-01 4.999631315071642601e-01
                    1.806761031013389140e-01 5.849149331180006905e-01 5.006736479698321585e-01
                    1.879334378838529440e-01 5.881246912104457492e-01 5.013896313522919757e-01
                    1.952924636272905801e-01 5.913080959200207598e-01 5.021159025525618880e-01
                    2.027443124327554802e-01 5.944650539898058694e-01 5.028546409645925364e-01
                    2.102776619158149840e-01 5.975956842983582984e-01 5.036102222376415138e-01
                    2.178835011739876371e-01 6.007001216855215597e-01 5.043855354959408954e-01
                    2.255521696872271886e-01 6.037786523436265984e-01 5.051841613954050070e-01
                    2.332749143253212143e-01 6.068316286127126702e-01 5.060092633925864503e-01
                    2.410430130087069522e-01 6.098595200479223211e-01 5.068641518721540562e-01
                    2.488489441844971561e-01 6.128628192186393875e-01 5.077515718650474907e-01
                    2.566840197606024554e-01 6.158422311486778655e-01 5.086750547256451149e-01
                    2.645425291505308918e-01 6.187983086880044503e-01 5.096365792929340444e-01
                    2.724162323640093031e-01 6.217319406932720893e-01 5.106395912364674050e-01
                    2.803003746086680792e-01 6.246437759477415641e-01 5.116857725517640620e-01
                    2.881881143962810587e-01 6.275347641630057982e-01 5.127779478270290126e-01
                    2.960748324121407205e-01 6.304057046767445049e-01 5.139178767894679867e-01
                    3.039555214161791530e-01 6.332575133378499643e-01 5.151075596956975478e-01
                    3.118254758119217707e-01 6.360911463304157465e-01 5.163488722475468862e-01
                    3.196814435511296515e-01 6.389074409416231060e-01 5.176430585409529384e-01
                    3.275188103490194735e-01 6.417074699067215615e-01 5.189919719540163623e-01
                    3.353355552299603914e-01 6.444920106677175520e-01 5.203963563323180663e-01
                    3.431277348645550562e-01 6.472621506436356809e-01 5.218577640245062321e-01
                    3.508936232465371674e-01 6.500187029053001719e-01 5.233768303957324619e-01
                    3.586308928051362699e-01 6.527625972609284455e-01 5.249544352005411918e-01
                    3.663368289249827603e-01 6.554948605761221625e-01 5.265915743478537525e-01
                    3.740119682194762429e-01 6.582160193662682790e-01 5.282880419647097980e-01
                    3.816510562742471135e-01 6.609275747235534570e-01 5.300456830101799577e-01
                    3.892580517875900981e-01 6.636294950347043642e-01 5.318631054022345817e-01
                    3.968261890078971788e-01 6.663236061846555813e-01 5.337425455110953454e-01
                    4.043593499172888350e-01 6.690098785171258999e-01 5.356826761355167887e-01
                    4.118554481683104340e-01 6.716893377862352965e-01 5.376841277217632165e-01
                    4.193119538935562440e-01 6.743631277083362852e-01 5.397475699351562684e-01
                    4.267325988888523436e-01 6.770311866202418649e-01 5.418718349209995511e-01
                    4.341129901223799714e-01 6.796950319734580415e-01 5.440580596655568701e-01
                    4.414547339722453279e-01 6.823550079810395408e-01 5.463056567639269501e-01
                    4.487598229028739172e-01 6.850113432670585922e-01 5.486140038956078824e-01
                    4.560240745098563253e-01 6.876655632831355502e-01 5.509839689679422170e-01
                    4.632498634530364812e-01 6.903178139102768007e-01 5.534147775766539157e-01
                    4.704389351347622594e-01 6.929683367461523247e-01 5.559058726658481220e-01
                    4.775890398856753039e-01 6.956182556045346077e-01 5.584575291026864230e-01
                    4.846997162615552246e-01 6.982683108816961637e-01 5.610695620964348818e-01
                    4.917741270603506742e-01 7.009183842407572529e-01 5.637411222277634026e-01
                    4.988124547076915882e-01 7.035690215604037956e-01 5.664719533702552434e-01
                    5.058115298480653221e-01 7.062215879970826782e-01 5.692622613027222833e-01
                    5.127738759690677606e-01 7.088760718122417703e-01 5.721113127309581659e-01
                    5.197010257328686933e-01 7.115326679134365007e-01 5.750185876848199484e-01
                    5.265932673521412921e-01 7.141918646847531527e-01 5.779837442967735717e-01
                    5.334494816269396145e-01 7.168545129510710545e-01 5.810065526719623286e-01
                    5.402690647877386176e-01 7.195213413684249382e-01 5.840866494502513495e-01
                    5.470548736207613283e-01 7.221921468817737999e-01 5.872233804977959881e-01
                    5.538072551957783363e-01 7.248673679611095100e-01 5.904163410524770894e-01
                    5.605265661533928023e-01 7.275474322700724583e-01 5.936651120515976654e-01
                    5.672131712764029166e-01 7.302327570205698892e-01 5.969692609600730782e-01
                    5.738649923321349489e-01 7.329244536181049874e-01 6.003283620189617809e-01
                    5.804845702488933279e-01 7.356223182652517067e-01 6.037418650719156288e-01
                    5.870727421002424062e-01 7.383266077479672118e-01 6.072092971993915400e-01
                    5.936298933201957784e-01 7.410376988478269977e-01 6.107301851240671819e-01
                    6.001564124519274124e-01 7.437559595200373685e-01 6.143040456187923715e-01
                    6.066526905401021796e-01 7.464817491552263595e-01 6.179303860708607044e-01
                    6.131191206125891080e-01 7.492154188267443615e-01 6.216087050224498034e-01
                    6.195560972423146406e-01 7.519573115237726535e-01 6.253384926911926822e-01
                    6.259635861936422296e-01 7.547078996712067722e-01 6.291191822939783407e-01
                    6.323414409254812796e-01 7.574676913865331374e-01 6.329501570682624090e-01
                    6.386911927456534466e-01 7.602366470910791874e-01 6.368310037415121361e-01
                    6.450132383997043695e-01 7.630150778921287458e-01 6.407611911364364810e-01
                    6.513079749051275957e-01 7.658032878094622742e-01 6.447401822997557153e-01
                    6.575757994324034073e-01 7.686015739346664377e-01 6.487674349657708284e-01
                    6.638171092147175933e-01 7.714102265818907345e-01 6.528424020163275943e-01
                    6.700323014813579503e-01 7.742295294309987641e-01 6.569645319377305226e-01
                    6.762217734101813038e-01 7.770597596640953508e-01 6.611332692747456941e-01
                    6.823859220949484161e-01 7.799011880964143995e-01 6.653480550814426797e-01
                    6.885251445237092760e-01 7.827540793025653532e-01 6.696083273682165160e-01
                    6.946398375648411561e-01 7.856186917391158042e-01 6.739135215439332471e-01
                    7.007303979576740005e-01 7.884952778644843674e-01 6.782630708517513041e-01
                    7.067972223051004477e-01 7.913840842570615264e-01 6.826564067967675342e-01
                    7.128407070658809852e-01 7.942853517324697243e-01 6.870929595632513376e-01
                    7.188612485448663270e-01 7.971993154607720511e-01 6.915721584188375681e-01
                    7.248592428796667431e-01 8.001262050844080154e-01 6.960934321026465144e-01
                    7.308350860228653989e-01 8.030662448375227580e-01 7.006562091939169123e-01
                    7.367891737192681090e-01 8.060196536672460388e-01 7.052599184573171698e-01
                    7.427219014782563411e-01 8.089866453573647531e-01 7.099039891607071828e-01
                    7.486327256130643759e-01 8.119678017125729896e-01 7.145873998348200029e-01
                    7.545223076817887398e-01 8.149632306197714948e-01 7.193096566367255251e-01
                    7.603915796412942241e-01 8.179729213555224643e-01 7.240704295566047222e-01
                    7.662409381414175824e-01 8.209970689783633313e-01 7.288691437189614986e-01
                    7.720707796269933310e-01 8.240358639329949941e-01 7.337052258319404219e-01
                    7.778815003305534770e-01 8.270894921786368092e-01 7.385781042557458820e-01
                    7.836734962813408645e-01 8.301581353160575327e-01 7.434872090029538416e-01
                    7.894471633387081244e-01 8.332419707110174656e-01 7.484319716624078245e-01
                    7.952028972601549173e-01 8.363411716110775718e-01 7.534118252377455249e-01
                    8.009401766621723207e-01 8.394563042330350777e-01 7.584255820860429376e-01
                    8.066584706735069332e-01 8.425879497641356464e-01 7.634719632765348818e-01
                    8.123599093553812711e-01 8.457355276105180675e-01 7.685515296147656938e-01
                    8.180448941186928558e-01 8.488991950022984900e-01 7.736637090147417961e-01
                    8.237138279483190439e-01 8.520791051833870311e-01 7.788079287668441264e-01
                    8.293671160911614271e-01 8.552754073408689317e-01 7.839836146965878383e-01
                    8.350048677827189847e-01 8.584883831868801440e-01 7.891899443847590900e-01
                    8.406246724049463159e-01 8.617194786977677712e-01 7.944239257301034529e-01
                    8.462299834685677036e-01 8.649674562084326279e-01 7.996873827745283325e-01
                    8.518212329614019973e-01 8.682324447773152043e-01 8.049797193798620132e-01
                    8.573988637074910768e-01 8.715145671024431273e-01 8.103003320898716222e-01
                    8.629633327895804840e-01 8.748139384163003962e-01 8.156486084917724533e-01
                    8.685112937853189941e-01 8.781325055194139084e-01 8.210202040569095638e-01
                    8.740469282893844616e-01 8.814686228402575097e-01 8.264178662740626624e-01
                    8.795708463666609411e-01 8.848223352046364898e-01 8.318410153090299852e-01
                    8.850836254841629724e-01 8.881937049832523412e-01 8.372889893254339411e-01
                    8.905836610629007666e-01 8.915838848352710677e-01 8.427587078675609078e-01
                    8.960718329388820402e-01 8.949928217018373600e-01 8.482495190363266158e-01
                    9.015509263802756745e-01 8.984195016922670307e-01 8.537628451897162352e-01
                    9.070218374764418279e-01 9.018638497314315217e-01 8.592980228813623667e-01
                    9.124832972967743538e-01 9.053269113370198129e-01 8.648517190341429295e-01
                    9.179365336943248188e-01 9.088084827838653901e-01 8.704232606428425889e-01
                    9.233851804191220980e-01 9.123071100087158936e-01 8.760148319228351355e-01
                    9.288308811012930821e-01 9.158223472192428272e-01 8.816263177427966502e-01
                    9.342716413123769437e-01 9.193556495997508016e-01 8.872531405900644375e-01
                    9.397124373542273812e-01 9.229047657151164819e-01 8.928998260600508052e-01
                    9.451563753217823161e-01 9.264682865796181055e-01 8.985697730985615639e-01
                    9.506035559047821826e-01 9.300461763422669392e-01 9.042646579991682199e-01
                    9.560531046363628382e-01 9.336385357088663461e-01 9.099886803187530182e-01
                    9.615066616129493982e-01 9.372434217616608665e-01 9.157560566822336989e-01
                    9.669573847273637002e-01 9.408624710997687268e-01 9.215795307640300971e-01
                    9.723870692594612786e-01 9.445024597948724621e-01 9.274670258562995873e-01
                    9.777785730890226068e-01 9.481687534308861354e-01 9.334364948680430318e-01
                    9.831050718338244510e-01 9.518727670560837018e-01 9.394860306213083101e-01
                    9.883417388454437402e-01 9.556282921109976458e-01 9.455836323325411685e-01
                    9.934918422996558141e-01 9.594375624216472387e-01 9.516983192548315040e-01
                    9.985763296811461798e-01 9.632965417140263442e-01 9.577895036430327247e-01
                    9.942114721489739848e-01 9.649414783718816002e-01 9.591713509300946461e-01
                    9.916915526798163460e-01 9.600677293546330260e-01 9.527406681900515428e-01
                    9.892073759214962125e-01 9.552017644060696311e-01 9.462702365737246657e-01
                    9.867719407557972167e-01 9.503380654950176476e-01 9.397586228881678050e-01
                    9.843739071729306067e-01 9.454788135288768602e-01 9.332265558186634280e-01
                    9.820182926871906526e-01 9.406217084851765664e-01 9.266730991029579201e-01
                    9.797019478013845317e-01 9.357670195072623764e-01 9.201050706729160256e-01
                    9.774207980730996725e-01 9.309154126689619391e-01 9.135293740478817037e-01
                    9.751815609868391688e-01 9.260643851639969171e-01 9.069399044850446900e-01
                    9.729704910473376822e-01 9.212176417076594070e-01 9.003536220502811327e-01
                    9.708034919455280631e-01 9.163699401531895106e-01 8.937532865598888376e-01
                    9.686635749176782939e-01 9.115260705110270756e-01 8.871590238431372732e-01
                    9.665594650012595546e-01 9.066829842770773862e-01 8.805615104307024099e-01
                    9.644872784946666444e-01 9.018415117599095643e-01 8.739656672915365743e-01
                    9.624420742140847862e-01 8.970028582796596428e-01 8.673774626543144795e-01
                    9.604345584612016262e-01 8.921632903093835720e-01 8.607852860310053478e-01
                    9.584498054598472594e-01 8.873272041947983801e-01 8.542061390733795001e-01
                    9.564989754176022041e-01 8.824906944948268661e-01 8.476279025304513937e-01
                    9.545749598026603833e-01 8.776557062005478915e-01 8.410587414202306267e-01
                    9.526745457712032517e-01 8.728229738617154787e-01 8.345024009748316374e-01
                    9.508087846256844111e-01 8.679885222040847337e-01 8.279470532549256800e-01
                    9.489631025192423186e-01 8.631568397645421609e-01 8.214088147188137734e-01
                    9.471457599180185261e-01 8.583248540920686009e-01 8.148789280451533834e-01
                    9.453550305765402451e-01 8.534927996091392632e-01 8.083595091982623826e-01
                    9.435830323389927665e-01 8.486630445127123501e-01 8.018591800395635794e-01
                    9.418425662160879730e-01 8.438308651791852633e-01 7.953645386207840451e-01
                    9.401224124677813876e-01 8.389997911592477209e-01 7.888877269212305476e-01
                    9.384205871436177571e-01 8.341702205035793627e-01 7.824309818729748844e-01
                    9.367501034134483318e-01 8.293372224294222050e-01 7.759809608693934990e-01
                    9.350962812812430025e-01 8.245056601392387607e-01 7.695531958095181979e-01
                    9.334610667718560295e-01 8.196745424298020888e-01 7.631458380069808811e-01
                    9.318539060108851357e-01 8.148400979617762552e-01 7.567494905851217535e-01
                    9.302622555140654947e-01 8.100065567310004155e-01 7.503772218901773039e-01
                    9.286876056094349741e-01 8.051730724177085241e-01 7.440277117015546837e-01
                    9.271395831502474705e-01 8.003357056772891776e-01 7.376916134646933632e-01
                    9.256059253327618697e-01 7.954987025858116789e-01 7.313814866564718464e-01
                    9.240863437601837260e-01 7.906618776065421628e-01 7.250978248201844778e-01
                    9.225925658427436282e-01 7.858203843156221780e-01 7.188294627149656169e-01
                    9.211126169365579930e-01 7.809784744974170856e-01 7.125884635347393692e-01
                    9.196455684792974594e-01 7.761362045971830215e-01 7.063759865046427278e-01
                    9.181979022265810420e-01 7.712906810954774928e-01 7.001861662941214481e-01
                    9.167672001571998130e-01 7.664424894801199484e-01 6.940217357022171463e-01
                    9.153481567945904729e-01 7.615934195093969628e-01 6.878880649909780987e-01
                    9.139404639768551331e-01 7.567432836340363123e-01 6.817857667786046960e-01
                    9.125539294301500126e-01 7.518876928783747582e-01 6.757062380855892725e-01
                    9.111780957235272593e-01 7.470305845086899765e-01 6.696595752391263368e-01
                    9.098122501649290594e-01 7.421719339914916169e-01 6.636468136699643638e-01
                    9.084563128768237128e-01 7.373114504136536462e-01 6.576684342801734084e-01
                    9.071185731273107011e-01 7.324452193667679856e-01 6.517176237314421527e-01
                    9.057893354554222842e-01 7.275770233727173464e-01 6.458034937673092779e-01
                    9.044682125236334080e-01 7.227067015228565428e-01 6.399268508324190696e-01
                    9.031548045364521382e-01 7.178340969900336432e-01 6.340885324865254136e-01
                    9.018536281154535539e-01 7.129568373829493488e-01 6.282852974032566706e-01
                    9.005618093543485969e-01 7.080758191596225881e-01 6.225202235043866272e-01
                    8.992759701121831872e-01 7.031922144373208283e-01 6.167967177972875081e-01
                    8.979956379512032960e-01 6.983059000713781606e-01 6.111157531581823399e-01
                    8.967203212860314077e-01 6.934167616666734313e-01 6.054783383481375791e-01
                    8.954526708976771054e-01 6.885231926405954717e-01 5.998830635410449252e-01
                    8.941915566016080952e-01 6.836253365962329243e-01 5.943315847613384051e-01
                    8.929334172678807802e-01 6.787245284319090022e-01 5.888273795844850556e-01
                    8.916776705242244194e-01 6.738207055230612808e-01 5.833715948784413685e-01
                    8.904237095992360018e-01 6.689138192653601989e-01 5.779654133729836829e-01
                    8.891709022410747565e-01 6.640038362746348843e-01 5.726100532447696567e-01
                    8.879187220887899690e-01 6.590906724758603952e-01 5.673066752652388134e-01
                    8.866700510547862457e-01 6.541724962989774461e-01 5.620541403257021118e-01
                    8.854200814778028228e-01 6.492513861023391231e-01 5.568566631651752363e-01
                    8.841680786914359880e-01 6.443273831269403784e-01 5.517155748042765762e-01
                    8.829132805112128723e-01 6.394005490891804255e-01 5.466322335097351104e-01
                    8.816548969914349554e-01 6.344709672901969189e-01 5.416080224154796730e-01
                    8.803921103664992254e-01 6.295387436822825755e-01 5.366443467601184070e-01
                    8.791240751896489680e-01 6.246040078804957485e-01 5.317426307330340718e-01
                    8.778499186810901911e-01 6.196669141071349252e-01 5.269043139249947050e-01
                    8.765687412960717628e-01 6.147276420564155019e-01 5.221308473833028430e-01
                    8.752796175219744734e-01 6.097863976665317542e-01 5.174236892760738504e-01
                    8.739815969115529715e-01 6.048434137862945814e-01 5.127843001752147023e-01
                    8.726737941588253999e-01 5.998988989392898263e-01 5.082140946168248741e-01
                    8.713557496267209102e-01 5.949528226638209905e-01 5.037142772058214035e-01
                    8.700256052408116281e-01 5.900059904309469250e-01 4.992866995883939452e-01
                    8.686823445912698061e-01 5.850587401058766623e-01 4.949327569364067592e-01
                    8.673249354213553586e-01 5.801114372307474287e-01 4.906538189940125583e-01
                    8.659523319452930856e-01 5.751644748762787529e-01 4.864512237593011101e-01
                    8.645634773771747605e-01 5.702182733234352208e-01 4.823262709991545383e-01
                    8.631573066573615671e-01 5.652732795695735168e-01 4.782802156484906031e-01
                    8.617327493598900823e-01 5.603299666552966629e-01 4.743142611496265482e-01
                    8.602887327613949475e-01 5.553888328102728478e-01 4.704295527914130193e-01
                    8.588244189007790963e-01 5.504502458314576296e-01 4.666271134383851993e-01
                    8.573388771464994784e-01 5.455146532392860514e-01 4.629079425847138496e-01
                    8.558307110590109845e-01 5.405828402802138610e-01 4.592730689337449212e-01
                    8.542988567054966564e-01 5.356554029582069054e-01 4.557233480219839428e-01
                    8.527422657992969057e-01 5.307329553519699594e-01 4.522595444303549317e-01
                    8.511599092771463537e-01 5.258161276783985816e-01 4.488823265086418490e-01
                    8.495507808257551918e-01 5.209055642285482790e-01 4.455922615760174454e-01
                    8.479139003246871642e-01 5.160019211927755478e-01 4.423898116530999292e-01
                    8.462483171732270160e-01 5.111058643932878676e-01 4.392753297751230135e-01
                    8.445531134701587117e-01 5.062180669436533442e-01 4.362490569289804720e-01
                    8.428274070171406507e-01 5.013392068558483183e-01 4.333111196492664408e-01
                    8.410703541185753362e-01 4.964699646160714575e-01 4.304615283001137493e-01
                    8.392811521535089581e-01 4.916110207509353791e-01 4.277001760608312164e-01
                    8.374593644075595256e-01 4.867627987429819503e-01 4.250269048992098009e-01
                    8.356043139975503076e-01 4.819259312622626301e-01 4.224414243859783702e-01
                    8.337148868611683472e-01 4.771014353948369036e-01 4.199431953989596344e-01
                    8.317904658009995789e-01 4.722899701137897588e-01 4.175316375108732991e-01
                    8.298304844983901418e-01 4.674921821164159108e-01 4.152060584430907197e-01
                    8.278344282915189867e-01 4.627087037890093568e-01 4.129656573155296440e-01
                    8.258018346451376779e-01 4.579401513029328075e-01 4.108095284678619508e-01
                    8.237322933175019735e-01 4.531871228543432051e-01 4.087366658020525345e-01
                    8.216256878784896633e-01 4.484499819155622347e-01 4.067461112065531847e-01
                    8.194815690724320811e-01 4.437294051177347876e-01 4.048366204140803060e-01
                    8.172995050178583076e-01 4.390260854893064946e-01 4.030068099238640067e-01
                    8.150792968857805132e-01 4.343405321271661679e-01 4.012553165972070901e-01
                    8.128207950280410543e-01 4.296732294003547947e-01 3.995807052170325946e-01
                    8.105238975108269850e-01 4.250246362540385792e-01 3.979814747698801614e-01
                    8.081885484630557670e-01 4.203951856851556590e-01 3.964560648094507256e-01
                    8.058148470586988799e-01 4.157851737769109879e-01 3.950029765120905423e-01
                    8.034026962411985329e-01 4.111951046689115152e-01 3.936204324142628108e-01
                    8.009521378981909745e-01 4.066253627118514569e-01 3.923066947182650144e-01
                    7.984632874676697023e-01 4.020762706733001512e-01 3.910600267933756480e-01
                    7.959362958998934534e-01 3.975481251157440554e-01 3.898786646276527490e-01
                    7.933713473226599033e-01 3.930411968130258504e-01 3.887608225579017307e-01
                    7.907686574783574507e-01 3.885557303926693296e-01 3.877046999689159890e-01
                    7.881284544156301752e-01 3.840919640902587529e-01 3.867084595587381712e-01
                    7.854510012103226302e-01 3.796501032751679605e-01 3.857702673246294345e-01
                    7.827365958532629397e-01 3.752303177293505598e-01 3.848883050999601374e-01
                    7.799855573279415033e-01 3.708327552306592834e-01 3.840607581802352732e-01
                    7.771982233463431422e-01 3.664575425267851405e-01 3.832858192234025463e-01
                    7.743749481448563010e-01 3.621047863617611329e-01 3.825616918198675998e-01
                    7.715160518337467188e-01 3.577746368683014655e-01 3.818864831434901075e-01
                    7.686219524815265380e-01 3.534671195520193154e-01 3.812584982442814296e-01
                    7.656930517948008497e-01 3.491822750013692245e-01 3.806760122706670524e-01
                    7.627297518608332494e-01 3.449201381780978570e-01 3.801373059981723590e-01
                    7.597324610418284552e-01 3.406807299395016586e-01 3.796406843015277532e-01
                    7.567015922742318379e-01 3.364640582281425152e-01 3.791844777860022275e-01
                    7.536375467881040180e-01 3.322701411099439062e-01 3.787669970190922220e-01
                    7.505407196227874556e-01 3.280990002753491619e-01 3.783865433211569540e-01
                    7.474115712082280982e-01 3.239505493622616417e-01 3.780416634360889150e-01
                    7.442505183508412170e-01 3.198247550352249502e-01 3.777308063416947026e-01
                    7.410579750279425726e-01 3.157215772518807695e-01 3.774524509829353947e-01
                    7.378343512384677449e-01 3.116409704467773545e-01 3.772051064266770948e-01
                    7.345800519344474200e-01 3.075828847103421748e-01 3.769873118279194468e-01
                    7.312954622460766663e-01 3.035472925968438762e-01 3.767975677519167510e-01
                    7.279809779377566237e-01 2.995341347142527755e-01 3.766344783074295766e-01
                    7.246370003939071047e-01 2.955433232127043786e-01 3.764967560209065978e-01
                    7.212639047983969709e-01 2.915748020657008555e-01 3.763830604225979481e-01
                    7.178620569041518351e-01 2.876285169108695472e-01 3.762920784118068407e-01
                    7.144318123927204667e-01 2.837044162122266955e-01 3.762225233422813453e-01
                    7.109735162895217675e-01 2.798024524280844916e-01 3.761731340017399616e-01
                    7.074875003639726767e-01 2.759225885803616163e-01 3.761426556123988463e-01
                    7.039740902088845731e-01 2.720647802144139371e-01 3.761299015349165442e-01
                    7.004335983291459788e-01 2.682289908573779469e-01 3.761337091569335600e-01
                    6.968663212773178461e-01 2.644152012227610760e-01 3.761529095527339495e-01
                    6.932725426766019883e-01 2.606234021780075572e-01 3.761863526116585033e-01
                    6.896525329109274294e-01 2.568535960113472738e-01 3.762329055321002591e-01
                    6.860065488492232966e-01 2.531057977294109973e-01 3.762914512548127810e-01
                    6.823348363761195801e-01 2.493800214422073891e-01 3.763609500820121467e-01
                    6.786376242649662105e-01 2.456763038128263466e-01 3.764403535110733556e-01
                    6.749151231751110425e-01 2.419947210999891518e-01 3.765285264384131692e-01
                    6.711675324616499516e-01 2.383353539955441192e-01 3.766243994612780699e-01
                    6.673950369964636309e-01 2.346983031289386346e-01 3.767269100578239382e-01
                    6.635978070649097837e-01 2.310836906390836831e-01 3.768350007100279009e-01
                    6.597759982863891093e-01 2.274916618065601082e-01 3.769476169820946132e-01
                    6.559297473209481089e-01 2.239223513030393353e-01 3.770639116649144307e-01
                    6.520591754729223588e-01 2.203759654266244650e-01 3.771828077558652681e-01
                    6.481643941443696599e-01 2.168527541906234979e-01 3.773031206968150975e-01
                    6.442454984979625321e-01 2.133529749202845993e-01 3.774237883397635329e-01
                    6.403025689859566105e-01 2.098769174022531714e-01 3.775437410683403772e-01
                    6.363356714241054091e-01 2.064249059420060206e-01 3.776618996147099727e-01
                    6.323448551853009247e-01 2.029972985879547887e-01 3.777771986896097389e-01
                    6.283301122985771592e-01 1.995944416686266654e-01 3.778890664877021521e-01
                    6.242914965761068302e-01 1.962168403108488501e-01 3.779958773221168133e-01
                    6.202290151835444521e-01 1.928649814733044143e-01 3.780964890449534654e-01
                    6.161426615427991749e-01 1.895393981082279800e-01 3.781897379927290359e-01
                    6.120324156302930918e-01 1.862406715933664081e-01 3.782744366000001524e-01
                    6.078982443238682976e-01 1.829694341941491276e-01 3.783493709986866516e-01
                    6.037399990506459035e-01 1.797263486622488471e-01 3.784140428470601503e-01
                    5.995576965897967403e-01 1.765121975913181707e-01 3.784665259488561584e-01
                    5.953512720943007208e-01 1.733277858871770660e-01 3.785054621982753553e-01
                    5.911206437729070728e-01 1.701739770998706713e-01 3.785295035265708874e-01
                    5.868657195674668037e-01 1.670516976788867514e-01 3.785372630540221883e-01
                    5.825863435035250060e-01 1.639619511220119508e-01 3.785276019365447775e-01
                    5.782823203682633251e-01 1.609058321765788335e-01 3.784994295515334839e-01
                    5.739536466014196758e-01 1.578844538430557720e-01 3.784505787767488694e-01
                    5.696001981012460691e-01 1.548990153349074916e-01 3.783794775494400686e-01
                    5.652218458222221242e-01 1.519507861175737884e-01 3.782845039525286057e-01
                    5.608184279263467298e-01 1.490411220014642157e-01 3.781641038386275300e-01
                    5.563896048179566289e-01 1.461715605343842650e-01 3.780173288805740439e-01
                    5.519354537512982661e-01 1.433434802585686063e-01 3.778414852793369194e-01
                    5.474558435278772395e-01 1.405584279152288785e-01 3.776347269993366451e-01
                    5.429506470730138812e-01 1.378180135121435668e-01 3.773951523619960002e-01
                    5.384196501144510316e-01 1.351239813932839096e-01 3.771211030874267456e-01
                    5.338626135422905872e-01 1.324781762459926737e-01 3.768109124380045194e-01
                    5.292796498492590151e-01 1.298822067930652524e-01 3.764618033826110377e-01
                    5.246706783487381509e-01 1.273378756852094063e-01 3.760716446399895441e-01
                    5.200356363922209457e-01 1.248470120733735367e-01 3.756382551974452033e-01
                    5.153742367301422656e-01 1.224117140591406694e-01 3.751600258018771838e-01
                    5.106866581816573714e-01 1.200336468026920456e-01 3.746341300161435961e-01
                    5.059729862980435477e-01 1.177145617740232852e-01 3.740580772636119544e-01
                    5.012332745984564575e-01 1.154562389937715539e-01 3.734295174784395543e-01
                    4.964674957217509177e-01 1.132605442683873864e-01 3.727463195168757570e-01
                    4.916758190194626121e-01 1.111290906599490258e-01 3.720059818290880060e-01
                    4.868585609140311798e-01 1.090632523851999269e-01 3.712058226357082269e-01
                    4.820159589014604840e-01 1.070644243738244350e-01 3.703434438399521023e-01
                    4.771482837998479165e-01 1.051338742417683159e-01 3.694164970204523168e-01
                    4.722559296975721854e-01 1.032725993184717139e-01 3.684225499867860298e-01
                    4.673393879073711177e-01 1.014813305146687883e-01 3.673591735272957459e-01
                    4.623991454541829804e-01 9.976065343840129218e-02 3.662241176536292220e-01
                    4.574358041489324234e-01 9.811082191626963045e-02 3.650151416090168244e-01
                    4.524502157463570762e-01 9.653152558968838837e-02 3.637298817526759542e-01
                    4.474429556959144683e-01 9.502266265951053725e-02 3.623665563347658880e-01
                    4.424148033564082039e-01 9.358361785401031474e-02 3.609233190012756665e-01
                    4.373666093617918915e-01 9.221344371517425920e-02 3.593984619211612608e-01
                    4.323000061476657274e-01 9.090969474799345806e-02 3.577897503619592023e-01
                    4.272154852828886629e-01 8.967151542383772211e-02 3.560963875772367726e-01
                    4.221140984792258188e-01 8.849693359303587026e-02 3.543172112905775273e-01
                    4.169969758462759302e-01 8.738362154375492463e-02 3.524512375858371849e-01
                    4.118658301965831270e-01 8.632803919827600203e-02 3.504973239811765007e-01
                    4.067221669675408213e-01 8.532676500724681312e-02 3.484548346345710534e-01
                    4.015666852977432533e-01 8.437757766440442952e-02 3.463238990387739746e-01
                    3.964007130297593218e-01 8.347701754004596686e-02 3.441044235724556866e-01
                    3.912256097683191602e-01 8.262141527259461715e-02 3.417965388380971303e-01
                    3.860430941874281596e-01 8.180633277581816909e-02 3.394004639140563162e-01
                    3.808555339288969832e-01 8.102605665725876039e-02 3.369164884312597086e-01
                    3.756631594887540060e-01 8.027857713191724476e-02 3.343459517649038371e-01
                    3.704673590324250032e-01 7.955977986332424257e-02 3.316898655204278401e-01
                    3.652695090362573227e-01 7.886552703397090025e-02 3.289494343067977944e-01
                    3.600709660701238435e-01 7.819169210462800779e-02 3.261260438198446687e-01
                    3.548730588141333908e-01 7.753419260275734581e-02 3.232212473818363851e-01
                    3.496770804198256477e-01 7.688902041627435069e-02 3.202367511798585031e-01
                    3.444853557545813350e-01 7.625034119537774102e-02 3.171743822649770728e-01
                    3.392982604990861795e-01 7.561586038621423422e-02 3.140362188799208365e-01
                    3.341166082970875029e-01 7.498254313274954619e-02 3.108243324765390669e-01
                    3.289414559917514524e-01 7.434697155739058982e-02 3.075408588184366798e-01
                    3.237737954171283072e-01 7.370590341662106026e-02 3.041879989569576392e-01
                    3.186145498482620964e-01 7.305628171756148315e-02 3.007680026557595920e-01
                    3.134645711821229530e-01 7.239524117563511663e-02 2.972831523530131137e-01
                    3.083246378402036969e-01 7.172011173017686647e-02 2.937357478389135412e-01
                    3.031954533707552635e-01 7.102841937235576664e-02 2.901280917978547591e-01
                    2.980776457172774063e-01 7.031788456468926474e-02 2.864624763353166292e-01
                    2.929717671102379239e-01 6.958641854495320467e-02 2.827411705802333475e-01
                    2.878782945311796904e-01 6.883211781082057557e-02 2.789664094251954607e-01
                    2.827976306923836725e-01 6.805325707674123037e-02 2.751403834400568682e-01
                    2.777301054710668571e-01 6.724828098296936618e-02 2.712652299698696257e-01
                    2.726759777346169922e-01 6.641579481976991883e-02 2.673430254060571443e-01
                    2.676354374924288515e-01 6.555455450916861104e-02 2.633757786005420098e-01
                    2.626097591395918363e-01 6.466149844681601255e-02 2.593660533804174051e-01
                    2.575992689989206608e-01 6.373522644476575794e-02 2.553159580518031269e-01
                    2.526028775756383737e-01 6.277674824557366584e-02 2.512267303676423147e-01
                    2.476205211412589313e-01 6.178532806314660647e-02 2.471000835041496368e-01
                    2.426520784951842757e-01 6.076032925135082391e-02 2.429376426307544024e-01
                    2.376973739685813714e-01 5.970120319443159712e-02 2.387409437758167274e-01
                    2.327561803770203386e-01 5.860747847771281133e-02 2.345114333677875140e-01
                    2.278285229664566147e-01 5.747825339035907838e-02 2.302506881609430733e-01
                    2.229180662380075839e-01 5.630663136815370479e-02 2.259630047830259447e-01
                    2.180204522985965676e-01 5.509893160286010588e-02 2.216467626153181270e-01
                    2.131352588343927157e-01 5.385491895770589538e-02 2.173030574920574165e-01
                    2.082620235717315138e-01 5.257438772048273617e-02 2.129328977906284059e-01
                    2.034002463374002811e-01 5.125715285178860520e-02 2.085372063771265272e-01]; 
               case 'phase' 
                    cmap = [6.583083928922510708e-01 4.699391690315133929e-01 4.941288203988051381e-02
                    6.643374189373471017e-01 4.662019008569991407e-01 5.766473450402211792e-02
                    6.702086925052345157e-01 4.624801381219734719e-01 6.534560309537773559e-02
                    6.760429905334627287e-01 4.586983759956768103e-01 7.273174322210870790e-02
                    6.817522846284524984e-01 4.549140651836585669e-01 7.979261956680192003e-02
                    6.874028047282801923e-01 4.510841669914616436e-01 8.667102950867147659e-02
                    6.929504980948593129e-01 4.472389296506211198e-01 9.335868960148416273e-02
                    6.984261912087648128e-01 4.433576785927097474e-01 9.992839268327721736e-02
                    7.038122981036579739e-01 4.394532762349419586e-01 1.063871045640915336e-01
                    7.091206923190102041e-01 4.355176525107516405e-01 1.127717449252612913e-01
                    7.143452449012380745e-01 4.315557562030358230e-01 1.190934789448820086e-01
                    7.194928861674689813e-01 4.275627171694253437e-01 1.253760590955056708e-01
                    7.245561927479047259e-01 4.235446971269447580e-01 1.316232516304018385e-01
                    7.295489482895208821e-01 4.194909849233816046e-01 1.378630483492892522e-01
                    7.344517242247444733e-01 4.154177405103107734e-01 1.440803933719045082e-01
                    7.392949550641365608e-01 4.112997327023164562e-01 1.503221736968965161e-01
                    7.440383351241327547e-01 4.071715772870146410e-01 1.565433461845777419e-01
                    7.487369523694534790e-01 4.029851907815384382e-01 1.628228161295872112e-01
                    7.533231937521204236e-01 3.988010690198531272e-01 1.690756638056752914e-01
                    7.578808294584169492e-01 3.945424511336589335e-01 1.754217924072464518e-01
                    7.623326022156785564e-01 3.902809567197497165e-01 1.817591538530497208e-01
                    7.667320478995347521e-01 3.859654943537603744e-01 1.881681875043557384e-01
                    7.710524721820763983e-01 3.816214148216601210e-01 1.946153193576832252e-01
                    7.752952778457107286e-01 3.772473246523442292e-01 2.011065208945331251e-01
                    7.794866593781552000e-01 3.728150850100374614e-01 2.076872953652798559e-01
                    7.835853362619955575e-01 3.683677247727679127e-01 2.142973641992331202e-01
                    7.876376313731141554e-01 3.638539950288945946e-01 2.210164757476915653e-01
                    7.916113385801130109e-01 3.593080404226161040e-01 2.277974012097270795e-01
                    7.955060552511121763e-01 3.547298979176854994e-01 2.346435260251332755e-01
                    7.993539838133075781e-01 3.500795929482780067e-01 2.416183193481767355e-01
                    8.031167084780931331e-01 3.454015163918872089e-01 2.486589162182153978e-01
                    8.068103261859892461e-01 3.406745225101673324e-01 2.558007514596662979e-01
                    8.104452043001210138e-01 3.358824841311982556e-01 2.630722168071251699e-01
                    8.139968010634672790e-01 3.310553831960054150e-01 2.704318260711638389e-01
                    8.174768947816700715e-01 3.261752637674419919e-01 2.779109619017897659e-01
                    8.208941481936247175e-01 3.212262889577221503e-01 2.855384604486080891e-01
                    8.242271261613036692e-01 3.162361985850312696e-01 2.932761677775464482e-01
                    8.274766137520660481e-01 3.112015433453460544e-01 3.011338788484150264e-01
                    8.306639856730879679e-01 3.060845940504859919e-01 3.091757851828964565e-01
                    8.337630658888238733e-01 3.009224359536462057e-01 3.173492086215453090e-01
                    8.367728587568070697e-01 2.957134612590129330e-01 3.256619937531243236e-01
                    8.396969301670653696e-01 2.904472328682123905e-01 3.341366497709652439e-01
                    8.425387334318170662e-01 2.851115076990556330e-01 3.427996193346179443e-01
                    8.452829693062083871e-01 2.797291658336165110e-01 3.516207809735176215e-01
                    8.479270414576970394e-01 2.743004518168249417e-01 3.606068061808916925e-01
                    8.504679257105771661e-01 2.688262350691851821e-01 3.697639476988339724e-01
                    8.529105574510703613e-01 2.632885874038991547e-01 3.791311611530527870e-01
                    8.552420044048544279e-01 2.577088839643039142e-01 3.886821696092926381e-01
                    8.574567263211506640e-01 2.520936658861107627e-01 3.984160108673813205e-01
                    8.595502320723124035e-01 2.464473708549377862e-01 4.083362533435646036e-01
                    8.615176695665305306e-01 2.407756349505418836e-01 4.184455712413106543e-01
                    8.633539245616504987e-01 2.350852138633325872e-01 4.287460573240959860e-01
                    8.650568452189218993e-01 2.293728768111193417e-01 4.392600781097817930e-01
                    8.666160631462880293e-01 2.236630759123471868e-01 4.499612660389528673e-01
                    8.680257787196800079e-01 2.179678524813165597e-01 4.608475800597766070e-01
                    8.692800281403405549e-01 2.123013227323138907e-01 4.719155444862356830e-01
                    8.703727414475561641e-01 2.066798807505955682e-01 4.831601540690036445e-01
                    8.712978068076887572e-01 2.011223984119477337e-01 4.945747939602412324e-01
                    8.720491402648608004e-01 1.956504128743852822e-01 5.061511790861612514e-01
                    8.726207597696474805e-01 1.902882889012951495e-01 5.178793172539123413e-01
                    8.730068619526127893e-01 1.850633392111822872e-01 5.297474998646168887e-01
                    8.732018998073337590e-01 1.800058813906028066e-01 5.417423233741458510e-01
                    8.732006592202997686e-01 1.751492049740893675e-01 5.538487436526324803e-01
                    8.729983321570818910e-01 1.705294177324829519e-01 5.660501641838993070e-01
                    8.725905843033689990e-01 1.661851370964102237e-01 5.783285576821829421e-01
                    8.719736126526835829e-01 1.621569796196943858e-01 5.906646606873059424e-01
                    8.711441430703839028e-01 1.584866680756561452e-01 6.030388085538075371e-01
                    8.700996606911847175e-01 1.552168662477456385e-01 6.154284378630663355e-01
                    8.688382315166471859e-01 1.523889228568965915e-01 6.278117548576560569e-01
                    8.673585803491290491e-01 1.500419870434277214e-01 6.401665059420836856e-01
                    8.656600953609360216e-01 1.482114935188816873e-01 6.524702171100530412e-01
                    8.637428203708292784e-01 1.469276207331668138e-01 6.647004334019646077e-01
                    8.616074355850228406e-01 1.462138565700954185e-01 6.768349518289993316e-01
                    8.592552279602999610e-01 1.460858181738819150e-01 6.888520417478360969e-01
                    8.566880526588371847e-01 1.465504601378146143e-01 7.007306474947510022e-01
                    8.539082940810565070e-01 1.476057634824317899e-01 7.124505416663439172e-01
                    8.509188132430276497e-01 1.492409439214724132e-01 7.239924974090546916e-01
                    8.477228732160004832e-01 1.514371672020923265e-01 7.353384896921373315e-01
                    8.443240935674275471e-01 1.541686506679828816e-01 7.464717393223647690e-01
                    8.407263939629991967e-01 1.574040314553000475e-01 7.573767829039294019e-01
                    8.369339386719339968e-01 1.611078555520196187e-01 7.680395184197363889e-01
                    8.329510832270686782e-01 1.652420527841307607e-01 7.784472273410004695e-01
                    8.287823240288615390e-01 1.697672910052525075e-01 7.885885755886976600e-01
                    8.244322514509230260e-01 1.746441391323516057e-01 7.984535959358352031e-01
                    8.199055067926670493e-01 1.798340046929527702e-01 8.080336545530442116e-01
                    8.152067432395495583e-01 1.852998414406232253e-01 8.173214043862717659e-01
                    8.103405908369269994e-01 1.910066437572490727e-01 8.263107279415592421e-01
                    8.053117579051806141e-01 1.969216010609510237e-01 8.349964504885767358e-01
                    8.001246695316600599e-01 2.030146524204510805e-01 8.433748621071933682e-01
                    7.947836730093801316e-01 2.092582614636481764e-01 8.514431956464665330e-01
                    7.892930180385833161e-01 2.156273737159832282e-01 8.591995655348122485e-01
                    7.836568069717744223e-01 2.220993550443281506e-01 8.666429444356050782e-01
                    7.778789796122015376e-01 2.286538551089048465e-01 8.737730828797295457e-01
                    7.719633002465959848e-01 2.352726496465096795e-01 8.805904299224873721e-01
                    7.659133466040326521e-01 2.419394735315277267e-01 8.870960555269299386e-01
                    7.597325004651533931e-01 2.486398530394249295e-01 8.932915751926430170e-01
                    7.534239396834068181e-01 2.553609429642716422e-01 8.991790771963389384e-01
                    7.469906314200031039e-01 2.620913721366889826e-01 9.047610526867282399e-01
                    7.404353264352245834e-01 2.688210993418901351e-01 9.100403287788533246e-01
                    7.337605543192734503e-01 2.755412805372748908e-01 9.150200047193097763e-01
                    7.269686195850667554e-01 2.822441475132592692e-01 9.197033911407090923e-01
                    7.200615985827201193e-01 2.889228976436723495e-01 9.240939523881464002e-01
                    7.130413372306516617e-01 2.955715940634778272e-01 9.281952518796100504e-01
                    7.059094495911268918e-01 3.021850754377659598e-01 9.320109004535763741e-01
                    6.986673173487639721e-01 3.087588744057577217e-01 9.355445076582962205e-01
                    6.913160902792203633e-01 3.152891437666022756e-01 9.387996359465969887e-01
                    6.838566878221142842e-01 3.217725894980000279e-01 9.417797577558794098e-01
                    6.762898018976056802e-01 3.282064097483991527e-01 9.444882154741288671e-01
                    6.686159011300095711e-01 3.345882390078672719e-01 9.469281843183131597e-01
                    6.608352366648267973e-01 3.409160967342457216e-01 9.491026381806494383e-01
                    6.529478497874137144e-01 3.471883397852293385e-01 9.510143185305323099e-01
                    6.449535815726188392e-01 3.534036180803885596e-01 9.526657064947888776e-01
                    6.368520848146350666e-01 3.595608329884161236e-01 9.540589982761437104e-01
                    6.286428385050449874e-01 3.656590980031973470e-01 9.551960841088614762e-01
                    6.203251651440072623e-01 3.716977013375789007e-01 9.560785309910512231e-01
                    6.118982511841425387e-01 3.776760701263490727e-01 9.567075694743777392e-01
                    6.033611709179768079e-01 3.835937359904362243e-01 9.570840848332117234e-01
                    5.947129141266446206e-01 3.894503017735960748e-01 9.572086129751365968e-01
                    5.859524178083830304e-01 3.952454093215529429e-01 9.570813414919564499e-01
                    5.770786022984990549e-01 4.009787082326203289e-01 9.567021162826109260e-01
                    5.680904120756430364e-01 4.066498255689147689e-01 9.560704542043091392e-01
                    5.589868615201727398e-01 4.122583365789312393e-01 9.551855622225308151e-01
                    5.497670858464057675e-01 4.178037365457583641e-01 9.540463635305430623e-01
                    5.404303973688802110e-01 4.232854139404694238e-01 9.526515310905910860e-01
                    5.309763471805271084e-01 4.287026251268353794e-01 9.509995290068892215e-01
                    5.214047922154096959e-01 4.340544709302417981e-01 9.490886620701871612e-01
                    5.117159675380545947e-01 4.393398754490472347e-01 9.469171337096092822e-01
                    5.019105635439048418e-01 4.445575675481795996e-01 9.444831124449268867e-01
                    4.919898075708288299e-01 4.497060655293146358e-01 9.417848067474301477e-01
                    4.819555492105473404e-01 4.547836655159373520e-01 9.388205479873396042e-01
                    4.718103483745911819e-01 4.597884341201938230e-01 9.355888808698555881e-01
                    4.615575649166388517e-01 4.647182059669867082e-01 9.320886604427106592e-01
                    4.511980082251949020e-01 4.695721776839077433e-01 9.283178636903177683e-01
                    4.407385246316099514e-01 4.743468819302843476e-01 9.242766906682317041e-01
                    4.301872191964954406e-01 4.790386406828319177e-01 9.199661970019400448e-01
                    4.195516590423077341e-01 4.836443962749994996e-01 9.153875906397229700e-01
                    4.088406331006777528e-01 4.881609386110861704e-01 9.105429315654867128e-01
                    3.980642083018061106e-01 4.925849423947825101e-01 9.054352272162426996e-01
                    3.872337717548847147e-01 4.969130108713233906e-01 9.000685170373873278e-01
                    3.763620564988861550e-01 5.011417254753421924e-01 8.944479427537497251e-01
                    3.654612663762428770e-01 5.052684034325813922e-01 8.885787674719856089e-01
                    3.545465439790672080e-01 5.092897968880314430e-01 8.824681923036010733e-01
                    3.436377931849474154e-01 5.132015795889635079e-01 8.761266360928694485e-01
                    3.327530939137686161e-01 5.170008198204416594e-01 8.695640864822369309e-01
                    3.219116583925260011e-01 5.206848733578030020e-01 8.627916585616799416e-01
                    3.111337198622898259e-01 5.242514414427805747e-01 8.558215206665023000e-01
                    3.004403986093719392e-01 5.276986238157412856e-01 8.486667948769511804e-01
                    2.898532566310526581e-01 5.310250549915590534e-01 8.413412305522999235e-01
                    2.793961594244663282e-01 5.342293148180886631e-01 8.338605094653681604e-01
                    2.690918119759744820e-01 5.373109871573702456e-01 8.262398418700288572e-01
                    2.589630009926549015e-01 5.402702010433201307e-01 8.184947475296583397e-01
                    2.490323931623479869e-01 5.431076290262466522e-01 8.106408952154365855e-01
                    2.393222899575573326e-01 5.458244839403564308e-01 8.026939223701164972e-01
                    2.298566351816503373e-01 5.484218890000336355e-01 7.946712158983202379e-01
                    2.206551002944838191e-01 5.509024141956824216e-01 7.865870626923798792e-01
                    2.117364094415158937e-01 5.532690129113547739e-01 7.784553296976853831e-01
                    2.031184306485467883e-01 5.555248911915373622e-01 7.702897347110843063e-01
                    1.948172015035146420e-01 5.576736467032974431e-01 7.621031776295352778e-01
                    1.868465966290397962e-01 5.597192229650957973e-01 7.539076339343501187e-01
                    1.792179884783634825e-01 5.616658613090964591e-01 7.457140662049784874e-01
                    1.719421959020442647e-01 5.635174702399022850e-01 7.375349843171632447e-01
                    1.650229512673002386e-01 5.652791533635259658e-01 7.293775439873284583e-01
                    1.584611602813638387e-01 5.669560017699032395e-01 7.212481886261494779e-01
                    1.522549918213680353e-01 5.685529672520610589e-01 7.131532051881548373e-01
                    1.463987618681510117e-01 5.700750569206646245e-01 7.050976866560473288e-01
                    1.408828405908224835e-01 5.715272914520931336e-01 6.970855353349069139e-01
                    1.356936625367466953e-01 5.729146658465781305e-01 6.891194758268863740e-01
                    1.308138525417579801e-01 5.742421127377883572e-01 6.812010763134157543e-01
                    1.262224734293325157e-01 5.755144682058210837e-01 6.733307769330967307e-01
                    1.218953938035638729e-01 5.767364399826563348e-01 6.655079242329828837e-01
                    1.178065439715534068e-01 5.779123523036159282e-01 6.577323260140637284e-01
                    1.139261298009009160e-01 5.790467957802227783e-01 6.499998366451625875e-01
                    1.102234782385779766e-01 5.801439808367224726e-01 6.423063681141680803e-01
                    1.066673235183970281e-01 5.812078215483343913e-01 6.346473314887408623e-01
                    1.032263067644274834e-01 5.822419757773840132e-01 6.270172879951815270e-01
                    9.986970116659238395e-02 5.832498248631478033e-01 6.194100074982554771e-01
                    9.656813269866712512e-02 5.842344549012296051e-01 6.118185305845139643e-01
                    9.329429206118067253e-02 5.851986396650793454e-01 6.042352341768219004e-01
                    9.002364276931534848e-02 5.861448252640329981e-01 5.966519005456781821e-01
                    8.673514013646366205e-02 5.870751166646858143e-01 5.890597894857422245e-01
                    8.341198598223345528e-02 5.879912662222103181e-01 5.814497133016774955e-01
                    8.004245400348602990e-02 5.888946643728962815e-01 5.738121141059205899e-01
                    7.662083060590932360e-02 5.897863326271246542e-01 5.661371427838456372e-01
                    7.314852485397654869e-02 5.906669189727291602e-01 5.584147388414824054e-01
                    6.963540714651897390e-02 5.915366957529474279e-01 5.506347102299997687e-01
                    6.610143501147561218e-02 5.923955600227664986e-01 5.427868121514007882e-01
                    6.257860760091599195e-02 5.932430363153763375e-01 5.348608238014407323e-01
                    5.911303759975024968e-02 5.940783314685343930e-01 5.268461376019544229e-01
                    5.576765010285392871e-02 5.949002975004614724e-01 5.187321979872163702e-01
                    5.262510565424809855e-02 5.957073200832398996e-01 5.105097808257538228e-01
                    4.978880680939656161e-02 5.964975038111670624e-01 5.021693595319391967e-01
                    4.738319269394732774e-02 5.972686215557058143e-01 4.937017361708459506e-01
                    4.555066847662456869e-02 5.980181252952141424e-01 4.850980867717663014e-01
                    4.444396189302631667e-02 5.987431566697054564e-01 4.763500004038906943e-01
                    4.421322948900235222e-02 5.994405566877633040e-01 4.674495124401222834e-01
                    4.498917867710291313e-02 6.001068740129568146e-01 4.583891327826852824e-01
                    4.686604485069371245e-02 6.007383712762072170e-01 4.491618701849053319e-01
                    4.988979024235475762e-02 6.013310288983726437e-01 4.397612541729841729e-01
                    5.405573006313083712e-02 6.018805543938667846e-01 4.301812015738862294e-01
                    5.932208540209162745e-02 6.023828876312081748e-01 4.204054325617814780e-01
                    6.560773880422715587e-02 6.028325831369402144e-01 4.104377156032376628e-01
                    7.281962363892094392e-02 6.032244155970721833e-01 4.002736262351211383e-01
                    8.086176781346332554e-02 6.035528335053932381e-01 3.899094149354802585e-01
                    8.964365767756551917e-02 6.038119404209635332e-01 3.793420777858104165e-01
                    9.908952486975769469e-02 6.039955435016706176e-01 3.685641184021282157e-01
                    1.091461686387718844e-01 6.040969460876761676e-01 3.575579876882719055e-01
                    1.197411868235576382e-01 6.041085844094219448e-01 3.463409605502589250e-01
                    1.308274635712054768e-01 6.040228036354146068e-01 3.349141605174003056e-01
                    1.423800347596975990e-01 6.038311914501962585e-01 3.232669997463486489e-01
                    1.543847008629261053e-01 6.035242532835117801e-01 3.113882295486853913e-01
                    1.667909278439370091e-01 6.030930063069259717e-01 2.993102878563658198e-01
                    1.795975741928860503e-01 6.025266800569213377e-01 2.870237044475537069e-01
                    1.927996570161342460e-01 6.018136381688595771e-01 2.745296424849442141e-01
                    2.063446461864067438e-01 6.009446614046377588e-01 2.618794002040416569e-01
                    2.202728729725795254e-01 5.999042993738278318e-01 2.490425136645111337e-01
                    2.344983329249494819e-01 5.986859114745567423e-01 2.361102156981361166e-01
                    2.490441577314199684e-01 5.972745984023483112e-01 2.230778046167814499e-01
                    2.638200588853058526e-01 5.956665601719514092e-01 2.100467332369221340e-01
                    2.788103972435367339e-01 5.938520960352462463e-01 1.970548426870546432e-01
                    2.939149436021079032e-01 5.918334783865150106e-01 1.842162055743405413e-01
                    3.090633958032663053e-01 5.896130197890373514e-01 1.716194213262129398e-01
                    3.241557701146778325e-01 5.872013194176191053e-01 1.593775348465770181e-01
                    3.391058987925449908e-01 5.846116417760637285e-01 1.475901239657965436e-01
                    3.537962403050050608e-01 5.818679304916246631e-01 1.363773405186440579e-01
                    3.681790539934771678e-01 5.789861015095073560e-01 1.258005415905540103e-01
                    3.821596573789931561e-01 5.759951191731197406e-01 1.159503977992311363e-01
                    3.957282445230750900e-01 5.729092809436007183e-01 1.068503820674790161e-01
                    4.088192648913888116e-01 5.697572697905517458e-01 9.855521202608327758e-02
                    4.214810649035653500e-01 5.665415862717396722e-01 9.104002246571415990e-02
                    4.336495334558000403e-01 5.632929561600926727e-01 8.434116239111955071e-02
                    4.453890783516751273e-01 5.600085906382493706e-01 7.841305443909030171e-02
                    4.567242130853677584e-01 5.566942951475172263e-01 7.322913012508130981e-02
                    4.676501707681009479e-01 5.533637269053324204e-01 6.876762134759795142e-02
                    4.781913781821610088e-01 5.500213008191600084e-01 6.498435721988443659e-02
                    4.883968612384750885e-01 5.466619506976900800e-01 6.182162837413415074e-02
                    4.982892398787829857e-01 5.432873953595327432e-01 5.922725775013665955e-02
                    5.078911374422868663e-01 5.398982668933685058e-01 5.714465998451528223e-02
                    5.172247456455653092e-01 5.364942872197172585e-01 5.551476288204022086e-02
                    5.263115025893324583e-01 5.330744300746458331e-01 5.427792556635833293e-02
                    5.351718645250624906e-01 5.296370658817842747e-01 5.337566853443905662e-02
                    5.438251536372012973e-01 5.261800887043684982e-01 5.275207547378853862e-02
                    5.522894663786476199e-01 5.227010255649303661e-01 5.235479075452074277e-02
                    5.605816294780567866e-01 5.191971290183284848e-01 5.213559944466891055e-02
                    5.687171933168381210e-01 5.156654540879882509e-01 5.205062439450063722e-02
                    5.767104547724090091e-01 5.121029206313181259e-01 5.206020204364426168e-02
                    5.845745037500035268e-01 5.085063619725217476e-01 5.212850674796504907e-02
                    5.923212894447653643e-01 5.048725602928031408e-01 5.222298774900729218e-02
                    5.999617038945115333e-01 5.011982688438658684e-01 5.231366903560310394e-02
                    6.075056816322426112e-01 4.974802205790658793e-01 5.237234455720112675e-02
                    6.149623152735116394e-01 4.937151222923196747e-01 5.237168182411356537e-02
                    6.223399877334323538e-01 4.898996328220020513e-01 5.228422613133864444e-02
                    6.296465225194588511e-01 4.860303233135512824e-01 5.208127381405584094e-02
                    6.368893542334475022e-01 4.821036169426133333e-01 5.173155232549151578e-02
                    6.440757220432393737e-01 4.781157049081125598e-01 5.119960076823739520e-02
                    6.512128893528150719e-01 4.740624350126631525e-01 5.044367478760234530e-02
                    6.583083928921535932e-01 4.699391690315524728e-01 4.941288204103298082e-02];
                otherwise
                    error('unknown colormap');
            end
            
            if nargin > 1
                % interpolating function
                map_ = @(n,cmap) interp1([0:size(cmap,1)-1]./(size(cmap,1)-1),cmap,linspace(0,1,n));
                % interpolate and return
                cmap = map_(n,cmap);
            end
            
            if nargout == 0; colormap(cmap); end
        end
        
        function [cmap] = clight_(cmap,alpha)
            brighten_ = @(x,alpha) x.*(1-alpha) + alpha;
            dim_      = @(x,alpha) x.*(1-alpha);
            if numel(alpha)==1; alpha = alpha*ones(size(cmap)); end
            ex_ = alpha(:)>=0; cmap(ex_) = brighten_(cmap(ex_),    alpha(ex_));
            ex_ = alpha(:)< 0; cmap(ex_) =      dim_(cmap(ex_),abs(alpha(ex_)));
        end
        
        function rgb = cmplx_color_(f,palette) % maps abs -> intensity, phase -> hue (palette)
            % [x,y]=meshgrid(linspace(-2,2,1000),linspace(-2,2,1000).*1i); z = x+y;
            % clist = am_lib.colormap_('red2blue',5);
            % clist = am_lib.rgb2hsl_(clist); clist(:,3) = 0.5; clist=am_lib.hsl2rgb_(clist);
            % C = am_lib.cmplx_color_(z,clist);
            % PP = surf(real(z),imag(z),abs(z),C,'edgecolor','none'); view([0 0 1]);

            % parse input
            if nargin<2; palette = am_lib.colormap('hsv',600); end
            % number of colors in palette encoding phase
            [m,n]=size(f); p = size(palette); p = p(1);
            % encode phase
            nphase = stepfct((angle(-f)+pi)/(2*pi),p);
            rgb(:,:,1) =reshape(palette(nphase,1),m,n);
            rgb(:,:,2) =reshape(palette(nphase,2),m,n);
            rgb(:,:,3) =reshape(palette(nphase,3),m,n);
            % encode modulous
%             bright = (log(abs(f))>=0).*((1-(1./(1+log(abs(f))))).^2)-(log(abs(f))<0).*(1-(1./(1-log(abs(f))))).^2;
%             bright = -((3.2/pi)*atan(abs(f))-0.2);
            bright = (3.2/pi)*atan(abs(f))-0.8;
            rgb = brightenRGB(rgb,bright);
            % deal with singularities
            rgb(isnan(rgb))=0.8;
            rgb(:,:,1) = rgb(:,:,1).*(abs(f)>0)+(1-rgb(:,:,1)).*(f==Inf);
            rgb(:,:,2) = rgb(:,:,2).*(abs(f)>0)+(1-rgb(:,:,2)).*(f==Inf);
            rgb(:,:,3) = rgb(:,:,3).*(abs(f)>0)+(1-rgb(:,:,3)).*(f==Inf);
            % reduce pure white and black to avoid printing problems
            rgb = 0.001 + 0.998*rgb;
            function y = stepfct(x,nmax)
                x = x-floor(x); y = floor(nmax*x)+1; y = int16(y); y = y+int16(y==0);
            end
            function RGB = brightenRGB(RGB,bright)
                if size(size(bright))==1; bright = bright * ones(size(RGB(:,:,1))); end
                RGB(:,:,1) = (bright>=0).* ((1-bright).*RGB(:,:,1) + bright.*ones(size(RGB(:,:,1)))) + (bright<0).*((1+bright).*RGB(:,:,1));
                RGB(:,:,2) = (bright>=0).* ((1-bright).*RGB(:,:,2) + bright.*ones(size(RGB(:,:,2)))) + (bright<0).*((1+bright).*RGB(:,:,2));
                RGB(:,:,3) = (bright>=0).* ((1-bright).*RGB(:,:,3) + bright.*ones(size(RGB(:,:,3)))) + (bright<0).*((1+bright).*RGB(:,:,3));
            end
        end
        
        function [th] = assign_cmap_(V)
            % assigns a number between [0,1] based on how close vectors V
            % in are to identity vectors in that same basis. It essentially
            % reduces the dimensionality of the data

            % set number of categories
            n = size(V,1);
            % get number of points
            m = size(V,2);

            if     n == 1
                c = ones(n,m);
            elseif n == 2
                p = [0 1]; 
                th = (p*abs(V)).';
            else
                % define categories around circle on complex plane
                p = exp(2*pi*1i*[1:n]/n);
                % get complex points
                c = (p*abs(V)).'; c = atan2d(imag(c),real(c));
                % 
                th = mod(mod(c/360,1)-0.5/n,1);
            end
        end
        
        
        % images
        
        function I = laplacian_interpolation_(I,mask)
            % laplacian interpolation of the masked region
            mask = mask==1;

            u = find(mask);
            w = find(~mask);

            M = size(mask,1);
            u_north = u - 1;
            u_east = u + M;
            u_south = u + 1;
            u_west = u - M;

            v = ones(size(u));
            ijv_mask = [...
                u  u         1.00*v
                u  u_north  -0.25*v
                u  u_east   -0.25*v
                u  u_south  -0.25*v
                u  u_west   -0.25*v ];

            ijv_nonmask = [w  w  1.00*ones(size(w))];

            ijv = [ijv_mask; ijv_nonmask];
            A = sparse(ijv(:,1),ijv(:,2),ijv(:,3));

            b = I(:);
            b(mask(:)) = 0;

            x = A\b;

            I = reshape(x,size(I));
        end
        

        % correlation and polynomial fitting
        
        function       plotcorr_(x,y)
            % plot correlation x vs y
            
            import am_lib.*
            
            % plot correlation for dft vs bvk forces on atoms
            hist2_(x,y); hold on; daspect([1 1 1]); box on;
            title(sprintf('R^2 = %f',pcorr_(x(:),y(:)).^2)); 
            
            % linear regression
            mv = minmax_([x(:);y(:)]);
            line( 0, 0, 'marker','.','markersize',20,'color','k');
            line( mv, mv, 'linewidth',1.5); 
            line( mv, pinterp_(mv,x(:),y(:),1) ,'linestyle','--','linewidth',1.5,'color','r'); 
            axis([mv,mv]);
            grid on;
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
        
        function [y] = peval_(C,x)
            n = numel(C)-1;
            % evaluate least squares fit polynomial of degree n to data
            y = (x(:).^[0:n])*C(:);
            y = reshape(y,size(x));
        end
        
        function [y] = pinterp_(x,x0,y0,n)
            % interpolate x0 vs y0 data at points x using a polynomial of degree n
            
            import am_lib.*
            
            y = peval_(x(:),pfit_(x0(:),y0(:),n));
        end

        function [C] = planefit_(varargin)
            % fit a plane to th date
            if     nargin == 3
                xx = x(:); yy = y(:); zz = z(:);
            elseif nargin == 1
                xx = varargin{1}(1,:); yy = varargin{1}(2,:);  zz = varargin{1}(3,:);
            end
            N = length(xx); O = ones(1,N); C = [xx;yy;O].'\zz.';
        end

        function Z   = plane_(X,Y,varargin)
            if numel(varargin) == 1
                C = varargin{1}(1:3);
            elseif numel(varargin) == 3
                [C(1),C(2),C(3)]=deal(varargin{:});
            else
                error('invalid input');
            end
            Z = X * C(1) + Y*C(2) + C(3);
%             normal = cross(p1 - p2, p1 - p3);
%             d = p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3); d = -d;
%             Z = (-d - (normal(1)*x) - (normal(2)*y))/normal(3);
        end


        % integral transforms related

        function [h]    = hilbert_(f) % hilbert transform (Kramers-Kronig transform)
            % Hilbert Transform is a 90 degree phase shift. Time domain to time domain.
            % t = [0:(N-1)]'/N; dt = t(2)-t(1); % time signal
            if any(abs(imag(f))>1E-8); error('hilbert only works on real-valued signals'); end
            N = numel(f); v = [0:(N-1)]'-floor(N/2);
            h = ifft(fftshift(sign(v(:))).*fft(f(:)))+f(:);
            h = reshape(h,size(f));
        end
        
        function [G,f,t]= stft_(t,y,w_,flag) % short-time fourier transform, y = signal, t = time vector, w_ = window function
            % options for windows:
            % 
            %	*) w = 0.03; w_ = @(x) am_lib.gauss_(x./w)*w; % gaussian window for gabor
            %   *) 
            % 
            % example:
            % 
            % N = 1001; t = [0:N-1]/N;
            % y = chirp(t,100,1,400)+chirp(t,-100,1,200); % define signal
            % w = 0.03; w_ = @(x) am_lib.gauss_(x./w)*w; % define window
            % am_lib.stft_(t,y,w);
            
% TO DO: ADD stft for 2D
% TO DO: ADD stft for 2D
% TO DO: ADD stft for 2D
%
% [wn,wm]=deal(21,19);
% w = am_lib.hannw_(wn).*am_lib.hannw_(wm).';
% 
% %%
% [n,m]=size(d); A = zeros(wn,wm,n-wn,m-wm,'single');
% for in = 1:(n-wn)
% for im = 1:(m-wm)
%     A(:,:,in,im) = fftshift(fftn(d([1:wn]+in-1,[1:wm]+im-1).*w));
% end
% end
% %%
% idx = kmeans(reshape(abs(A).^2,wn*wm,in*im).',2);
% %%
% am_lib.imagesc_(reshape(idx,in,im))

            if nargin<4; flag=''; end
            if nargout==0
                if contains(flag,'detail')
                    % define plotting grid
                    a=6; b=6; g=reshape(1:a*b,b,a).';
                    % plot curve
                    subplot(a,b,g(1,1:end-1)); plot(t,y); axis tight; set(gca,'xtick',[]);
                    % plot periodogram
                    [f,yf] = am_lib.fft_(t,y);
                    subplot(a,b,g(2:end,end)); loglog(abs(yf).^2,(1:numel(f))); axis tight; set(gca,'ytick',[]);
                    % plot short time FT
                    subplot(a,b,am_lib.flatten_(g(2:end,1:end-1))); am_lib.stft_(t,y,w_,'half'); axis tight; set(gca,'yscale','log'); box on; 
                else
                    [G,f,t] = am_lib.stft_(t,y,w_,flag);
                    surf(t,f,log(abs(G)),'edgecolor','none'); view([0 0 1]);
                end
                return;
            end
            % main part:
            N  = numel(t);
            fs = 1./(t(2)-t(1));% define sampling frequency
            f  = fftshift(fs/N*([0:N-1]-ceil(N/2))); % define frequencies
            G  = fft(w_(t.'-t).*y(:),[],1); % apply fft
            [f,t] = ndgrid(f,t);
            if contains(flag,'half')
                f = f(1:floor(end/2),:);
                t = t(1:floor(end/2),:);
                G = G(1:floor(end/2),:);
            end
        end
        
        function [g]    = dct_(f) % discrete cosine transform
            % Equivalent to this transform kernel:
            % N = 1000;
            % k = [0:(N-1)]'/N; r = [0:(N-1)]'-floor(N/2);
            % K = cos(2*pi*r*k.'); % Define cosine transform kernel.
            % f = K * g; % Evaluate transform.
            g = real(fftshift(fft(f)));
        end
        
        function [g]    = dst_(f) % discrete sine transform 
            % Equivalent to this transform kernel:
            % N = 1000;
            % k = [0:(N-1)]'/N; r = [0:(N-1)]'-floor(N/2);
            % K = sin(2*pi*r*k.'); % Define sine transform kernel.
            % f = K * g; % Evaluate transform.
            g = -imag(fftshift(fft(f)));
        end

        function [a]    = fft_autocorrelation_(x)
            n = numel(x);
            %FFT method based on zero padding
            f = fft([x; zeros(n,1)]); % zero pad and FFT
            a  = ifft(f.*conj(f)); % abs()^2 and IFFT
            % circulate to get the peak in the middle and drop one
            % excess zero to get to 2*n-1 samples
            a = [a(n+2:end); a(1:n)];
        end
        
        function [y,x]  = fftdn_(y) % down-sample y by half, x = [0,1)
            n = numel(y); n = floor(n/2); 
            x = [0:n-1]/n; y = fft(y); 
            y = real(ifft(y(1:n))); 
        end

        function [y,x]  = fftup_(y) % up-sample y two fold, x = [0,1)
            n = numel(y); x = [0:2*n-1]/(2*n); y = fft(y); 
            y = [y(1:floor(n/2)),zeros(1,n),y(floor(n/2)+1:end)]; 
            y = 2*real(ifft(y)); 
        end
        
        function [y,a,b]= fftup2_(y) % up-sample y two fold (2D), x = [0,1)
            [n,m] = size(y); a = [0:2*n-1]/(2*n); b = [0:2*m-1]/(2*m); t = zeros(2*n,2*m); 
            t([1:floor(n/2),n+[floor(n/2)+1:n]],[1:floor(m/2),m+[floor(m/2)+1:m]]) = fftn(y); y = 2*real(ifftn(t));
        end

        function [r,gr] = fft_(k,gk,flag)
            
            if nargin < 3; flag='half'; end
            
            % interpolate gk on equidistant points
            kp= linspace(min(k),max(k),numel(k));
            g = interp1(k,gk,kp); k=kp; clear ke;
            % apply fft
            N = numel(k);
            r = fftshift(([0:(N-1)]'-floor(N/2))/(k(end)-k(1))); gr = fft(g); 
            switch flag
                case 'half'
                    gr = gr(1:floor(N/2)); r = r(1:floor(N/2));
                case 'whole'
                    % do nothing
                otherwise; error('fft_: Unknown flag.');
            end
        end
        
        function [g]    = lopass_(y,w) % low pass filter
            
            g = real(ifft(fft(y) .* reshape(fftshift(am_lib.gaussw_(numel(y),w)),size(y)) ));
        end
        
        function [g]    = lopass2_(y,w) % low pass filter (2D)
            w = fftshift(am_lib.gaussw_(size(y,1),w)*am_lib.gaussw_(size(y,2),w).');
            g = real(ifftn(fftn(y) .* w ));
        end
        
        function [y,bg] = fft_filter_(y,m,npasses,flag)
            % m = # of fft components to keep
            n = size(y); d = ndims(y); bg = zeros(size(y));
            if numel(m)==1; m=repmat(m,1,d); end
            
            switch flag
                case {'lo','low','low-pass','LP'}            
                    for i = 1:npasses
                        yf = fftn(y); % yf = yf(1:max(floor(end/2),1),1:max(floor(end/2),1),1:max(floor(end/2),1));
                        yf(1:end>m(1),:,:) = 0;
                        yf(:,1:end>m(2),:) = 0; if d > 2
                        yf(:,:,1:end>m(3)) = 0; 
                        if d > 3; error('not yet implemented'); end; end
                        b = real(ifftn(yf,n)); bg = b + bg; y = y - b; 
                    end
                case {'hi','high','high-pass','HP'}            
                    for i = 1:npasses
                        yf = fftn(y); % yf = yf(1:max(floor(end/2),1),1:max(floor(end/2),1),1:max(floor(end/2),1));
                        yf(1:end<=m(1),:,:) = 0;
                        yf(:,1:end<=m(2),:) = 0; 
                        if d > 2
                        yf(:,:,1:end<=m(3)) = 0; 
                        if d > 3; error('not yet implemented'); end; end
                        b = real(ifftn(yf,n)); bg = b + bg; y = y - b; 
                    end
                otherwise
                    error('unknown filter type');
            end
        end
        
        function          plot_power_spectrum_(t,y,leg)
            % demo: am_lib.plot_power_spectrum_(x,sin(x/5*2*pi))
            
            % count number of scans
            nys=size(y,2);
            % 
            if isempty(t); for i = 1:nys; t(:,i) = [1:size(y,1)]; end; end
            % apply fft
            for i = 1:nys; [f(:,i),yf(:,i)] = am_lib.fft_(t(:,i).',y(:,i).'); end
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

        function [x,f]  = get_statistical_function(D,flag,scanaxis)
            % [x,hhcf] = get_statistical_function(D,'rHHCF',1);
            % [~,acf ] = get_statistical_function(D,'rACF',1);
            % % convert acf to hhcf:
            % rms = std(F.F(:)).^2; hhcf_from_acf = 2*(rms-acf);
            % % compare
            % loglog(x,hhcf,'-',x,hhcf_from_acf.','.')
            
            % move scan axis to first dimension
            D = permute(D,circshift([1,2,3],1-scanaxis)); n = size(D,1);
            
            switch flag
                % "radial" HHCF (as defined in gywddion)
                case {'rHHCF','height-height correlation'}
                    x = ([1:n]-1)./n; f = zeros(1,n);
                    for i = 1:n
                        L = (n-i+1); ex1_=i:n; ex2_=1:L;
                        f(i) = am_lib.sum_( (D(ex1_,:)-D(ex2_,:)).^2 ,[1,2] ) ./ (L.*n);
                    end
                % "radial" ACF (as defined in gywddion)
                case {'rACF','autocorrelation'}
                    x = ([1:n]-1)./n; f = zeros(1,n);
                    for i = 1:n
                        L = (n-i+1); ex1_=i:n; ex2_=1:L;
                        f(i) = am_lib.sum_( (D(ex1_,:).*D(ex2_,:))   ,[1,2] ) ./ (L.*n);
                    end
                otherwise
                    error('method unknown');
            end
        end


        % image processing
        
        function [cluster,neighbor] = floodfill_(F,G,P,i,maxclustersize)
            % F = scalar field
            % G(2,n) = edge list; edge n connects G(1,n) to G(2,n) with weight G(3,:)
            % P, probability that a point will be incorporated into the cluster:
            %       p = 1-exp(-2/kT(k))     for Wolff
            %       p = 1                   for flood fill
            % i = index of seed
            % 
            
            % no input is passed
            if nargin == 0
                % define dimensions, seed, and target cluster
                rng(1); n=[4,4]; j = 5; F = zeros(n); F(j+[1:2])=1;
                % compute dimensions
                p = cumprod(n);
                % define neighbors
                n_ = @(i) [mod(i   -1-1,p(end))+1, ...  % up
                           mod(i   +1-1,p(end))+1, ...  % down
                           mod(i-p(1)-1,p(end))+1, ...  % left
                           mod(i+p(1)-1,p(end))+1];     % right       
                % build edge list
                G = zeros(2,4*p(end)); G(1,:) = repelem(1:p(end),4);
                for i = 1:p(end); G(2,4*(i-1)+[1:4]) = n_(i); end
                % perform floodfill
                [ind,nbr] = am_lib.floodfill_( F , G , j , 1 , 2 )  % 1-exp(-2*beta)
                % for j=5;
                %   ind should equal [5,9]
                %   nbr should equal [4,6,1,8,10,13]
                return;
            end
            
            if isempty(i); error('seed must not be empty'); end
            if nargin < 5; maxclustersize = Inf; end 
            % initialize queue and cluster
            cluster = zeros(1,numel(F)); ic = 1;  cluster(1) = i;  nn = sum(G(1,:)==i);
            queue = zeros(1,numel(F)); nq=1; queue(1:nn+1) = [i, G(2,G(1,:)==i)]; iq=sum(queue~=0); 
            while nq~=iq
                % cycle queue
                nq=nq+1; q=queue(nq:iq); nq=iq;
                % get aligned spins and add it to cluster with probability 1-exp(-2/kT)s
                ex_ = F(q)==F(i); ex_(ex_) = rand(1,sum(ex_)) <= P; ncs = sum(ex_); 
                % limit maximum cluster size
                if ~isinf(maxclustersize)
                    m = min(maxclustersize-ic,ncs); v = [true(1,m),false(1,ncs-m)];
                    ex_(ex_) = v(randperm(ncs)); ncs = m;
                end
                % build cluster
                if ncs~=0
                    % add new queue points to cluster
                    cluster(ic+[1:ncs]) = q(ex_); ic=ic+ncs;
                    % loop over neighbors which have never been considered
                    n=G(2,any(G(1,:)==q(ex_).',1)); n=unique(n(n~=0)); ex_=am_lib.setdiffi_(n,queue(1:iq)); 
                    nns=sum(ex_); queue(iq+[1:nns]) = n(ex_); iq=iq+nns;
                end
                % % animate floodfill (slows everything down)
                % sp_DEBUG__ = zeros(size(F)); figure(1);
                % hold on;
                % sp_DEBUG__(queue(queue~=0))=1;     spy(sp_DEBUG__,'r'); sp_DEBUG__(:) = 0;
                % sp_DEBUG__(cluster(cluster~=0))=1; spy(sp_DEBUG__,'k'); sp_DEBUG__(:) = 0;
                % hold off;
                % drawnow;
            end
            cluster = cluster(1:ic); neighbor = queue(am_lib.setdiffi_(queue(1:iq),cluster)); 
        end

        function [D,R] = DT(img)
            % Two-dimensional generalized distance transform
            %
            % Input: f - the sampled function
            %
            % Output: D - distance transform
            %         R - power diagram
            %
            % Based on the paper:
            % P. Felzenszwalb, D. Huttenlocher
            % Distance Transforms of Sampled Functions
            % Cornell Computing and Information Science Technical Report TR2004-1963,
            % September 2004
            %
            %
            % This is a simple MATLAB implmentation of the generalized distance
            % transform algorithm. The function DT() gives the distance transform 
            % of a 2D image by calling DT1() for each dimension. By using DT1(), 
            % this could be easily extended to higher dimensions. It seems to have 
            % problems with inf values, so for points in the image with "no" parabola 
            % centered there, they should instead be given a large numeric value 
            % (such as 1e10). I also modified the algorithm so that the second argument 
            % returns the power diagram of the input. The power diagram is a diagram 
            % where each point is assigned to the point that is closest to it with 
            % respect to the distance transform. If all input points have the same 
            % value, this function reduces to giving the standard distance transform 
            % and the Voronoi diagram.
            % 
            % % EXAMPLE:
            % % 
            % % Create some random points
            % X = randi(100,50,2);
            % % Create an image with rando values at these points
            % img = sparse(X(:,2), X(:,1), rand(50,1)*20,100,100);
            % % Set all other values to a high number
            % img(img==0) = 1e10;
            % % Call the function
            % [D R] = DT(img);
            % % Plot the results
            % figure;
            % subplot(1,2,1);
            % imagesc(D);
            % title('Generalized Distance transform');
            % axis image;
            % subplot(1,2,2);
            % imagesc(R);
            % title('Power diagram');
            % axis image;

            D = zeros(size(img));
            R = zeros(size(img));

            for i = 1:size(img,1)
                [d r] = DT1(img(i,:));
                R(i,:) = r;
                D(i,:) = d;
            end
            for j = 1:size(img,2)
                [d r] = DT1(D(:,j));
                D(:,j) = d;
                R(:,j) = sub2ind(size(img), r, R(r,j));
            end
            function [D R] = DT1(f)
                % One-dimensional generalized distance transform
                %
                % Input: f - the sampled function
                %
                % Output: D - distance transform
                %         R - power diagram
                %
                % Based on the paper:
                % P. Felzenszwalb, D. Huttenlocher
                % Distance Transforms of Sampled Functions
                % Cornell Computing and Information Science Technical Report TR2004-1963, September 2004

                n = numel(f);
                D = zeros(n,1);
                R = zeros(n,1);

                k = 1; % Index of the rightmost parabola in the lower envelope
                v = ones(n,1); % Locations of the parabolas in the lower envelope
                z = ones(n,1); %Locations of boundaries between parabolas
                z(1) = -inf;
                z(2) = inf;

                for q = 2:n
                    s = ((f(q) + q^2) - (f(v(k)) + v(k)^2))/(2*q - 2*v(k));
                    while s <= z(k)
                        k = k - 1;
                        s = ((f(q) + q^2) - (f(v(k)) + v(k)^2))/(2*q - 2*v(k));
                    end
                    k = k + 1;
                    v(k) = q;
                    z(k) = s;
                    z(k+1) = inf;
                end
                k = 1;
                for q = 1:n
                    while z(k+1) < q
                        k = k+1;
                    end
                    D(q) = (q-v(k))^2 + f(v(k));
                    R(q) = v(k);
                end
            end
        end

        
        % vector calculous
        
        function [C]    = curl_(f)
            % f = [(x,y,z),x,y,z) -> [(dx,dy,dz),x,y,z)
            % curl of vector field
            
            [m,n(1),n(2),n(3)] = size(f); if m~=2&&m~=3; error('curl_ requires (x,y) or (x,y,z) in first dimension'); end

            % generate Fourier mesh
            fftmesh = @(N) fftshift([0:(N-1)]'-floor(N/2)); 
            [r(1,:,:,:),r(2,:,:,:),r(3,:,:,:)]=ndgrid(fftmesh(n(1)),fftmesh(n(2)),fftmesh(n(3)));

            % compute curl
            C = zeros([m,n]);
            for i = 1:m; C(i,:,:,:) =  fftn(f(i,:,:,:)); end
            C = cross( 2i*pi*r, C, 1);
            for i = 1:m; C(i,:,:,:) = ifftn(C(i,:,:,:)); end
            
            if all(isreal(f(:))); C = real(C); end
        end

        function [G]    = grad_(f)
            % f = [(x,y,z),x,y,z] -- > [(x,y,z),(dx,dy,dz),x,y,z]
            % f = [x,y,z]         -- > [      1,(dx,dy,dz),x,y,z]
            % gradient of vector or scalar field
            
            n = size(f); 
            
            switch numel(n)
                % vector field
                % 2D [(x,y)  ,x,y,z]
                % 3D [(x,y,z),x,y,z]
                case 4
                    
                    m=n(1); n=n(2:4);                    
                    if m~=2&&m~=3; error('grad_ requires (x,y) or (x,y,z) in first dimension'); end

                    % generate Fourier mesh
                    fftmesh = @(N) fftshift([0:(N-1)]'-floor(N/2)); 
                    [r(1,:,:,:),r(2,:,:,:),r(3,:,:,:)]=meshgrid(fftmesh(n(1)),fftmesh(n(2)),fftmesh(n(3)));

                    % compute gradient d/dx + d/dy + d/dz
                    G = zeros([m,m,n]);
                    for i = 1:m; fi = fftn(f(i,:,:,:));
                    for j = 1:m; G(i,j,:,:,:) = ifftn( dot( 1i*2*pi*r(j,:,:,:), fi, 1) );
                    end; end
                
                % scalar field
                % 3D [(x,y,z),x,y,z]
                case 3
                    
                    % generate Fourier mesh
                    fftmesh = @(N) fftshift([0:(N-1)]'-floor(N/2)); 
                    [r(1,:,:,:),r(2,:,:,:),r(3,:,:,:)]=meshgrid(fftmesh(n(1)),fftmesh(n(2)),fftmesh(n(3)));

                    % compute gradient d/dx, d/dy, d/dz
                    G = zeros([1,3,n]); fi(1,:,:,:) = fftn(f);
                    for j = 1:3; G(1,j,:,:,:) = ifftn( dot( 1i*2*pi*r(j,:,:,:), fi(1,:,:,:), 1) );
                    end
                    
                otherwise
                    error('invalid input size');
            end
            
            if all(isreal(f(:))); G = real(G); end
        end
        
        function [D]    = div_(f)
            % f = [(x,y,z),x,y,z) -- > d/dx+d/dy+d/dz [x,y,z]
            % divergence of vector field (trace over tensor of vector field gradient)
            
            
            [m,n(1),n(2),n(3)] = size(f); if m~=2&&m~=3; error('div_ requires (x,y) or (x,y,z) in first dimension'); end

            % generate Fourier mesh
            fftmesh = @(N) fftshift([0:(N-1)]'-floor(N/2)); 
            [r(1,:,:,:),r(2,:,:,:),r(3,:,:,:)]=meshgrid(fftmesh(n(1)),fftmesh(n(2)),fftmesh(n(3)));

            % compute gradient d/dx + d/dy + d/dz
            D = zeros([1,n]);
            for i = 1:m; D = D + ifftn(dot( 1i*2*pi*r(i,:,:,:), fftn(f(i,:,:,:)), 1)); end
            D = squeeze(D);
            
            if all(isreal(f(:))); D = real(D); end
        end
        
        
        % interpolation

        function [fq,S,iS,c] = fftinterp_(f,q,n,algo,R)
            % R is only required for algo = 'star'
            % R must be input in real space!
            %
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
                    if nargin < 3
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
                    % evaluate kernels
                    K = cos_kernel_(k,r,n); iK = icos_kernel_(q,r,n);
                    % get dft expansion coefficients: c [ nrs, n ]
                    c = K * f.';
                    % interpolate f on q
                    fq = ( iK * c ).';
                case 'star'
                    % get the correct orientation
                    if size(f,1)~=prod(n(1:3)); f = f.'; end

                    % define function to get kernel
                     fft_kernel_ = @(r,k,n) exp(-2i*pi*am_lib.matmul_(reshape(r,1,3,[],1),reshape(k,3,1,1,[])))./sqrt(prod(n(1:3)));
                    ifft_kernel_ = @(k,r,n) exp(+2i*pi*am_lib.matmul_(reshape(k,1,3,[],1),reshape(r,3,1,1,[])))./sqrt(prod(n(1:3)));

                    % generate direct and reciprocal meshes r and k
                    m_ = @(i) [0:(n(i)-1)]-floor(n(i)/2); [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); r = [Y{1}(:),Y{2}(:),Y{3}(:)].';
                    m_ = @(i) [0:(n(i)-1)]./n(i);         [Y{1:3}]=ndgrid(m_(1),m_(2),m_(3)); k = [Y{1}(:),Y{2}(:),Y{3}(:)].';

                    % get irreducible maps in direct and reciprocal space (must do k and r space explicitly!)
                    %     (cannot just use f2i from fbz/ibz because R ~= R.' in fractional coordinates)
                    Rk = permute(R,[2,1,3]);
                    PM = am_lib.member_(am_lib.matmul_(R ,r),r);              A = am_lib.get_connectivity(PM); r_w=sum(A,2).'; r_i2f = round(am_lib.findrow_(A)).'; r_f2i = round(([1:size(A,1)]*A)); 
                    PM = am_lib.member_(am_lib.mod_(am_lib.matmul_(Rk,k)),k); A = am_lib.get_connectivity(PM); k_w=sum(A,2).'; k_i2f = round(am_lib.findrow_(A)).'; k_f2i = round(([1:size(A,1)]*A));

                    % get star function (reciprocal to real) and inverse start function (real to reciprocal space)
                    a = max(r_f2i); b = max(k_f2i); S = zeros(a,b); iS = zeros(size(q,2),a);
                    for j = 1:b;  S(:,j) = sum( fft_kernel_(r(:,r_i2f),k(:,k_f2i==j),n),2);  end
                    for i = 1:a; iS(:,i) = sum(ifft_kernel_(q,r(:,r_f2i==i),n),2);           end

                    % get expansion coefficients: c [ nrs, n ]
                    c = S * f(k_i2f,:);
                    % interpolate
                    fq = real(iS*c).';
                    
                otherwise
                    error('unknown algorithm');
            end
            
        end          


        % linear interpolation
        
        function y    = linspacen_(v1, v2, n)
            % linearly interpolate vectors
            n  = double(n); v1 = squeeze(v1); v2 = squeeze(v2);
            NDim = ndims(v1);
            if NDim==2 && any(size(v1)==1)
                NDim = NDim-1;
                if all(size(v1)==1)
                    NDim = 0;
                end
            end
            pp      = (0:n-2)./(floor(n)-1);
            Sum1    = tensor_product(v1, ones(1,n-1));
            Sum2    = tensor_product((v2-v1), pp);
            y = cat(NDim+1, Sum1  + Sum2, shiftdim(v2, size(v1, 1)==1 ));

            function Z = tensor_product(X,Y)
                sX=size(X);sY=size(Y); ndim1=ndims(X);ndim2=ndims(Y); indperm=[ndim2+1:ndim1+ndim2,1:ndim2];
                Z=squeeze(repmat(X,[ones(1,ndims(X)),sY]).*permute(repmat(Y,[ones(1,ndims(Y)),sX]),indperm));
            end
        end
        

        % utilities
        function [ha, pos] = subplot_(Nh, Nw, gap, marg_h, marg_w)

        % tight_subplot creates "subplot" axes with adjustable gaps and margins
        %
        % [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
        %
        %   in:  Nh      number of axes in hight (vertical direction)
        %        Nw      number of axes in width (horizontaldirection)
        %        gap     gaps between the axes in normalized units (0...1)
        %                   or [gap_h gap_w] for different gaps in height and width 
        %        marg_h  margins in height in normalized units (0...1)
        %                   or [lower upper] for different lower and upper margins 
        %        marg_w  margins in width in normalized units (0...1)
        %                   or [left right] for different left and right margins 
        %
        %  out:  ha     array of handles of the axes objects
        %                   starting from upper left corner, going row-wise as in
        %                   subplot
        %        pos    positions of the axes objects
        %
        %  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
        %           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
        %           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

        % Pekka Kumpulainen 21.5.2012   @tut.fi
        % Tampere University of Technology / Automation Science and Engineering


            if nargin<3; gap = .02; end
            if nargin<4 || isempty(marg_h); marg_h = .05; end
            if nargin<5; marg_w = .05; end

            if numel(gap)==1; 
                gap = [gap gap];
            end
            if numel(marg_w)==1; 
                marg_w = [marg_w marg_w];
            end
            if numel(marg_h)==1; 
                marg_h = [marg_h marg_h];
            end

            axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
            axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

            py = 1-marg_h(2)-axh; 

            % ha = zeros(Nh*Nw,1);
            ii = 0;
            for ih = 1:Nh
                px = marg_w(1);

                for ix = 1:Nw
                    ii = ii+1;
                    ha(ii) = axes('Units','normalized', ...
                        'Position',[px py axw axh], ...
                        'XTickLabel','', ...
                        'YTickLabel','');
                    px = px+axw+gap(2);
                end
                py = py-axh-gap(1);
            end
            if nargout > 1
                pos = get(ha,'Position');
            end
            ha = ha(:);
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
                
        function [X] = call_system_(fnc,fname,fpattern)
            % calls fnc, waits until fpattern is found in fname
            system(fnc);

            switch nargin
                case 1
                    % do nothing
                case 2
                    while true
                        if exist(fname,'file'); break; end
                        pause(5); % Wait 5 second.
                    end
                case 3
                    while true
                        if system(sprintf('grep "%s" %s',fpattern,fname))==0; break; end
                        pause(5); % Wait 5 second.
                    end
            end            
        end
        
        function [n] = count_lines_(fname)
            if ispc
                [~,a] = system(sprintf('type %s | find /c /v ""',fname));
            elseif or(ismac,isunix)
                [~,a] = system(sprintf('wc -l %s',fname));
            end
            n = sscanf(a,'%i'); n = max(n(:));
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








