classdef am_field

    properties (Constant)
        tiny = 1E-8; 
        mu0  = 1; % vacuum permeability
    end
    
    properties
        T = []; % type (scalar/vector)
        d = []; % dimensions (2 or 3)
        s = []; % scheme{1:d} ('fourier/chebyshev/legendre/cdiff') for each dimension
        n = []; % grid points / dimension          3D: [n(1),n(2),n(3)]     2D: [n(1),n(2)] 
        a = []; % lattice/grid spacing             3D: [a(1),a(2),a(3)]     2D: [n(1),n(2)] 
        R = []; % cartesian coordinates                [   x , y , z  ]
        F = []; % field
                % vector field
                %     [ (             F(x)             )              ]
                %     [ (             F(y)             ) , x , y , z  ]
                %     [ (             F(z)             )              ]
                % scalar field
                %     [ (               1              ) , x , y , z  ]
        J = []; % jacobian
                % J = jacobian (for vector field) [3,3,n(1),n(2),n(3)]
                %     [ ( dF(x)/dx  dF(x)/dy  dF(x)/dz )              ]
                %     [ ( dF(y)/dx  dF(y)/dy  dF(y)/dz ) , x , y , z  ]
                %     [ ( dF(z)/dx  dF(z)/dy  dF(z)/dz )              ]
                % J = jacobian (for scalar field) [3,3,n(1),n(2),n(3)]
                %     [ (  dF/dx        0        0     )              ]
                %     [ (     0       dF/dy      0     ) , x , y , z  ]
                %     [ (     0         0      dF/dz   )              ]
        H = []; % hessian
                % H = jacobian (for scalar field) [3,3,n(1),n(2),n(3)]
                %     [ ( d2F/dxdx  d2F/dxdy  d2F/dxdz )              ]
                %     [ ( d2F/dydx  d2F/dydy  d2F/dydz ) , x , y , z  ]
                %     [ ( d2F/dzdx  d2F/dzdy  d2F/dzdz )              ]
        D = []; % divergence
                % trace of jacobian = 
                %     [         dF/dx + dF/dy + dF/z     , x , y , z  ]
        L = []; % laplacian
                % trace of hessian = 
                %     [d^2F/dx^2 + d^2F/dy^2 + d^2F/dz^2 , x , y , z  ]
        C = []; % curl
                % skew  of jacobian = 
                %     [ (     dFz/dy    -    dFy/dz    ) 
                %     [ (     dFx/dz    -    dFz/dx    ) , x , y , z  ]
                %     [ (     dFy/dx    -    dFx/dy    ) 
    end

    methods (Static)
        
        function F = demo()
            
            import am_field.*
            
            F   = define_field([2,2].^[7,8],[2*pi,2*pi],{'cdiff','fourier'});
            F.R = get_collocation_points(F);
            F.F = cat(1,sin(F.R(1,:,:,:)), sin(F.R(2,:,:,:))); 
            F.F = sum(F.F,1);
            F.T = get_field_type(F);
            F.J = get_jacobian(F);
            F.D = get_divergence(F);
            F.C = get_curl(F);
            F.H = get_hessian(F);
            
            figure(1); plot_jacobian(F);
            figure(2); plot_hessian(F);
            
        end
        
        function F = demo_biot_savart()
            
            import am_field.*
            
            F   = define_field([2,2,2].^[5,5,5],[1,1,1],{'chebyshev','chebyshev','chebyshev'});
            F.R = get_collocation_points(F);

            D_ = @(N)      sparse([1:N],1+[1:N],1,N,N+1)-sparse([1:N],[1:N],1,N,N+1); % Define forward-difference and forward-mean transforms.
            M_ = @(N) 0.5*(sparse([1:N],1+[1:N],1,N,N+1)+sparse([1:N],[1:N],1,N,N+1));

            I_ = {@(d,th,M) [d*sin(2*pi*th) ,d*cos(2*pi*th),zeros(M+1,1)]; % Define a library of current paths: circle, solenoid, and straight wire.
                  @(d,th,M) [d*sin(20*pi*th),d*cos(20*pi*th),linspace(min(F.R(3,:)),max(F.R(3,:)),M+1)'];
                  @(d,th,M) [zeros(M+1,1),zeros(M+1,1),mpgrid(M+1)];};

            iI = 1; r = 0.5; M = 500; th = [0:M]'/M; % Select and construct a current path.
            dI = (D_(M)*I_{iI}(r,th,M)).'; 
             I = (M_(M)*I_{iI}(r,th,M)).';

            F.F = am_field.get_vector_potential(F.R,dI,I); % get magnetic vector potential 
            F = F.get_derivatives(); % compute derivatives

            subplot(1,2,1); F.plot_field('F'); axis tight; title('vector potential A'); % plot magnetic vector potential
            line(I(1,:),I(2,:),I(3,:),'linewidth',2,'color',[1 1 1]*0.5); % draw wire
            subplot(1,2,2); F.plot_field('C'); axis tight; title('magnetic field B'); % plot magnetic field (curl of vector potential)
            line(I(1,:),I(2,:),I(3,:),'linewidth',2,'color',[1 1 1]*0.5); % draw wire
            
        end

        function F = demo_hilliard_cahn()
            
            import am_field.* 

            F = am_field.define_field([2,2].^[6,6],[2,2].^[6,6],{'pdiff','pdiff'}); % initialize object
            F.F = rand([1,F.n]); F.F = F.F - mean(F.F(:)); % initialize random field
            F = F.solve_differential_equation('CH58',{@(x)(x^2-1)^2/4,1},'explicit',0.01,10000); % solve
        end
        
        function [F] = define_field(n,a,s)
            m = numel(n);
            if numel(s)~=m; error('s dimensions mismatch'); end
            if numel(a)~=m; error('a dimensions mismatch'); end
            F = am_field(); F.a=a; F.n=n; F.s=s; F.d=numel(n); F.R = get_collocation_points(F);
        end
        
    end 
    
    methods 

        function [F] = extract_slice(F,i,j)
            % i = dimension extracting
            % j = position along i
            import am_field.*
                  
            % trim dimension
            F.d   = F.d-1;
            % trim scheme, grid size, grid spacing
            for f = {'s','n','a'}; if ~isempty(F.(f{:}))
                F.(f{:})(i) = []; 
            end; end
            % trim coordinates, field, jacobian, hexxian, divergence, curl
            for f = {'R','F','J','H','D','C'}; if ~isempty(F.(f{:}))
                switch f{:}
                    case 'R'; o=1; case 'F'; o=1; case 'J'; o=2;
                    case 'H'; o=2; case 'D'; o=1; case 'C'; o=1;
                end
                % permute
                p=[1:ndims(F.(f{:}))]; p([1,i+o])=p([i+o,1]);     
                % trim and permute back to orginal
                F.(f{:}) = permute(F.(f{:}),p); F.(f{:}) = permute(F.(f{:})(j,:,:,:,:,:,:,:),p);
                % permute trimmed to end
                p=[1:ndims(F.(f{:}))+1]; p([i+o,end])=p([end,i+o]); F.(f{:}) = permute(F.(f{:}),p);
                switch f{:}
                    case 'R'; F.(f{:})(i,:,:,:,:,:) = []; 
                    case 'F'; F.(f{:})(i,:,:,:,:,:) = []; 
                end
            end; end
        end
        
        function [F] = get_derivatives(F)
            import am_field.*
            F.T = get_field_type(F);
            F.J = get_jacobian(F);
            F.D = get_divergence(F);
            F.C = get_curl(F);
            if strcmp(F.T,'scalar')
                F.H = get_hessian(F); 
                F.L = get_laplacian(F);
            end
        end
        
        function [F] = solve_differential_equation(F,equation,x,solver,dt,M)
            % examples
            % F = am_field.define_field([2,2].^[6,6],[2,2].^[6,6],{'pdiff','pdiff'}); % initialize field
            % F.F = rand([1,F.n]); F.F = F.F-mean(F.F(:)); % initialize field
            % F = F.solve_differential_equation('SW76',{0.1,1.0},'implicit',2,500);         % hexagons
            % F = F.solve_differential_equation('SW76',{0.3,0.0},'implicit',2,500);         % finger prints
            % F = F.solve_differential_equation('SW76',{0.8,1.85},'explicit',0.03,5000);    % maze
            % F = F.solve_differential_equation('CH58',{@(x)(x^2-1)^2/4,1},'explicit',0.01,10000);
            % F = F.solve_differential_equation('CH58',{@(x)(x^2-1)^2/4,1},'implicit',1,1000);
            % F = F.solve_differential_equation('GL50',{@(x)(x^2-1)^2/4,1},'explicit',0.01,5000);
            %
            % rho = am_lib.gauss_(am_lib.normc_(F.R-[30;30]));                              % Poisson solver
            % F = F.solve_differential_equation('poisson',{-rho},'explicit',0.01,5000);
            %
            import am_field.*
            
            [~,L] = F.get_flattened_differentiation_matrices(); % Get laplacian.
            
            N=prod(F.n); % Simplify notation: get total number of grid points.

			switch solver
            case {'explicit','implicit','crank-nicolson','jacobi','gauss-seidel'}
            	if any(strcmp(solver,{'explicit','crank-nicolson'}))
	                % equations of the form F(n+1) = F(n) + dt * LHS
                    switch equation
                        case 'SW76' % Swift-Hohenberg (PRA 1976), x = {eps, g1}
                            LHSe_ = @(U,x) x{1}*U(:) - (L+speye(N))^2*U(:) + x{2}*U(:).^2 - U(:).^3; nargs=2;
                        case 'LP97' % Lifshitz-Petrich (PRL 1997), x = {e, g1, q}
                            LHSe_ = @(U,x) x{1}*U(:) - (L+speye(N))^2*(L+x{3}.^2*speye(N))^2*U(:) + x{2}*U(:).^2 - U(:).^3; nargs=3;
                        case 'CH58' % Cahn-Hilliard (J. Chem. Phys. 1958), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(diff(x{1}(z),z)); % P.E. derivative
                            LHSe_ = @(U,x) L*( x{1}(U(:)) - x{2}*L*U(:) ); nargs=2;
                        case 'GL50' % Ginzburg-Landau (Zh. Eksp. Teor. Fiz. 1950), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(diff(x{1}(z),z)); % P.E. derivative
                            LHSe_ = @(U,x)   ( x{1}(U(:)) + x{2}*L*U(:) ); nargs=2;
                        case 'poisson' % Poisson equation
                            LHSe_ = @(U,x) L * U(:) - x{1}(:); nargs=1;
                        otherwise
                            error('unknown propagator');
                    end
                end
                if any(strcmp(solver,{'implicit','crank-nicolson'}))
                % equations of the form F(n+1) = (1 - dt*LHS)\F(n)
                    switch equation
                        case 'SW76' % Swift-Hohenberg (PRA 1976), x = {eps, g1}
                            LHSi_ = @(U,x) x{1}*speye(N) - (L+speye(N))^2 + x{2}*spdiags(U(:),0,N,N) - spdiags(U(:).^2,0,N,N); nargs=2;
                        case 'LP97' % Lifshitz-Petrich (PRL 1997), x = {e, g1, q}
                            LHSi_ = @(U,x) x{1}*speye(N) - (L+speye(N))^2*(L+x{3}.^2*speye(N))^2 + x{2}*spdiags(U(:),0,N,N) - spdiags(U(:).^2,0,N,N); nargs=3;
                        case 'CH58' % Cahn-Hilliard (J. Chem. Phys. 1958), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(expand(diff(x{1}(z),z)/z)); % P.E. derivative with field factored out
                            LHSi_ = @(U,x) L*( spdiags(x{1}(U(:)),0,N,N) - x{2}*L ); nargs=2;
                        case 'GL50' % Ginzburg-Landau (Zh. Eksp. Teor. Fiz. 1950), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(expand(diff(x{1}(z),z)/z)); % P.E. derivative with field factored out
                            LHSi_ = @(U,x)   ( spdiags(x{1}(U(:)),0,N,N) + x{2}*L ); nargs=2;
                        case 'poisson' % Poisson equation
                            LHSi_ = @(U,x) L - spdiags(x{1}(:),0,N,N); nargs=1;
                        otherwise
                            error('unknown propagator');
                    end
                end
                if any(strcmp(solver,{'jacobi','gauss-seidel'}))
                % equations of the form F(n+1) = (LHS*F(n)-x{1}), i.e. Ax-b = 0 
                    switch equation
                        case 'poisson' % Poisson equation
                            LHSx_ = @(x) L; nargs=1; % b is implicitly assumed
                        otherwise
                            error('unknown propagator');
                    end
                end
            otherwise
                error('unknown solver');
            end
            
            if ~iscell(x); x=num2cell(x); end % convert array to cells if necessary (cells allow passing functions)
            if numel(x)~=nargs; error('incorrect number of parameters'); end % check parameters
            if isempty(F.F); error('field has not been initialized'); end % check parameters

            % Propagate in time 
            switch solver
                case 'explicit' % unstable for large steps
                    for i = [1:M]
                        F.F(1:N) = F.F(:) + dt*LHSe_(F.F,x);
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0; F.plot_field('F'); title(num2str(i)); drawnow; end
                    end
                case 'implicit' % more stable
                    for i = [1:M]
                        F.F(1:N) = (speye(N)-dt*LHSi_(F.F,x))\F.F(:);
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0; F.plot_field('F'); title(num2str(i)); drawnow; end
                    end
                case 'crank-nicolson'
                    for i = [1:M]
                        F.F(1:N) = ( (speye(N)-dt*LHSi_(F.F,x))\F.F(:) + F.F(:)+dt*LHSe_(F.F,x) )/2;
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0; F.plot_field('F'); title(num2str(i)); drawnow; end
                    end
                case 'jacobi'
                    % not working
                    % not working
                    % not working
                    LDU_ = @(A) deal( -(tril(A)-spdiags(diag(A),0,N,N)), spdiags(diag(A),0,N,N), -(triu(A)-spdiags(diag(A),0,N,N)) ); % A = D - L - U;
                    w = 2/3; [L,D,U] = LDU_(LHSx_(x)); Rj = D\(L+U); Rw = (1-w)*speye(N) + w*Rj; 
                    for i = [1:M]
                        F.F(1:N) = Rw * F.F(:);
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0; F.plot_field('F'); title(num2str(i)); drawnow; end
                    end
                case 'gauss-seidel'
                    % not working
                    % not working
                    % not working
                    % not working
                    LDU_ = @(A) deal( -(tril(A)-spdiags(diag(A),0,N,N)), spdiags(diag(A),0,N,N), -(triu(A)-spdiags(diag(A),0,N,N)) ); % A = D - L - U;
                    [L,D,U] = LDU_(LHSx_(x)); Rg = (D-L)\U; 
                    for i = [1:M]
                        F.F(1:N) = Rg*F.F(:);
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0; F.plot_field('F'); title(num2str(i)); drawnow; end
                    end
                otherwise
                    error('unknown solver');
            end
        end
        
        function [h] = plot_field(F,field)
            
            sl_ = @(field,i)   squeeze(F.(field)(i,:,:,:,:,:));
            
            switch size(F.(field),1)
                case {1} % scalar
                    switch F.d
                        case 2 % 2D
                            set(gcf,'color','w');
                            h = surf(sl_('R',1), sl_('R',2), squeeze(F.(field))); 
                            h.EdgeColor= 'none'; h.LineWidth = 1; 
                            view([0 0 1]); daspect([1 1 1]); axis tight;
                        case 3 % 3D
                            error('not yet implemented');
                        otherwise; error('invalid field dimension');
                    end
                case {2,3} % vector field
                    switch F.d
                        case 2 % 2D
                            set(gcf,'color','w');
                            h = quiver(sl_('R',1), sl_('R',2), sl_(field,1), sl_(field,2) ); 
                            h.LineWidth = 1; view([0 0 1]); daspect([1 1 1]); axis tight;
                        case 3 % 3D
                            quiver3(...
                                    sl_('R',2)  , sl_('R',1)  , sl_('R',3),...
                                    sl_(field,2), sl_(field,1), sl_(field,3),...
                                   'AutoScaleFactor',3,'ShowArrowHead','off','linewidth',0.5);
                            streamslice(...
                                    sl_('R',2)  , sl_('R',1)  , sl_('R',3),...
                                    sl_(field,2), sl_(field,1), sl_(field,3),...
                                   [min(F.R(2,:))],[min(F.R(1,:))],[min(F.R(2,:))]);
                            set(gca,'DataAspectRatio',[1,1,1],'CameraPosition',[2,1,1],'Box','on');
                        otherwise; error('invalid field dimension');
                    end
            end

        end

        function [ax]= plot_jacobian(F)
            
            sl_ = @(R,i)   squeeze(R(i,:,:,:,:,:));
            sl2_= @(R,i,j) squeeze(R(i,j,:,:,:,:));
            
            mm = am_lib.minmax_(F.J(:)); xyz=['x','y','z']; set(gcf,'color','w');
            for j = 1:3
            for i = 1:3
                n=j+3*(i-1);
                ax{n} = subplot(3,3,n); h = surf(sl_(F.R,1), sl_(F.R,2), sl2_(F.J,i,j) ); 
                h.EdgeColor= 'none'; h.LineWidth = 1; view([0 0 1]); daspect([1 1 1]); axis tight;
                title(sprintf('d_%sf_%s',xyz(i),xyz(j))); caxis(mm);
            end
            end
            linkaxes([ax{:}],'xy');
            % colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);
        end

        function [ax]= plot_hessian(F)
            
            sl_ = @(R,i)   squeeze(R(i,:,:,:,:,:));
            sl2_= @(R,i,j) squeeze(R(i,j,:,:,:,:));
            
            mm = am_lib.minmax_(F.H(:)); xyz=['x','y','z']; set(gcf,'color','w');
            for j = 1:3
            for i = 1:3
                n=j+3*(i-1);
                ax{n} = subplot(3,3,n); h = surf(sl_(F.R,1), sl_(F.R,2), sl2_(F.H,i,j) ); 
                h.EdgeColor= 'none'; h.LineWidth = 1; view([0 0 1]); daspect([1 1 1]); axis tight;
                title(sprintf('d_%s_%sf',xyz(i),xyz(j))); caxis(mm);
            end
            end
            linkaxes([ax{:}],'xy');
            % colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);
        end
        
        function [h] = plot_divergence(F)
            set(gcf,'color','w');
            
            sl_ = @(R,i) squeeze(R(i,:,:,:,:,:));
            
            switch F.d
                case 2
                    mm = am_lib.minmax_(F.D(:)); 
                    h = surf( sl_(F.R,1), sl_(F.R,2), sl_(F.D,1) ); 
                    h.EdgeColor= 'none'; h.LineWidth = 1; view([0 0 1]); daspect([1 1 1]); axis tight;
                    caxis(mm); % colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);
                case 3
                    isosurface( sl_(F.R,1), sl_(F.R,2), sl_(F.R,3), sl_(F.D,1) );
                    daspect([1 1 1]); axis tight;
            end
        end
        
    end
   
    methods % internal stuff
        
        function [R] = get_collocation_points(F)
            for i = 1:F.d % loop over dimensions
                n = F.n;
                switch F.s{i}
                    case 'chebyshev'; R{i} = am_field.chebyshevUr_(n(i),'edge');
                    case 'legendre';  R{i} = am_field.legendrer_(n(i));
                    case 'fourier';   R{i} = am_field.fourierr_(n(i));
                    case 'cdiff';     R{i} = am_field.cdiff_(n(i));
                    case 'pdiff';     R{i} = am_field.pdiff_(n(i));
                    otherwise; error('unknown s');
                end
                n(i) = 1; R{i} = repmat(permute(F.a(i)*R{i},circshift([1,2,3],i-1)),n);
            end
            R = permute(cat(4,R{:}),[4,1,2,3]);
        end

        function [T] = get_field_type(F)
            switch size(F.F,1)
                case {1}  ; T = 'scalar';
                case {2,3}; T = 'vector';
                otherwise; error('unknown field type');
            end
        end

        function [J] = get_jacobian(F)
            % define matmul which supports sparse matrices
            matmul_ = @(A,B) reshape(A*reshape(B,size(B,1),[]),size(A,1),size(B,2),size(B,3),size(B,4),size(B,5));
            % allocate space
            J = zeros([3,3,F.n]);
            for i = 1:F.d % loop over dimensions
                switch F.s{i}
                    case 'chebyshev'; [~,D] = am_field.chebyshevUr_(F.n(i),'edge');
                    case 'legendre';  [~,D] = am_field.legendrer_(F.n(i));
                    case 'fourier';   [~,D] = am_field.fourierr_(F.n(i));
                    case 'cdiff';     [~,D] = am_field.cdiff_(F.n(i));
                    case 'pdiff';     [~,D] = am_field.pdiff_(F.n(i));
                    otherwise; error('unknown s');
                end 
                D = D(:,:,1)/F.a(i); % keep only first derivative
                if strcmp(F.s{i},'cdiff'); D=sparse(D); end % speed up finite difference with sparse matrices
                p = [1:F.d+1]; p([1,i+1])=p([i+1,1]);
                switch F.T % evaluate derivatives
                    case 'scalar'; J(i,i,:,:,:)     = permute(matmul_(D,permute(F.F,p)),p);
                    case 'vector'; J(i,1:F.d,:,:,:) = permute(matmul_(D,permute(F.F,p)),p);
                    otherwise; error('unknown field type');
                end
            end
        end

        function [H] = get_hessian(F)
            if ~strcmp(F.T,'scalar'); error('hessian is only defined for scalar fields'); end
            % define matmul which supports sparse matrices
            matmul_ = @(A,B) reshape(A*reshape(B,size(B,1),[]),size(A,1),size(B,2),size(B,3),size(B,4),size(B,5));
            % allocate space
            H = zeros([3,3,F.n]);
            for i = 1:F.d % loop over dimensions
                switch F.s{i}
                    case 'chebyshev'; [~,D] = am_field.chebyshevUr_(F.n(i),'edge');
                    case 'legendre';  [~,D] = am_field.legendrer_(F.n(i));
                    case 'fourier';   [~,D] = am_field.fourierr_(F.n(i));
                    case 'cdiff';     [~,D] = am_field.cdiff_(F.n(i));
                    case 'pdiff';     [~,D] = am_field.pdiff_(F.n(i));
                    otherwise; error('unknown s');
                end 
                D = D(:,:,1)/F.a(i); % keep only first derivative
                if strcmp(F.s{i},'cdiff'); D=sparse(D); end % speed up finite difference with sparse matrices
                p = [1:F.d+1]; p([1,i+1])=p([i+1,1]); % evaluate hessian from jacobian
                switch F.T % evaluate derivatives
                    case 'scalar'; Ji = am_lib.diag_(F.J);
                        H(i,:,:,:,:) = permute(matmul_(D,permute(Ji,p)),p);
                    otherwise; error('hessian is only defined for scalar fields');
                end
            end
        end

        function [D] = get_divergence(F)
            D = am_lib.trace_(F.J);
        end
        
        function [L] = get_laplacian(F)
            L = am_lib.trace_(F.H);
        end

        function [C] = get_curl(F)
            C = cat(1, F.J(3,2,:,:,:)-F.J(2,3,:,:,:), ...
                       F.J(1,3,:,:,:)-F.J(3,1,:,:,:), ...
                       F.J(2,1,:,:,:)-F.J(1,2,:,:,:));
            C = permute(C,[1,3,4,5,6,7,8,2]);
        end
        
        function [D,L] = get_flattened_differentiation_matrices(F)
            % checked with this code:
            % % 2D difference matrix
            % F = am_field.define_field([500,500],[1,1]*2*pi,{'chebyshev','chebyshev'});
            % F.F = sum(sin(F.R),1); F = F.get_derivatives; 
            % [D,L] = get_flattened_differentiation_matrices(F);
            % L = reshape(L*F.F(:),F.n);
            % max(abs(L(:)-F.L(:)))
            % subplot(1,2,1); surf(L,'edgecolor','none'); view([0 0 1]); daspect([1 1 1]); axis tight;
            % subplot(1,2,2); surf(squeeze(F.L),'edgecolor','none'); view([0 0 1]); daspect([1 1 1]); axis tight;
            
            for i = 1:F.d % loop over dimensions
                switch F.s{i}
                    case 'chebyshev'; [~,Q{i}] = am_field.chebyshevUr_(F.n(i),'edge'); 
                    case 'legendre';  [~,Q{i}] = am_field.legendrer_(F.n(i));
                    case 'fourier';   [~,Q{i}] = am_field.fourierr_(F.n(i));
                    case 'cdiff';     [~,Q{i}] = am_field.cdiff_(F.n(i));
                    case 'pdiff';     [~,Q{i}] = am_field.pdiff_(F.n(i));
                    otherwise; error('unknown s');
                end 
                Q{i} = Q{i}./reshape(F.a(i).^[1:2],1,1,2);
            end
            % note: although the specral differentiation matrix may be full, 
            %       when spanning multiple dimensions it will become sparse.
            % divergence
            D = cellfun(@(x)sparse(x(:,:,1)),Q,'UniformOutput',false);
            D = am_field.get_flattened_divergence(D{:});
            % laplacian
            L = cellfun(@(x)sparse(x(:,:,2)),Q,'UniformOutput',false);
            L = am_field.get_flattened_divergence(L{:});
        end
        
    end
     
    methods (Static) % polynomials and spectral methods

        function [x]     = canonicalr_(n) % roots of canonical all-one polynomial
           x(:,1) = zeros(n,1);
        end
        
        function [x,D,w] = clenshawcurtisr_(n) % Clenshaw-Curtis collocation [-1,+1]
            if     n==1; w=2; x=0;
            elseif n==2; w=[1;1]; x=[-1;1];
            else
                N=n-1; c=zeros(n,2);
                c(1:2:n,1)=(2./[1 1-(2:2:N).^2 ])'; c(2,2)=1;
                f=real(ifft([c(1:n,:);c(N:-1:2,:)]));
                w=2*([f(1,1); 2*f(2:N,1); f(n,1)])/2;
                x=N*f(n:-1:1,2);
            end
            D = am_field.get_differentiation_matrix(x);
        end
        
        function [x,D,w] = chebyshevTr_(n,flag) % roots of Chebyshev T (1st kind) (-1,+1)
            if nargin<2; flag=''; end
            f_ = @(n) -cos(pi*(2*[1:n]-1)/(2*n)).';
            if contains(flag,'edge')
                n=n-2; x(:,1)=[-1;f_(n);1]; 
            else
                x(:,1)=f_(n); 
            end
            if nargout < 2; return; end
            D = am_field.get_differentiation_matrix(x); D = cat(3,D,D^2);
            if nargout < 3; return; end
            w(:,1) = am_field.get_integration_weights(x);
        end
        
        function [x,D,w] = chebyshevUr_(n,flag) % roots of Chebyshev U (2nd kind) [-1,+1] 
            if nargin<2; flag=''; end
            f_ = @(n) -cos(pi*[1:n]/(n+1)).';
            if contains(flag,'edge')
                n=n-2; x(:,1)=[-1;f_(n);1]; 
            else
                x(:,1)=f_(n); 
            end
            if nargout < 2; return; end
            D = am_field.get_differentiation_matrix(x); D = cat(3,D,D^2);
            if nargout < 3; return; end
            w(:,1) = am_field.get_integration_weights(x);
        end
        
        function [x,D,w] = legendrer_(n) % roots of Legendre [-1,+1] 
            % Integrates functions between [-1,+1] without integration weight function:
            %   _
            %  /   + 1              __ n
            %  |       f(x) dx  =  \      F  w
            % _/   - 1             /__ i   i  i
            %
            % asciiTeX "\int_{-1}^{+1} f(x) dx = \sum_i^n  F_i w_i"

            % A. N. Lowan, N. Davids, and A. Levenson, Bulletin of the American Mathematical Society 48, 739 (1942).
            A = (1:n-1)./sqrt(4*(1:n-1).^2-1); J = diag(A,-1)+diag(A,1); [V,x] = eig(J,'vector'); [x,fwd]=sort(x(:));
            if nargout < 2; return; end
            D = am_field.get_differentiation_matrix(x); D = cat(3,D,D^2);
            if nargout < 3; return; end
            w(:,1) = 2*V(1,fwd)'.^2;
        end
        
        function [x,D,w] = laguerrer_(n,flag) % roots of Laguerre [0,+Inf]
            % Integrates functions between [0,+Inf] weighed by exp(-x):
            %   _
            %  /  oo   - x              __
            %  |     e     f(x) dx  =  \      F  . w
            % _/  0                    /__ i   i    i
            % 
            % asciiTeX "\int_0^\infty e^{-x} f(x) dx = \sum_i  F_i . w_i"
            
            if nargin<2; flag=''; end
            J_ = @(n) diag(1:2:2*n-1)-diag(1:n-1,1)-diag(1:n-1,-1);
            if contains(flag,'edge')
                n=n-1; [V,D]=eig(J_(n),'vector'); [x,fwd]=sort(D(:)); x=[0;x];
            else
                       [V,D]=eig(J_(n),'vector'); [x,fwd]=sort(D(:));
            end
            if nargout < 2; return; end
            alpha = exp(-x/2); m = 2; % order of differentiation
            beta  = (-0.5).^[1:m].'*ones(1,numel(x));
            D = am_field.get_differentiation_matrix(x, alpha, beta);
            if nargout < 3; return; end
            w(:,1) = V(1,fwd).^2;
            if contains(flag,'edge')
                w = [0;w];
            end
        end
        
        function [x,D,w] = hermiter_(n) % roots of Hermite polynomial (harmonic oscilator solution) [-Inf,+Inf] 
            % Integrates functions between [-Inf,+Inf] weighed by exp(-x^2):
            %   _            2
            %  /  oo      - x               __ n
            %  |        e      f(x) dx  =  \      F  . w
            % _/   - oo                    /__ i   i    i
            % 
            %  asciiTeX "\int_{-\infty}^{\infty} e^{-x^2} f(x) dx = \sum_i^n  F_i . w_i"

            % R. E. Greenwood and J. J. Miller, Bulletin of the American Mathematical Society 54, 765 (1948).
            % George Papazafeiropoulos
            A = sqrt((1:n-1)/2); J = diag(A,1)+diag(A,-1); [V,x]=eig(J,'vector'); [x,fwd]=sort(x(:));
            if nargout < 2; return; end
            alpha = exp(-x.^2/2); beta = [ones(size(x'));-x']; m = 2; % order of differentiation
            for l = 3:m+1; beta(l,:) = -x'.*beta(l-1,:)-(l-2)*beta(l-2,:); end; beta(1,:) = [];
            D = am_field.get_differentiation_matrix(x, alpha, beta);
            if nargout < 3; return; end
            w = sqrt(pi)*V(1,fwd)'.^2;
        end

        function [x,D,w] = fourierr_(n) % roots of fourier function [0,1) 
            x(1:n,1) = [0:(n-1)]/n;
            if nargout < 2; return; end
            D = get_fourier_differentiation_matrix(n);
            if nargout < 3; return; end
            w(1:n,1) = 1/n;
            function [D] = get_fourier_differentiation_matrix(n)
                re = [0,0.5*(-1).^(1:n-1).*cot((1:n-1)*pi/n)]; 
                im = (-1).^(1:n)*sqrt(-1)/2;
                D  = real(2*pi*toeplitz(re+im,-re+im)); 
                % D = D./(2*pi); % use this if x = [0,1); if x = [0,2pi) comment it out.
            end
        end
        
        function [x,D,w] = cdiff_(n) % evenly spaced central difference [0,1)
            x(1:n,1) = [0:n-1]/n;
            if nargout < 2; return; end
            % get first and second derivative
            D = zeros(n,n,2);
            for i = 1:2
                [c,v] = am_field.get_differentiation_weights([-1,0,1],i); nvs = numel(v); m = ceil(nvs/2);
                D(:,:,i) = toeplitz([c(m:-1:1),zeros(1,n-m)],[c(m:end),zeros(1,n-m)])*n.^(i);
            end
            if nargout < 3; return; end
            w(1:n,1) = 1;
        end        

        function [x,D,w] = pdiff_(n) % evenly spaced  periodic central difference [0,1)
            x(1:n,1) = [0:n-1]/n;
            if nargout < 2; return; end
            % get first and second derivative
            D = zeros(n,n,2);
            for i = 1:2
                [c,v] = am_field.get_differentiation_weights([-1,0,1],i); nvs = numel(v); m = ceil(nvs/2);
                D(:,:,i) = am_lib.circulant_(circshift([c,zeros(1,n-nvs)],m-nvs))*n.^(i);
            end
            if nargout < 3; return; end
            w(1:n,1) = 1;
        end        
    
        function [D]     = get_flattened_divergence(Dx,Dy,Dz) 
            switch nargin
                case 1
                    D = Dx;
                case 2
                    n(1) = size(Dx,1); n(2) = size(Dy,1); 
                    D = kron(eye(n(2)),Dx) + ... 
                        kron(Dy,eye(n(1)));
                case 3
                    n(1) = size(Dx,1); n(2) = size(Dy,1); n(3) = size(Dz,1);
                    D = kron(eye(n(3)),kron(eye(n(2)),Dx)) + ... 
                        kron(eye(n(3)),kron(Dy,eye(n(1)))) + ...
                        kron(Dz,kron(eye(n(2)),eye(n(1))));
                otherwise
                    error('not yet implemented');
            end
        end
        
        function [L]     = get_flattened_laplacian(varargin) 
            L = am_field.get_divergence_matrix(varargin{:})^2;
        end

        function [c,x]   = get_differentiation_weights(x,n) 
            % x = collocation points or number of collocation points
            % n = order of differentiation
            % for example, 
            %   x = [-2,0,1,2];
            %   n = 3;
            % for centered space methods:
            %   x = [-1,1]; % x is even
            %   n = 1;
            %
            nxs = numel(x);
            if nxs == 1; N = x; x = [0:(N-1)]-(N-1)/2; nxs = numel(x); end
            if nxs <= n; error('n is not bigger than the number of elements in x'); end
            algo = 1; 
            switch algo
                case 1
                    % This algorithm is numerically more stable, based on the recursion formula in:
                    % B. Fornberg, "Calculation of weights in finite difference formulas", SIAM Review 40 (1998), pp. 685-691.
                    c1 = 1; c4 = x(1); C = zeros(nxs-1,n+1); C(1,1) = 1;
                    for i=1:nxs-1
                        i1 = i+1; mn = min(i,n); c2 = 1; c5 = c4; c4 = x(i1);
                        for j=0:i-1
                            j1 = j+1; c3 = x(i1) - x(j1); c2 = c2*c3;
                            if j==i-1
                                for s=mn:-1:1
                                    s1 = s+1; C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
                                end
                                C(i1,1) = -c1*c5*C(i1-1,1)/c2;
                            end
                            for s=mn:-1:1
                                s1 = s+1; C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
                            end
                            C(j1,1) = c4*C(j1,1)/c3;
                        end
                        c1 = c2;
                    end
                    c = C(:,end).';
                case 2
                    % explicit methodology
                    % Obtained by taylor expanding at each stencil point and setting the sum of coefficients equal to 
                    % all derivatives equal to zero, except to the those derivative for which the sum equals factorial(n).
                    A = [x.^([1:nxs].'-1)]; B = zeros(nxs,1); B(n+1) = factorial(n); c = ( A \ B ).';
            end
        end
        
        function [D]     = get_differentiation_matrix(varargin) 
            % get_differentiation_matrix(x)
            % get_differentiation_matrix(x,order)
            % get_differentiation_matrix(x,alpha,beta)
            
            % select inputs
            switch numel(varargin)
                case 1; x = varargin{1}(:);
                case 2; x = varargin{1}(:); N = numel(x); alpha = ones(N,1);      M = varargin{2}; B = zeros(M,N);
                case 3; x = varargin{1}(:); N = numel(x); alpha = varargin{2}(:); B = varargin{3}; M = size(B,1); 
            end
            
            % select the algorithm
            switch numel(varargin)
                case 1
                    % based on the function by Greg von Winckel
                    % only first derivative is output
                    x=sort(x(:));                       N=length(x); N1=N+1; N2=N*N;
                    X=repmat(x,1,N);                    Xdiff=X-X'+eye(N);
                    W=repmat(1./prod(Xdiff,2),1,N);     D=W./(W'.*Xdiff); 
                    D(1:N1:N2)=1-sum(D);                D=-D';
                case {2,3}
                    % based on the function by Weideman and Reddy
                    % all derivatives up to M-th is output, weight functions B can be utilized
                        L = logical(eye(N));
                       XX = x(:,ones(1,N)); DX = XX-XX'; DX(L) = ones(N,1);
                        c = alpha.*prod(DX,2); C = c(:,ones(1,N)); C = C./C';
                        Z = 1./DX; Z(L) = zeros(N,1); X = Z'; X(L) = [];
                        X = reshape(X,N-1,N); Y = ones(N-1,N); DM = eye(N);

                    for i = 1:M
                        Y = cumsum([B(i,:); i*Y(1:N-1,:).*X]); 
                        DM = i*Z.*(C.*repmat(diag(DM),1,N) - DM);  DM(L) = Y(N,:); 
                        D(:,:,i) = DM;
                    end
            end
        end
  
        function [w]     = get_integration_weights(x) 
            % works only for polynomials integrating between -1 and 1
            % Greg von Winckel
            [x,fwd]=sort(x(:));
            [xl,wl]=lgwt(numel(x)+4,min(x),max(x));
            w(fwd)=lagrange(x,xl)*wl;
            function [x,w,L]=lgwt(N,a,b)
                N=N-1; N1=N+1; N2=N+2;
                y1=cos((2*(0:N)'+1)*pi/(2*N+2)); y=y1; L=zeros(N1,N2); Lp=zeros(N1,N2); y0=2;
                % Iterate until new points are uniformly within epsilon of old points
                while max(abs(y-y0))>eps
                    L(:,1)=1;    Lp(:,1)=0;
                    L(:,2)=y;    Lp(:,2)=1;    
                    for k=2:N1
                        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
                    end
                    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
                    y0=y;
                    y=y0-L(:,N2)./Lp;
                end
                % Linear map from[-1,1] to [a,b]
                x=(a*(1-y)+b*(1+y))/2;      
                % Compute the weights
                w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
            end
            function L=lagrange(x,xl)
                n=length(x); N=length(xl); x=x(:); X=repmat(x,1,n);
                % Compute the weights
                wg=1./prod(X-X'+eye(n),2); xdiff=repmat(xl.',n,1)-repmat(x,1,N);
                % find all the points where the difference is zero
                zerodex=(xdiff==0); 
                % See eq. 3.1 in Ref (1)
                lfun=prod(xdiff,1);
                % kill zeros
                xdiff(zerodex)=eps;
                % Compute lebesgue function
                L=diag(wg)*repmat(lfun,n,1)./xdiff;
            end
        end
        
    end
    
    % electricity/magnetism
    
    methods (Static)
        
        function [A] = get_vector_potential(R,dI,I)
            % get the magnetic vector potential A [3,x,y,z] at positons R [x,y,z] given current flow dI [3,(x,y,z)] at positions I [3,(x,y,z)]
            I = reshape(I,3,[]); dI = reshape(dI,3,[]); M = size(I,2);
            if size(I,2)~=size(dI,2); error('dI and I dimension mismatch'); end
            A = am_field.mu0/(4*pi)*sum(reshape(dI,3,1,1,1,M)./(am_field.tiny+am_lib.normc_(R-reshape(I,3,1,1,1,M))),5);
        end
        
    end
end

