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

            subplot(1,2,1); F.plot_vector_field('F'); axis tight; title('vector potential A'); % plot magnetic vector potential
            line(I(1,:),I(2,:),I(3,:),'linewidth',2,'color',[1 1 1]*0.5); % draw wire
            subplot(1,2,2); F.plot_vector_field('C'); axis tight; title('magnetic field B'); % plot magnetic field (curl of vector potential)
            line(I(1,:),I(2,:),I(3,:),'linewidth',2,'color',[1 1 1]*0.5); % draw wire
            
        end

        function [F] = define_field(n,a,s)
            m = numel(n);
            if numel(s)~=m; error('s dimensions mismatch'); end
            if numel(a)~=m; error('a dimensions mismatch'); end
            F = am_field(); F.a=a; F.n=n; F.s=s; F.d=numel(n);
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
            if strcmp(F.T,'scalar'); F.H = get_hessian(F); end
        end
        
        function [T] = get_field_type(F)
            switch size(F.F,1)
                case {1}  ; T = 'scalar';
                case {2,3}; T = 'vector';
                otherwise; error('unknown field type');
            end
        end

        function [R] = get_collocation_points(F)
            for i = 1:F.d % loop over dimensions
                n = F.n;
                switch F.s{i}
                    case 'chebyshev'; R{i} = am_lib.chebyshevUr_(n(i),'edge');
                    case 'legendre';  R{i} = am_lib.legendrer_(n(i));
                    case 'fourier';   R{i} = am_lib.fourierr_(n(i));
                    case 'cdiff';     R{i} = am_lib.cdiff_(n(i));
                    otherwise; error('unknown s');
                end
                n(i) = 1; R{i} = repmat(permute(F.a(i)*R{i},circshift([1,2,3],i-1)),n);
            end
            R = permute(cat(4,R{:}),[4,1,2,3]);
        end
        
        function [J] = get_jacobian(F)
            % define matmul which supports sparse matrices
            matmul_ = @(A,B) reshape(A*reshape(B,size(B,1),[]),size(A,1),size(B,2),size(B,3),size(B,4),size(B,5));
            % allocate space
            J = zeros([3,3,F.n]);
            for i = 1:F.d % loop over dimensions
                switch F.s{i}
                    case 'chebyshev'; [~,D] = am_lib.chebyshevUr_(F.n(i),'edge');
                    case 'legendre';  [~,D] = am_lib.legendrer_(F.n(i));
                    case 'fourier';   [~,D] = am_lib.fourierr_(F.n(i));
                    case 'cdiff';     [~,D] = am_lib.cdiff_(F.n(i));
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
                    case 'chebyshev'; [~,D] = am_lib.chebyshevUr_(F.n(i),'edge');
                    case 'legendre';  [~,D] = am_lib.legendrer_(F.n(i));
                    case 'fourier';   [~,D] = am_lib.fourierr_(F.n(i));
                    case 'cdiff';     [~,D] = am_lib.cdiff_(F.n(i));
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
        
        function [D] = get_laplacian(F)
            D = am_lib.trace_(F.H);
        end

        function [C] = get_curl(F)
            C = cat(1, F.J(3,2,:,:,:)-F.J(2,3,:,:,:), ...
                       F.J(1,3,:,:,:)-F.J(3,1,:,:,:), ...
                       F.J(2,1,:,:,:)-F.J(1,2,:,:,:));
            C = permute(C,[1,3,4,5,6,7,8,2]);
        end
        
        function [h] = plot_field(F)
            
            sl_ = @(R,i)   squeeze(R(i,:,:,:,:,:));
            
            switch F.T
                case 'scalar'
                    switch F.d
                        case 2
                            set(gcf,'color','w');
                            h = surf(sl_(F.R,1), sl_(F.R,2), squeeze(F.F)); 
                            h.EdgeColor= 'none'; h.LineWidth = 1; view([0 0 1]); daspect([1 1 1]); axis tight;
                        case 3
                            error('not yet implemented');
                        otherwise; error('invalid field dimension');
                    end
                case 'vector'
                    switch F.d
                        case 2
                            set(gcf,'color','w');
                            h = quiver(sl_(F.R,1), sl_(F.R,2), sl_(F.F,1), sl_(F.F,2) ); 
                            h.LineWidth = 1; view([0 0 1]); daspect([1 1 1]); axis tight;
                        case 3
                            error('not yet implemented');
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
        
        function [h] = plot_vector_field(F,field)
            figure(1); set(gcf,'color','w');
            quiver3(F.R(1,:),F.R(2,:),F.R(3,:),...
                    F.(field)(1,:),F.(field)(2,:),F.(field)(3,:),'AutoScaleFactor',6,'ShowArrowHead','off','linewidth',0.5);
            % line(I(1,:),I(2,:),I(3,:),'linewidth',2,'color',[1 1 1]*0.5); axis([-1 +1 -1 +1 -1 +1]);
            streamslice(squeeze(F.R(2,:,:,:)),squeeze(F.R(1,:,:,:)),squeeze(F.R(3,:,:,:)),...
                        squeeze(F.(field)(2,:,:,:)),squeeze(F.(field)(1,:,:,:)),squeeze(F.(field)(3,:,:,:)),...
                        [min(F.R(2,:))],[min(F.R(1,:))],[min(F.R(2,:))]);
            set(gca,'DataAspectRatio',[1,1,1],'CameraPosition',[2,1,1],'Box','on');
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

