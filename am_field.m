classdef am_field

    properties (Constant)
        tiny = 1E-8; 
        mu0  = 1; % vacuum permeability
    end

    properties
        T = []; % field type (scalar/vector)
        Y = []; % coordinate type (cartesian/polar/cylindrical/spherical)
        d = []; % dimensions (2 or 3)
        s = []; % scheme{1:d} ('fourier/chebyshev/legendre/cdiff/discrete') for each dimension
        n = []; % grid points / cart. dimension    3D: [n(1),n(2),n(3)]     2D: [n(1),n(2)] 
        a = []; % lattice/grid spacing             3D: [a(1),a(2),a(3)]     2D: [n(1),n(2)] 
        R = []; % coordinates                          [  x , y , z ] , [ r , th , z ] , [  r , phi , chi ]
        F = []; % field
                % vector field
                %     [ (             F(x)             )                 ]
                %     [ (             F(y)             ) , x , y , z , t ]
                %     [ (             F(z)             )                 ]
                % scalar field
                %     [ (               1              ) , x , y , z , t ]
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
        Q = []; % topological charge
                %     [ (               1              ) , x , y , z  ]
        % AFM
       rACF=[]; % radial autocorrelation function
      rHHCF=[]; % radial height-height function
        % TEM
        E    = []; % (GPA) strain field
        W    = []; % (GPA) rotation field
        Ar_g = []; % (GPA) amplitude
        Pr_g = []; % (GPA) phase
        dPr_g= []; % (GPA) phase derivative
    end

    methods (Static)
        
        function [F] = define(n,a,s) 
            % F = define(n,a,s)
            ndims = sum(n~=1); ndiscretes=sum(strcmp(s,'discrete'));
            if numel(s)~=ndims; error('s dimensions mismatch'); end
            if numel(a)~=ndims; error('a dimensions mismatch'); end
            if ndiscretes>1; error('only one discrete dimension is supported at this time'); end
            F = am_field(); F.a=a; F.n=n; F.s=s; F.d=ndims-ndiscretes; F.R = F.get_collocation_points(); F.Y = 'cartesian';
        end
        
        function [F] = demo()
            
            F   = am_field.define([2,2].^[7,8],[2*pi,2*pi],{'cdiff','fourier'});
            F.F = cat(1,sin(F.R(1,:,:,:)), sin(F.R(2,:,:,:))); 
            F.F = sum(F.F,1);
            F.T = F.get_field_type();
            F.J = F.get_jacobian();
            F.D = F.get_divergence();
            F.C = F.get_curl();
            F.H = F.get_hessian();
            
            figure(1); plot_jacobian(F);
            figure(2); plot_hessian(F);
            
        end

        function [F] = demo_biot_savart()
            
            F   = am_field.define([2,2,2].^[5,5,5],[1,1,1],{'chebyshev','chebyshev','chebyshev'});
            F.R = get_collocation_points(F);

            D_ = @(N)      sparse([1:N],1+[1:N],1,N,N+1)-sparse([1:N],[1:N],1,N,N+1); % Define forward-difference and forward-mean transforms.
            M_ = @(N) 0.5*(sparse([1:N],1+[1:N],1,N,N+1)+sparse([1:N],[1:N],1,N,N+1));

            I_ = {@(d,th,M) [d*sin(2*pi*th) ,d*cos(2*pi*th),zeros(M+1,1)]; % Define a library of current paths: circle, solenoid, and straight wire.
                  @(d,th,M) [d*sin(20*pi*th),d*cos(20*pi*th),linspace(min(F.R(3,:)),max(F.R(3,:)),M+1)'];
                  @(d,th,M) [zeros(M+1,1),zeros(M+1,1),mpgrid(M+1)];};

            iI = 2; r = 0.5; M = 1000; th = [0:M]'/M; % Select and construct a current path.
            dI = (D_(M)*I_{iI}(r,th,M)).'; 
             I = (M_(M)*I_{iI}(r,th,M)).';

            F.F = am_field.get_vector_potential(F.R,dI,I); % get magnetic vector potential 
            F = F.get_derivatives(); % compute derivatives

            subplot(1,2,1); F.plot_field('F'); axis tight; title('vector potential A'); % plot magnetic vector potential
            line(I(1,:),I(2,:),I(3,:),'linewidth',2,'color',[1 1 1]*0.5); % draw wire
            subplot(1,2,2); F.plot_field('C'); axis tight; title('magnetic field B'); % plot magnetic field (curl of vector potential)
            line(I(1,:),I(2,:),I(3,:),'linewidth',2,'color',[1 1 1]*0.5); % draw wire
            
        end

        function [F] = demo_hilliard_cahn()
            F = am_field.define([2,2].^[6,6],[2,2].^[6,6],{'pdiff','pdiff'}); % initialize object
            F.F = rand([1,F.n]); F.F = F.F - mean(F.F(:)); % initialize random field
            F = F.evolve_scalar_field('CH58',{@(x)(x^2-1)^2/4,1},'explicit',0.01,10000); % solve
        end
        
        function [F] = demo_qin_pablo()
            % Cahn-Hilliard Polymer (mean = -0.45: hexagons, mean = 0: lines)
            F = am_field.define([2,2].^7,[2,2].^7,{'pdiff','pdiff'});
            F.F = rand([1,F.n]); F.F = F.F-mean(F.F(:)); F.F=F.F-0.000; % initialize field, lines    (mean = 0)
            F.F = rand([1,F.n]); F.F = F.F-mean(F.F(:)); F.F=F.F-0.455; % initialize field, hexagons (mean = -0.455)
            F = F.evolve_scalar_field('QP13',{@(x)(x^2-1)^2/4,0.5,0.1},'explicit',0.05,10000);
            %
        end
        
        function [F] = demo_complex_ginzburg_landau()
            F = am_field.define([2,2].^[7,7],[2,2].^[7,7],{'pdiff','pdiff'});
            F.F = rand([1,F.n]); F.F = F.F-mean(F.F(:)); % initialize field
            F = F.evolve_scalar_field('GLXX',{1},'explicit',0.1,1000);
        end

        function [F] = demo_poisson_equation()
            F = am_field.define([2,2].^[7,7],[2,2].^[7,7],{'pdiff','pdiff'});
            F.F = zeros([1,F.n]);
            rho = am_lib.gauss_(am_lib.normc_(F.R-[30;30])); % put a charge at (30,30)
            F = F.evolve_scalar_field('poisson',{-rho},'explicit',0.05,5000);
        end
        
        function [F] = demo_parallel_plate_capacitor()
            clear;clc;
            
            F = am_field.define([2,2].^[7,7],[2,2].^[7,7],{'cdiff','cdiff'});
            F.F = zeros([1,F.n]);
            % create boundary conditions
            dirichlet = [];
            bc = find(F.R(1,:,:)==54  & F.R(2,:,:)<100 & F.R(2,:,:)>20 ) ; dirichlet = [dirichlet;[bc,-ones(size(bc))]];
            bc = find(F.R(1,:,:)==74  & F.R(2,:,:)<100 & F.R(2,:,:)>20 ) ; dirichlet = [dirichlet;[bc,+ones(size(bc))]];
            dirichlet = dirichlet.';
            % find steadys tate field
            F = F.evolve_scalar_field('laplace',{},'steadystate',0,0,'dirichlet',dirichlet);
        end
        
        function [F] = demo_metal_inside_field()
            clear;clc
            % initialize (can also use {'cdiff','pdiff'} for infinately long capacitor plates)
            F = am_field.define([2,2].^[7,7],[2 2].^7,{'cdiff','pdiff'});
            F.F = zeros([1,F.n]);
            % create capactior plates
            dirichlet = [];
            bc = find(F.R(1,:,:)==min(F.R(1,:))) ; dirichlet = [dirichlet;[bc,-ones(size(bc))]];
            bc = find(F.R(1,:,:)==max(F.R(1,:))) ; dirichlet = [dirichlet;[bc,+ones(size(bc))]];
            dirichlet = dirichlet.'; plates = dirichlet;
            % find steady state field
            F = F.evolve_scalar_field('laplace',{},'steadystate',0,0,'dirichlet',plates);
            % create conductor [inside conductor: 1) equipotential, 2) electric field, 3) no charge]
            ex_=F.define_mask('square',F.n/2,F.n/8);
            dirichlet = []; mean_ = @(f,ex_) repmat(mean(f(ex_(:))), sum(ex_(:)), 1);
            bc = find(ex_); dirichlet = [dirichlet;[bc,mean_(F.F,ex_)]]; 
            dirichlet = dirichlet.'; metal = dirichlet;
            % find steady state field with the metal
            F = F.evolve_scalar_field('laplace',{},'steadystate',0,0,'dirichlet',[plates,metal]);
        end

        function [F] = demo_dielectric_inside_field()
            clear;clc
            % initialize (can also use {'cdiff','pdiff'} for infinately long capacitor plates)
            F = am_field.define([2,2].^7,[2 2].^7,{'cdiff','pdiff'});
            F.F = zeros([1,F.n]);
            % create capactior plates
            dirichlet = [];
            bc = find(F.R(1,:,:)==min(F.R(1,:))) ; dirichlet = [dirichlet;[bc,-ones(size(bc))]];
            bc = find(F.R(1,:,:)==max(F.R(1,:))) ; dirichlet = [dirichlet;[bc,+ones(size(bc))]];
            dirichlet = dirichlet.'; plates = dirichlet;
            % define dielectric constant
            ex_ = F.define_mask('circle',F.n/2,F.n(1)/8); epsilon = 3*ex_ + 1*(~ex_); 
            % find steady state field
            F = F.evolve_scalar_field('nlaplace',{epsilon},'steadystate',0,0,'dirichlet',plates);
            % F = F.evolve_scalar_field('nlaplace',{epsilon},'explicit',0.05,100000,'dirichlet',plates);
            % get derivatives
            F = F.get_derivatives();
            % plot field
            clf; am_lib.set_plot_defaults_(); hold on; clist=am_lib.colormap_('magma',100); colormap(clist); n=2;
            J = am_lib.diag_(F.J); % J = J./am_lib.normc_(J); J(isnan(J))=0;
            hc = contour(squeeze(F.R(1,:,:)),squeeze(F.R(2,:,:)),squeeze(F.F(1,:,:)),10,'linewidth',2);
            hq = am_lib.quiverc_( ...
                   squeeze(F.R(1,1:n:end,1:n:end)),squeeze(F.R(2,1:n:end,1:n:end)),...
                   squeeze(J(1,1:n:end,1:n:end)),squeeze(J(2,1:n:end,1:n:end)),squeeze(F.F(1,1:n:end,1:n:end)),clist);
            hq.LineWidth=1.2;
            daspect([1 1 1]); axis tight; box on; xticks([]); yticks([]);
        end

        function [F] = demo_dielectric_half_space()
            clear;clc
            % initialize (can also use {'cdiff','pdiff'} for infinately long capacitor plates)
            F = am_field.define([2,2].^7,[2,2].^7,{'cdiff','cdiff'});
            F.F = zeros([1,F.n]);
            % create capactior plates
            dirichlet = [];
            bc = find(F.R(1,:,:)==F.n(1)/2 & F.R(2,:,:)==F.n(2)/2) ; dirichlet = [dirichlet;[bc,1]];
            dirichlet = dirichlet.'; charge = dirichlet;
            % define dielectric constant
            ex_ = F.define_mask('halfspace',F.n(1)/3); epsilon = 30*ex_ + 1*(~ex_); 
            % find steady state field
            F = F.evolve_scalar_field('nlaplace',{epsilon},'steadystate',0,0,'dirichlet',charge);
            % F = F.evolve_scalar_field('nlaplace',{epsilon},'explicit',0.05,100000,'dirichlet',plates);
            % get derivatives
            F = F.get_derivatives();
            % plot field
            clf; am_lib.set_plot_defaults_(); hold on; clist=am_lib.colormap_('magma',100); colormap(clist); n=2;
            J = am_lib.diag_(F.J); % J = J./am_lib.normc_(J); J(isnan(J))=0;
            hc = contour(squeeze(F.R(1,:,:)),squeeze(F.R(2,:,:)),squeeze(F.F(1,:,:)),10,'linewidth',2);
            hq = am_lib.quiverc_( ...
                   squeeze(F.R(1,1:n:end,1:n:end)),squeeze(F.R(2,1:n:end,1:n:end)),...
                   squeeze(J(1,1:n:end,1:n:end)),squeeze(J(2,1:n:end,1:n:end)),squeeze(F.F(1,1:n:end,1:n:end)),clist);
            hq.LineWidth=1.2;
            daspect([1 1 1]); axis tight; box on; xticks([]); yticks([]); 
        end
        
        function [F] = demo_charge_in_dielectric()
            clear;clc
            % initialize (can also use {'cdiff','pdiff'} for infinately long capacitor plates)
            F = am_field.define([2,2].^6,[2,2].^6,{'cdiff','cdiff'}); 
            F.F = zeros([1,F.n]);
            % define dielectric constant
            epsilon = F.define_mask('annulus',F.n/2,[F.n(1)/8,F.n(1)/4]); epsilon = 2.5*epsilon + 1*(~epsilon); 
            % define charge_
            charge  = F.define_mask('point',F.n/2); charge = -10*charge;
            % find steady state field
            F = F.evolve_scalar_field('npoisson',{charge,epsilon},'explicit',0.1,1E5);
            % F = F.evolve_scalar_field('nlaplace',{epsilon},'explicit',0.05,100000,'dirichlet',plates);
            axis tight;
            % get derivatives
            F = F.get_derivatives();
            % plot field
            clf; am_lib.set_plot_defaults_(); hold on; clist=am_lib.colormap_('magma',100); colormap(clist); n=2;
            J = am_lib.diag_(F.J); % J = J./am_lib.normc_(J); J(isnan(J))=0;
            hc = contour(squeeze(F.R(1,:,:)),squeeze(F.R(2,:,:)),squeeze(F.F(1,:,:)),10,'linewidth',2);
            hq = am_lib.quiverc_( ...
                   squeeze(F.R(1,1:n:end,1:n:end)),squeeze(F.R(2,1:n:end,1:n:end)),...
                   squeeze(J(1,1:n:end,1:n:end)),squeeze(J(2,1:n:end,1:n:end)),squeeze(F.F(1,1:n:end,1:n:end)),clist);
            hq.LineWidth=1.2;
            daspect([1 1 1]); axis tight; box on; xticks([]); yticks([]); 
            %
            plot(abs(diff(squeeze(F.F(1,end/2,:)))))
        end
        
        function [F] = demo_heat_diffusion()

            clear;clc
            
            F = am_field.define([2].^[9],1,{'chebyshev'});
            F.F = zeros([1,F.n]);
            % create boundary conditions
            dirichlet = [];
            dirichlet = [dirichlet;[1,10]];
            dirichlet = dirichlet.';
            neumann   = [];
            neumann   = [neumann;[F.n,0]];
            neumann   = neumann.';

            F = F.evolve_scalar_field('dissipative_diffusion',{8,5},'steadystate',0,0,'dirichlet',dirichlet,'neumann',neumann);

        end
        
        function [F] = demo_vortex_unbinding()
            % make animation of vortex pair separation
            % m = winding number
            % n = phase offset
            % s = separation
            
            clear;clc;

            coneplot2_ = @(x,y,u,v,c) ...
                coneplot(cat(3,x,x),cat(3,y,y),3*cat(3,-ones(size(x)),ones(size(x))), ...
                         cat(3,u,u),cat(3,v,v),3*cat(3,-ones(size(x)),ones(size(x))), ...
                         (x-u/2)*0.90,(y-v/2)*0.90,rescale(c,-1,1)*0.5, cat(3,c,c) );

            n = 4;
            F = am_field.define([2,2].^[n,n+1],[2,2].^[n,n+1],{'cdiff','cdiff'});
            F.R = F.R - [2;2].^[n-1;n];

            for i = [1:200]

                m = +1; n = pi/2; s = (i-1)/10;
                th_ = @(s) atan2(F.R(2,:,:)-s+0.5,F.R(1,:,:)+0.5);
                F.F = cat(1,sin(m*th_(s)-m*th_(-s)+n),cos(m*th_(s)-m*th_(-s)+n));
                F = F.get_topological_charge();

                % plot
                clist = am_lib.clight_(am_lib.colormap_('hsv',256),-0.2);
                u = squeeze(F.F(1,:,:)); x = squeeze(F.R(1,:,:));
                v = squeeze(F.F(2,:,:)); y = squeeze(F.R(2,:,:));
                q = squeeze(F.Q); w = am_lib.gaussw_(21,5)*am_lib.gaussw_(21,5).';
                q = conv2( q, w, 'same'); th = mod(atan2(u,v)+pi/2,2*pi); th(1) = 0; th(end)=2*pi;

                switch 'fancy'
                    case 'simple'
                        if i == 1
                            am_lib.set_plot_defaults_(); hold on; offset_ = 100;
                            hs = surf(x,y,q-offset_,'edgecolor','none');
                            hq = am_lib.quiverc_(x-u/2,y-v/2,u,v,th,clist);
                            hq.AutoScaleFactor=0.5; hq.LineWidth = 0.8;
                            shading interp; box on; axis tight; view([0 0 1]); daspect([1 1 1]);
                            yticks([]); xticks([]); axis([-F.n(1)/2 F.n(1)/2 -F.n(2)/2 F.n(2)/2]-1/2);
                            caxis([-1,1]*3-offset_); colormap(am_lib.clight_(flipud(am_lib.colormap_('red2blue',100)),0.5));
                            set(gca, 'XLimMode', 'manual', 'YLimMode', 'manual');
                        else
                            am_lib.quiverc_(hq,th,clist);
                            hq.XData = x-u/2; hq.UData = u;
                            hq.YData = y-v/2; hq.VData = v;
                            hs.CData = q-offset_;
                        end
                case 'fancy'
                        cla;
                        am_lib.set_plot_defaults_(); offset_ = 100;
                        hc = coneplot2_(x.',y.',u.',v.',q.');
                        hc.EdgeColor='none'; hc.DiffuseStrength = 0.1; hc.AmbientStrength=0.6;
                        shading interp; box on; axis tight; view([0 0 1]); daspect([1 1 1]);
                        yticks([]); xticks([]); axis([-F.n(1)/2 F.n(1)/2 -F.n(2)/2 F.n(2)/2]*0.92-1/2);
                        camlight(45,-5); lighting gouraud;
                        caxis([-2 2]/2); colormap(flipud(am_lib.colormap_('red2blue',100)));
                        box on; % set(gcf,'renderer','painters');
                end
                
                drawnow();

                % get frame
                f(i) = getframe();

            end

            % save movie going backwards and forwards
            mov=f([1:200,199:-1:1]);
            
            % write movie
            h = VideoWriter('vortex_fancy.mp4','MPEG-4'); h.FrameRate=60; h.Quality=100; 
            open(h); for i = 1:numel(mov); writeVideo(h,mov(i)); end; close(h);
 
        end
        
        function [F] = demo_ising_metropolis_periodic_bc()
            F = am_field.define([2,2].^7,[2,2].^7,{'pdiff','pdiff'});
            % F.F = 2*round(rand([1,F.n]))-1; % initialize binary field
            F.F = ones([1,F.n]); % initialize the field

            M = 50000; % monte carlo samples
            
            % compute energy difference of flipping ith spin
            % first neighbor
            en1_ = @(F,x)  2 * F.F(1,x(1),x(2)) .* ( ...
                              F.F(1,mod(x(1)-1-1,F.n(1))+1,x(2)) + ...
                              F.F(1,mod(x(1)+1-1,F.n(1))+1,x(2)) + ...
                              F.F(1,x(1),mod(x(2)-1-1,F.n(2))+1) + ...
                              F.F(1,x(1),mod(x(2)+1-1,F.n(2))+1) );
            % second neighbor
            en2_ = @(F,x) 2 * F.F(1,x(1),x(2)) .* ( ...
                              F.F(1,mod(x(1)-1-1,F.n(1))+1,mod(x(2)-1-1,F.n(2))+1) + ...
                              F.F(1,mod(x(1)+1-1,F.n(1))+1,mod(x(2)-1-1,F.n(2))+1) + ...
                              F.F(1,mod(x(1)-1-1,F.n(1))+1,mod(x(2)+1-1,F.n(2))+1) + ...
                              F.F(1,mod(x(1)+1-1,F.n(1))+1,mod(x(2)+1-1,F.n(2))+1) );

            J1 = 1;
            J2 = 0;
            
            kT_list = [1.3:0.1:3.5]; % meV
            for k = 1:numel(kT_list)
                kT = kT_list(k);
            % monte carlo samples
            for j = 1:M
                % pick a random coordinate to flip
                for i = 1:F.d; x(i) = randi(F.n(i)); end
                % calculate energy of flipping
                E = J1*en1_(F,x) + J2*en2_(F,x);
                % flip?
                if E <= 0 || rand(1) <= exp(-E/kT)
                    F.F(1,x(1),x(2)) = - F.F(1,x(1),x(2));
                end
                % plot?
                if mod(j,2000)==0; F.plot_field('F'); drawnow; end
            end
                m(k) = mean(F.F(:));
                s(k) = std(F.F(:));
            end
            
            clf; h = semilogx(kT_list,s,'s-',kT_list,m,'o-'); 
            line([1 1]*2/log(1+sqrt(2)),get(gca,'YLim'),'Color','k','linewidth',0.5); % ideal Tc
            
        end

        function [F] = demo_ising_metropolis_helical_bc()
            %isng model with helical boundary conditions
            clear;clc;

            F = am_field.define([2,2].^7,[2,2].^7,{'pdiff','pdiff'});
            % F.F = 2*round(rand([1,F.n]))-1; % initialize binary field
            F.F = ones([1,F.n]); % initialize the field

            M = 100000; % monte carlo samples

            % neighbors (helical boundary condition)
            p = cumprod(F.n);
            n_ = @(i) [mod(i   -1-1,p(end))+1;  % up
                       mod(i   +1-1,p(end))+1;  % down
                       mod(i-p(1)-1,p(end))+1;  % left
                       mod(i+p(1)-1,p(end))+1]; % right
            en_= @(F,i) 2 * F.F(1,i) * sum(F.F(1,n_(i)));

            k=0; kT_list = [1:0.3:3.5]; % initalize and loop over temperatures
            for kT = kT_list
                % run monte carlo
                for j = 1:M
                    % pick a random coordinate to flip
                    i = randi(p(end));
                    % calculate energy of flipping
                    E = en_(F,i);
                    % flip?
                    if E <= 0 || rand(1) <= exp(-E/kT)
                        F.F(1,i) = - F.F(1,i);
                    end
                    % plot?
                    if mod(j,2000)==0; F.plot_field('F'); drawnow; end
                end
                % record statistics
                k=k+1;
                m(k) = mean(F.F(:));
                s(k) = std(F.F(:));
            end

            clf; h = semilogx(kT_list,s,'s-',kT_list,m,'o-'); 
            line([1 1]*2/log(1+sqrt(2)),get(gca,'YLim'),'Color','k','linewidth',0.5); % ideal Tc 
        end

        function [F] = demo_nudge_elastic_band()
            clear;clc;
            % define central-difference mesh between [0,1) 
            F = am_field.define([2,2].^8,[1,1],{'cdiff','cdiff'}); 
            % scale x to [-1.5,1.2) and y to [-0.2,2.0)
            s = [ -1.5,1.2; -0.2,2.0 ]; F.R = F.R.*(s(:,2)-s(:,1))+s(:,1);
            % load Mueller potential
            F = F.load_model('Mueller');
            F = F.get_derivatives();
            % get minimum energy path using nudge elastic band
            vi = [-1 ;0.5]; vf = [0.7;0.5]; nnodes = 10; ediff = 1E-3; flag = '';
            get_minimum_energy_path(F,vi,vf,nnodes,ediff,flag)
        end

        function [F] = demo_gpa()
            fname = '11_SnO_FS20nm_ADF_lowI30mrad_stack.dm3'; % good
            % fname = '10_SnO_FS10nm_ADF_lowI30mrad_stack.dm3'; % defect
            % fname = '08_SnO_FS10nm_ADF_lowI30mrad_stack.dm3'; % defect
            [d] = DM3Import(fname); % high mag
            d=d.image_data;
            % correct drift and average over stack
            [d,~] = am_field.correct_stack_drift(d);
            % 
            [~,th] = am_lib.imhorizonlevel_(d,'hough');
            am_lib.imagesc_(flipud(dd));colormap('gray');
            [E,W,Ar_g,Pr_g,dPr_g] = am_field.gpa(d,'rot:field',th);
            [h,ax] = am_lib.overlay_(d, squeeze(E(2,2,:,:)), 5);%[-0.4:0.05:0.4]+0.025)
            % am_lib.imagesc_(squeeze(Pr_g{1})) 
        end

        function [F] = load_image(A,dx,dy)
            % load image as scalar field
            if nargin == 1; dx=1;dy=1; end
            switch ndims(A)
                case 3 % stack
                    error('not yet implemented');
                    F = am_field.define(size(A),[dx,dy,1],{'fourier','fourier'});
                case 2 % image
                    F = am_field.define(size(A),[dx,dy],{'fourier','fourier'});
            end
            F.F = permute(A,[9,1:8]);
        end

    end

    methods 

        function [F] = set(F,varargin) % set field properties 
            % parse remaining input
            nvs = numel(varargin);
            if ~am_lib.iseven_(nvs); error('{key, value} inputs required'); end; i=1; 
            while true; F.(varargin{i}) = varargin{i+1}; i=i+2; if i>nvs; break; end; end
        end

        function [F] = load_model(F,model)
            switch model
                case 'Mueller'
                    % parameters in Mueller potential
                    a = [-1 -1 -6.5 0.7]; b = [0 0 11 0.6]; c = [-10 -10 -6.5 0.7]; 
                    A = [-200 -100 -170 15]; X = [1 0 -0.5 -1]; Y = [0 0.5 1.5 1];
                    F.F = sum( A(:).*exp(a(:).*(F.R(1,:)-X(:)).^2 + ...
                                         b(:).*(F.R(1,:)-X(:)).*(F.R(2,:)-Y(:)) + ...
                                         c(:).*(F.R(2,:)-Y(:)).^2), 1);
                    F.F = min(F.F,200); 
                    F.F = reshape(F.F,[1,F.n]);
                case 'otherwise'
                    error('unknown model');
            end
        end
        
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
            % trim coordinates, field, jacobian, hessian, divergence, curl
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

        function ex_ = define_mask(F,shape,varargin)
            switch shape
                case 'square'
                    % (center, width)
                    [c,w]=deal(varargin{:});
                    if F.d==1; ex_ = true(1,F.n); else; ex_ = true(F.n); end
                    for i = 1:F.d
                        ex_(ex_) = F.R(i,ex_)<(c(i)+w(i)) & F.R(i,ex_)>(c(i)-w(i));
                    end
                case 'circle'
                    % (center, radius)
                    [c,r]=deal(varargin{:});
                    ex_ = am_lib.normc_(F.R-c(:)) < r;
                case {'annulus','ring'}
                    % (center, inner and outer radii)
                    [c,r]=deal(varargin{:});
                    R = am_lib.normc_(F.R-c(:));
                    ex_ = R < r(2) & R > r(1);
                case {'point'}
                    % (r)
                    r = deal(varargin{:});
                    if F.d==1; ex_ = true(1,F.n); else; ex_ = true(F.n); end
                    for i = 1:F.d; ex_(ex_) = F.R(i,ex_)==r(i); end
                case 'halfspace'
                    % (edge)
                    [e]=deal(varargin{:});
                    ex_ = F.R(1,:,:) < e;
                otherwise 
                    error('unknown shape');
            end
        end

        
%         function [F] = crop(F,x,dx)
%             % x = coordinates to begin croping
%             % dx= crop size
%                   
%             % trim dimension
%             F.d   = F.d-1;
%             % trim scheme, grid size, grid spacing
%             for f = {'s','n','a'}; if ~isempty(F.(f{:}))
%                 F.(f{:})(i) = []; 
%             end; end
%             % trim coordinates, field, jacobian, hessian, divergence, curl
%             for f = {'R','F','J','H','D','C'}; if ~isempty(F.(f{:}))
%                 switch f{:}
%                     case 'R'; o=1; case 'F'; o=1; case 'J'; o=2;
%                     case 'H'; o=2; case 'D'; o=1; case 'C'; o=1;
%                 end
%                 % permute
%                 p=[1:ndims(F.(f{:}))]; p([1,i+o])=p([i+o,1]);     
%                 % trim and permute back to orginal
%                 F.(f{:}) = permute(F.(f{:}),p); F.(f{:}) = permute(F.(f{:})(j,:,:,:,:,:,:,:),p);
%                 % permute trimmed to end
%                 p=[1:ndims(F.(f{:}))+1]; p([i+o,end])=p([end,i+o]); F.(f{:}) = permute(F.(f{:}),p);
%                 switch f{:}
%                     case 'R'; F.(f{:})(i,:,:,:,:,:) = []; 
%                     case 'F'; F.(f{:})(i,:,:,:,:,:) = []; 
%                 end
%             end; end
%         end

        % differentiation
        
        function [F] = get_derivatives(F)
            F.T = F.get_field_type();
            F.J = F.get_jacobian();
            F.Q = F.get_topological_charge();
            F.D = F.get_divergence();
            F.C = F.get_curl();
            if contains(F.T,'scalar')
                F.H = F.get_hessian(); 
                F.L = F.get_laplacian();
            end
        end
        
        function [F] = evolve_scalar_field(F,equation,x,algorithm,dt,M,varargin)
            % examples
            % F = am_field.define([2,2].^[6,6],[2,2].^[6,6],{'pdiff','pdiff'}); % initialize field
            % F.F = rand([1,F.n]); F.F = F.F-mean(F.F(:)); % initialize field
            % F = F.evolve_scalar_field('SW76',{0.1,1.0},'implicit',2,500);         % hexagons
            % F = F.evolve_scalar_field('SW76',{0.3,0.0},'implicit',2,500);         % finger prints
            % F = F.evolve_scalar_field('SW76',{0.8,1.85},'explicit',0.03,5000);    % maze
            % F = F.evolve_scalar_field('CH58',{@(x)(x^2-1)^2/4,1},'explicit',0.01,10000);
            % F = F.evolve_scalar_field('CH58',{@(x)(x^2-1)^2/4,1},'implicit',1,1000);
            % F = F.evolve_scalar_field('GL50',{@(x)(x^2-1)^2/4,1},'explicit',0.01,5000);
            %
            
            import am_field.*
            
            % parse boundary conditions if present
                p = inputParser; 
                checkdbc_ = @(x) isempty(x) || (isnumeric(x) && size(x,1) == 2); 
                checknbc_ = @(x) isempty(x) || (isnumeric(x) && size(x,1) == F.d+1);
                addParameter(p,'dirichlet',[],checkdbc_);
                addParameter(p,'neumann'  ,[],checknbc_);
                parse(p,varargin{:});
                dirichlet = p.Results.dirichlet; % dirichlet( (index, value), n)
                neumann   = p.Results.neumann;   %   neumann( (index, value), n)
            
            N=prod(F.n); % Simplify notation: get total number of grid points.

            [D,L,G] = F.get_flattened_differentiation_matrices(); I = speye(N); % Get flattened divergence, laplacian, gradient operators
            
			switch algorithm
            case {'explicit','implicit','steadystate','crank-nicolson','jacobi','gauss-seidel'}
            	if any(strcmp(algorithm,{'explicit','crank-nicolson'}))
	                % equations of the form F(n+1) = F(n) + dt * LHS
                    %                       LHS = ( F(n+1) - F(n)) / dt = dF/dt
                    switch equation
                        case 'SW76' % Swift-Hohenberg (PRA 1976), x = {eps, g1}
                            LHSe_ = @(U,x) x{1}*U(:) - (L+I)^2*U(:) + x{2}*U(:).^2 - U(:).^3; nargs=2;
                        case 'LP97' % Lifshitz-Petrich (PRL 1997), x = {e, g1, q}
                            LHSe_ = @(U,x) x{1}*U(:) - (L+I)^2*(L+x{3}.^2*speye(N))^2*U(:) + x{2}*U(:).^2 - U(:).^3; nargs=3;
                        case 'GLXX' % Complex Ginzburg-Landau (Kramer & Aranson, Rev Mod Phys 2002; Arason & Tang PRL 1998)
                            LHSe_ = @(U,x) ( I*(1-1i*x{1}) + L - spdiags((1-1i*x{1})*abs(U(:)).^2,0,N,N) )*U(:); nargs=1;
                        case 'GL50' % Ginzburg-Landau (Zh. Eksp. Teor. Fiz. 1950), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(diff(x{1}(z),z)); % P.E. derivative
                            LHSe_ = @(U,x)   ( x{1}(U(:)) + x{2}*L*U(:) ); nargs=2;
                        case 'CH58' % Cahn-Hilliard (J. Chem. Phys. 1958), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(diff(x{1}(z),z)); % P.E. derivative
                            LHSe_ = @(U,x) L*( x{1}(U(:)) - x{2}*L*U(:) ); nargs=2;
                        case 'QP13' % Qin-Pablo (Soft Matter, 2013, 9, 11467)
                            syms z; x{1} = matlabFunction(diff(x{1}(z),z)); % P.E. derivative
                            LHSe_ = @(U,x) L*( x{1}(U(:)) - x{2}*L*U(:) ) - x{3}*(U(:) - mean(U(:))); nargs=3;
                        case 'LS91' % Lai-das Sarma (PRL 1991)
                            % NEED TO TEST
                            LHSe_ = @(U,x) -x{1}*L^2*U(:) + x{2}*L*(D*U(:))^2 + x{3}*randn(size(U)); nargs=2;
                        case 'dissipative_diffusion' % Diffusion equation with dissipation, x = { diffusivity, dissipation }
                            OP = ( x{1}*L - x{2}*I );
                            LHSe_ = @(U,x) OP * U(:); nargs=2;
                        case 'poisson' % Poisson equation, x = { charge density }
                            LHSe_ = @(U,x) L * U(:) - x{1}(:); nargs=1;
                        case 'laplace' % Laplace equation, x = { }
                            LHSe_ = @(U,x) L * U(:); nargs=0;
                        case 'npoisson' % Poisson equation with a spatial-dependent dielectric constant, x = { charge density, dielectric constant }
                            % seems ok
                            GE = cellfun(@(G) spdiags(G*x{2}(:),0,N,N), G, 'UniformOutput', false); 
                            OP = ( spdiags(x{2}(:),0,N,N) * L  +  cat(2,GE{:}) * cat(1,G{:}) );
                            LHSe_ = @(U,x) OP * U(:) - x{1}(:); nargs=2;
                        case 'nlaplace' % Laplace equation with a spatial-dependent dielectric constant, x = { dielectric constant }
                            GE = cellfun(@(G) spdiags(G*x{1}(:),0,N,N), G, 'UniformOutput', false); 
                            OP = ( spdiags(x{1}(:),0,N,N) * L  +  cat(2,GE{:}) * cat(1,G{:}) );
                            LHSe_ = @(U,x) OP * U(:); nargs=1;
                        otherwise
                            % Ingredients available for designing custom equations:
                            % U (field), L (laplacian), D (divergence), G{i} (gradient along i), I (identity)
                            % example: F = F.evolve_scalar_field('8*L-5*eye(N)',{},'steadystate',0,0,'dirichlet',dirichlet,'neumann',neumann);
                            eval(['LHSe_ = @(U,x) ',equation,';']); nargs=0;
                    end
                end
                if any(strcmp(algorithm,{'implicit','crank-nicolson','steadystate'}))
                    % equations of the form F(n+1) = (1 - dt*LHS)\F(n)
                    %                       LHS * F(n+1) = ( F(n+1) - F(n) ) / dt
                    %                       i.e., LHS has a factor of F(n+1) removed
                    switch equation
                        case 'SW76' % Swift-Hohenberg (PRA 1976), x = {eps, g1}
                            LHSi_ = @(U,x) x{1}*I - (L+I)^2 + x{2}*spdiags(U(:),0,N,N) - spdiags(U(:).^2,0,N,N); nargs=2;
                        case 'LP97' % Lifshitz-Petrich (PRL 1997), x = {e, g1, q}
                            LHSi_ = @(U,x) x{1}*I - (L+I)^2*(L+x{3}.^2*speye(N))^2 + x{2}*spdiags(U(:),0,N,N) - spdiags(U(:).^2,0,N,N); nargs=3;
                        case 'GLXX' % Complex Ginzburg-Landau (Kramer & Aranson, Rev Mod Phys 2002)
                            LHSe_ = @(U,x) ( speye(N)*(1-1i*x{1}) + L - spdiags((1-1i*x{1})*abs(U(:)).^2,0,N,N) ); nargs=1;
                        case 'GL50' % Ginzburg-Landau (Zh. Eksp. Teor. Fiz. 1950), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(expand(diff(x{1}(z),z)/z)); % P.E. derivative with field factored out
                            LHSi_ = @(U,x)   ( spdiags(x{1}(U(:)),0,N,N) + x{2}*L ); nargs=2;
                        case 'CH58' % Cahn-Hilliard (J. Chem. Phys. 1958), x = {P.E., gamma^2}
                            syms z; x{1} = matlabFunction(expand(diff(x{1}(z),z)/z)); % P.E. derivative with field factored out
                            LHSi_ = @(U,x) L*( spdiags(x{1}(U(:)),0,N,N) - x{2}*L ); nargs=2;
                        case 'QP13' % Qin-Pablo (Soft Matter, 2013, 9, 11467)
                            % Something is wrong about this implicit implementation; results do no match explicit version. Compare:
                            % F = F.evolve_scalar_field('QP13',{@(x)(x^2-1)^2/4,1,0.05},'implicit',1,100);
                            % F = F.evolve_scalar_field('QP13',{@(x)(x^2-1)^2/4,1,0.05},'explicit',0.01,10000);
                            error('something is wrong with this implementation');
                            syms z; x{1} = matlabFunction(expand(diff(x{1}(z),z)/z)); % P.E. derivative with field factored out
                            LHSi_ = @(U,x) L*( spdiags(x{1}(U(:)),0,N,N) - x{2}*L ) - x{3}*spdiags(U(:) - mean(U(:)),0,N,N); nargs=3;
                        case 'poisson' % Poisson equation, x = { charge density }
                            error('The Poisson equation cannot be relaxed using an implicit method.');
                            % The Poisson equation cannot be relaxed using an implicit method because
                            % the field needs to be factored out from the charge.
                            LHSi_ = @(U,x) L - spdiags(x{1}(:)./U(:),0,N,N); nargs=1;
                        case 'dissipative_diffusion' % x = {diffusivity, dissipation}
                            LHSi_ = @(U,x) x{1}*L - x{2}*speye(N); nargs=2;
                        case 'laplace'  % Laplace equation, x = { }
                            LHSi_ = @(U,x) L; nargs=0;
                        case 'nlaplace' % Laplace equation with a spatial-dependent dielectric constant, x = { dielectric constant }
                            GE = cellfun(@(G) spdiags(G*x{1}(:),0,N,N), G, 'UniformOutput', false); 
                            OP = ( spdiags(x{1}(:),0,N,N) * L  +  cat(2,GE{:}) * cat(1,G{:}) );
                            LHSi_ = @(U,x) OP; nargs=1;
                        otherwise
                            % Ingredients available for designing custom equations:
                            % U (field), L (laplacian), D (divergence), G{i} (gradient along i), I (identity)
                            % example: F = F.evolve_scalar_field('8*L-5*eye(N)',{},'steadystate',0,0,'dirichlet',dirichlet,'neumann',neumann);
                            eval(['LHSi_ = @(U,x) ',equation,';']); nargs = 0;
                    end
                end
            otherwise
                error('unknown solver');
            end

            if ~iscell(x); x=num2cell(x); end % convert array to cells if necessary (cells allow passing functions)
            if numel(x)~=nargs; error('incorrect number of parameters'); end % check parameters
            if isempty(F.F); error('field has not been initialized'); end % check parameters

            % Propagate in time 
            switch algorithm
                case 'steadystate' % dF/dt = 0
                    % initialize system of equations
                    n = prod(F.n); V = zeros(n,1); O = LHSi_(F.F,x);
                    % add dirichlet boundary conditions
                    if ~isempty(dirichlet)
                        % make the b.c. robust by removing coupling to everything else
                        O(dirichlet(1,:),:) = [];
                        V(dirichlet(1,:),:) = [];
                        % add boundary conditions
                        Op = speye(n);
                        Op = Op(dirichlet(1,:),:);
                        Vp = dirichlet(2,:).';
                        % augment
                        O = [O;Op]; V = [V;Vp];
                    end
                    % add neuamnn boundary conditions for each dimension
                    if ~isempty(neumann); for j = 1:F.d
                        Op = G{j}(neumann(1,:),:);
                        Vp = neumann(j+1,:).';
                        % augment
                        O = [O;Op]; V = [V;Vp];
                    end; end
                    % evaluate
                    F.F(1:N) = O\V;
                    % plot
                    F.plot_field('F');
                case 'explicit' % unstable for large steps
                    for i = [1:M]
                        UP = F.F(:); F.F(1:N) = F.F(:) + dt*LHSe_(F.F,x);
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0
                            F.plot_field('F'); ediff = norm(UP(:)-F.F(:));
                            title(sprintf('%i, %0.5g',i,ediff)); drawnow; 
                            if ediff<F.tiny; break; end
                        end
                        % b.c.
                        if ~isempty(dirichlet); F.F(:,dirichlet(1,:)) = dirichlet(2:end,:); end
                        if ~isempty(neumann); error('not yet implemented'); end
                    end
                case 'implicit' % more stable
                    for i = [1:M]
                        UP = F.F(:); F.F(1:N) = (speye(N) - dt*LHSi_(F.F,x))\F.F(:);
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0
                            F.plot_field('F'); ediff = norm(UP(:)-F.F(:));
                            title(sprintf('%i, %0.5g',i,ediff)); drawnow; 
                            if ediff<F.tiny; break; end
                        end
                        % b.c.
                        if ~isempty(dirichlet); F.F(:,dirichlet(1,:)) = dirichlet(2:end,:); end
                        if ~isempty(neumann); error('not yet implemented'); end
                    end
                case 'crank-nicolson'
                    for i = [1:M]
                        UP = F.F(:); F.F(1:N) = ( (speye(N)-dt*LHSi_(F.F,x))\F.F(:) + F.F(:)+dt*LHSe_(F.F,x) )/2;
                        if any(isnan(F.F(:))); warning('NaN'); break; end
                        if mod(i,round(M/100))==0
                            F.plot_field('F'); ediff = norm(UP(:)-F.F(:));
                            title(sprintf('%i, %0.5g',i,ediff)); drawnow; 
                            if ediff<F.tiny; break; end
                        end
                        % b.c.
                        if ~isempty(dirichlet); F.F(:,dirichlet(1,:)) = dirichlet(2:end,:); end
                        if ~isempty(neumann); error('not yet implemented'); end
                    end
                otherwise
                    error('unknown solver');
            end
        end

        function [F] = get_minimum_energy_path(F,vi,vf,nnodes,ediff,flag)
            % nudge_elastic_band
            % nnodes = 10; 
            % % starting and ending positions
            % vi = [-1 ;0.5]; vf = [0.7;0.5];
            % vi = [-.52 ;1.5]; vf = [0.6;0.5];
            
            % build nodes
            r = am_lib.linspacen_(vi, vf, nnodes);

            % define spring constant, time step
            k = 10; dt = 0.01;

            % store force (negative gradient) in memory
            J = -am_lib.diag_(F.J); 
            
            % get local effective mass
            if contains(flag,'effectivemass'); M = 1./am_lib.normc_(am_lib.diag_(F.H)); M = max(M./max(M(:)),0.1); end

            % initialize
            if nargout==0 || contains(flag,'draw'); clf; F.plot_field('F'); hold on; end
            nsteps=1000;  mass(1:nnodes)=1; v=zeros(F.d,nnodes,2); a=zeros(F.d,nnodes,2); PE=Inf(2,1);
            % loop over md
            for j = 1:nsteps
                if contains(flag,'dropedge') % drop points outside boundary?
                    ex_ = true(1,nnodes);
                    ex_(ex_) = (min(F.R(1,:))<r(1,ex_) & max(F.R(1,:))>r(1,ex_));
                    ex_(ex_) = (min(F.R(2,:))<r(2,ex_) & max(F.R(2,:))>r(2,ex_));
                    nnodes=sum(ex_); r = r(:,ex_); mass=mass(ex_); v = v(:,ex_,:); a=a(:,ex_,:);
                end
                if contains(flag,'droptop') % drop points ontop of each other?
                    ex_ = false(1,nnodes); [~,i]=am_lib.uniquec_(r); ex_(i) = true;
                    nnodes=sum(ex_); r = r(:,ex_); mass=mass(ex_); v = v(:,ex_,:); a=a(:,ex_,:);
                end
                % find index of closest mesh point to each node
                i = knnsearch(F.R(:,:).',r(:,:,1).').'; r(1:F.d,:,1) = F.R(:,i);
                % get gradient
                G = J(1:F.d,i);
                % update mass based on local curvature?
                if contains(flag,'effectivemass'); mass(1:nnodes) = M(i); end
                % get potential and kinetic energies per node
                PE(2) = sum(F.F(i)-min(F.F(:)))./nnodes;
                KE(2) = am_lib.sum_(mass.*v(:,:,2).^2,[1,2])/2./nnodes;

                % get tangents
                t = conv2(r(:,:,1),-[-1,0,1],'same'); t=t./am_lib.normc_(t); t(:,[1,end])=0;
                % get perpendicular forces
                f_perp = G - sum(G.*t,1).*t; f_perp(:,[1,end]) = 0;
                % get parallel (spring) forces: f_para = k ( |r_{i+1}?ri|?|ri?r_{i?1}| ) . t
                f_para = k*[0,diff(abs(am_lib.normc_(diff(r(:,:,1),1,2))),1,2),0].*t;
                % get total force
                f_total= f_para + f_perp;
                % relax endpoints?
                if contains(flag,'relax_ends')
                    f_total(:,[1,nnodes]) = G(1:F.d,[1,nnodes]);
                end
                % get force = ma on each point as the negative of the gradient
                a(:,:,2) = f_total./mass;
                % integrate
                integrator = 'sd';
                switch integrator
                    case {'vvhs','velocity-verlet-half-step'}
                        v(:,:,2) = v(:,:,1) + a(:,:,1)/2 .* dt;
                        r(:,:,2) = r(:,:,1) + v(:,:,2) .* dt;
                        v(:,:,2) = v(:,:,1) + a(:,:,2)/2 .* dt;
                    case {'vvd','velocity-verlet-damped'}
                        r(:,:,2) = r(:,:,1) + v(:,:,1) .* dt + a(:,:,1)/2 .* dt.^2;
                        v(:,:,2) = v(:,:,1) + (a(:,:,1)+a(:,:,2))/2 .* dt;
                        ex_=sum(v(:,:,2).*f_total,1)<0; v(:,ex_,2)=0;
                    case {'vv','velocity-verlet'}
                        r(:,:,2) = r(:,:,1) + v(:,:,1) .* dt + a(:,:,1)/2 .* dt.^2;
                        v(:,:,2) = v(:,:,1) + (a(:,:,1)+a(:,:,2))/2 .* dt;
                    case {'sd','steepest-descent'}
                        r(:,:,2) = r(:,:,1) + v(:,:,2) .* dt + a(:,:,2)/2 .* dt.^2;
                        v(:,:,2) = 0;
                    otherwise
                        error('unknown integrator');
                end
                % update
                r = circshift(r,-1,3); v = circshift(v,-1,3); PE = circshift(PE,-1,1); KE = circshift(KE,-1,1);
                % draw
                if nargout==0 || contains(flag,'draw')
                    % print 
                    dPE = abs(PE(2)-PE(1)); dKE = abs(KE(2)-KE(1));
                    if j==1; fprintf('%10s %10s %10s %10s %10s\n','TE','PE','dPE','KE','dKE'); 
                    else;    fprintf('%10.5g%10.5g %10.5g %10.5g %10.5g\n',PE(1)+KE(1),PE(1),dPE,KE(1),dKE); end
                    % plot
                    hold on; h = plot3(r(1,:,1),r(2,:,1),ones(size(r(2,:,1)))*1E10,'Marker','o','Color','k');
                    % h=plot3(r(1,:,1),r(2,:,1),ones(size(r(2,:,1)))*100,'.','Color','r');
                    drawnow; delete(h);
                end
                % check energy difference and break if smaller
                if j>2; if dPE<ediff; break; end; end
            end
            hold on; h = plot3(r(1,:,1),r(2,:,1),ones(size(r(2,:,1)))*1E10,'Marker','o','Color','k');
        end

        function [r] = relax_points(F,r,ediff,flag)
            % relax the coordinates of r [(x,y),n] points
            
            % get number of nodes
            nnodes = size(r,2);
            
            % define time step
            dt = 0.01;

            % store force (negative gradient) in memory
            J = -am_lib.diag_(F.J); 
            
            % get local effective mass
            if contains(flag,'effectivemass'); M = 1./am_lib.normc_(am_lib.diag_(F.H)); M = max(M./max(M(:)),0.1); end

            % initialize
            if nargout==0 || contains(flag,'draw'); clf; F.plot_field('F'); hold on; end
            nsteps=1000;  mass(1:nnodes)=1; v=zeros(F.d,nnodes,2); a=zeros(F.d,nnodes,2); PE=Inf(2,1);
            % loop over md
            for j = 1:nsteps
                if contains(flag,'dropedge') % drop points outside boundary?
                    ex_ = true(1,nnodes);
                    ex_(ex_) = (min(F.R(1,:))<r(1,ex_) & max(F.R(1,:))>r(1,ex_));
                    ex_(ex_) = (min(F.R(2,:))<r(2,ex_) & max(F.R(2,:))>r(2,ex_));
                    nnodes=sum(ex_); r = r(:,ex_); mass=mass(ex_); v = v(:,ex_,:); a=a(:,ex_,:);
                end
                if contains(flag,'droptop') % drop points ontop of each other?
                    ex_ = false(1,nnodes); [~,i]=am_lib.uniquec_(r); ex_(i) = true;
                    nnodes=sum(ex_); r = r(:,ex_); mass=mass(ex_); v = v(:,ex_,:); a=a(:,ex_,:);
                end
                % find index of closest mesh point to each node
                i = knnsearch(F.R(:,:).',r(:,:,1).').'; r(1:F.d,:,1) = F.R(:,i);
                % get gradient
                G = J(1:F.d,i);
                % update mass based on local curvature?
                if contains(flag,'effectivemass'); mass(1:nnodes) = M(i); end
                % get potential and kinetic energies per node
                PE(2) = sum(F.F(i)-min(F.F(:)))./nnodes;
                KE(2) = am_lib.sum_(mass.*v(:,:,2).^2,[1,2])/2./nnodes;

                % get total force
                f_total= G;
                
                % get force = ma on each point as the negative of the gradient
                a(:,:,2) = f_total./mass;
                % integrate
                integrator = 'sd';
                switch integrator
                    case {'vvhs','velocity-verlet-half-step'}
                        v(:,:,2) = v(:,:,1) + a(:,:,1)/2 .* dt;
                        r(:,:,2) = r(:,:,1) + v(:,:,2) .* dt;
                        v(:,:,2) = v(:,:,1) + a(:,:,2)/2 .* dt;
                    case {'vvd','velocity-verlet-damped'}
                        r(:,:,2) = r(:,:,1) + v(:,:,1) .* dt + a(:,:,1)/2 .* dt.^2;
                        v(:,:,2) = v(:,:,1) + (a(:,:,1)+a(:,:,2))/2 .* dt;
                        ex_=sum(v(:,:,2).*f_total,1)<0; v(:,ex_,2)=0;
                    case {'vv','velocity-verlet'}
                        r(:,:,2) = r(:,:,1) + v(:,:,1) .* dt + a(:,:,1)/2 .* dt.^2;
                        v(:,:,2) = v(:,:,1) + (a(:,:,1)+a(:,:,2))/2 .* dt;
                    case {'sd','steepest-descent'}
                        r(:,:,2) = r(:,:,1) + v(:,:,2) .* dt + a(:,:,2)/2 .* dt.^2;
                        v(:,:,2) = 0;
                    otherwise
                        error('unknown integrator');
                end
                % update
                r = circshift(r,-1,3); v = circshift(v,-1,3); PE = circshift(PE,-1,1); KE = circshift(KE,-1,1);
                dPE = abs(PE(2)-PE(1)); dKE = abs(KE(2)-KE(1));
                % draw
                if nargout==0 || contains(flag,'draw')
                    % print 
                    if j==1; fprintf('%10s %10s %10s %10s %10s\n','TE','PE','dPE','KE','dKE'); 
                    else;    fprintf('%10.5g%10.5g %10.5g %10.5g %10.5g\n',PE(1)+KE(1),PE(1),dPE,KE(1),dKE); end
                    % plot
                    hold on; h = plot3(r(1,:,1),r(2,:,1),ones(size(r(2,:,1)))*1E10,'Marker','o','LineStyle','-','Color','k');
                    drawnow; delete(h);
                end
                % check energy difference and break if smaller
                if j>2; if dPE<ediff; break; end; end
            end
            hold on; h = plot3(r(1,:,1),r(2,:,1),ones(size(r(2,:,1)))*1E10,'Marker','o','LineStyle','-','Color','k');
        end
        
        function [F] = simulate_ising_(F,kT,M,algorithm,boundary)
            % 2D ising model
            % F = am_field.define([2,2].^7,[2,2].^7,{'pdiff','pdiff'});
            % F.F = 2*round(rand([1,F.n]))-1; % initialize binary field
            % F = F.simulate_ising_(2,20000,'MT53','pbc')

            % multiplication is faster than division
            beta = 1./kT;
            
            % compute dimensions
            p = cumprod(F.n);
            
            % define boundary conditions
            switch boundary
                case {'pbc','periodic'}
                    x = reshape(1:p(end),F.n);
                    [y{1:2}] = ndgrid(1:F.n(1),1:F.n(2));
                    % get neighbors
                    n_ = @(i) [x(mod(y{1}(i)+1-1,F.n(1))+1,y{2}(i)), ...
                               x(mod(y{1}(i)-1-1,F.n(1))+1,y{2}(i)), ...
                               x(y{1}(i),mod(y{2}(i)+1-1,F.n(2))+1), ...
                               x(y{1}(i),mod(y{2}(i)-1-1,F.n(2))+1)];
                    % get energy gain from flipping spin
                    en_ = @(F,i) 2 * F.F(1,i) * sum( F.F(1,n_(i)) );
                case {'hbc','helical'}
                    n_ = @(i) [mod(i   -1-1,p(end))+1, ...  % up
                               mod(i   +1-1,p(end))+1, ...  % down
                               mod(i-p(1)-1,p(end))+1, ...  % left
                               mod(i+p(1)-1,p(end))+1];     % right
                    en_= @(F,i) 2 * F.F(1,i) * sum(F.F(1,n_(i)));
                otherwise
                    error('unknown boundary condition');
            end
            
            % run monte carlo 
            switch algorithm
                case 'MT53' % Metropolis (J. Chem. Phys. 21, 1087 (1953).)
                    for j = 1:M
                        % pick a random coordinate to flip
                        i = randi(p(end));
                        % calculate energy of flipping
                        E = en_(F,i);
                        % flip?
                        if E <= 0 || rand(1) <= exp(-E*beta); F.F(1,i) = -F.F(1,i); end
                        % plot?
                        if mod(i,round(M*0.01))==0; F.plot_field('F'); title(num2str(j)); drawnow; end
                    end
                case 'WL89' % Wolff (Phys. Rev. Lett. 62, 361 (1989).)
                    % build edge list
                    G = zeros(2,4*p(end)); G(1,:) = repelem(1:p(end),4);
                    for i = 1:p(end); G(2,4*(i-1)+[1:4]) = n_(i); end
                    % run algorithm
                    for j = 1:M
                        % pick a random coordinate to flip
                        i = randi(p(end));
                        % flood fill
                        cluster = am_lib.floodfill_( F.F(:), G, 1-exp(-2*beta), i );
                        % flip cluster
                        F.F(1,cluster) = - F.F(1,cluster);
                        % plot?
                        if mod(i,round(M*0.01))==0; F.plot_field('F'); title(num2str(j)); drawnow; end
                    end
                case 'SW87' % Swendsen-Wang (Phys. Rev. Lett. 58, 86 (1987).)
                    % build edge list
                    G = zeros(2,4*p(end)); G(1,:) = repelem(1:p(end),4);
                    for i = 1:p(end); G(2,4*(i-1)+[1:4]) = n_(i); end
                    % run algorithm
                    for j = 1:M
                        % divide field into clusters by linking parallel spins with probability p = 1-exp(-2*beta)
                        i=0; C = zeros([1,F.n]); 
                        while any(C(:)==0)
                            inds = am_lib.floodfill_( F.F(:), G, 1-exp(-2*beta), find(~C(:),1) ); i=i+1; C(inds)=i;
                        end
                        nclusters=i;
                        % flip each cluster with probability 1/2
                        for i = 1:nclusters
                            if rand()>.5; F.F(1,C(:)==i) = -F.F(1,C(:)==i); end
                        end
                        % plot?
                        if mod(i,round(M*0.01))==0; F.plot_field('F'); title(num2str(j)); drawnow; end
                    end
                case {'KL87','MULTI'} % Kandel-Loh Multigrid (Phys. Rev. Lett. 60, 1591 (1988)) 
                    % NOT WORKING PROPERLY. TO DO:
                    % D. Kandel et al., PRL 60, 1591 (1988) uses weights to describe bonds between high level clusters
                    % weights N{l}(3,:) would need to equal the neighbors shared between clusters, not just set to 1 as it is now.
                    %
                    % Data structures:
                    % I maps clusters between levels
                    % N maps clusters within levels (weights are ignored).
                    
                    % set maximum cluster dimension
                    cdim = 5;
                    % define multigrid pattern
                    updn=[1,1,1,1,1,1,1,1,-1,0,1,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0]; % updn=repmat(updn,1,10);
                    % updn=[1,1,1,1,0];
                    % initialize: field S(i) for i clusters
                    S = [F.F(1,:),NaN(1,100)];
                    % initialize: level counter
                    l=1;
                    % initialize: edge list G(1:3,i) for i clusters (cluster-to-cluster connections within same lvl)
                    N{l} = ones(3,4*p(end)); N{l}(1,:) = repelem(1:p(end),4);
                    for i = 1:p(end); N{l}(2,4*(i-1)+[1:4]) = n_(i); end
                    % create edge list P(1:3,i) for i clusters ( cluster(lvl-1)-to-cluster(lvl) )
                    I{l} = [1:p(end);zeros(1,p(end));ones(1,p(end))];

                    % perform monte carlo run
                    for k = 1:numel(updn); l=l+updn(k);
                        switch updn(k)
                            case  0 % do M metropolis sweeps
                                x = am_lib.minmax_(I{l}(1,:)); dx=x(2)-x(1)+1;
                                for j = 1:M*dx
                                    % pick a random cluster to flip
                                    i = x(1)+randi(dx)-1;
                                    % calculate the energy of flipping
                                    E = 2*S(i)*sum(S(N{l}(2,N{l}(1,:)==i)));
                                    % flip?
                                    if E <= 0 || rand(1) <= exp(-E/kT); S(i) = -S(i); end
                                end
                            case -1 % refine
                                % distribute course field values on fine clusters 
                                S(I{l+1}(2,:))=S(I{l+1}(1,:));
                                % reset l+1 level
                                x=am_lib.minmax_(I{l+1}(1,:)); S(x(1):x(2))=NaN; I(l+1)=[]; N(l+1)=[]; 
                            case +1 % corsen
                                % Sj marks which spins have already been considered consider, j is the next position to sweep, iC counts and labels clusters
                                i=0; nIs=0; nNs=0; N{l}=zeros(3,4*p(end)); I{l}=zeros(3,4*p(end)); 
                                x=am_lib.minmax_(I{l-1}(1,:)); j=x(1); iC=x(2); Sj=S; Sj(1:x(1)-1)=NaN; 
                                % divide field into clusters by linking parallel spins with probability p = 1-exp(-2/kT)
                                while ~isempty(j)
                                    % perform flood fill
                                    [I_,N_] = am_lib.floodfill_( Sj , N{l-1} , 1-exp(-2/kT) , j , cdim );
                                    % mark sites as visited and find next node which has not been assigned a cluster
                                    i=i+1; iC=iC+1; Sj(I_)=NaN; j=find(~isnan(Sj(:)),1);
                                    % save clusters and neighbors (in terms of last lvl's cluster) 
                                    nI = numel(I_); I{l}(:,nIs+[1:nI]) = [repmat(iC,1,nI);I_;ones(1,nI)]; nIs=nIs+nI; %
                                    nN = numel(N_); N{l}(:,nNs+[1:nN]) = [repmat(iC,1,nN);N_;ones(1,nN)]; nNs=nNs+nN; % 
                                end
                                % trim un-used space
                                N{l} = N{l}(:,1:nNs); I{l} = I{l}(:,1:nIs);
                                % rewrite same-level cluster-cluster edge list in terms of clusters on this lvl
                                aux=[]; aux(I{l}(2,:))=I{l}(1,:); N{l}(2,:)=aux(N{l}(2,:));
                                % pass on field values on to new level. S(1,2,3,4,...,nclusters)
                                aux=[]; aux(I{l}(1,:))=I{l}(2,:); S(aux~=0)=S(aux(aux~=0));
                                % keep only one instance of each edge (not consistent with Kendal's original algorithm).
                                N{l}=am_lib.uniquec_(N{l});

                                % plot graph
                                subplot(1,2,2);
                                X=cat(2,I{2:end}); X = unique(sort(X(1:2,:)',2),'rows')'; 
                                plot(graph(X(1,:),X(2,:)),'Layout','layered','Sinks',unique(I{1}(1,:)),'Sources',unique(I{end}(1,:)))
                                ylim([0, max(cumsum(updn))+2]); drawnow; 
                        end
                        % distribute course field values on fine clusters and update F
                        for i = l-1:-1:1; S(I{i+1}(2,:))=S(I{i+1}(1,:)); end
                        F.F = reshape(S(1:p(end)),[1,F.n]);
                        % plot
                        subplot(1,2,1);
                        F.plot_field('F'); drawnow;
                    end                    
                otherwise
                    error('unknown algorithm');
            end
           
        end
        
        function [h] = plot_field(F,field)
            
            sl_ = @(field,i)   squeeze(F.(field)(i,:,:,:,:,:));
            
            switch size(F.(field),1)
                case {1} % scalar
                    switch F.d
                        case 1 % 1D
                            plot(F.R,F.(field));
                        case 2 % 2D
                            if isreal(F.(field))
                                set(gcf,'color','w');
%                                 imagesc(flipud(squeeze(F.(field)).'));
                                h = surf(sl_('R',1), sl_('R',2), squeeze(F.(field)));  
                                h.EdgeColor = 'none'; h.LineWidth = 1; 
                                view([0 0 1]); daspect([1 1 1]); axis tight;
                            else
                                switch 1
                                    case 1 
                                        cmap  = am_lib.colormap_('spectral',200); n = size(cmap,1); 
                                        phase = angle(F.(field))/(2*pi)+1/2; % [0,1]
                                        phase = 2*(phase-1/2); % [-1,1]
                                        amp   = log(abs(F.(field))); amp = (amp-min(amp(:)))./(max(amp(:))-min(amp(:)));
                                        cmap  = reshape( am_lib.clight_(cmap(ceil((n-1)*amp+1),:),phase) , [F.n,3]);
                                    case 2 % phase
                                        cmap  = am_lib.colormap_('jet',200); n = size(cmap,1); 
                                        phase = angle(F.(field))/(2*pi)+1/2; % [0,1]
                                        cmap  = reshape( cmap(ceil(n*phase),:) , [F.n,3]);
                                    case 3 % amplitude
                                        cmap  = am_lib.colormap_('jet',200); n = size(cmap,1); 
                                        amp   = log(abs(F.(field))); amp = (amp-min(amp(:)))./(max(amp(:))-min(amp(:))); 
                                        cmap  = reshape( cmap(ceil((n-1)*amp+1),:) , [F.n,3]);
                                end
                                set(gcf,'color','w');
                                h = surf(sl_('R',1), sl_('R',2), squeeze(abs(F.(field))), cmap); h.EdgeColor= 'none'; h.LineWidth = 1; 
                                view([0 0 1]); daspect([1 1 1]); axis tight;
                            end
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
                otherwise
                    error('field is not scalar or vector');
            end

        end

        function [ax]= plot_jacobian(F)
            
            if isempty(F.J); F.J = F.get_jacobian(); end
            
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
        
        function [h,f] = plot_statistical_function(F,varargin)
            
            if F.d~=2; error('plot_statistical_function is only implemented for 2d field'); end
            if ~all(contains(F.s,'diff')); error('statistical functions are only implemented for finite differences'); end

            [x,f] = am_lib.get_statistical_function(squeeze(F.F),varargin{:});
            
            h = am_lib.overlay_(squeeze(F.F),x,f);
        end
        
        % statistical quantities


        % TO DO: drop particles (nodes) on a converged potential and move MD based on Lagrangian mechanics
        % evolve each particle using different algorithms, including verlet. 
        % this is a precursor to the nudge elastic band method for finding minimum energy path.
        % can also add nose hoover and make the simulation an NVT ensamble instead of NVE
        

%                 % 6) current temperature
%                 Tj = 2/3*KE(j)/uc.natoms/k_boltz;
% 
%                 % 5) compute Nose-Hoover drag: p_eta = KE - TE
%                 nosehoover = v(:,:,j)/Q * ( Tj - T ) / uc.natoms;
% 
%                 % 6) get acceleration
%                 acc = f(:,:,j) ./ uc.mass(uc.species);
% 
%                 % ***) update md [frac]: x' = x + v * dt; v' = v + a * dt; Nose-Hoover dv/dt becomes a - p_eta / Q * v;
%                 if j ~= nsteps
%                     u(:,:,j+1) = u(:,:,j) + dt * v(:,:,j);
%                     v(:,:,j+1) = v(:,:,j) + dt * (acc - nosehoover);
%                 end
    end
   
    methods %(Access = protected) % internal stuff
        
        function [R] = get_collocation_points(F)
            for i = 1:F.d % loop over dimensions
                n = F.n(1:F.d);
                switch F.s{i}
                    case 'chebyshev'; R{i} = am_field.chebyshevUr_(n(i),'edge');
                    case 'legendre';  R{i} = am_field.legendrer_(n(i));
                    case 'fourier';   R{i} = am_field.fourierr_(n(i));
                    case 'cdiff';     R{i} = am_field.cdiff_(n(i));
                    case 'pdiff';     R{i} = am_field.pdiff_(n(i));
                    otherwise; error('unknown s');
                end
                n(i) = 1; R{i} = repmat(permute(F.a(i)*R{i},circshift([1,2,3],i-1)),n(1:F.d));
            end
            R = permute(cat(4,R{:}),[4,1,2,3]);
        end

        function [T] = get_field_type(F)
            if ~isempty(F.T); T = F.T; return; end
            switch size(F.F,1)
                case {1}  ; T = 'scalar';
                case {2,3}; T = 'vector';
                otherwise; error('unknown field type');
            end
        end

        function [J] = get_jacobian(F)
            if isempty(F.T); F.T = get_field_type(F); end
            % define matmul which supports sparse matrices
            matmul_ = @(A,B) reshape(A*reshape(B,size(B,1),[]),size(A,1),size(B,2),size(B,3),size(B,4),size(B,5));
            % get differentiation matrices
            Q = F.get_differentiation_matrices();
            % allocate space
            J = zeros([3,3,F.n]);
            for i = 1:F.d
                D = Q{i}(:,:,1); % keep only first derivative
                if strcmp(F.s{i},'cdiff'); D=sparse(D); end % speed up finite difference with sparse matrices
                p = [1:F.d+1]; p([1,i+1])=p([i+1,1]); % permute accordingly
                switch F.T % evaluate derivatives
                    case 'scalar'; J(i,i,:,:,:)     = permute(matmul_(D,permute(F.F,p)),p);
                    case 'vector'; T                = permute(matmul_(D,permute(F.F,p)),p);
                                   J(i,1:F.d,:,:,:) = reshape(T(1:F.d,:),size(J(i,1:F.d,:,:,:)));
                    otherwise; error('unknown field type');
                end
            end
            % undefine boundaries
            if F.d>0 && strcmp(F.s{1},'cdiff'); J(1,1,[1,end],:,:)=0; end
            if F.d>1 && strcmp(F.s{2},'cdiff'); J(2,2,:,[1,end],:)=0; end
            if F.d>2 && strcmp(F.s{3},'cdiff'); J(3,3,:,:,[1,end])=0; end
        end

        function [H] = get_hessian(F)
            if ~contains(F.T,'scalar'); error('hessian is only defined for scalar fields'); end
            % define matmul which supports sparse matrices
            matmul_ = @(A,B) reshape(A*reshape(B,size(B,1),[]),size(A,1),size(B,2),size(B,3),size(B,4),size(B,5));
            % get differentiation matrices
            Q = F.get_differentiation_matrices();
            % allocate space
            H = zeros([3,3,F.n]);
            for i = 1:F.d % loop over dimensions
                D = Q{i}(:,:,1); % keep only first derivative
                if strcmp(F.s{i},'cdiff'); D=sparse(D); end % speed up finite difference with sparse matrices
                p = [1:F.d+1]; p([1,i+1])=p([i+1,1]); % evaluate hessian from jacobian
                switch F.T % evaluate derivatives
                    case 'scalar'
                        Ji = am_lib.diag_(F.J); H(i,:,:,:,:) = permute(matmul_(D,permute(Ji,p)),p);
                    otherwise; error('hessian is only defined for scalar fields');
                end
            end
            % undefine boundaries
            if F.d>0 && strcmp(F.s{1},'cdiff'); H(1,1,[1,end],:,:)=0; end
            if F.d>1 && strcmp(F.s{2},'cdiff'); H(2,2,:,[1,end],:)=0; end
            if F.d>2 && strcmp(F.s{3},'cdiff'); H(3,3,:,:,[1,end])=0; end
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
        
        function [Q] = get_topological_charge(F)
            % copy
            N = F;
            % normalize field
            N.F = N.F./am_lib.normc_(N.F);
            % recompute jacobian
            N.J = N.get_jacobian();
            % compute charge
            switch F.d
                case 3
                    Q = cross(F.J(:,1,:,:,:),F.J(:,2,:,:,:),1); 
                    Q = dot( F.F, permute(Q,[1,3,4,5,2]), 1);
                case 2
                    Q = cross(N.J(:,1,:,:),N.J(:,2,:,:),1); 
                    Q = permute(Q(3,:,:,:),[1,3,4,2]);
                otherwise
                    error('ERROR [get_topological_charge]: invalid dimension.');
            end
        end
        
        function [D,L,G] = get_flattened_differentiation_matrices(F)
            % checked with this code:
            % % 2D difference matrix
            % F = am_field.define([500,500],[1,1]*2*pi,{'chebyshev','chebyshev'});
            % F.F = sum(sin(F.R),1); F = F.get_derivatives; 
            % [D,L] = get_flattened_differentiation_matrices(F);
            % L = reshape(L*F.F(:),F.n);
            % max(abs(L(:)-F.L(:)))
            % subplot(1,2,1); surf(L,'edgecolor','none'); view([0 0 1]); daspect([1 1 1]); axis tight;
            % subplot(1,2,2); surf(squeeze(F.L),'edgecolor','none'); view([0 0 1]); daspect([1 1 1]); axis tight;
            
            Q = F.get_differentiation_matrices();
            % note: although the specral differentiation matrix may be full, 
            %       when spanning multiple dimensions it will become sparse.
            % divergence
            D = cellfun(@(x)sparse(x(:,:,1)),Q,'UniformOutput',false);
            D = am_field.get_flattened_divergence_(D{:});
            % laplacian
            if nargout < 2; return; end
            L = cellfun(@(x)sparse(x(:,:,2)),Q,'UniformOutput',false);
            L = am_field.get_flattened_divergence_(L{:});
            % gradient
            if nargout < 3; return; end
            G = cellfun(@(x)sparse(x(:,:,1)),Q,'UniformOutput',false);
            G = am_field.get_flattened_gradients_(G{:});
            
        end
        
        function [Q]     = get_differentiation_matrices(F)
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
            for i = 1:F.d
                switch F.Y % convert to cylindrical or spherical coordinates
                    case {'cartesian'} % do nothing
                    case {'polar'};       if i == 1; elseif i == 2; Q{i}(:,:,1) = Q{i}(:,:,1)./permute(F.R(1,:,:),[2,3,1]); else; error('invalid dimension'); end
                    case {'cylindrical'}; if i == 1; elseif i == 2; Q{i}(:,:,1) = Q{i}(:,:,1)./permute(F.R(1,:,:),[2,3,1]); elseif i == 3; else; error('invalid dimension'); end
                    case {'spherical'};   if i == 1; elseif i == 2; Q{i}(:,:,1) = Q{i}(:,:,1)./permute(F.R(1,:,:),[2,3,1]); elseif i == 3; Q{i}(:,:,1) = Q{i}(:,:,1)./( permute(F.R(1,:,:),[2,3,1]).*sin(permute(F.R(2,:,:),[2,3,1])) ); else; error('invalid dimension'); end
                    otherwise; error('unknown coordinate type');
                end
            end
        end
        
    end
     
    methods (Static, Access = protected) % polynomials and spectral methods

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
            D = am_field.get_differentiation_matrix_(x);
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
            D = am_field.get_differentiation_matrix_(x); D = cat(3,D,D^2);
            if nargout < 3; return; end
            w(:,1) = am_field.get_integration_weights_(x);
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
            D = am_field.get_differentiation_matrix_(x); D = cat(3,D,D^2);
            if nargout < 3; return; end
            w(:,1) = am_field.get_integration_weights_(x);
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
            D = am_field.get_differentiation_matrix_(x); D = cat(3,D,D^2);
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
            D = am_field.get_differentiation_matrix_(x, alpha, beta);
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
            D = am_field.get_differentiation_matrix_(x, alpha, beta);
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
                [c,v] = am_field.get_differentiation_weights_([-1,0,1],i); nvs = numel(v); m = ceil(nvs/2);
                D(:,:,i) = toeplitz([c(m:-1:1),zeros(1,n-m)],[c(m:end),zeros(1,n-m)])*n.^(i);
            end
            if nargout < 3; return; end
            w(1:n,1) = 1;
        end        

        function [x,D,w] = pdiff_(n) % evenly spaced periodic central difference [0,1)
            x(1:n,1) = [0:n-1]/n;
            if nargout < 2; return; end
            % get first and second derivative
            D = zeros(n,n,2);
            for i = 1:2
                [c,v] = am_field.get_differentiation_weights_([-1,0,1],i); nvs = numel(v); m = ceil(nvs/2);
                D(:,:,i) = am_lib.circulant_(circshift([c,zeros(1,n-nvs)],m-nvs))*n.^(i);
            end
            if nargout < 3; return; end
            w(1:n,1) = 1;
        end
    
        function [D]     = get_flattened_gradients_(Dx,Dy,Dz) 
            switch nargin
                case 1
                    D{1} = Dx;
                case 2
                    n(1) = size(Dx,1); n(2) = size(Dy,1); 
                    D{1} = kron(eye(n(2)),Dx);
                    D{2} = kron(Dy,eye(n(1)));
                case 3
                    n(1) = size(Dx,1); n(2) = size(Dy,1); n(3) = size(Dz,1);
                    D{1} = kron(eye(n(3)),kron(eye(n(2)),Dx));
                    D{2} = kron(eye(n(3)),kron(Dy,eye(n(1))));
                    D{3} = kron(Dz,kron(eye(n(2)),eye(n(1))));
                otherwise
                    error('not yet implemented');
            end
        end
        
        function [D]     = get_flattened_divergence_(Dx,Dy,Dz)
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

        function [L]     = get_flattened_laplacian_(varargin) 
            L = am_field.get_divergence_matrix(varargin{:})^2;
        end

        function [c,x]   = get_differentiation_weights_(x,n) 
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
        
        function [D]     = get_differentiation_matrix_(varargin) 
            % get_differentiation_matrix_(x)
            % get_differentiation_matrix_(x,order)
            % get_differentiation_matrix_(x,alpha,beta)
            
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
  
        function [w]     = get_integration_weights_(x) 
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
    
    methods (Static) % electricity/magnetism
        
        function [A] = get_vector_potential(R,dI,I)
            % get the magnetic vector potential A [3,x,y,z] at positons R [x,y,z] given current flow dI [3,(x,y,z)] at positions I [3,(x,y,z)]
            I = reshape(I,3,[]); dI = reshape(dI,3,[]); M = size(I,2);
            if size(I,2)~=size(dI,2); error('dI and I dimension mismatch'); end
            A = am_field.mu0/(4*pi)*sum(reshape(dI,3,1,1,1,M)./(am_field.tiny+am_lib.normc_(R-reshape(I,3,1,1,1,M))),5);
        end
        
    end
    
    methods (Static) % TEM analysis
        
        function [E,W,Ar_g,Pr_g,dPr_g] = gpa(Ir,flag,th,ex_r_,ex_g_)
            %
            % Geometric phase analysis. 
            %
            % Required inputs: Ir is the real-space image. Flag options are 'rot:image/rot:field, interactive, activate'.
            % Optional inputs: ex_r_ and ex_g_ are real-space mask of reference region and reciprocal-space mask of bragg reflections. 
            %                  If not provided, the function runs interactively. 
            % Outputs:         E and W are strain and rotational fields. 
            %                  Ar_g, Pr_g, and dPr_g are real-space amplitude, phase, and phase gradient maps.
            %
            % Antonio Mei April 2018
            %

            if nargin==1; flag=''; end
            if nargin>=2 || isempty(flag); flag=[flag,',interactive']; end
            if nargin<3; th=[]; end

            % rotate image?
            if contains(flag,{'rot:image','rot:field'})
                if contains(flag,'interactive') && isempty(th)
                    imagesc_(Ir); title('Select two points to level horizon'); p = ginput(2); 
                    if isempty(p); th=0; else; dp = diff(p.',1,2); th = atan2d(dp(2),dp(1))+90; end
                end
                if contains(flag,{'rot:image'}); Ir = imrotate(Ir,th,'bicubic');  end
            end

            % get real space image and mask reference ROI
            if contains(flag,'interactive') || isempty(ex_r_)
                figure(1); clf; hr = imagesc_(Ir); title('Select reference ROI');
                er = imellipse(gca); ex_r_ = createMask(er,hr); 
            end

            % get reciprocal space image and mask off reflection of interest
            nbraggs=2; Ik = fftshift(fftn(Ir));
            if contains(flag,'interactive') || isempty(ex_g_)
                figure(2); clf; hk = imagesc_(log(abs(Ik).^2)); title('Zoom in on the image'); 
                input('Press enter when done'); title(sprintf('Select %i reflections',nbraggs));
                ex_g_ = cell(1,nbraggs); for q = 1:nbraggs; hg = imellipse(gca); ex_g_{q} = createMask(hg,hk); end
            end

            % compute grids
            [n,m] = size(Ir); fft_k_ = @(n,Dx) ([0:(n-1)]-floor(n/2))./Dx; % do not apply fftshift here.
            [r{1:2}] = ndgrid([1:n], [1:m]);
            [k{1:2}] = ndgrid(fft_k_(n,r{1}(end)-r{1}(1)), fft_k_(m,r{2}(end)-r{2}(1))); 

            % allocate space for Bragg vectors K((x,y),(1,2)) and phase gradient G((d/dx,d/dy),(1,2),(x,y))
            ndims=2; K = zeros(ndims,nbraggs); G = zeros(ndims,nbraggs,n,m); 
            for q = 1:nbraggs % loop over Bragg vectors
                % get position of bragg reflection for reference ROI
                Rk_r = fftshift(fftn(Ir.*ex_r_)).*ex_g_{q}; [~,i]=max(abs(Rk_r(:))); [i,j]=ind2sub([n,m],i); % ind2subs checked.
                Ik_g = Ik.*ex_g_{q}; Ik_g = circshift(circshift(Ik_g,1-i,1),1-j,2); % moves bragg reference: Ik_g(i,j) -> Ik_g(1,1) % gives the passive strain

                % get amplitude, phase, and phase gradient images
                Ir_g = ifftn(fftshift(Ik_g)); Ar_g{q} = abs(Ir_g); Pr_g{q} = angle(Ir_g); 
                % Pr_g = Pr_g - 2*pi*(r{1}*k{1}(i,j)+r{2}*k{2}(i,j)); %% add shift?
                [T{2:-1:1}] = gradient(exp(1i.*Pr_g{q})); dPr_g{q} = cat(3,T{:}); dPr_g{q} = imag(exp(-1i.*Pr_g{q}).*dPr_g{q});

                % save referenece recirpocal lattice point and gradients
                K(:,q)     = [k{1}(i,j);k{2}(i,j)];
                G(:,q,:,:) = permute(dPr_g{q},[3,4,1,2]);
            end  

            % get direct lattice vectors [ Eq 36  M.J. Hytch et al. Ultramicroscopy 74 (1998) 131?146 ]
            A = inv(K).';

            % get strain field [ Eq 42  M.J. Hytch et al. Ultramicroscopy 74 (1998) 131?146 ]
            E = am_lib.matmul_(A,permute(G,[2,1,3,4]))/(-2*pi);

            % exclude meaningless region that was cut off by rotating the image
            if     contains(flag,'rot:image')
                ex_rotated_ = abs(Ir)<eps; ex_rotated_ = permute(ex_rotated_,[3,4,1,2]);
                ex_rotated_ = repmat(ex_rotated_,[2,2,1,1]); E(ex_rotated_) = NaN;
            elseif contains(flag,'rot:field')
                R = rotz(th); E = am_lib.matmul_(am_lib.matmul_(R(1:2,1:2),E),R(1:2,1:2).');
            end

            % (anti-)symmeterize strain E and rotational W fields
            W = (E - permute(E,[2,1,3,4]))/2;
            E = (E + permute(E,[2,1,3,4]))/2;

            % convert strain from passive (space-stretching) to active (lattice-stretching) representation
            if contains(flag,'activate'); activate_ = @(E) 1./(E+1)-1; else; activate_ = @(E) E; end
            E = activate_(E); W = activate_(W);

            function h=imagesc_(A)
                h=imagesc(A); axis tight; daspect([1 1 1]); axis off; 
            end
        end

        function [s_averaged, s_shifted, xs, ys] = correct_stack_drift(stack) % correct drift 
            %  CORRECT_DRIFT aligns stack to first image in stack, by cross-correlation
            % passes out cropped image and aligned stack (not cropped)
            % correlation to first image of the stack
            % stack should be a 3 dimensional matrix, [x, y, images]
            % Megan Holtz 2015 based off off correct_adf_drift by
            % David Muller 2003
            %  Added bandpass filter, 2/21/2005   DM

            [xsize,ysize,numIm] = size(stack);

            s_shifted=stack;

            if numIm == 1
                xs=0; ys=0; s_averaged = stack;
            else
                xs=zeros(1,numIm); ys=xs;
                for j=2:size(stack,3)
                    [s_shifted(:,:,j), xs(j), ys(j)] = correct_adf_drift(stack(:,:,j),stack(:,:,1));
                end
                s_averaged = sum(s_shifted,3);
                % crop the resulting image
                xdw = [1, xs(xs<=xsize/2)];
                ydw = [1, ys(ys<=ysize/2)];
                xup = [xsize, xs(xs>xsize/2)];
                yup = [ysize, ys(ys>ysize/2)];
                yup = min(yup); ydw = max(ydw);  
                xdw = max(xdw); xup = min(xup);
                s_averaged = s_averaged(xdw+2:xup, ydw+2:yup);
            end

            function [shifted_adf1,xshf,yshf] = correct_adf_drift(adf1,adf2)
                %  CORRECT_DRIFT aligns array adf1 to adf2, by cross-correlation 
                %  2015 edits ESP and MEH
                %  David Muller 2003
                %  Added bandpass filter, 2/21/2005   DM

                n=4;  % number of pixels to average over
                adf1=double(adf1); adf2=double(adf2);

                [xr,xc] = size(adf1);

                ix=1:xr; iy=1:xc; [ry,rx]=meshgrid(iy,ix); % initialize grid
                k   = sqrt( ((rx-xr/2-0.5)/xr).^2+((ry-xc/2-0.5)/xc).^2 ); % k space mesh
                hpk = fftshift( (k<1/2/n).*sin( 2*n*pi*k).^2 );  % bandpass filter for the images
                % remove discontinuities at the edges of the images
                fil = (sin( pi*rx/xr ) .* sin(pi*ry/xc) ).^2;  %damp edges out to 0
                ref1=( adf1-mean(mean(adf1)) ).*fil;  
                ref2=( adf2-mean(mean(adf2)) ).*fil;
                % cross correlate
                crrd = abs(ifft2( hpk.*fft2(ref2) .* conj(fft2(ref1)) ));
                %find peak of cross correlation
                mm=max(max(crrd));
                [xshf,yshf] = find(crrd==mm);  % amount to shift by
                xshf=xshf(1); yshf=yshf(1);
                xshf = xshf - 1; yshf = yshf - 1;
                w=-exp( -(2*pi*1i)*(xshf*rx ./xr + yshf * ry ./xc ));
                shifted_adf1=abs(ifft2( fft2( adf1 ) .* w )); % shift IF ADF
            end
        end

    end

end

