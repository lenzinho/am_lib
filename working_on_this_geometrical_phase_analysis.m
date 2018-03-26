clear;clc;cd('/Users/lenzinho/Downloads/SnO_ABM');
% [dm3struct] = DM3Import('SnO_0084.dm3'); % high mag
[dm3struct] = DM3Import('SnO_0017.dm3'); % low mag
%%

    flag = 'dm3,rot:field';
    if     contains(flag,'dm3')
        Ir = dm3struct.image_data;  % real space image
        dx = dm3struct.xaxis.scale;
        dy = dm3struct.yaxis.scale;
    elseif contains(flag,'test')
        % load test image
        [n,m]=deal(1000,1200); dx=1;dy=1; [a,c]=deal(50*pi/(2*n),70*pi/(2*m));
        [r{1:2}] = ndgrid(([0:n-1]-floor(n/2)).*dx,([0:m-1]-floor(m/2)).*dy);

        Sx = 0.35*am_lib.fermi_dirac_(r{1}./50); % define strain passively: expanding the coordinate
        Sy = 0.20*am_lib.fermi_dirac_(r{2}./80); % is equivalent to shrinking the lattice! 
        Ix = sin(a*r{1}.*(1 + Sx));
        Iy = sin(c*r{2}.*(1 + Sy));
        Ir = (Ix.*Iy).^2;
    end

    % rotate image?
    if contains(flag,{'rot:image','rot:field'})
        imagesc_(Ir); title('Select two points to level horizon'); p = ginput(2); 
        if isempty(p); th=0; else
            dp = diff(p.',1,2); th = atan2d(dp(2),dp(1))+90; 
            if contains(flag,{'rot:image'}); Ir = imrotate(Ir,th,'bicubic');  end
        end
    end
    
    % get real space image and mask reference ROI
    figure(1); clf; hr = imagesc_(Ir); title('Select reference ROI');
    er = imellipse(gca); ex_r_ = createMask(er,hr); 

    % get reciprocal space image and mask off reflection of interest
    figure(2); clf; nbraggs=2; ex_g_ = cell(1,nbraggs); Ik = fftshift(fftn(Ir)); 
    hk = imagesc_(log(abs(Ik).^2)); title('Zoom in on the image'); input('Press enter when done'); title(sprintf('Select %i reflections',nbraggs));
    for q = 1:nbraggs; hg = imellipse(gca); ex_g_{q} = createMask(hg,hk); end
    
    % compute grids
    [n,m] = size(Ir); fft_k_ = @(n,Dx) ([0:(n-1)]-floor(n/2))./Dx; % do not apply fftshift here.
    [r{1:2}] = ndgrid([1:n], [1:m]);
    [k{1:2}] = ndgrid(fft_k_(n,r{1}(end)-r{1}(1)), fft_k_(m,r{2}(end)-r{2}(1))); 
%%
    % allocate space for Bragg vectors K((x,y),(1,2)) and phase gradient G((d/dx,d/dy),(1,2),(x,y))
    ndims=2; K = zeros(ndims,nbraggs); G = zeros(ndims,nbraggs,n,m); 
    for q = 1:nbraggs % loop over Bragg vectors
        % get position of bragg reflection for reference ROI
        Rk_r = fftshift(fftn(Ir.*ex_r_)).*ex_g_{q}; [~,i]=max(abs(Rk_r(:))); [i,j]=ind2sub([n,m],i); % ind2subs checked.
        Ik_g = Ik.*ex_g_{q}; 
        Ik_g = circshift(circshift(Ik_g,1-i,1),1-j,2); % moves bragg reference: Ik_g(i,j) -> Ik_g(1,1) % gives the passive strain
        
        % get amplitude, phase, and phase gradient images
        Ir_g = ifftn(fftshift(Ik_g)); Ar_g = abs(Ir_g); Pr_g = angle(Ir_g); 
%         Pr_g = Pr_g - 2*pi*(r{1}*k{1}(i,j)+r{2}*k{2}(i,j)); %% add shift?
        [T{2:-1:1}] = gradient(exp(1i.*Pr_g)); dPr_g = cat(3,T{:}); dPr_g = imag(exp(-1i.*Pr_g).*dPr_g);
        
        % save referenece recirpocal lattice point and gradients
        K(:,q)     = [k{1}(i,j);k{2}(i,j)];
        G(:,q,:,:) = permute(dPr_g,[3,4,1,2]);
    end  
    
%%
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
    
    % convert strain from passive (fabric-stretching) to active (lattice-strteching) representation
    activate_ = @(E) 1./(E+1)-1;
    % activate_ = @(E) E;
    
    % plot fields
    figure(2); clf;
    ax(1) = subplot(3,2,1); imagesc_(flipud(Ir.'));
    ax(2) = subplot(3,2,3); imagesc_(activate_(flipud(permute(E(1,1,:,:),[4,3,1,2])))); title('exx');
    ax(3) = subplot(3,2,4); imagesc_(activate_(flipud(permute(E(2,2,:,:),[4,3,1,2])))); title('eyy');
    ax(4) = subplot(3,2,5); imagesc_(activate_(flipud(permute(E(1,2,:,:),[4,3,1,2])))); title('exy');
    ax(5) = subplot(3,2,6); imagesc_(activate_(flipud(permute(W(1,2,:,:),[4,3,1,2])))); title('wxy');
    linkaxes(ax);
    
    if contains(flag,'test')
        ax(6) = subplot(3,2,2); imagesc_(activate_(flipud(Sy.')));
    end

    colormap(ax(1),colormap('gray')); 
    for i = 2:numel(ax)
        colormap(ax(i),flipud(am_lib.colormap_('red2blue'))); 
        % caxis(ax(i),3*[-1,1]*abs(activate_(nanstd(abs(E(:)))))); 
        caxis(ax(i),[-0.4 0.4]); 
    end

    %% plot averaged profile along x or y direction
    subplot(1,2,1); 
    set(gcf,'color','w');
    am_lib.set_plot_defaults_; 
    smooth_ = @(y) activate_(conv(y,ones(1,15)/15,'same'));
    z = zeros(m,1); y = zeros(m,3); s = zeros(m,3);
    x(:,1) = r{2}(1,:)*dy - 44;
    y(:,1) = smooth_(squeeze(nanmean(E(1,1,:,:),3)));
    y(:,2) = smooth_(squeeze(nanmean(E(2,2,:,:),3)));
    y(:,3) = smooth_(squeeze(nanmean(E(1,2,:,:),3)));
    s(:,1) = smooth_(squeeze(nanstd(E(1,1,:,:),[],3))); % std error
    s(:,2) = smooth_(squeeze(nanstd(E(2,2,:,:),[],3)));
    s(:,3) = smooth_(squeeze(nanstd(E(1,2,:,:),[],3)));
    plot(y,x); axis tight;
    % errorbar(repmat(x,1,3),y,s); axis tight;
    legend('exx','eyy','exy'); ylabel('d [nm]'); xlabel('strain');
%     xlim([-0.1 0.5]); ylim([-15 45]); % active strain
    % xlim([-0.4 0.1]); ylim([-15 45]); % passive strain
    
  
    %% export figures individually
    figure(1); imagesc_(flipud(Ir.')); colormap(colormap('gray')); caxis([200 1000]); daspect([1 1 1]); box on; axis off;
%     print('SnO_0084_img','-dpdf','-r300');
    %%
    V = [-0.25:0.05:0.25]; set(gcf,'Renderer','painters');
    figure(1); clf; imcontour_(activate_(flipud(permute(E(1,1,:,:),[4,3,1,2]))),V); caxis([-1 1]*0.4); colormap(flipud(am_lib.colormap_('red2blue'))); 

%     print('SnO_0084_yy_','-dpdf','-r300');
    
    
%%
    % over lay the strain on the image
    figure(1); clf; imagesc_(flipud(Ir.')); colormap(colormap('gray')); caxis([200 1000]); daspect([1 1 1]); box on; axis off;
    % overlay statistical function on field
    ax(1) = gca; ax(2) = axes('position',get(gca,'position')); linkaxes(ax);
    % plot isostrain
    V = [-0.4:0.05:0.4]+0.025;
    [~,h]=imcontour_(activate_(flipud(permute(E(1,1,:,:),[4,3,1,2]))),V); caxis(gca,[-1 1]*0.4); colormap(gca,flipud(am_lib.colormap_('red2blue'))); 
    % remove white background
    set(gca,'color','none'); 
    set(h,'linewidth',2)
    %%
    
function h=imagesc_(A)
    h=imagesc(A); axis tight; daspect([1 1 1]); axis off; 
end
function [varargout]=imcontour_(A,varargin)
    [varargout{1:2}]=imcontour(A,varargin{:}); axis tight; daspect([1 1 1]); axis off; 
end