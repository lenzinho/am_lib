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
            % Antonio Mei May 2018
            %

            if nargin>=2; flag=[flag,'interactive']; end
            if contains(flag,{'rot:image','rot:field'})
                if contains(flag,'interactive') || isempty(th)
                    imagesc_(Ir); title('Select two points to level horizon'); p = ginput(2); 
                    if isempty(p); th=0; else; dp = diff(p.',1,2); th = atan2d(dp(2),dp(1))+90; end
                end
                if contains(flag,{'rot:image'}); Ir = imrotate(Ir,th,'bicubic');  end
            end
            if contains(flag,'interactive') || isempty(ex_r_)
                figure(1); clf; hr = imagesc_(Ir); title('Select reference ROI');
                er = imellipse(gca); ex_r_ = createMask(er,hr); 
            end
            nbraggs=2; Ik = fftshift(fftn(Ir));
            if contains(flag,'interactive') || isempty(ex_g_)
                figure(2); clf; hk = imagesc_(log(abs(Ik).^2)); title('Zoom in on the image'); 
                input('Press enter when done'); title(sprintf('Select %i reflections',nbraggs));
                ex_g_ = cell(1,nbraggs); for q = 1:nbraggs; hg = imellipse(gca); ex_g_{q} = createMask(hg,hk); end
            end
            [n,m] = size(Ir); fft_k_ = @(n,Dx) ([0:(n-1)]-floor(n/2))./Dx;
            [r{1:2}] = ndgrid([1:n], [1:m]);
            [k{1:2}] = ndgrid(fft_k_(n,r{1}(end)-r{1}(1)), fft_k_(m,r{2}(end)-r{2}(1))); 
            ndims=2; K = zeros(ndims,nbraggs); G = zeros(ndims,nbraggs,n,m); 
            for q = 1:nbraggs
                Rk_r = fftshift(fftn(Ir.*ex_r_)).*ex_g_{q}; [~,i]=max(abs(Rk_r(:))); [i,j]=ind2sub([n,m],i);
                Ik_g = Ik.*ex_g_{q}; Ik_g = circshift(circshift(Ik_g,1-i,1),1-j,2);
                Ir_g = ifftn(fftshift(Ik_g)); Ar_g = abs(Ir_g); Pr_g = angle(Ir_g); 
                [T{2:-1:1}] = gradient(exp(1i.*Pr_g)); dPr_g = cat(3,T{:}); dPr_g = imag(exp(-1i.*Pr_g).*dPr_g);
                K(:,q)     = [k{1}(i,j);k{2}(i,j)];
                G(:,q,:,:) = permute(dPr_g,[3,4,1,2]);
            end  
            A = inv(K).';
            E = am_lib.matmul_(A,permute(G,[2,1,3,4]))/(-2*pi);
            if     contains(flag,'rot:image')
                ex_rotated_ = abs(Ir)<eps; ex_rotated_ = permute(ex_rotated_,[3,4,1,2]);
                ex_rotated_ = repmat(ex_rotated_,[2,2,1,1]); E(ex_rotated_) = NaN;
            elseif contains(flag,'rot:field')
                R = rotz(th); E = am_lib.matmul_(am_lib.matmul_(R(1:2,1:2),E),R(1:2,1:2).');
            end
            W = (E - permute(E,[2,1,3,4]))/2;
            E = (E + permute(E,[2,1,3,4]))/2;
            if contains(flag,'activate'); activate_ = @(E) 1./(E+1)-1; else; activate_ = @(E) E; end
            E = activate_(E); W = activate_(W);
        end