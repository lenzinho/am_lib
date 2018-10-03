classdef am_group < matlab.mixin.Copyable % everything is modified implicitly by reference (this is done so that am_group can host am_group subgroups)

    properties (Constant)
        tiny = 1E-8; 
    end

    properties
        T  = []; % group type (space/point/translational/[empty])
        nSs= []; % order of the group
        % representation
        GR = []; % ground representation (for subgroups)
        S  = []; % input representation
        % properties
        MT = []; % multiplication table
        CT = []; % character table
        % indices of important symmetries
        E  = []; % index of identity
        N  = []; % symmetry generators N(:,ngens)
        I  = []; % symmetry inverses
        C  = []; % symmetry conjugates
        % subgroup properties
        isproper = [];
        isnormal = [];
        nHs= []; % number of subgroups
        H  = []; % subgroup
        Rc = []; % subgroup coset representatives
        Hc = []; % subgroup conjugates
        Sc = []; % subgroup stabilizers
        Oc = []; % subgroup orbit (Rc x s2g)
        Qc = []; % subgroup quotient group index (for normal subgroups only; 0 if group is not a normal group)
        s2g= []; % subgroup-to-group symmetry indexing
    end

    methods (Static) % demo

        function demo_get_crystal_field()
            clear;clc;
            G = am_group();
            G.define('pg','o_h');
            G.get_group_properties();
            A=G.expand_crystal_field([1:6]);
        end
        
    end

    methods % functions that operate on and return the class
        
        function         define(G,P,S) % define(S), define('MT',MT), define('pg',pg_id), define('sg',sg_id), define('dg',dg_id)
            switch upper(P)
                case 'MT'
                    MT = S; S = [];
                case 'PG'
                    if ischar(S) || numel(S)==1;     S = am_group.generate_pg(S); end
                    MT = am_group.get_multiplication_table(S,am_group.tiny);
                case 'DG'
                    if ischar(S) || numel(S)==1; [~,S] = am_group.generate_pg(S); end
                    MT = am_group.get_multiplication_table(S,am_group.tiny);
                case 'SG'
                    if ischar(S) || numel(S)==1;     S = am_group.generate_sg(S); end
                    MT = am_group.get_multiplication_table(S,am_group.tiny);
                case 'CG'
                    MT = am_group.get_multiplication_table(S,am_group.tiny);
                otherwise
                    error('invalid input');
            end
            G.S = S; G.T = P; G.MT = MT; G.get_order(); G.get_identity(); G.get_symmetry_inverses();

        end
        
        function         relabel_symmetries(G,fwd)
            % sort
            [~,fwd] = sort(fwd);
            % relabel multiplication table
            rev(fwd) = [1:G.nSs]; 
            % two dimensional (symmetry indicies)
            for i_ = {'MT'}
                if ~isempty(G.(i_{:}))
                    G.(i_{:})(:,:) = rev(G.(i_{:})(fwd,fwd)); 
                end
            end
            % one dimensional (symmetry indicies)
            for i_ = {'I'}
                if ~isempty(G.(i_{:}))
                    G.(i_{:})(:) = rev(G.(i_{:})(fwd)); 
                end
            end
            % three dimensional (symmetry indicies)
            for i_ = {'S'}
                if ~isempty(G.(i_{:}))
                    G.(i_{:})(:,:,:) = rev(G.(i_{:})(:,:,fwd)); 
                end
            end
            % one dimensional (not symmetry indices)
            for i_ = {'C'}
                if ~isempty(G.(i_{:}))
                    G.(i_{:})(:) = G.(i_{:})(fwd); 
                end
            end
            % position irrelevant
            for i_ = {'N','E'}
                if ~isempty(G.(i_{:}))
                    G.(i_{:})(:) = rev(G.(i_{:})(:));
                end
            end
            % subgroups (position irrelevant)
            for i  = 1:G.nHs
            for i_ = {'s2g'}
                if ~isempty(G.H(i).(i_{:}))
                    G.H(i).(i_{:})(:) = rev(G.H(i).(i_{:})(:)); 
                end
            end
            end
        end

        function         relabel_subgroups(G,fwd)
            % sort
            [~,fwd] = sort(fwd);
            % relabel multiplication table
            rev(fwd) = [1:G.nHs]; 
            % one dimensional (not symmetry indices)
            if ~isempty(G.H); G.H = G.H(fwd); end
            % one dimensional (symmetry indices)
            for i = 1:G.nHs
            for i_ = {'Hc'}
                if ~isempty(G.H(i).(i_{:}))
                    G.H(i).(i_{:})(:) = rev(G.H(i).(i_{:}));
                end
            end
            end
        end

        function         get_group_properties(G)
            G.get_symmetry_conjugates();
            G.get_symmetry_generators();
            G.get_subgroups();
            G.get_subgroup_conjugates();
            G.identify_normal_subgroups();
            G.identify_proper_subgroups();
            G.get_subgroup_quotients();
        end
        
        function         get_order(G)
            G.nSs = size(G.MT,1);
        end
        
        function         get_identity(G)
            G.E = find(all(G.MT==[1:G.nSs].',1));
        end
        
        function         get_symmetry_inverses(G)
            if isempty(G.E); G = G.get_identity(); end
            G.I = [G.MT==G.E]*[1:G.nSs].';
        end

        function         get_symmetry_conjugates(G)
            %
            % for AX = XB, if elements A and B are conjugate pairs for some other element X in the group,  they are in the same class
            %
            if isempty(G.E); G.get_identity(); end
            if isempty(G.I); G.get_inverses(); end
            
            % initialize
            C_ = zeros(G.nSs,1); C = zeros(G.nSs,1); k = 0;
            % determine conjugacy classes
            for i = 1:G.nSs
            if C(i)==0
                k=k+1;
                % conjugate each element with all other group elements
                % A = X(j) * B * X(j)^-1
                for j = 1:G.nSs
                    C_(j) = G.MT(j,G.MT(i,G.I(j)));
                end
                % for each subgroup element created by conjugation find the corresponding index of the element in the group
                % in order to save the class class_id number
                for j = 1:G.nSs
                    C( C_(j) ) = k;
                end
            end
            end
            % relabel classes based on how many elements each class has
            G.C = am_lib.reindex_using_occurances(C);
        end

        function         get_symmetry_generators(G)
            % returns all possible sets of generators with ngens elements
            n = 0; N = [];
            for ngens = 1:G.nSs
                % generate a list of all possible generators (incldues order of generator for completness)
                P = am_lib.perm_norep_(G.nSs,ngens); nPs = size(P,2);
                % loop over the list one set of generators at a time, checking whether the entire multiplication table can be generated
                for i = 1:nPs
                    if numel( G.expand_multiplication_table(P(:,i)) ) == G.nSs
                        n=n+1; N(:,n) = P(:,i);
                    end
                end
                % smallest generating sets found 
                if n~=0; break; end
            end
            % save
            G.N = N;
        end

        function         get_subgroups(G)
            % G. Butler, Fundamental Algorithms for Permutation Groups (2014), p 46.
            G.H = am_group(); h = 1; level = 1;
            % initialize identity group
            G.H(1).define('MT',1); G.H(1).s2g=G.E; 
            % loop over layers
            for i = 2:numel(factor(G.nSs))+1
                for j = find(level==i-1) % loop over subgroups in level i-1
                   for g = find( ~G.contains(G.H(j),[1:G.nSs]) ) % loop over elements not in this subgroup
                       if G.stabilizes(G.H(j),g)  % select only elements that stabilize the group
                           s2g = G.expand_multiplication_table([g;G.H(j).s2g]);
                           ex_ = [G.H.nSs]==numel(s2g);
                           if ~any(ex_) || ~any(all([G.H(ex_).s2g]==s2g(:),1))
                               h = h + 1; level(h) = i; 
                               % add group
                               fwd(s2g)=[1:numel(s2g)]; G.H(h).s2g=s2g;
                               G.H(h).define('MT', fwd(G.MT(s2g,s2g))); 
                           end
                       end
                   end
                end
            end
            G.nHs = h;
            
            % check for uniqueness
            % X = arrayfun(@(x) any(x.s2g.'==[1:G.nSs].',2), G.H,'uniformoutput',false); size([X{:}]), size(am_lib.uniquec_(double([X{:}])));
        end
        
        function         get_subgroup_coset_representatives(G,H)
            switch nargin
                case 1
                    for i = 1:G.nHs
                        G.get_coset_representatives(G.H(i));
                    end
                case 2
                    % get coset representative
                    ex_ = true(1,G.nSs); Rc = zeros(1,G.nSs); k=0;
                    for i = [G.E,1:G.nSs]; if ex_(i)
                        k=k+1; Rc(k) = i;
                        ex_(G.MT(i,H.s2g)) = false; 
                    end; end
                    H.Rc = Rc(1:k);
                otherwise
                    error('invalid inputs');
            end
        end

        function         get_subgroup_stabilizers(G,H)
            switch nargin
                case 1
                    for i = 1:G.nHs
                        G.get_subgroup_stabilizers(G.H(i));
                    end
                case 2
                    % loop over coset representatives
                    s2g = false(1,G.nSs);
                    for j = 1:G.nSs
                        % stabilizer?
                        if all(sort(G.MT(G.MT(j,H.s2g),G.I(j)))==H.s2g(:)); s2g(j) = true; end
                    end
                    % find matching subgroup
                    s2g = find(s2g); o_ = find([G.H.nSs]==numel(s2g)); o_ = o_(all([G.H(o_).s2g]==s2g(:),1));
                    % save
                    if ~isempty(o_)
                        H.Sc=o_;
                    else
                        error('failed to find stabilizer subgroup'); 
                    end
                otherwise
                    error('invalid input');
            end
        end

        function         get_subgroup_conjugates(G,H)
            switch nargin
                case 1
                    for i = 1:G.nHs
                        G.get_subgroup_conjugates(G.H(i));
                    end
                case 2
                    % get coset representatives
                    if isempty(H.Rc); G.get_subgroup_coset_representatives(H); end
                    %
                    Hc = false(1,H.nHs);
                    % loop over coset representatives
                    for j = [G.E,H.Rc]
                        % operate on the subgroup with the coset rep
                        s2g = sort(G.MT(G.MT(j,H.s2g),G.I(j)));
                        % search over subgroups with the same order ...
                        o_ = find([G.H.nSs]==numel(s2g));
                        % to find matching subgroups
                        o_ = o_(all([G.H(o_).s2g]==s2g,1));
                        if ~isempty(o_)
                            Hc(o_)=true;
                        else
                            error('failed to find conjguate subgroup'); 
                        end
                    end
                    % update
                    H.Hc = find(Hc);
                otherwise
                    error('invalid input');
            end
        end
        
        function         get_subgroup_quotients(G,H)
            % J. S. Lomont, Applications of Finite Groups (1959), p 21.
            switch nargin
                case 1
                    for i = find([G.H.isnormal])
                        G.get_subgroup_quotients(G.H(i)); 
                    end
                    % confirm
                    if ~all( [G.H([G.H([G.H.isnormal]).Qc]).nSs].*[G.H([G.H.isnormal]).nSs] == G.nSs )
                        error('failed to find quotient subgroup');
                    end
                case 2
                    % get quotient group: cosets are elevated to the status of symmetries
                    cosets = sort(G.MT(H.Rc,H.s2g).'); nCs = size(cosets,2); MT = zeros(nCs,nCs);
                    for i = 1:nCs; for j = 1:nCs
                        MT(i,j) = find(all(sort(G.MT(cosets(:,i),cosets(1,j)))==cosets,1));
                    end; end
                    % create factor group
                    Q = am_group(); Q.define('MT',MT);
                    % find a subgroup that is isomorphic to this multiplication table
                    for o_ = find([G.H.nSs]==nCs)
                        if ~isempty(get_bijections(Q,G.H(o_))); break; end
                    end
                    % save
                    if isempty(o_) || H.nSs * G.H(o_).nSs ~= G.nSs
                        error('failed to find quotient subgroup'); 
                    else
                        H.Qc = o_;
                    end
            end
        end

        function         get_subgroup_groundreps(G,H)
            switch nargin
                case 1
                    for i = 1:G.nHs
                        G.get_subgroup_groundreps(G.H(i));
                    end
                case 2
                    n = numel(H.Rc); GR = zeros(n,n,G.nSs); s2g = permute(H.s2g,[3,2,1]);
                    for i = 1:G.nSs
                        GR(:,:,i) = sum( ( G.MT(G.multiply(G.I(H.Rc),i),H.Rc) == s2g ) .* s2g , 3);
                    end
                    H.GR = GR;
            end
        end
        
        function         identify_normal_subgroups(G,H) % identify subgroup that are normal (invariant, self-conjugate)
            switch nargin
                case 1
                    G.identify_normal_subgroups(G.H);
                case 2
                    for i = 1:numel(H)
                        if isempty(G.H(i).Hc)
                            G.get_subgroup_conjugates(H(i))
                        end
                        H(i).isnormal = numel(G.H(i).Hc)==1;
                    end
            end
        end
        
        function         identify_proper_subgroups(G,H)
            switch nargin
                case 1
                    ex_ = ~( [G.H.nSs]==1 |  [G.H.nSs]==G.nSs);
                    ex_ = num2cell(ex_); [G.H.isproper] = ex_{:};
                case 2
                    ex_ = ~(   [H.nSs]==1 |    [H.nSs]==G.nSs);
                    ex_ = num2cell(ex_); [H.isproper] = ex_{:};
            end
        end
        
        function         plot_cayley_graph(G)
            if isempty(G.N); G = G.get_generators_and_subgroups(); end
            i=[]; j=[]; x = G.MT(1,:); flatten_ = @(x) x(:);
            for k = 1:numel(G.N(:,1))
                i = [i;flatten_(x)];
                j = [j;flatten_(G.MT(G.N(k,1),x))];
            end
            plot(digraph(i,j,'OmitSelfLoops'),'layout','subspace3');
        end

    end
   
    methods %(Access = protected) % internal stuff; functions for which the class is input, but not output

        function [x]   = multiply(G,u,v)
            % x = u * v
            [n] = size(u); [m] = size(v);
            if numel(u) == 1
                u = repmat(u,size(v));
            elseif numel(v) == 1
                v = repmat(v,size(u));
            elseif ~all(n == m)
                error('outer products not yet implemented');
            end
            x = G.MT(sub2ind([G.nSs,G.nSs], u(:),v(:)));
        end
        
        function [x]   = conjugate(G,u,v)
            % x = u * v * u^-1
            [n] = size(u); [m] = size(v);
            if numel(u) == 1
                u = repmat(u,size(v));
            elseif numel(v) == 1
                v = repmat(v,size(u));
            elseif ~all(n == m)
                error('outer products not yet implemented');
            end
            x = G.multiply( G.multiply(u,v), G.I(u) );
        end
        
        function [u]   = expand_multiplication_table(G,u)
			% initialize
            x = false(G.nSs,1); x([G.E;u(:)]) = true; u = false(G.nSs,1);
            % expand multiplication table
            while true
                % expand
                u(G.MT(x,x)) = true;
                % procede ...
                switch sum(u)
                    case G.nSs;  break; % whole group, generators found!
                    case sum(x); break; % subgroup found!
                    otherwise; x = u;   % update
                end
            end
            % convert back to indices
            u = find(u);
        end

        function [L]   = contains(G,H,u)
            L = any(u(:).'==H.s2g(:),1);
        end
        
        function [L]   = stabilizes(G,H,u)
            n = numel(u); L = false(1,n);
            for j = 1:n
                L(j) = all( sort(G.MT(G.MT(u(j),H.s2g),G.I(u(j)))) == H.s2g(:) );
            end
        end
        
        function [Hcc] = get_subgroup_conjugacy_classes(G)
            Hcc = zeros(1,G.nHs); k=0;
            for i = 1:G.nHs; if all(Hcc(G.H(i).Hc)==0); k=k+1; Hcc(G.H(i).Hc) = k; end; end
        end

        function [IR]  = get_irreducible_representations(G)
            % MATHEMATICS OF COMPUTATION, VOLUME 24, NUMBER 111, JULY, 1970
            % Computing Irreducible Representations of Groups
            % By John D. Dixon

            % get regular rep G by putting identity along diagonal of multiplciation table
            RR = G.get_regular_representations();
            
            % initialize decomposition loop
            U = eye(G.nSs); inds = ones(G.nSs,1); ninds = 1;

            % loop until irreps are fully decomposed
            while true
                % loop over cycle structures
                for j = 1:max(inds)
                    ex_ = inds==j; 
                    H = dixon_decomposition_( RR(ex_,ex_,:) ); [Vp,~] = eig(H,'vector');
                    for ig = 1:G.nSs
                        RR(ex_,ex_,ig) = Vp \ RR(ex_,ex_,ig) * Vp;
                    end
                    U(:,ex_) = (U(:,ex_)*Vp);
                end
                A = am_lib.merge_(sum(abs(RR),3)>am_lib.eps); inds = [1:size(A,1)]*A;
                if ninds == max(inds); break; else; ninds = max(inds); end
            end

            % get irreducible representations
            nirreps = max(inds); IR = cell(1,nirreps);
            for i = 1:nirreps
                IR{i} = RR(inds==i,inds==i,1:G.nSs);
            end
            
            % sort the irreps by size
            IR = IR(am_lib.rankc_(cellfun(@(x)size(x,1),IR)));
            
            function H = dixon_decomposition_(rr)
                nbases=size(rr,1);
                nsyms =size(rr,3);

                for r = 1:nbases
                for s = 1:nbases
                    Hrs = zeros(nbases);
                    if     r==s
                        Hrs(r,s) = 1;
                    elseif r>s
                        Hrs(r,s) = 1;
                        Hrs(s,r) = 1;
                    elseif r<s
                        Hrs(r,s) = sqrt(-1);
                        Hrs(s,r) =-sqrt(-1);
                    end

                    H(1:nbases,1:nbases) = 0;
                    for q = 1:nsyms
                        H = H + rr(:,:,q) \ Hrs * rr(:,:,q);
                    end
                    H = H / nsyms;

                    if any(abs( H(1,1)*eye(nbases)-H ) >am_lib.eps); return; end
                end
                end
                H = eye(nbases);
            end
        end

        function [RR]  = get_regular_representations(G)
            % get regular rep G by putting identity along diagonal of multiplciation table
            RR = double(am_lib.accessc_(G.MT,G.I)==permute([1:G.nSs],[1,3,2]));
        end

        function [N]   = expand_generators(G,N)
            % expands generators, recording the order in which the elements are generated.
            while true
                % expand
                N = unique([am_lib.flatten_(G.MT(N,N));N],'stable');
                % status?
                if numel(N) == G.nSs; break; end
            end
        end

        function [B]   = get_bijections(G,H)
            % for i = 1:size(B,2)
            %     G.MT - relabel_multiplication_table(H.MT,B(:,i))
            % end

            % compare symmetries
            if G.nSs ~= H.nSs; B = []; return; end
            % get generators
            if isempty(G.N); G.get_symmetry_generators(); end
            if isempty(H.N); H.get_symmetry_generators(); end
            % compare number of generators
            if any(size(G.N)~=size(H.N)); B = []; return; end

            % sorting if desired (changes order of bjc but does not affect final result)
            Gfwd = [1:G.nSs]; Hfwd = [1:H.nSs];

            % compare multiplication tables, looking for bijections between symmetries g and h, i.e. isomorphism between G and H
            j=0; B=[]; [Gg] = G.expand_generators(G.N(:,1));
            for i = 1:size(H.N,2)
                % expand generators
                Hg = H.expand_generators(H.N(:,i));
                % map Hg onto Gg
                fwd(Gg) = Hg;
                % check for isomorphism
                if am_lib.all_( G.MT == relabel_multiplication_table(H.MT,fwd) )
                    % save bijection which takes symmetries of H onto G in original order
                    j=j+1; B(Gfwd,j) = Hfwd(fwd);
                end
            end

            function MT = relabel_multiplication_table(MT,fwd)
                rev(fwd) = [1:size(MT,1)]; MT = rev(MT(fwd,fwd));
            end
        end

        function [eq]  = expand_crystal_field(G,order)
                % G.expand_crystal_field([1:4]) expands crystal field potential up to 4th order
                % eq = G.expand_crystal_field([1:4])
                % subs(eq,'z',0)
                syms x y z
                eq = 0;
                for n = order
                    % get unique indices (because multiplication operations commutes)
                    M = am_lib.kronpow_([x;y;z],n); 
                    % cannot take unique here
                    [M,~,j_]=unique(M,'stable'); J = (j_==[1:max(j_)]).'; 
                    % m_ = true(size(M)); % [~,m_] = unique(M); 
                    % loop over idntity and generators
                    W = sum( am_lib.kronpow_(G.S(:,:,G.N(:,1)),n) - eye(3^n), 3);
                    % solve
                    N = J*null(W,'r'); N = am_lib.frref_(N.').';
                    % expand
                    c = sym(sprintf('c%02i%s',n,repmat('_%d',1,1)),[1,size(N,2)], 'real').';
                    % collect
                    eq = eq + collect( M.'*N*c , c);
                end

        end

    end
     
    methods (Static) %, Access = protected) % static functions for which the class is neither input nor output

        
        function [CT,cc_id,ir_id] = get_character_table(IR)
            % R=generate_pg(32,true);
            % IR = get_irreducible_representations(R);
            % [CT,cc_id,ir_id] = get_character_table(IR);
            % print_character_table(R,CT,cc_id)
            
            import am_lib.*
            
            %
            nirreps = numel(IR);
            nsyms   = size(IR{1},3);
            
            % get character table
            CT = zeros(nirreps,nsyms);
            for i = 1:nsyms; for j = 1:nirreps
                CT(j,i) = trace(IR{j}(:,:,i));
            end; end
        
            % correct numerical error
            CT = wdv_(CT);

            % get irreducible irreps
            [CT,~,cc_id]=uniquec_(CT);

            % get irreducible classes
            [CT,~,ir_id]=uniquec_(CT.'); CT=CT.';
            
            % check 
            if ~isdiag_(CT'*CT)
                sym(CT)
                error('Character table is not orthogonal')
            end
        end

        function                print_character_table(R,CT,cc_id)
            % clear;clc
            % R=generate_pg(32,true);
            % IR = get_irreducible_representations(R);
            % [CT,cc_id,ir_id] = get_character_table(IR);
            % print_character_table(R,CT,cc_id)

            
            %     clear;clc; import am_dft.* am_lib.*
            %     [R,W]=generate_pg('o_h',false);
            %     % d = permute(det_(R),[3,1,2]);
            %     % j=1/2; [W] = get_wigner(j,d,'spherical');
            % 
            %     % (1i*eq_(d,-1) + eq_(d,1));
            %     % W = cat(3,1i*W,-1i*W);
            %     % W = cat(3,kron_(R,eye(2)),kron_(R,[0 1; 1 0]));
            %     % W = complete_group(W);
            %     % MT = am_dft.get_multiplication_table(W);
            %     % MTk = am_dft.get_multiplication_table(K);
            %     IR = get_irreducible_representations(W);
            %     [CT,cc_id,ir_id] = get_character_table( IR );
            %     % sym(CT)
            %     print_character_table(W,CT,cc_id)
            
            import am_dft.* am_lib.*
            
            % identify the prototypical symmetries
            [~,cc_rep]=unique(cc_id,'stable');

            % name of prototypical symmetries
            ps_name_long = get_long_ps_name(R);
                
            % get number of classes and number of elements in each class
            nclasses = numel(cc_rep); nelems = sum(cc_id(cc_rep)==cc_id',2).'; nsyms = sum(nelems);
            
            % calculate orbital characters
            [chi,J] = get_orbital_characters(R(:,:,cc_rep));

            % print character table% class label 
                print_table_(nclasses,wdv_(CT),chi,ps_name_long,nelems,cc_id,cc_rep)
            
            % DECOMPOSITIONS
                decomp_ = @(chi) wdv_( chi * (CT.*nelems).' / nsyms );
                
            % IRREP x IRREP decomposition
                decomp = decomp_( CT.^2 ); decomp = wdv_(decomp); nirreps = size(CT,1);
                for i = 1:nirreps; row_name{i} = sprintf('irr # %g',i); end
                for i = 1:nirreps; col_name{i} = sprintf('irr # %g',i); end
                print_table_decomposition_('irrep x irrep', row_name, col_name, decomp);
                
            % ORBITAL decomposition
                decomp = decomp_( chi ); decomp = wdv_(decomp);
                for i = 1:nirreps; row_name{i} = sprintf('irr # %g',i); end
                clear row_name; for i = 1:numel(J); row_name{i} = sprintf('%s',sym(J(i))); end
                clear col_name; for i = 1:nirreps;  col_name{i} = sprintf('irr # %g',i); end
                print_table_decomposition_('orbitals', row_name, col_name, decomp);
                
                
            function [chi,J] = get_orbital_characters(W)
                % identify double groups symmetries (negative ps_id = double valued)
                [ps_id,~,d] = am_dft.identify_point_symmetries(W);
                % get rotation angle
                alpha = am_lib.R_angle_(permute(d,[3,1,2]).*W);
                % compute character
                character_ = @(l,alpha) sin((l+1/2).*(alpha+1E-12))./sin((alpha+1E-12)./2); 
                chi_l = [0.0:5.0]; chi_O3 = character_(chi_l(:), alpha) .* (am_lib.eq_(d,-1).*1i + 1*am_lib.eq_(d,1));% .* (d).^(chi_l(:));% .* sign(ps_id);
                chi_j = [0.5:5.5]; chi_U2 = character_(chi_j(:), alpha) .* (d).^(chi_j(:)) .* sign(ps_id);
                J = [chi_l(:);chi_j(:)]; chi = [chi_O3;chi_U2]; chi = am_lib.wdv_(chi);
            end
            
            function print_table_(nclasses,CT,chi,ps_name_long,nelems,cc_id,cc_rep)
                fprintf('      '); for j = 1:nclasses; fprintf('%12s',['#',num2str(j)]); end; fprintf('\n');
                print_bar(nclasses);
                % print class name 
                fprintf('      '); for j = 1:nclasses; fprintf('%12s',ps_name_long{cc_rep(j)}); end; fprintf('\n');
                % print class elements
                fprintf('      '); for j = 1:nclasses; fprintf('%12i',nelems(j)); end; fprintf('\n');
                print_bar(nclasses);
                % print character table
                for l = 1:size(CT,1); fprintf('      '); fprintf('%12s',sym(CT(l,:))); fprintf('\n'); end
                % print SO(3) and SU(2) characters single and double groups
                print_bar(nclasses);
                for l = 1:size(chi,1); fprintf('      '); fprintf('%12s',sym(chi(l,:))); fprintf('\n'); end
                print_bar(nclasses);

                % symmetries in each class
                    fprintf('      '); fprintf('Symmetries in each classes:\n')
                    for j = 1:nclasses
                        fprintf('%12s:%s\n',['#',num2str(j)], sprintf(' %-12s',ps_name_long{cc_id==j}));
                    end
            end
            
            function print_table_decomposition_(table_name,row_name,col_name,decomp)
                [nrows,~] = size(decomp); 
                fprintf('\n');
                fprintf('      '); fprintf('%12s\n',table_name); 
                fprintf('      '); fprintf('%12s\n','==========='); 

                fprintf('      '); fprintf('%12s',''); fprintf('%12s',col_name{:}); fprintf('\n');
                print_bar(numel(col_name)+1)
                % fprintf('      '); fprintf('%12s',repmat(' -----------',1,ncols+1)); fprintf('\n');
                for k = 1:nrows
                    fprintf('        '); fprintf('%-10s',row_name{k}); fprintf('%12g',decomp(k,:)); fprintf('\n');
                end
            end
            
            function print_bar(nclasses)
                fprintf('      '); fprintf('%12s',repmat(' -----------',1,nclasses)); fprintf('\n');
            end
        end


        function [fwd,M]      = find_pointgroup_transformation(g,h,algo)
            % finds the permutation fwd and rotation M which transforms point group g into h:
            % matmul_(matmul_(M(:,:,k),g(:,:,:)),inv(M(:,:,k))) == h(:,:,bjc(:,k)) for any k

            import am_dft.* am_lib.*
            
            % algo = 1, returns only 1 bijection of many possibilities
            if nargin<3; algo=0; end

            % find bijection bjc: h -> g and generators for h
            G=get_multiplication_table(g);
            H=get_multiplication_table(h);
            [bjc,gen] = find_isomorphic_bijection(G,H);
            srt_ = rankc_(bjc); bjc = bjc(:,srt_); gen = gen(:,srt_);

            % check diagonalize symmetry values to make sure the symmetries are of the correct type
            % this removes the possibility that multiplication table may be obeyed with symmetries of differnet types
            nbjcs=size(bjc,2); ex_ = true(1,nbjcs); 
            g_ps_id = identify_point_symmetries(g); 
            h_ps_id = identify_point_symmetries(h);
            for j = 1:nbjcs;   ex_ = all(h_ps_id(bjc(:,j))==g_ps_id,1); end
            bjc = bjc(:,ex_); gen = gen(:,ex_); nbjcs=size(bjc,2);

            % find integer transformation matrix (try all possible bijections and
            % integer transformations with nonzero-determinants); in reality, not all
            % matrices are tried, instead: 
            %       first try the identity
            %       then try all matrices with elements in [-1,0,1] with determinant 1 
            %       then try all matrices with elements in [-2,-1,0,1,2] with determinant 1 
            q = 0; ngens = size(gen,1); 
            for m = 0:1
                % generate the possible transformation matrices
                switch m
                    case 0
                        % first try the identity
                        X = eye(3); nXs = 1;
                    otherwise
                        % then try other integer matrices with nonzero determinants
                        N=9; Q=[-m:m]; nQs=numel(Q);[Y{N:-1:1}]=ndgrid(1:nQs); X=reshape(Q(reshape(cat(N+1,Y{:}),[],N)).',3,3,[]);
                        ex_ = false(1,size(X,3)); for i = 1:size(X,3); ex_(i) = abs(det(X(:,:,i)))<am_lib.eps; end
                        X = X(:,:,~ex_); nXs = size(X,3);
                        % sort it to make it nice (put positive values first)
                        X=X(:,:,rankc_( -[max(reshape(X,9,[]));min(reshape(X,9,[]));sum(reshape(X,9,[]),1)] ));
                end
                % get inverse elements
                invX = zeros(3,3,nXs); for i = 1:nXs; invX(:,:,i) = inv(X(:,:,i)); end

                % find similarity transform which converts all symmetries from g to h by
                % checking each matrix generated above
                for k = 1:nXs; for j = 1:nbjcs
                    gencheck_ = true;
                    % check generators
                    for i = 1:ngens
                        if any(any(abs( X(:,:,k) * g(:,:,gen(i,j)) * inv(X(:,:,k)) - h(:,:,bjc(gen(i,j),j)) )>am_lib.eps ))
                            gencheck_ = false; break;
                        end
                    end
                    % check whether a match has been found
                    if gencheck_
                        % at this point a bijection permutation and transformation matrix have been
                        % found such that: matmul_(matmul_(X(:,:,k),g(:,:,:)),invX(:,:,k)) == h(:,:,bjc(:,j))
                        % save it 
                        q = q+1;
                        M(:,:,q) = X(:,:,k); 
                        fwd(:,q) = bjc(:,j);
                        % only one is needed
                        if algo==1
                            return; 
                        end
                    end
                end;end
            end
        end

        
    end
    
    % basic functions needed to create group object
    
    methods (Static)

        function S            = generate_sg(sg_code,from_memory)

            import am_lib.* am_dft.*

            if nargin<2; from_memory=false; end

            if from_memory
                % ~ 500x faster for large space groups
                S = load_sg_symmetries(sg_code);
                S = reshape(S,4,4,[]);
            else
                
                % load recipe
                recipe = get_recipe(sg_code);
                
                % allocate space
                S = zeros(4,4,sscanf(recipe(2),'%i')+double(recipe(1)=='1'));

                % mix ingredients
                nsyms = 0; k = 0;
                nsyms = nsyms+1; S(1:4,1:4,nsyms) = eye(4);
                %
                k = k+1;
                if recipe(1) == '1'
                    nsyms = nsyms+1; S(1:3,1:3,nsyms) = -eye(3); S(4,4,nsyms)=1;
                end
                k = k+1; ngens = sscanf(recipe(2),'%i');
                for i = 1:ngens
                    nsyms = nsyms+1;
                for j = 1:4
                    k=k+1;
                    switch j
                        case 1
                            S(1:3,1:3,nsyms) = decode_R(recipe(k)); S(4,4,nsyms)=1;
                        otherwise
                            S(j-1,4,nsyms) = mod(decode_T(recipe(k)),1);
                    end
                end
                end

                S = am_cell.complete_group(S);
            end

            function recipe  = get_recipe(sg_id)
                sg_recipe_database = {...
                    '000                                     ','100                                     ','01cOOO0                                 ', ...
                    '01cODO0                                 ','02aDDOcOOO0                             ','01jOOO0                                 ', ...
                    '01jOOD0                                 ','02aDDOjOOO0                             ','02aDDOjOOD0                             ', ...
                    '11cOOO0                                 ','11cODO0                                 ','12aDDOcOOO0                             ', ...
                    '11cOOD0                                 ','11cODD0                                 ','12aDDOcOOD0                             ', ...
                    '02bOOOcOOO0                             ','02bOODcOOD0                             ','02bOOOcDDO0                             ', ...
                    '02bDODcODD0                             ','03aDDObOODcOOD0                         ','03aDDObOOOcOOO0                         ', ...
                    '04aODDaDODbOOOcOOO0                     ','03aDDDbOOOcOOO0                         ','03aDDDbDODcODD0                         ', ...
                    '02bOOOjOOO0                             ','02bOODjOOD0                             ','02bOOOjOOD0                             ', ...
                    '02bOOOjDOO0                             ','02bOODjDOO0                             ','02bOOOjODD0                             ', ...
                    '02bDODjDOD0                             ','02bOOOjDDO0                             ','02bOODjDDO0                             ', ...
                    '02bOOOjDDD0                             ','03aDDObOOOjOOO0                         ','03aDDObOODjOOD0                         ', ...
                    '03aDDObOOOjOOD0                         ','03aODDbOOOjOOO0                         ','03aODDbOOOjODO0                         ', ...
                    '03aODDbOOOjDOO0                         ','03aODDbOOOjDDO0                         ','04aODDaDODbOOOjOOO0                     ', ...
                    '04aODDaDODbOOOjBBB0                     ','03aDDDbOOOjOOO0                         ','03aDDDbOOOjDDO0                         ', ...
                    '03aDDDbOOOjDOO0                         ','12bOOOcOOO0                             ','03bOOOcOOOhDDD1BBB                      ', ...
                    '12bOOOcOOD0                             ','03bOOOcOOOhDDO1BBO                      ','12bDOOcOOO0                             ', ...
                    '12bDOOcDDD0                             ','12bDODcDOD0                             ','12bDOOcOOD0                             ', ...
                    '12bOOOcDDO0                             ','12bDDOcODD0                             ','12bOODcODD0                             ', ...
                    '12bOOOcDDD0                             ','03bOOOcDDOhDDO1BBO                      ','12bDDDcOOD0                             ', ...
                    '12bDODcODD0                             ','12bDODcODO0                             ','13aDDObOODcOOD0                         ', ...
                    '13aDDObODDcODD0                         ','13aDDObOOOcOOO0                         ','13aDDObOOOcOOD0                         ', ...
                    '13aDDObODOcODO0                         ','04aDDObDDOcOOOhODD1OBB                  ','14aODDaDODbOOOcOOO0                     ', ...
                    '05aODDaDODbOOOcOOOhBBB1ZZZ              ','13aDDDbOOOcOOO0                         ','13aDDDbOOOcDDO0                         ', ...
                    '13aDDDbDODcODD0                         ','13aDDDbODOcODO0                         ','02bOOOgOOO0                             ', ...
                    '02bOODgOOB0                             ','02bOOOgOOD0                             ','02bOODgOOF0                             ', ...
                    '03aDDDbOOOgOOO0                         ','03aDDDbDDDgODB0                         ','02bOOOmOOO0                             ', ...
                    '03aDDDbOOOmOOO0                         ','12bOOOgOOO0                             ','12bOOOgOOD0                             ', ...
                    '03bOOOgDDOhDDO1YBO                      ','03bOOOgDDDhDDD1YYY                      ','13aDDDbOOOgOOO0                         ', ...
                    '04aDDDbDDDgODBhODB1OYZ                  ','03bOOOgOOOcOOO0                         ','03bOOOgDDOcDDO0                         ', ...
                    '03bOODgOOBcOOO0                         ','03bOODgDDBcDDB0                         ','03bOOOgOODcOOO0                         ', ...
                    '03bOOOgDDDcDDD0                         ','03bOODgOOFcOOO0                         ','03bOODgDDFcDDF0                         ', ...
                    '04aDDDbOOOgOOOcOOO0                     ','04aDDDbDDDgODBcDOF0                     ','03bOOOgOOOjOOO0                         ', ...
                    '03bOOOgOOOjDDO0                         ','03bOOOgOODjOOD0                         ','03bOOOgDDDjDDD0                         ', ...
                    '03bOOOgOOOjOOD0                         ','03bOOOgOOOjDDD0                         ','03bOOOgOODjOOO0                         ', ...
                    '03bOOOgOODjDDO0                         ','04aDDDbOOOgOOOjOOO0                     ','04aDDDbOOOgOOOjOOD0                     ', ...
                    '04aDDDbDDDgODBjOOO0                     ','04aDDDbDDDgODBjOOD0                     ','03bOOOmOOOcOOO0                         ', ...
                    '03bOOOmOOOcOOD0                         ','03bOOOmOOOcDDO0                         ','03bOOOmOOOcDDD0                         ', ...
                    '03bOOOmOOOjOOO0                         ','03bOOOmOOOjOOD0                         ','03bOOOmOOOjDDO0                         ', ...
                    '03bOOOmOOOjDDD0                         ','04aDDDbOOOmOOOjOOO0                     ','04aDDDbOOOmOOOjOOD0                     ', ...
                    '04aDDDbOOOmOOOcOOO0                     ','04aDDDbOOOmOOOcDOF0                     ','13bOOOgOOOcOOO0                         ', ...
                    '13bOOOgOOOcOOD0                         ','04bOOOgOOOcOOOhDDO1YYO                  ','04bOOOgOOOcOOOhDDD1YYY                  ', ...
                    '13bOOOgOOOcDDO0                         ','13bOOOgOOOcDDD0                         ','04bOOOgDDOcDDOhDDO1YBO                  ', ...
                    '04bOOOgDDOcDDDhDDO1YBO                  ','13bOOOgOODcOOO0                         ','13bOOOgOODcOOD0                         ', ...
                    '04bOOOgDDDcOODhDDD1YBY                  ','04bOOOgDDDcOOOhDDD1YBY                  ','13bOOOgOODcDDO0                         ', ...
                    '13bOOOgDDDcDDD0                         ','04bOOOgDDDcDDDhDDD1YBY                  ','04bOOOgDDDcDDOhDDD1YBY                  ', ...
                    '14aDDDbOOOgOOOcOOO0                     ','14aDDDbOOOgOOOcOOD0                     ','05aDDDbDDDgODBcDOFhODB1OBZ              ', ...
                    '05aDDDbDDDgODBcDOBhODB1OBZ              ','01nOOO0                                 ','01nOOC0                                 ', ...
                    '01nOOE0                                 ','02aECCnOOO0                             ','11nOOO0                                 ', ...
                    '12aECCnOOO0                             ','02nOOOfOOO0                             ','02nOOOeOOO0                             ', ...
                    '02nOOCfOOE0                             ','02nOOCeOOO0                             ','02nOOEfOOC0                             ', ...
                    '02nOOEeOOO0                             ','03aECCnOOOeOOO0                         ','02nOOOkOOO0                             ', ...
                    '02nOOOlOOO0                             ','02nOOOkOOD0                             ','02nOOOlOOD0                             ', ...
                    '03aECCnOOOkOOO0                         ','03aECCnOOOkOOD0                         ','12nOOOfOOO0                             ', ...
                    '12nOOOfOOD0                             ','12nOOOeOOO0                             ','12nOOOeOOD0                             ', ...
                    '13aECCnOOOeOOO0                         ','13aECCnOOOeOOD0                         ','02nOOObOOO0                             ', ...
                    '02nOOCbOOD0                             ','02nOOEbOOD0                             ','02nOOEbOOO0                             ', ...
                    '02nOOCbOOO0                             ','02nOOObOOD0                             ','02nOOOiOOO0                             ', ...
                    '12nOOObOOO0                             ','12nOOObOOD0                             ','03nOOObOOOeOOO0                         ', ...
                    '03nOOCbOODeOOC0                         ','03nOOEbOODeOOE0                         ','03nOOEbOOOeOOE0                         ', ...
                    '03nOOCbOOOeOOC0                         ','03nOOObOODeOOO0                         ','03nOOObOOOkOOO0                         ', ...
                    '03nOOObOOOkOOD0                         ','03nOOObOODkOOD0                         ','03nOOObOODkOOO0                         ', ...
                    '03nOOOiOOOkOOO0                         ','03nOOOiOODkOOD0                         ','03nOOOiOOOeOOO0                         ', ...
                    '03nOOOiOODeOOO0                         ','13nOOObOOOeOOO0                         ','13nOOObOOOeOOD0                         ', ...
                    '13nOOObOODeOOD0                         ','13nOOObOODeOOO0                         ','03bOOOcOOOdOOO0                         ', ...
                    '05aODDaDODbOOOcOOOdOOO0                 ','04aDDDbOOOcOOOdOOO0                     ','03bDODcODDdOOO0                         ', ...
                    '04aDDDbDODcODDdOOO0                     ','13bOOOcOOOdOOO0                         ','04bOOOcOOOdOOOhDDD1YYY                  ', ...
                    '15aODDaDODbOOOcOOOdOOO0                 ','06aODDaDODbOOOcOOOdOOOhBBB1ZZZ          ','14aDDDbOOOcOOOdOOO0                     ', ...
                    '13bDODcODDdOOO0                         ','14aDDDbDODcODDdOOO0                     ','04bOOOcOOOdOOOeOOO0                     ', ...
                    '04bOOOcOOOdOOOeDDD0                     ','06aODDaDODbOOOcOOOdOOOeOOO0             ','06aODDaDODbODDcDDOdOOOeFBF0             ', ...
                    '05aDDDbOOOcOOOdOOOeOOO0                 ','04bDODcODDdOOOeBFF0                     ','04bDODcODDdOOOeFBB0                     ', ...
                    '05aDDDbDODcODDdOOOeFBB0                 ','04bOOOcOOOdOOOlOOO0                     ','06aODDaDODbOOOcOOOdOOOlOOO0             ', ...
                    '05aDDDbOOOcOOOdOOOlOOO0                 ','04bOOOcOOOdOOOlDDD0                     ','06aODDaDODbOOOcOOOdOOOlDDD0             ', ...
                    '05aDDDbDODcODDdOOOlBBB0                 ','14bOOOcOOOdOOOeOOO0                     ','05bOOOcOOOdOOOeOOOhDDD1YYY              ', ...
                    '14bOOOcOOOdOOOeDDD0                     ','05bOOOcOOOdOOOeDDDhDDD1YYY              ','16aODDaDODbOOOcOOOdOOOeOOO0             ', ...
                    '16aODDaDODbOOOcOOOdOOOeDDD0             ','07aODDaDODbODDcDDOdOOOeFBFhBBB1ZZZ      ','07aODDaDODbODDcDDOdOOOeFBFhFFF1XXX      ', ...
                    '15aDDDbOOOcOOOdOOOeOOO0                 ','15aDDDbDODcODDdOOOeFBB0                 ','01dOOO0                                 ', ...
                    '11dOOO0                                 ','02dOOOfOOO0                             ','02dOOOlOOO0                             ', ...
                    '02dOOOlDDD0                             ','12dOOOfOOO0                             ','12dOOOfDDD0                             '};
                recipe = sg_recipe_database{sg_id};
            end
            function R       = decode_R(str_r)
                switch str_r
                case 'a'; R = [  1  0  0;  0  1  0;  0  0  1]; % a
                case 'b'; R = [ -1  0  0;  0 -1  0;  0  0  1]; % b
                case 'c'; R = [ -1  0  0;  0  1  0;  0  0 -1]; % c
                case 'd'; R = [  0  0  1;  1  0  0;  0  1  0]; % d
                case 'e'; R = [  0  1  0;  1  0  0;  0  0 -1]; % e
                case 'f'; R = [  0 -1  0; -1  0  0;  0  0 -1]; % f
                case 'g'; R = [  0 -1  0;  1  0  0;  0  0  1]; % g
                case 'h'; R = [ -1  0  0;  0 -1  0;  0  0 -1]; % h
                case 'i'; R = [  1  0  0;  0  1  0;  0  0 -1]; % i
                case 'j'; R = [  1  0  0;  0 -1  0;  0  0  1]; % j
                case 'k'; R = [  0 -1  0; -1  0  0;  0  0  1]; % k
                case 'l'; R = [  0  1  0;  1  0  0;  0  0  1]; % l
                case 'm'; R = [  0  1  0; -1  0  0;  0  0 -1]; % m
                case 'n'; R = [  0 -1  0;  1 -1  0;  0  0  1]; % n
                end
            end
            function T       = decode_T(str_t)
                switch str_t
                    case 'A'; T = 1/6; % A
                    case 'B'; T = 1/4; % B
                    case 'C'; T = 1/3; % C
                    case 'D'; T = 1/2; % D
                    case 'E'; T = 2/3; % E
                    case 'F'; T = 3/4; % F
                    case 'G'; T = 5/6; % G
                    case 'O'; T =   0; % O
                    case 'X'; T =-3/8; % X
                    case 'Y'; T =-1/4; % Y
                    case 'Z'; T =-1/8; % Z
                end
            end
        end

        function [R,W]        = generate_pg(pg_code,from_memory)
            %     1.   c_1        9.   c_3        17. d_4       25. c_6v
            %     2.   s_2        10.  s_6        18. c_4v      26. d_3h
            %     3.   c_2        11.  d_3        19. d_2d      27. d_6h
            %     4.   c_1h       12.  c_3v       20. d_4h      28. t
            %     5.   c_2h       13.  d_3d       21. c_6       29. t_h
            %     6.   d_2        14.  c_4        22. c_3h      30. o
            %     7.   c_2v       15.  s_4        23. c_6h      31. t_d
            %     8.   d_2h       16.  c_4h       24. d_6       32. o_h
            
            import am_lib.* am_dft.*
            
            if nargin<2; from_memory=true; end
            
            if ischar(pg_code)
                pg={'c_1' ,'s_2' ,'c_2' ,'c_1h','c_2h','d_2' ,'c_2v','d_2h', ...
                    'c_3' ,'s_6' ,'d_3' ,'c_3v','d_3d','c_4' ,'s_4' ,'c_4h', ...
                    'd_4' ,'c_4v','d_2d','d_4h','c_6' ,'c_3h','c_6h','d_6' , ...
                    'c_6v','d_3h','d_6h','t'   ,'t_h' ,'o'   ,'t_d' ,'o_h'};
                pg_code = find(string(pg)==pg_code); %#ok<STRCLQT>
            end

            if from_memory
                % generated with this code:
                %
                % clear;clc;import am_dft.*
                % for i = 1:32
                %     R{i} = generate_pg(i,false);
                %     x(i) = identify_pointgroup(R{i});
                %     d{i} = decode_pg(x(i));
                %     fprintf('case %i; R = [',i); fprintf('%s',sym(R{i}(1))); fprintf(',%s',sym(R{i}(2:end))); fprintf(']; %% %s\n',d{i});
                % end
                % for i = 1:32
                %     [~,W{i}] = generate_pg(i,false);
                %     fprintf('case %i; W = [',i); fprintf('%s',sym(W{i}(1))); fprintf(',%s',sym(W{i}(2:end))); fprintf('];\n');
                % end
                %
                
                % single-valued representation
                switch pg_code
                case 1;  R = [1,0,0,0,1,0,0,0,1]; % c_1
                case 2;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1]; % s_2
                case 3;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1]; % c_2
                case 4;  R = [1,0,0,0,1,0,0,0,1,1,0,0,0,-1,0,0,0,1]; % c_1h
                case 5;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,1]; % c_2h
                case 6;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1]; % d_2
                case 7;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % c_2v
                case 8;  R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % d_2h
                case 9;  R = [1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,1,0,0,1,0,0,0,1,1,0,0]; % c_3
                case 10; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,1,0,0,0,1,1,0,0,0,-1,0,0,0,-1,-1,0,0]; % s_6
                case 11; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,0,-1,0,-1,0,0,0,1,1,0,0,0,1,0,0,-1,0,-1,0,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,0,1,0,0,0,1,1,0,0]; % d_3
                case 12; R = [1,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0]; % c_3v
                case 13; R = [1,0,0,0,1,0,0,0,1,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,0,0,-1,-1,0,0,0,-1,0,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,0,1,0,0,0,1,1,0,0,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-1/3,2/3,2/3,0,-1,0,0,0,-1,-1,0,0]; % d_3d
                case 14; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,1]; % c_4
                case 15; R = [1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,1,0,-1,0,0,0,0,-1]; % s_4
                case 16; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1]; % c_4h
                case 17; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,1,1,0,0,0,-1,0,0,0,-1]; % d_4
                case 18; R = [1,0,0,0,1,0,0,0,1,0,1,0,-1,0,0,0,0,1,1,0,0,0,-1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,1,0,0,0,0,1,0,-1,0,-1,0,0,0,0,1,0,-1,0,1,0,0,0,0,1,-1,0,0,0,1,0,0,0,1]; % c_4v
                case 19; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,-1,0,1,0,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,1,0,1,0,-1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,-1]; % d_2d
                case 20; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,1,0,1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,1,0,0,0,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,0,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1]; % d_4h
                case 21; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,0,0,1,1,0,0,0,1,0,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,1,0,0,0,1,1,0,0,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3]; % c_6
                case 22; R = [1,0,0,0,1,0,0,0,1,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,0,0,1,1,0,0,0,1,0,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,0,1,0,0,0,1,1,0,0,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % c_3h
                case 23; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,0,-1,-1,0,0,0,-1,0,0,1,0,0,0,1,1,0,0,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,0,-1,0,0,0,-1,-1,0,0,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % c_6h
                case 24; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1,0,0,0,0,-1,0,-1,0,0,0,1,1,0,0,0,1,0,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,-1,0,-1,0,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,0,1,0,0,0,1,1,0,0,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3]; % d_6
                case 25; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3]; % c_6v
                case 26; R = [1,0,0,0,1,0,0,0,1,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % d_3h
                case 27; R = [1,0,0,0,1,0,0,0,1,-1/3,2/3,2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1,0,0,0,0,-1,0,-1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,-1,0,0,0,-1,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,2/3,-1/3,2/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,0,-1,0,-1,0,0,0,0,-1,1,0,0,0,0,1,0,1,0,0,0,-1,0,-1,0,-1,0,0,0,1,0,0,0,1,1,0,0,0,0,-1,-1,0,0,0,-1,0,-2/3,1/3,-2/3,1/3,-2/3,-2/3,-2/3,-2/3,1/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-2/3,-2/3,1/3,-2/3,1/3,-2/3,1/3,-2/3,-2/3,2/3,-1/3,2/3,2/3,2/3,-1/3,-1/3,2/3,2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3,-2/3,1/3,-2/3,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,-1,0,0,0,-1,-1,0,0,2/3,-1/3,2/3,-1/3,2/3,2/3,2/3,2/3,-1/3,2/3,2/3,-1/3,2/3,-1/3,2/3,-1/3,2/3,2/3,-2/3,1/3,-2/3,-2/3,-2/3,1/3,1/3,-2/3,-2/3]; % d_6h
                case 28; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,0,0,1,1,0,0,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,1,0,0,-1,0,0,0,-1,1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,0,-1,0,0,0,-1]; % t
                case 29; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,0,0,1,1,0,0,-1,0,0,0,-1,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,1,0,0,0,-1,0,0,0,1,0,-1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,1,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,0,-1,1,0,0,0,0,-1,-1,0,0,0,1,0,0,-1,0,0,0,1,1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,0,0,-1,1,0,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,-1,0,0,0,-1,0,0,1,0,0,0,1,-1,0,0,0,0,1,1,0,0,0,-1,0,-1,0,0,0,-1,0,0,0,1,0,0,-1,1,0,0,0,1,0,1,0,0,0,-1,0,0,0,-1,0,0,1,-1,0,0,0,1,0,1,0,0,0,1,0,0,0,-1,-1,0,0,0,1,0,0,0,1]; % t_h
                case 30; R = [1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,-1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,0,0,-1,0,1,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,0,-1,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,1,0,0,0,-1,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,0,0,-1,0,-1,0,-1,0,0,-1,0,0,0,0,-1,0,-1,0,0,-1,0,-1,0,0,0,0,-1]; % o
                case 31; R = [1,0,0,0,1,0,0,0,1,-1,0,0,0,-1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,1,0,0,0,-1,0,0,0,-1,0,-1,0,-1,0,0,0,0,1,0,-1,0,0,0,1,-1,0,0,0,1,0,-1,0,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,0,-1,0,1,0,0,0,0,-1,1,0,0,0,0,1,0,1,0,0,-1,0,0,0,-1,1,0,0,0,0,1,0,1,0,1,0,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,-1,0,0,0,0,-1,0,1,0,0,0,-1,0,1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,1,0,0,0,0,-1,0,-1,0,0,0,-1,0,-1,0,1,0,0,0,0,-1,1,0,0,0,-1,0]; % t_d
                case 32; R = [1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,1,0,-1,0,0,0,0,1,-1,0,0,0,-1,0,0,0,-1,0,0,1,1,0,0,0,1,0,0,0,1,0,-1,0,1,0,0,0,-1,0,0,0,-1,-1,0,0,-1,0,0,0,0,1,0,1,0,-1,0,0,0,-1,0,0,0,1,0,-1,0,1,0,0,0,0,-1,1,0,0,0,0,-1,0,1,0,0,0,-1,-1,0,0,0,-1,0,0,-1,0,1,0,0,0,0,1,0,-1,0,0,0,-1,1,0,0,0,0,-1,0,1,0,-1,0,0,0,0,1,0,1,0,-1,0,0,1,0,0,0,0,-1,0,-1,0,0,-1,0,0,0,1,-1,0,0,1,0,0,0,1,0,0,0,-1,0,0,-1,0,1,0,1,0,0,0,0,-1,-1,0,0,0,1,0,-1,0,0,0,0,1,0,-1,0,1,0,0,0,0,1,0,-1,0,0,1,0,-1,0,0,0,0,-1,0,0,-1,1,0,0,0,-1,0,0,1,0,0,0,1,-1,0,0,0,1,0,0,0,-1,-1,0,0,0,0,-1,0,-1,0,1,0,0,0,0,1,-1,0,0,0,-1,0,0,1,0,0,0,-1,1,0,0,0,1,0,1,0,0,0,0,-1,0,0,1,0,-1,0,-1,0,0,-1,0,0,0,1,0,0,0,-1,0,0,1,1,0,0,0,-1,0,-1,0,0,0,0,-1,0,1,0,1,0,0,0,-1,0,0,0,-1,0,0,1,-1,0,0,0,1,0,0,0,-1,0,-1,0,-1,0,0,0,-1,0,0,0,1,1,0,0,-1,0,0,0,0,-1,0,-1,0,0,0,-1,1,0,0,0,1,0,0,-1,0,-1,0,0,0,0,1,1,0,0,0,-1,0,0,0,1,0,-1,0,-1,0,0,0,0,-1,-1,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,1]; % o_h
                end
                R = reshape(R,3,3,[]);
                
                % double-valued representation
                switch pg_code
                case 1; W = [1,0,0,1,-1,0,0,-1];
                case 2; W = [1,0,0,1,1i,0,0,1i,-1,0,0,-1,-1i,0,0,-1i];
                case 3; W = [1,0,0,1,0,-1,1,0,-1,0,0,-1,0,1,-1,0];
                case 4; W = [1,0,0,1,0,-1i,1i,0,-1,0,0,-1,0,1i,-1i,0];
                case 5; W = [1,0,0,1,0,-1,1,0,1i,0,0,1i,0,-1i,1i,0,-1,0,0,-1,0,1,-1,0,-1i,0,0,-1i,0,1i,-1i,0];
                case 6; W = [1,0,0,1,1i,0,0,-1i,0,-1,1,0,0,1i,1i,0,-1,0,0,-1,-1i,0,0,1i,0,1,-1,0,0,-1i,-1i,0];
                case 7; W = [1,0,0,1,1i,0,0,-1i,0,-1i,1i,0,0,-1,-1,0,-1,0,0,-1,-1i,0,0,1i,0,1i,-1i,0,0,1,1,0];
                case 8; W = [1,0,0,1,1i,0,0,-1i,0,-1,1,0,1i,0,0,1i,0,1i,1i,0,-1,0,0,1,0,-1i,1i,0,0,-1,-1,0,-1,0,0,-1,-1i,0,0,1i,0,1,-1,0,-1i,0,0,-1i,0,-1i,-1i,0,1,0,0,-1,0,1i,-1i,0,0,1,1,0];
                case 9; W = [1,0,0,1,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-1,0,0,-1,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2];
                case 10; W = [1,0,0,1,1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,-1,0,0,-1,-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2];
                case 11; W = [1,0,0,1,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-1,0,0,-1,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2];
                case 12; W = [1,0,0,1,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-1,0,0,-1,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2];
                case 13; W = [1,0,0,1,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,-6^(1/2)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,6^(1/2)/6,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*3^(1/2))/3,6^(1/2)*(- 1/6 - 1i/6),6^(1/2)*(- 1/6 + 1i/6),-(2^(1/2)*3^(1/2))/3,6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,-6^(1/2)/6,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,-1,0,0,-1,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,6^(1/2)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*3^(1/2))/3,6^(1/2)*(1/6 + 1i/6),6^(1/2)*(1/6 - 1i/6),(2^(1/2)*3^(1/2))/3,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2];
                case 14; W = [1,0,0,1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-1,0,0,-1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2)];
                case 15; W = [1,0,0,1,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),1i,0,0,-1i,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),-1,0,0,-1,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-1i,0,0,1i,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2)];
                case 16; W = [1,0,0,1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,1i,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-1,0,0,1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),-1,0,0,-1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,-1i,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),1,0,0,-1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2)];
                case 17; W = [1,0,0,1,0,-1,1,0,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,1i,1i,0,-1,0,0,-1,0,1,-1,0,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,-1i,-1i,0];
                case 18; W = [1,0,0,1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,-1i,1i,0,1i,0,0,-1i,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,-1,-1,0,-1,0,0,-1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,1i,-1i,0,-1i,0,0,1i,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,1,1,0];
                case 19; W = [1,0,0,1,0,-1,1,0,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,1i,0,0,-1i,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,1i,1i,0,-1,0,0,-1,0,1,-1,0,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,-1i,0,0,1i,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,-1i,-1i,0];
                case 20; W = [1,0,0,1,0,-1,1,0,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,1i,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,-1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),0,1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-1,0,0,1,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),0,-1,-1,0,-1,0,0,-1,0,1,-1,0,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,-1i,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),0,-1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,1,0,0,-1,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),0,1,1,0];
                case 21; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6)];
                case 22; W = [1,0,0,1,3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 23; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6),- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 24; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6)];
                case 25; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,-6^(1/2)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,6^(1/2)/6,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*3^(1/2))/3,6^(1/2)*(- 1/6 - 1i/6),6^(1/2)*(- 1/6 + 1i/6),-(2^(1/2)*3^(1/2))/3,6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,-6^(1/2)/6,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,6^(1/2)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,-6^(1/2)/6,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*3^(1/2))/3,6^(1/2)*(1/6 + 1i/6),6^(1/2)*(1/6 - 1i/6),(2^(1/2)*3^(1/2))/3,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6)];
                case 26; W = [1,0,0,1,3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 27; W = [1,0,0,1,-(3^(1/2)*1i)/3,3^(1/2)*(1/3 - 1i/3),3^(1/2)*(- 1/3 - 1i/3),(3^(1/2)*1i)/3,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1i,0,0,1i,(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,-(6^(1/2)*1i)/6,3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(1/2 + 1i/6),3^(1/2)/3,3^(1/2)*(1/3 + 1i/3),3^(1/2)*(1/3 - 1i/3),-3^(1/2)/3,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(- 1/6 + 1i/6),6^(1/2)*(1/6 + 1i/6),(2^(1/2)*3^(1/2)*1i)/3,-6^(1/2)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,6^(1/2)/6,-(6^(1/2)*1i)/6,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,(6^(1/2)*1i)/6,3^(1/2)*(1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/2 - 1i/6),3^(1/2)*(1/6 + 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/2),0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,(2^(1/2)*3^(1/2))/3,6^(1/2)*(- 1/6 - 1i/6),6^(1/2)*(- 1/6 + 1i/6),-(2^(1/2)*3^(1/2))/3,6^(1/2)/6,6^(1/2)/6 - (2^(1/2)*3^(1/2)*1i)/3,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,-6^(1/2)/6,3^(1/2)*(- 1/6 + 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/2),-1,0,0,-1,(3^(1/2)*1i)/3,3^(1/2)*(- 1/3 + 1i/3),3^(1/2)*(1/3 + 1i/3),-(3^(1/2)*1i)/3,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-1i,0,0,-1i,-(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2)*1i)/3 + 6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(- 1/2 - 1i/6),-3^(1/2)/3,3^(1/2)*(- 1/3 - 1i/3),3^(1/2)*(- 1/3 + 1i/3),3^(1/2)/3,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)*(1/6 - 1i/6),6^(1/2)*(- 1/6 - 1i/6),-(2^(1/2)*3^(1/2)*1i)/3,6^(1/2)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,- (2^(1/2)*3^(1/2))/3 - (6^(1/2)*1i)/6,-6^(1/2)/6,(6^(1/2)*1i)/6,(2^(1/2)*3^(1/2))/3 + (6^(1/2)*1i)/6,(6^(1/2)*1i)/6 - (2^(1/2)*3^(1/2))/3,-(6^(1/2)*1i)/6,3^(1/2)*(- 1/2 - 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/2 + 1i/6),3^(1/2)*(- 1/6 - 1i/2),3^(1/2)*(- 1/6 - 1i/6),3^(1/2)*(- 1/6 + 1i/6),3^(1/2)*(1/6 - 1i/2),0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-(2^(1/2)*3^(1/2))/3,6^(1/2)*(1/6 + 1i/6),6^(1/2)*(1/6 - 1i/6),(2^(1/2)*3^(1/2))/3,-6^(1/2)/6,(2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,- (2^(1/2)*3^(1/2)*1i)/3 - 6^(1/2)/6,6^(1/2)/6,3^(1/2)*(1/6 - 1i/2),3^(1/2)*(1/6 + 1i/6),3^(1/2)*(1/6 - 1i/6),3^(1/2)*(- 1/6 - 1i/2)];
                case 28; W = [1,0,0,1,0,-1,1,0,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1i,0,0,-1i,0,1i,1i,0,-1,0,0,-1,0,1,-1,0,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,-1i,0,0,1i,0,-1i,-1i,0];
                case 29; W = [1,0,0,1,0,-1,1,0,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1i,0,0,1i,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,0,-1i,1i,0,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1i,0,0,-1i,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,0,1i,1i,0,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,-1,0,0,1,0,-1,-1,0,-1,0,0,-1,0,1,-1,0,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,-1i,0,0,-1i,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,0,1i,-1i,0,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,-1i,0,0,1i,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,0,-1i,-1i,0,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1,0,0,-1,0,1,1,0];
                case 30; W = [1,0,0,1,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,1i,0,0,-1i,2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,0,-1,1,0,0,1i,1i,0,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,-1,0,0,-1,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,-1i,0,0,1i,-2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,-2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,0,1,-1,0,0,-1i,-1i,0,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0];
                case 31; W = [1,0,0,1,1i,0,0,-1i,0,-1,1,0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,0,1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,-1,0,0,-1,-1i,0,0,1i,0,1,-1,0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,0,-1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2];
                case 32; W = [1,0,0,1,1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1i,0,0,1i,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,1i,0,0,-1i,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,2^(1/2)*(1/2 + 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,2^(1/2)/2,1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 - 1i/2,-1,0,0,1,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(1/2 + 1i/2),1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 + 1i/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(1/2 + 1i/2),0,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,0,-1,1,0,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,0,1i,1i,0,- 1/2 + 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,1/2 + 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,-1i,1i,0,0,2^(1/2)*(- 1/2 - 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,-1,-1,0,-2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,2^(1/2)/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(1/2 + 1i/2),0,-1,0,0,-1,- 1/2 + 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 + 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),-1i,0,0,-1i,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,-(2^(1/2)*1i)/2,2^(1/2)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,-1i,0,0,1i,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(1/2 - 1i/2),-2^(1/2)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,2^(1/2)*(- 1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 + 1i/2),- 1/2 - 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,2^(1/2)/2,2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,- 1/2 - 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 + 1i/2,1,0,0,-1,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-2^(1/2)/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)*(1/2 - 1i/2),0,0,2^(1/2)*(- 1/2 - 1i/2),- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,- 1/2 + 1i/2,1/2 + 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,- 1/2 - 1i/2,0,2^(1/2)*(1/2 - 1i/2),2^(1/2)*(- 1/2 - 1i/2),0,-(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-(2^(1/2)*1i)/2,0,1,-1,0,- 1/2 - 1i/2,1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,-2^(1/2)/2,-(2^(1/2)*1i)/2,0,-1i,-1i,0,1/2 - 1i/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,1/2 + 1i/2,1/2 - 1i/2,1/2 - 1i/2,-(2^(1/2)*1i)/2,-2^(1/2)/2,2^(1/2)/2,(2^(1/2)*1i)/2,- 1/2 - 1i/2,- 1/2 + 1i/2,- 1/2 - 1i/2,1/2 - 1i/2,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(1/2 - 1i/2),0,0,1i,-1i,0,0,2^(1/2)*(1/2 + 1i/2),2^(1/2)*(- 1/2 + 1i/2),0,0,1,1,0,2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,-2^(1/2)/2,2^(1/2)/2,-(2^(1/2)*1i)/2,(2^(1/2)*1i)/2,-2^(1/2)/2,0,2^(1/2)*(- 1/2 + 1i/2),2^(1/2)*(- 1/2 - 1i/2),0];
                end
                W = reshape(W,2,2,[]);
                
            else
                % load recipe
                recipe = get_recipe(pg_code);
                % initialize
                R = zeros(3,3,numel(recipe)+1);
                % make generators
                nsyms = 0;
                nsyms = nsyms+1; R(:,:,nsyms) = eye(3);
                for k = 1:numel(recipe)
                    nsyms=nsyms+1; R(:,:,nsyms) = decode_R(recipe(k)); 
                end
                
                % add elements until no new elements are generated
                R = am_cell.complete_group(R);
                
                % [IMPORTANT] use rhombohedral rather than hexagonal setting for hexagonal point groups
                if any(pg_code==[9:13,21:27])
                    T = [-1 1 1; 2 1 1; -1 -2 1].'/3;
                    R = matmul_(matmul_(inv(T),R),(T));
                end
                
                % generate double group
                if nargout > 1
                    % get double group
                    j=1/2; [W] = get_wigner(j,R,'spherical');                     
                    % Add plus and minus in accordance with 
                    % V. Heine, Group Theory in Quantum Mechanics (Elsevier, 2014), p 62, eq 8.24.
                    W = cat(3,W,-W); % W = am_cell.complete_group(reshape(uniquec_(reshape(wdv_(W),4,[])),2,2,[]));
                    % remove numerical noise
                    W = wdv_(W);
                end
            end

            function recipe  = get_recipe(sg_id)
                pg_recipe_database = {...
                     '','h','c','j','ch','bc','bj','bch','n','hn','en','kn','fhn','g','m','gh','cg',...
                    'gj','cm','cgh','bn','in','bhn','ben','bkn','ikn','benh','cd','cdh','dg','bcld','dgh'};
                recipe = pg_recipe_database{sg_id};
            end
            function R       = decode_R(str_r)
                switch str_r
                case 'a'; R = [  1  0  0;  0  1  0;  0  0  1]; % a
                case 'b'; R = [ -1  0  0;  0 -1  0;  0  0  1]; % b
                case 'c'; R = [ -1  0  0;  0  1  0;  0  0 -1]; % c
                case 'd'; R = [  0  0  1;  1  0  0;  0  1  0]; % d
                case 'e'; R = [  0  1  0;  1  0  0;  0  0 -1]; % e
                case 'f'; R = [  0 -1  0; -1  0  0;  0  0 -1]; % f
                case 'g'; R = [  0 -1  0;  1  0  0;  0  0  1]; % g
                case 'h'; R = [ -1  0  0;  0 -1  0;  0  0 -1]; % h
                case 'i'; R = [  1  0  0;  0  1  0;  0  0 -1]; % i
                case 'j'; R = [  1  0  0;  0 -1  0;  0  0  1]; % j
                case 'k'; R = [  0 -1  0; -1  0  0;  0  0  1]; % k
                case 'l'; R = [  0  1  0;  1  0  0;  0  0  1]; % l
                case 'm'; R = [  0  1  0; -1  0  0;  0  0 -1]; % m
                case 'n'; R = [  0 -1  0;  1 -1  0;  0  0  1]; % n
                end
            end 
        end

        function [S,MT]       = complete_group(S,tol)

            import am_lib.*
            
            % set tolernece
            if nargin<2; tol=am_dft.tiny; end

            % get size
            s=size(S); d=s(1)*s(2);
            
            % switch
            if s(1)==4 && s(2)==4 && all_(eq_(S(4,1:4,:), [0,0,0,1], tol))
                algo = 2; % seitz symmetry
            else
                algo = 1; % point symmetry
            end

            % exclude empty symmetries
            S = S(:,:,any(reshape(S,d,[])~=0,1)); 
            
            % add identity if not present 
            E = eye(s(1),s(2)); if ~any(all(all(S==E))); S = cat(3,E,S); end
            
            % get number of symmetries at this point
            nsyms=size(S,3);

            % allocate space
            S=cat(3,S,zeros(size(S,1),size(S,1),2*192-nsyms));

            % add symmetry until no new symmetry is generated
            nsyms_last=0; MT = zeros(2*192,2*192); 
            while nsyms ~= nsyms_last
                nsyms_last=nsyms;
                for i = 1:nsyms_last
                for j = 1:nsyms_last
                    if MT(i,j)==0
                        A = symmul_(S(:,:,i),S(:,:,j),algo,tol);
                        A_id = member_(A(:),reshape(S(:,:,1:nsyms),d,[]),tol);
                        if A_id == 0
                            nsyms = nsyms+1; S(:,:,nsyms) = A; MT(i,j) = nsyms;
                        else
                            MT(i,j) = A_id;
                        end
                    end
                end
                end
            end
            
            % trim output
            S  = S(:,:,1:nsyms);
            MT = MT(1:nsyms,1:nsyms);
            
            function C = symmul_(A,B,algo,tol)
                switch algo
                    case 1
                        % point symmetry
                        C = A*B;
                    case 2
                        % seitz symmetry
                        C = A*B; C(1:3,4)=mod(C(1:3,4)+tol,1)-tol;
                end
            end
        end
             
        function [MT]         = get_multiplication_table(S,tol)
            % [MT,E,I] = get_multiplication_table(S)
            % get multiplication table: S(:,:,i)*S(:,:,j) = S(:,:,MT(i,j))

            import am_lib.* am_dft.*

            if nargin<2; tol=am_dft.tiny; end

            s = size(S); if numel(s)<3; s(3) = 1; end

            if iscell(S)  % seitz operator combined with permutation (represented as a two-part cell)
                s = size(S{1}); if numel(s)<3; s(3) = 1; end
                if     s(1)==4 && s(2)==4 && all_(eq_(S{1}(4,1:4,:), [0,0,0,1], tol)) % seitz operator (applies mod to translational components)
                    md_ = @(X) [X(1:12,:);mod_(X(13:15,:),tol);X(16:end,:)];
                    rs_ = @(X) md_(reshape(X,s(1)*s(2),[])); nsyms = s(3);

                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref, tol), s(3), s(3)); 
                elseif s(1)==3 && s(2)==3                                          % point operator
                    rs_ = @(X) reshape(X,s(1)*s(2),[]); nsyms = s(3);

                    ref = [rs_(S{1});reshape(S{2},size(S{2},1),[])];
                    opr = [rs_(matmulp_(S{1},permute(S{1},[1,2,4,3]))); reshape(operm_(S{2},S{2}),size(S{2},1),[])];

                    MT = reshape( member_( opr , ref, tol), s(3), s(3));
                end
            elseif s(1)==4 && s(2)==4 && all_(eq_(S(4,1:4,:), [0,0,0,1], tol))     % seitz operator (applies mod to translational components)
                md_ = @(X) [X(1:12,:);mod_(X(13:15,:),tol);X(16:end,:)];
                rs_ = @(X) md_(reshape(X,s(1)*s(2),[])); nsyms = s(3);

                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))), rs_(S), tol) , nsyms, nsyms);
            else                                                                   % point operator
                rs_ = @(X) reshape(X,s(1)*s(2),[]); nsyms = s(3); 

                MT = reshape( member_( rs_(matmulp_(S,permute(S,[1,2,4,3]))), rs_(S), tol ) , nsyms, nsyms);
            end

            if any(sum(MT,1)~=sum([1:nsyms])) || any(sum(MT,2)~=sum([1:nsyms]))
                error('MT is incorrect. Check for mistakes in the symmetry and ensure that the symmetry is in the primitive basis');
            end
        end

    end

end

