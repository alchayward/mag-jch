function new_configs = Get_Config_Data(configs)
%Take in a set of Pfaffian configs and sort them into sets of different
%lattice configurations. Create the lattice object then send off the
%calculation.
pf = Pfaffians();


%sort into lattice configs.
%First number of particles.
lattice_configs = unique_lattice_configs(configs);
system_cache_file = 'gjch_pfaffian_system_cache.mat';
field_compare_list = {'lattice_dims','nParticles'};
get_from_cachef = @(c)get_from_cache(c,system_cache_file,field_compare_list);
lattice_configs = cellfun(get_from_cachef,lattice_configs,'UniformOutput',false);
new_configs={};
for ii=1:length(lattice_configs)
    lc = lattice_configs{ii};
    if ~lc.done
    lc.s = pf.new_jch_system(lc.lattice_dims,lc.nParticles,...
       0,sqrt(2)); %Make our lattice system
    lc.done = true;
    put_in_cache(system_cache_file,lc,field_compare_list);
    end
    
    
    new2_configs = pfaffian_overlap_and_energy(lc.s,system_filter(lc,configs));
    for jj =1:length(new2_configs)
        new_configs{end+1} = new2_configs{jj}; %#ok<AGROW>
    end
end


end

function lattice_configs = unique_lattice_configs(configs)
    lattice_configs = {};
    for ii = 1:length(configs)
        c = configs{ii};
        new = true;
        for jj = 1:length(lattice_configs)
            l = lattice_configs{jj};
            if (l.nParticles == c.nParticles && all(l.lattice_dims == c.lattice_dims))
                new=false;
            end
        end
        if new
            lattice_configs{end+1} = struct('nParticles',c.nParticles,...
               'lattice_dims',c.lattice_dims,'done',false); %#ok<AGROW>
        end
    end
end

function new_configs = system_filter(lc,configs)
new_configs = {};
    for ii = 1:length(configs)
        c = configs{ii};
        if (lc.nParticles == c.nParticles && all(lc.lattice_dims == c.lattice_dims))
            new_configs{end+1} = c; %#ok<AGROW>
        end
    end
end



function new_configs = pfaffian_overlap_and_energy(s,configs,varargin)

if nargin == 3
    opts = varargin{1};
else
    opts = struct();
end

if ~isfield(opts,'cache')
    opts.cache = true;
end
cache = opts.cache;

if ~isfield(opts,'cache_file')
  opts.cache_file = './pfaff_g_and_ol_save.mat';
  
end

cache_file = opts.cache_file;
%%% Create cache file if it doesnt' exist %%%


%libraries
hub = HubbardLibrary();
pf = Pfaffians();
wf = Wavefunctions();

field_compare_list = {'lattice_dims','nParticles',...
    'beta','delta','beta_2','kappa'};

if cache
    get_from_cachef = @(c)get_from_cache(c,cache_file,field_compare_list);
    new_configs = cellfun(get_from_cachef,configs,'UniformOutput',false);
else
    new_configs = configs;
end

   nLevels=6;
   GSopts = struct('nLevels',nLevels,'V0',rand(s.hilb_dim,nLevels),...
       'method','eigs');
   jchpfaff = s.PfaffianState(struct('delta',new_configs{1}.delta));
   GSopts.V0(:,1:3) = jchpfaff; %Use pfaffian as inital guess.
   for ii = 1:length(new_configs)
       c = new_configs{ii};
       if ~isfield(c,'done')
           c.done = false;
       end
       
       if ~c.done
           disp('finding groundstate of:')
           disp(c);
           h = pf.make_jch_ham_params(s,c.kappa,c.delta,c.beta_2);
           jchpfaff = s.PfaffianState(c);
            
           [V,D,~] = hub.Groundstate(h,GSopts); %find_wavefunction(parameters)
           GSopts.V0 = V; %Use solution as next inital guess.
           %This step takes the longest it seems. Do need to construct a new
           %state for each 

           ol = wf.wavefunction_overlap(V, jchpfaff); %overlap matrix
           c.D = D;
           c.ol = ol;
           c.done = true;
           new_configs{ii}=c;
           %save these results each iteration, cause it would be a shame to
           %lose them.
           put_in_cache(cache_file,c,field_compare_list);
       end
   end
   
end

function new_c = get_from_cache(c,cache_file,field_list)
cache = load_cache(cache_file);
ind = in_cache(c,cache,field_list);
if ind ~= 0
    new_c = cache{ind};
else
    new_c = c;
end
end

function ind = in_cache(c,cache,field_list)

ind = 0;
for ii=1:length(cache)
    if all(cellfun(@(field)(all(cache{ii}.(field)==c.(field))),...
            field_list)) && cache{ii}.done == true
        ind = ii;
    end
end

end

function put_in_cache(cache_file,c,field_list)
cache = load_cache(cache_file);

if ~in_cache(c,cache,field_list)
    cache{end+1} = c; %#ok<NASGU>
end
save(cache_file,'cache');

end

function cache = load_cache(cache_file)
    if ~exist(cache_file,'file')
      cache = {};
      save(cache_file,'cache');
    else
    f = load(cache_file);
    cache=f.cache;
    end
end
