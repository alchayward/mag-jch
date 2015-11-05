classdef Cache

properties
    file
    field_list
end %properties

methods
    function c = Cache(file,field_list)
       c.file = file;
       c.field_list=field_list;
    end
    
function items = get_from_cache(c,items)
%want to make it so we can pass a cell array. make code so it can take
%a cell array, then jsut wrap a single struct in a cell.
struct_flag=0;
if isstruct(items)
    struct_flag=1;
    items = {items};
end

cache = c.load_cache(c.file);
for ii = 1:length(items)
    item = items{ii};
    ind = c.in_cache(item,cache,c.field_list);
    if ind ~= 0
        items{ii} = cache{ind};
    end
end
%if it was a struct, then jsut return the struct:
if struct_flag
    items = items{1};
end
end


function put_in_cache(c,items)
    %%% pul
if isstruct(items)
    items = {items};
end
cache = c.load_cache(c.file);
for ii = 1:length(items)
item = items{ii};

ind = c.in_cache(item,cache,c.field_list);
if ind == 0
    cache{end+1} = item; %#ok<AGROW,NASGU>
else
    cache{ind} = item; %#ok<NASGU> REPLACE 
end
end
save(c.file,'cache');
end %put_in_cache
end %methods

methods (Static = true)

function ind = in_cache(item,cache,field_list)
ind = 0;
comp_field = @(field,x)(isequal(cache{x}.(field),item.(field)));
comp_all_fields = @(x)all(cellfun(@(f)comp_field(f,x),...
            field_list));
for ii=1:length(cache)
    if comp_all_fields(ii)
        ind = ii;
        break; %Quit if we find the value
    end %if
end %for
end %in_cache

function cache = load_cache(cache_file)
    if ~exist(cache_file,'file')
      cache = {};
      save(cache_file,'cache');
    else
    f = load(cache_file);
    cache=f.cache;
    end
end %load_cache
end %static methods
end %classdef
