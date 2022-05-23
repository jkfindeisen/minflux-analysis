classdef Configuration < specmx.internal.Configuration
    
    properties ( Dependent )
        
        % read only properties
        Name
        NumberStacks
    end
    
    methods
        function obj = Configuration( varargin )
            obj = obj@specmx.internal.Configuration( varargin{:} );
        end
        
        % read only properties
        
        function name = get.Name( obj )
            name = obj.name();
        end
        
        function number_stacks = get.NumberStacks( obj )
            number_stacks = obj.number_of_stacks();
        end
    end
    
end