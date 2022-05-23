classdef Measurement < specmx.internal.Measurement
    
    properties ( Dependent )
        
        % read only properties
        Name
        NumberConfigurations
        NumberStacks
    end
    
    methods
        function obj = Measurement( varargin )
            obj = obj@specmx.internal.Measurement( varargin{:} );
        end
        
        % read only properties
        
        function name = get.Name( obj )
            name = obj.name();
        end
        
        function number_configurations = get.NumberConfigurations( obj )
            number_configurations = obj.number_of_configurations();
        end
        
        function number_stacks = get.NumberStacks( obj )
            number_stacks = obj.number_of_stacks();
        end
    end
    
end