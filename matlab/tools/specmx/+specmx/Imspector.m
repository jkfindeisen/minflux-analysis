classdef Imspector < specmx.internal.Imspector
    
    properties ( Dependent )
        
        % read only properties
        Version
        Host
    end
    
    methods
        function obj = Imspector( varargin )
            obj = obj@specmx.internal.Imspector( varargin{:} );
        end
        
        % read only properties
        
        function version = get.Version( obj )
            version = obj.version();
        end
        
        function host = get.Host( obj )
            host = obj.host();
        end
    end
    
end