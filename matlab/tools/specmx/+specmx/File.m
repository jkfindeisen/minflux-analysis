classdef File < specmx.internal.File
    
    properties ( Dependent )
        
        % read and writeable properties
        Description
        
        % read only properties
        NumberStacks
    end
    
    methods
        function obj = File( varargin )
            
            if ~(nargin == 1 && isa(varargin{1}, 'SwigRef'))
                % any default parameters handling can come here
                
                % default access mode = read
                if nargin == 1
 	                varargin = [varargin, specmx.File.Read];
                end
            end
            obj = obj@specmx.internal.File( varargin{:} );
        end
        
        function write( obj, stacks, compression )
            % Enhanced write functions takes single stacks or cell arrays
            % of stacks and stacks can be Matlab arrays or specmx.Stacks
            assert( nargin >= 2, 'Not enough parameters' );
            if nargin == 2
                compression = true;
            end
            
            % wrap in cell array
            if ~iscell( stacks )
                stacks = { stacks };
            end
            
            % iterate over cell array
            for i = 1 : numel( stacks )
                stack = stacks{ i };
                
                % wrap in stack
                if ~isa( stack, 'specmx.Stack' )
                    % convert to uint8 if logical
                    if isa( stack, 'logical' )
                        stack = uint8( stack );
                    end
                    stack = specmx.Stack( stack );
                end
                
                % write using super class
                write@specmx.internal.File( obj, stack, compression );
            end
        end
        
        % read and writeable properties
        
        function description = get.Description( obj )
            description = obj.description();
        end
        
        function obj = set.Description( obj, value )
            obj.set_description( value );
        end
        
        % read only properties
        
        function number_stacks = get.NumberStacks( obj )
            number_stacks = obj.number_of_stacks();
        end
        
    end
    
    methods ( Static )
        
        function stacks = read_file( path )
            % Convenience function, opens a file for reading, reads all
            % stacks in the file and closes the file again.
            
            % open file
            f = specmx.File( path, specmx.File.Read );
            n = f.NumberStacks;
            
            % read all stacks in file into cell array
            stacks = cell(n, 1);
            for i = 1 : n
                stacks{ i } = f.read( i - 1 );
            end
        end
        
        function write_file( path, stacks, compression )
            % Convenience function, opens a file for writing, writes some
            % stacks (overwrites existing file and stacks) and closes the
            % file again.
            assert( nargin >= 2, 'Not enough parameters' );
            if nargin == 2
                compression = true;
            end
            
            % open file
            f = specmx.File( path, specmx.File.Write );
            
            % write
            f.write( stacks, compression );
        end
    end
end