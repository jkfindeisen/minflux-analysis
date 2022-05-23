classdef Stack < specmx.internal.Stack
    
    properties ( Dependent )
        
        % read and writeable properties
        Name
        Description
        Labels
        Lengths
        Offsets
        PixelLengths
        
        % read only properties
        Type
        Size
        NumberDimensions
    end
    
    methods
        function obj = Stack( varargin )
            obj = obj@specmx.internal.Stack( varargin{:} );
        end
        
        % read and writeable properties
        
        function name = get.Name( obj )
            name = obj.name();
        end
        
        function obj = set.Name( obj, value )
            obj.set_name( value );
        end
        
        function description = get.Description( obj )
            description = obj.description;
        end
        
        function obj = set.Description( obj, value )
            obj.set_description( value );
        end
        
        function labels = get.Labels( obj )
            % output is row vector
            labels = obj.labels().';
        end
        
        function obj = set.Labels( obj, value )
            % are converted to column vector
            obj.set_labels( value(:) );
        end
        
        function lengths = get.Lengths( obj )
            % output is row vector
            lengths = obj.lengths().';
        end
        
        function obj = set.Lengths( obj, value )
            % are converted to column vector
            obj.set_lengths( value(:) );
        end
        
        function offsets = get.Offsets( obj )
            % output is row vector
            offsets = obj.offsets().';
        end
        
        function obj = set.Offsets( obj, value )
            % are converted to column vector
            obj.set_offsets( value(:) );
        end
        
        function pixel_lengths = get.PixelLengths( obj )
            % this uses Lengths and Sizes and divides both
            pixel_lengths = obj.Lengths ./ obj.Size;
        end
        
        function obj = set.PixelLengths( obj, value )
            % use set Lengths internally, no need to use both, just one or
            % the other
            obj.Lengths = value .* obj.Size;
        end
        
        % read only properties
        
        function type = get.Type( obj )
            type = obj.type();
        end
        
        function size = get.Size( obj )
            % are converted to row double vector
            % internally always given as column
            size = double( obj.sizes() ).';
        end
        
        function number_dimensions = get.NumberDimensions( obj )
            % without trailing zero dimensions
            size = double( obj.sizes() ).';
            
            number_dimensions = numel(size);                        
            while ( size(number_dimensions) <= 1  && number_dimensions > 1 )
                number_dimensions = number_dimensions - 1;
            end
        end
        
    end
    
    methods ( Static )
        
        function new_stack = reorder( stack, order )
            % given a stack and a new dimension order, returns a reordered
            % stack, the new stack is a copy
            %
            % similar to matlab function permute, but also permutes
            % meta-data
            
            assert( nargin == 2 );
            
            data = stack.data();
            dims = size( data );
            
            assert( length( order ) == length( dims ) );
            
            % create new stack
            new_data = permute( data, order);
            new_stack = specmx.Stack( new_data );
            
            % copy name, description
            new_stack.Name = [ stack.Name, '.REORDERED' ];
            new_stack.Description = stack.Description;            
            
            % copy and reorder labels, lengths, offsets
            labels = stack.Labels;
            new_stack.Labels = labels( order );
            
            lengths = stack.Lengths;
            new_stack.Lengths = lengths( order );
            
            offsets = stack.Offsets;
            new_stack.Offsets = offsets( order );
            
        end
        
        function new_stack = resize( stack, left_top, right_bottom )
            % given a left, top corner and a right, bottom corner the sub
            % data stack is copied and a new stack is created, name,
            % description are copied, the length is adjusted, the corners
            % are included.
            %
            % current limitations
            % - only cut out, no extend
            % - ignores offset
            % - trailing singleton dimensions are cut off (general Matlab issue)
            
            % sanity check
            assert( nargin == 3 );
            assert( all( right_bottom >= left_top ) );
            
            data = stack.data();
            dims = size( data );
            
            % currently you can only make it smaller, not larger
            assert( all( left_top >= 1 ) );
            assert( all( right_bottom <= dims ) );
            
            % cut out
            switch length(dims)
                case 1
                    data = data( left_top(1) : right_bottom(1) );
                case 2
                    data = data( left_top(1) : right_bottom(1), left_top(2) : right_bottom(2) );
                case 3
                    data = data( left_top(1) : right_bottom(1), left_top(2) : right_bottom(2), left_top(3) : right_bottom(3) );
                case 4
                    data = data( left_top(1) : right_bottom(1), left_top(2) : right_bottom(2), left_top(3) : right_bottom(3), left_top(4) : right_bottom(4) );
                otherwise
                    error('Unsupported number of dimension!');
            end
            dims = size( data );
            
            % create new stack
            new_stack = specmx.Stack( data );
            
            % copy name, description
            new_stack.Name = [ stack.Name, '.RESIZE' ];
            new_stack.Description = stack.Description;
            new_stack.Labels = stack.Labels;
            
            % adjust length, ignore offset
            new_stack.Lengths = dims .* stack.PixelLengths( 1 : length(dims) );
            
        end
        
        function stack = findByName( stacks, name )
            % given a cell array of stacks, finds the one whose name
            % matches name and returns it, will throw an error if no stack
            % with this name exists
            
            assert( nargin == 2 );
            
            index = find( cellfun( @(x) strcmp(x.Name, name), stacks ) );
            
            assert( length( index ) == 1 );
            
            stack = stacks{ index };
            
        end
        
        function cloned_stack = clone( stack, type )
            % creates a new (empty) stack of the same type and size and
            % copies the properties Name, Description, Labels, Lengths, Offsets
            % optionally a different type can be given
            
            assert( nargin >= 1 );
            if nargin == 1 || isempty(type)
                type = stack.type();
            end
            
            cloned_stack = specmx.Stack( type, stack.sizes() );
            
            % copies properties
            properties_to_copy = {'Name', 'Description', 'Labels', 'Lengths', 'Offsets'};
            for i = 1 : length(properties_to_copy)
                property = properties_to_copy{i};
                cloned_stack.(property) = stack.(property);
            end
            
        end
        
    end
    
end