classdef File < SwigRef
    %Usage: File ()
    %
  methods
    function this = swig_this(self)
      this = specmx(3, self);
    end
    function self = File(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = specmx(38, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function delete(self)
      if self.swigPtr
        specmx(39, self);
        self.swigPtr=[];
      end
    end
    function varargout = description(self,varargin)
    %Usage: retval = description ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(40, self, varargin{:});
    end
    function varargout = set_description(self,varargin)
    %Usage: set_description (to)
    %
    %to is of type std::string const &. 
      [varargout{1:nargout}] = specmx(41, self, varargin{:});
    end
    function varargout = number_of_stacks(self,varargin)
    %Usage: retval = number_of_stacks ()
    %
    %retval is of type std::size_t. 
      [varargout{1:nargout}] = specmx(42, self, varargin{:});
    end
    function varargout = read(self,varargin)
    %Usage: retval = read (position)
    %
    %position is of type std::size_t. position is of type std::size_t. retval is of type Remote::Stack::Shared. 
      [varargout{1:nargout}] = specmx(43, self, varargin{:});
    end
    function varargout = write(self,varargin)
    %Usage: write (stack)
    %
    %stack is of type Remote::Stack::Shared. 
      [varargout{1:nargout}] = specmx(44, self, varargin{:});
    end
  end
  methods(Static)
    function v = Read()
      persistent vInitialized;
      if isempty(vInitialized)
        vInitialized = specmx(0, 1);
      end
      v = vInitialized;
    end
    function v = Write()
      persistent vInitialized;
      if isempty(vInitialized)
        vInitialized = specmx(0, 2);
      end
      v = vInitialized;
    end
    function v = Append()
      persistent vInitialized;
      if isempty(vInitialized)
        vInitialized = specmx(0, 3);
      end
      v = vInitialized;
    end
  end
end
