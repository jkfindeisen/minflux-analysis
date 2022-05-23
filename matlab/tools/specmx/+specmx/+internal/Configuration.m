classdef Configuration < SwigRef
    %Usage: Configuration ()
    %
  methods
    function this = swig_this(self)
      this = specmx(3, self);
    end
    function delete(self)
      if self.swigPtr
        specmx(45, self);
        self.swigPtr=[];
      end
    end
    function varargout = name(self,varargin)
    %Usage: retval = name ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(46, self, varargin{:});
    end
    function varargout = parameters(self,varargin)
    %Usage: retval = parameters (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type Json::Value. 
      [varargout{1:nargout}] = specmx(47, self, varargin{:});
    end
    function varargout = set_parameters(self,varargin)
    %Usage: set_parameters (at, from)
    %
    %at is of type std::string const &. from is of type Json::Value const &. 
      [varargout{1:nargout}] = specmx(48, self, varargin{:});
    end
    function varargout = parameters_(self,varargin)
    %Usage: retval = parameters_ (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type OProp. 
      [varargout{1:nargout}] = specmx(49, self, varargin{:});
    end
    function varargout = set_parameters_(self,varargin)
    %Usage: set_parameters_ (at, from)
    %
    %at is of type std::string const &. from is of type OProp const &. 
      [varargout{1:nargout}] = specmx(50, self, varargin{:});
    end
    function varargout = number_of_stacks(self,varargin)
    %Usage: retval = number_of_stacks ()
    %
    %retval is of type std::size_t. 
      [varargout{1:nargout}] = specmx(51, self, varargin{:});
    end
    function varargout = stack_names(self,varargin)
    %Usage: retval = stack_names ()
    %
    %retval is of type std::vector< std::string >. 
      [varargout{1:nargout}] = specmx(52, self, varargin{:});
    end
    function varargout = stack(self,varargin)
    %Usage: retval = stack (name)
    %
    %name is of type std::string const &. name is of type std::string const &. retval is of type Remote::Stack::Shared. 
      [varargout{1:nargout}] = specmx(53, self, varargin{:});
    end
    function self = Configuration(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        error('No matching constructor');
      end
    end
  end
  methods(Static)
  end
end
