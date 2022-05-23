classdef StackReference < SwigRef
    %Usage: StackReference ()
    %
  methods
    function this = swig_this(self)
      this = specmx(3, self);
    end
    function delete(self)
      if self.swigPtr
        specmx(5, self);
        self.swigPtr=[];
      end
    end
    function varargout = parameters(self,varargin)
    %Usage: retval = parameters (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type Json::Value. 
      [varargout{1:nargout}] = specmx(6, self, varargin{:});
    end
    function varargout = parameters_(self,varargin)
    %Usage: retval = parameters_ (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type OProp. 
      [varargout{1:nargout}] = specmx(7, self, varargin{:});
    end
    function self = StackReference(varargin)
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
