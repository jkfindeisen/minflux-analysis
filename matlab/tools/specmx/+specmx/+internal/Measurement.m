classdef Measurement < SwigRef
    %Usage: Measurement ()
    %
  methods
    function this = swig_this(self)
      this = specmx(3, self);
    end
    function delete(self)
      if self.swigPtr
        specmx(54, self);
        self.swigPtr=[];
      end
    end
    function varargout = name(self,varargin)
    %Usage: retval = name ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(55, self, varargin{:});
    end
    function varargout = parameters(self,varargin)
    %Usage: retval = parameters (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type Json::Value. 
      [varargout{1:nargout}] = specmx(56, self, varargin{:});
    end
    function varargout = set_parameters(self,varargin)
    %Usage: set_parameters (at, from)
    %
    %at is of type std::string const &. from is of type Json::Value const &. 
      [varargout{1:nargout}] = specmx(57, self, varargin{:});
    end
    function varargout = parameters_(self,varargin)
    %Usage: retval = parameters_ (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type OProp. 
      [varargout{1:nargout}] = specmx(58, self, varargin{:});
    end
    function varargout = set_parameters_(self,varargin)
    %Usage: set_parameters_ (at, from)
    %
    %at is of type std::string const &. from is of type OProp const &. 
      [varargout{1:nargout}] = specmx(59, self, varargin{:});
    end
    function varargout = number_of_configurations(self,varargin)
    %Usage: retval = number_of_configurations ()
    %
    %retval is of type int. 
      [varargout{1:nargout}] = specmx(60, self, varargin{:});
    end
    function varargout = configuration_names(self,varargin)
    %Usage: retval = configuration_names ()
    %
    %retval is of type std::vector< std::string >. 
      [varargout{1:nargout}] = specmx(61, self, varargin{:});
    end
    function varargout = active_configuration(self,varargin)
    %Usage: retval = active_configuration ()
    %
    %retval is of type Remote::Configuration::Shared. 
      [varargout{1:nargout}] = specmx(62, self, varargin{:});
    end
    function varargout = configuration(self,varargin)
    %Usage: retval = configuration (name)
    %
    %name is of type std::string const &. name is of type std::string const &. retval is of type Remote::Configuration::Shared. 
      [varargout{1:nargout}] = specmx(63, self, varargin{:});
    end
    function varargout = activate(self,varargin)
    %Usage: activate (that)
    %
    %that is of type Remote::Configuration::Shared. 
      [varargout{1:nargout}] = specmx(64, self, varargin{:});
    end
    function varargout = clone(self,varargin)
    %Usage: retval = clone (that)
    %
    %that is of type Remote::Configuration::Shared. that is of type Remote::Configuration::Shared. retval is of type Remote::Configuration::Shared. 
      [varargout{1:nargout}] = specmx(65, self, varargin{:});
    end
    function varargout = remove(self,varargin)
    %Usage: remove (that)
    %
    %that is of type Remote::Configuration::Shared. 
      [varargout{1:nargout}] = specmx(66, self, varargin{:});
    end
    function varargout = number_of_stacks(self,varargin)
    %Usage: retval = number_of_stacks ()
    %
    %retval is of type std::size_t. 
      [varargout{1:nargout}] = specmx(67, self, varargin{:});
    end
    function varargout = stack_names(self,varargin)
    %Usage: retval = stack_names ()
    %
    %retval is of type std::vector< std::string >. 
      [varargout{1:nargout}] = specmx(68, self, varargin{:});
    end
    function varargout = stack(self,varargin)
    %Usage: retval = stack (name)
    %
    %name is of type std::string const &. name is of type std::string const &. retval is of type Remote::Stack::Shared. 
      [varargout{1:nargout}] = specmx(69, self, varargin{:});
    end
    function varargout = create_stack(self,varargin)
    %Usage: retval = create_stack (type, sizes)
    %
    %type is of type Remote::Type. sizes is of type std::vector< std::size_t > const &. type is of type Remote::Type. sizes is of type std::vector< std::size_t > const &. retval is of type Remote::Stack::Shared. 
      [varargout{1:nargout}] = specmx(70, self, varargin{:});
    end
    function varargout = update(self,varargin)
    %Usage: update ()
    %
      [varargout{1:nargout}] = specmx(71, self, varargin{:});
    end
    function varargout = save_as(self,varargin)
    %Usage: save_as (path)
    %
    %path is of type std::string const &. 
      [varargout{1:nargout}] = specmx(72, self, varargin{:});
    end
    function self = Measurement(varargin)
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
