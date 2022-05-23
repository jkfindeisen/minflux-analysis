classdef Imspector < SwigRef
    %Usage: Imspector ()
    %
  methods
    function this = swig_this(self)
      this = specmx(3, self);
    end
    function delete(self)
      if self.swigPtr
        specmx(73, self);
        self.swigPtr=[];
      end
    end
    function varargout = host(self,varargin)
    %Usage: retval = host ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(74, self, varargin{:});
    end
    function varargout = set_host(self,varargin)
    %Usage: set_host (host)
    %
    %host is of type std::string const &. 
      [varargout{1:nargout}] = specmx(75, self, varargin{:});
    end
    function varargout = version(self,varargin)
    %Usage: retval = version ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(76, self, varargin{:});
    end
    function varargout = device_drivers(self,varargin)
    %Usage: retval = device_drivers ()
    %
    %retval is of type Json::Object. 
      [varargout{1:nargout}] = specmx(77, self, varargin{:});
    end
    function varargout = parameters(self,varargin)
    %Usage: retval = parameters (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type Json::Value. 
      [varargout{1:nargout}] = specmx(78, self, varargin{:});
    end
    function varargout = set_parameters(self,varargin)
    %Usage: set_parameters (at, from)
    %
    %at is of type std::string const &. from is of type Json::Value const &. 
      [varargout{1:nargout}] = specmx(79, self, varargin{:});
    end
    function varargout = parameters_(self,varargin)
    %Usage: retval = parameters_ (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type OProp. 
      [varargout{1:nargout}] = specmx(80, self, varargin{:});
    end
    function varargout = set_parameters_(self,varargin)
    %Usage: set_parameters_ (at, from)
    %
    %at is of type std::string const &. from is of type OProp const &. 
      [varargout{1:nargout}] = specmx(81, self, varargin{:});
    end
    function varargout = state(self,varargin)
    %Usage: retval = state (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type Json::Value. 
      [varargout{1:nargout}] = specmx(82, self, varargin{:});
    end
    function varargout = set_state(self,varargin)
    %Usage: set_state (at, from)
    %
    %at is of type std::string const &. from is of type Json::Value const &. 
      [varargout{1:nargout}] = specmx(83, self, varargin{:});
    end
    function varargout = state_(self,varargin)
    %Usage: retval = state_ (at)
    %
    %at is of type std::string const &. at is of type std::string const &. retval is of type OProp. 
      [varargout{1:nargout}] = specmx(84, self, varargin{:});
    end
    function varargout = set_state_(self,varargin)
    %Usage: set_state_ (at, from)
    %
    %at is of type std::string const &. from is of type OProp const &. 
      [varargout{1:nargout}] = specmx(85, self, varargin{:});
    end
    function varargout = measurement_names(self,varargin)
    %Usage: retval = measurement_names ()
    %
    %retval is of type std::vector< std::string >. 
      [varargout{1:nargout}] = specmx(86, self, varargin{:});
    end
    function varargout = active_measurement(self,varargin)
    %Usage: retval = active_measurement ()
    %
    %retval is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(87, self, varargin{:});
    end
    function varargout = measurement(self,varargin)
    %Usage: retval = measurement (name)
    %
    %name is of type std::string const &. name is of type std::string const &. retval is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(88, self, varargin{:});
    end
    function varargout = create_measurement(self,varargin)
    %Usage: retval = create_measurement ()
    %
    %retval is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(89, self, varargin{:});
    end
    function varargout = open(self,varargin)
    %Usage: retval = open (path)
    %
    %path is of type std::string const &. path is of type std::string const &. retval is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(90, self, varargin{:});
    end
    function varargout = activate(self,varargin)
    %Usage: activate (that)
    %
    %that is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(91, self, varargin{:});
    end
    function varargout = start(self,varargin)
    %Usage: start (that)
    %
    %that is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(92, self, varargin{:});
    end
    function varargout = pause(self,varargin)
    %Usage: pause (that)
    %
    %that is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(93, self, varargin{:});
    end
    function varargout = stop(self,varargin)
    %Usage: stop (that)
    %
    %that is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(94, self, varargin{:});
    end
    function varargout = close(self,varargin)
    %Usage: close (that)
    %
    %that is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(95, self, varargin{:});
    end
    function varargout = active_stack(self,varargin)
    %Usage: retval = active_stack ()
    %
    %retval is of type Remote::Stack::Shared. 
      [varargout{1:nargout}] = specmx(96, self, varargin{:});
    end
    function varargout = connect_begin(self,varargin)
    %Usage: connect_begin (to, flag)
    %
    %to is of type Omas::Shared_Slot< void () >. flag is of type int. 
      [varargout{1:nargout}] = specmx(97, self, varargin{:});
    end
    function varargout = disconnect_begin(self,varargin)
    %Usage: disconnect_begin (from, flag)
    %
    %from is of type Omas::Shared_Slot< void () >. flag is of type int. 
      [varargout{1:nargout}] = specmx(98, self, varargin{:});
    end
    function varargout = connect_end(self,varargin)
    %Usage: connect_end (to, flag)
    %
    %to is of type Omas::Shared_Slot< void () >. flag is of type int. 
      [varargout{1:nargout}] = specmx(99, self, varargin{:});
    end
    function varargout = disconnect_end(self,varargin)
    %Usage: disconnect_end (from, flag)
    %
    %from is of type Omas::Shared_Slot< void () >. flag is of type int. 
      [varargout{1:nargout}] = specmx(100, self, varargin{:});
    end
    function self = Imspector(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = specmx(101, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
    function varargout = run(self,varargin)
    %Usage: run (that)
    %
    %that is of type Remote::Measurement::Shared. 
      [varargout{1:nargout}] = specmx(102, self, varargin{:});
    end
  end
  methods(Static)
  end
end
