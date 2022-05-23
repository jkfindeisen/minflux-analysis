classdef Stack < SwigRef
    %Usage: Stack ()
    %
  methods
    function this = swig_this(self)
      this = specmx(3, self);
    end
    function delete(self)
      if self.swigPtr
        specmx(8, self);
        self.swigPtr=[];
      end
    end
    function varargout = name(self,varargin)
    %Usage: retval = name ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(9, self, varargin{:});
    end
    function varargout = set_name(self,varargin)
    %Usage: set_name (to)
    %
    %to is of type std::string const &. 
      [varargout{1:nargout}] = specmx(10, self, varargin{:});
    end
    function varargout = description(self,varargin)
    %Usage: retval = description ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(11, self, varargin{:});
    end
    function varargout = set_description(self,varargin)
    %Usage: set_description (to)
    %
    %to is of type std::string const &. 
      [varargout{1:nargout}] = specmx(12, self, varargin{:});
    end
    function varargout = type(self,varargin)
    %Usage: retval = type ()
    %
    %retval is of type Remote::Type. 
      [varargout{1:nargout}] = specmx(13, self, varargin{:});
    end
    function varargout = number_of_elements(self,varargin)
    %Usage: retval = number_of_elements ()
    %
    %retval is of type std::size_t. 
      [varargout{1:nargout}] = specmx(14, self, varargin{:});
    end
    function varargout = number_of_dimensions(self,varargin)
    %Usage: retval = number_of_dimensions ()
    %
    %retval is of type std::size_t. 
      [varargout{1:nargout}] = specmx(15, self, varargin{:});
    end
    function varargout = size(self,varargin)
    %Usage: retval = size (dimension)
    %
    %dimension is of type int. dimension is of type int. retval is of type std::size_t. 
      [varargout{1:nargout}] = specmx(16, self, varargin{:});
    end
    function varargout = sizes(self,varargin)
    %Usage: retval = sizes ()
    %
    %retval is of type std::vector< std::size_t >. 
      [varargout{1:nargout}] = specmx(17, self, varargin{:});
    end
    function varargout = label(self,varargin)
    %Usage: retval = label (dimension)
    %
    %dimension is of type int. dimension is of type int. retval is of type std::string. 
      [varargout{1:nargout}] = specmx(18, self, varargin{:});
    end
    function varargout = set_label(self,varargin)
    %Usage: set_label (dimension, to)
    %
    %dimension is of type int. to is of type std::string const &. 
      [varargout{1:nargout}] = specmx(19, self, varargin{:});
    end
    function varargout = labels(self,varargin)
    %Usage: retval = labels ()
    %
    %retval is of type std::vector< std::string >. 
      [varargout{1:nargout}] = specmx(20, self, varargin{:});
    end
    function varargout = set_labels(self,varargin)
    %Usage: set_labels (to)
    %
    %to is of type std::vector< std::string > const &. 
      [varargout{1:nargout}] = specmx(21, self, varargin{:});
    end
    function varargout = length(self,varargin)
    %Usage: retval = length (dimension)
    %
    %dimension is of type int. dimension is of type int. retval is of type double. 
      [varargout{1:nargout}] = specmx(22, self, varargin{:});
    end
    function varargout = set_length(self,varargin)
    %Usage: set_length (dimension, to)
    %
    %dimension is of type int. to is of type double. 
      [varargout{1:nargout}] = specmx(23, self, varargin{:});
    end
    function varargout = lengths(self,varargin)
    %Usage: retval = lengths ()
    %
    %retval is of type std::vector< double >. 
      [varargout{1:nargout}] = specmx(24, self, varargin{:});
    end
    function varargout = set_lengths(self,varargin)
    %Usage: set_lengths (to)
    %
    %to is of type std::vector< double > const &. 
      [varargout{1:nargout}] = specmx(25, self, varargin{:});
    end
    function varargout = offset(self,varargin)
    %Usage: retval = offset (dimension)
    %
    %dimension is of type int. dimension is of type int. retval is of type double. 
      [varargout{1:nargout}] = specmx(26, self, varargin{:});
    end
    function varargout = set_offset(self,varargin)
    %Usage: set_offset (dimension, to)
    %
    %dimension is of type int. to is of type double. 
      [varargout{1:nargout}] = specmx(27, self, varargin{:});
    end
    function varargout = offsets(self,varargin)
    %Usage: retval = offsets ()
    %
    %retval is of type std::vector< double >. 
      [varargout{1:nargout}] = specmx(28, self, varargin{:});
    end
    function varargout = set_offsets(self,varargin)
    %Usage: set_offsets (to)
    %
    %to is of type std::vector< double > const &. 
      [varargout{1:nargout}] = specmx(29, self, varargin{:});
    end
    function varargout = data(self,varargin)
    %Usage: retval = data ()
    %
    %retval is of type Remote::Data. 
      [varargout{1:nargout}] = specmx(30, self, varargin{:});
    end
    function varargout = set_data(self,varargin)
    %Usage: set_data (data)
    %
    %data is of type Remote::Data const. 
      [varargout{1:nargout}] = specmx(31, self, varargin{:});
    end
    function varargout = signal_data_changed(self,varargin)
    %Usage: signal_data_changed ()
    %
      [varargout{1:nargout}] = specmx(32, self, varargin{:});
    end
    function varargout = meta_data(self,varargin)
    %Usage: retval = meta_data ()
    %
    %retval is of type Json::Object. 
      [varargout{1:nargout}] = specmx(33, self, varargin{:});
    end
    function varargout = share(self,varargin)
    %Usage: retval = share ()
    %
    %retval is of type std::string. 
      [varargout{1:nargout}] = specmx(34, self, varargin{:});
    end
    function varargout = update_from_ref(self,varargin)
    %Usage: update_from_ref ()
    %
      [varargout{1:nargout}] = specmx(35, self, varargin{:});
    end
    function varargout = reference(self,varargin)
    %Usage: retval = reference ()
    %
    %retval is of type Remote::StackReference::Shared. 
      [varargout{1:nargout}] = specmx(36, self, varargin{:});
    end
    function self = Stack(varargin)
      if nargin==1 && strcmp(class(varargin{1}),'SwigRef')
        if ~isnull(varargin{1})
          self.swigPtr = varargin{1}.swigPtr;
        end
      else
        tmp = specmx(37, varargin{:});
        self.swigPtr = tmp.swigPtr;
        tmp.swigPtr = [];
      end
    end
  end
  methods(Static)
  end
end
