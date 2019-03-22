%EXAMPLE_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef example_interface < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = example_interface(varargin)
            this.objectHandle = example_mex('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            example_mex('delete', this.objectHandle);
        end
        
        %%Calculate new input
        function varargout = Evaluate(this, varargin)
            [varargout{1:nargout}] = example_mex('Evaluate', this.objectHandle, varargin{:});
        end
    
        %%Set current state
        function varargout = set_current_state(this, varargin)
            [varargout{1:nargout}] = example_mex('set_current_state', this.objectHandle, varargin{:});
        end
        
        %%Set desire state
        function varargout = set_desire_state(this, varargin)
            [varargout{1:nargout}] = example_mex('set_desire_state', this.objectHandle, varargin{:});
        end
        
        %%Set desire stacionary input
        function varargout = set_desire_stacionary_input(this, varargin)
            [varargout{1:nargout}] = example_mex('set_desire_stacionary_input', this.objectHandle, varargin{:});
        end
    end
end