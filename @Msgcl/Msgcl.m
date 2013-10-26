classdef Msgcl
  %MSGCL This class allows provides a logging utility for Matlab
  %   Allows the user to provide log level sensitive reporting for Matlab
  %   functions
  
  properties
    logLevel = 2;
  end
  
  properties (SetAccess = private)
    ALL  = -1;
    ERR  =  2;
    WARN = 10;
    PED  = 99;
  end
  
  properties (SetAccess = private, GetAccess = private)
    fileName = '';
    fid = 1;
    logSource = 'Default ';
  end

  
  methods
    %----------------------------------------------------------------------
    % Name:
    % Description:
    % Input:
    %
    % Output:
    %
    function obj = Msgcl(loglevel, fileName, logSource)
      if nargin == 3
        obj.logSource = obj.setName8(logSource);
      end
      
      if nargin >= 2
        obj.logLevel = loglevel;
        obj.fileName = fileName;
        obj.fid = fopen(fileName,'a');
      else
        obj.logLevel = loglevel;
      end
    end
    
    %----------------------------------------------------------------------
    % Name:
    % Description:
    % Input:
    %
    % Output:
    %
    function pmsg(obj, msgtype, varargin)
      if obj.logLevel >= msgtype
        outstr = [datestr(now,13) ' |' obj.loglvlName(msgtype) '| ' varargin{1} '\n'];
        fprintf(outstr, varargin{2:nargin-2});
        if obj.fid ~= 1
          fprintf(obj.fid, outstr, varargin{2:nargin-2});
        end
      end
    end
    
    %----------------------------------------------------------------------
    % Name:
    % Description:
    % Input:
    %
    % Output:
    %
    function ck = lvlck(obj, msgtype)
      if obj.logLevel >= msgtype
        ck = true;
      else
        ck = false;
      end
     end
    
    
    %----------------------------------------------------------------------
    % Name:
    % Description:
    % Input:
    %
    % Output:
    %
    function delete(obj)
      fclose(obj.fid);
    end 
  
  end % Public Methods
  
  methods (Access = private)
    %----------------------------------------------------------------------
    % Name:
    % Description:
    % Input:
    %
    % Output:
    %
    function str = setName8(obj,in_str)
      str = [in_str, '        '];
      str = str(1:8);
    end

        %----------------------------------------------------------------------
    % Name:
    % Description:
    % Input:
    %
    % Output:
    %
    function str = loglvlName(obj,loglvl)
      switch loglvl
        case obj.ALL
          str = 'ALL ';
        case obj.ERR
          str = 'ERR ';
        case obj.WARN
          str = 'WARN';
        case obj.PED
          str = 'PED ';
      end
    end

  end % Private Methods
  
end

