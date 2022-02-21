classdef TreeNode %< handle
    
    
    properties (Constant)
        % global constants see Par file
    end
    
    properties
        id = 0;
        type = 0;
        x=0;
        y=0;
        z=0;
        label = '';
        dist_from_parent = 0;
        childrenIndex = [];
        parentIndex = 0;
        parentid = 0;
    end
    methods
        function obj = TreeNode(swcrow)
             obj.id = str2double(swcrow{1});
             obj.type = str2double(swcrow{2});
             
             obj.x = str2double(swcrow{3});
             obj.y = str2double(swcrow{4});
             obj.z = str2double(swcrow{5});
             
             obj.dist_from_parent = str2double(swcrow{7});
             obj.label = lower(swcrow{8});
             obj.parentid = str2double(swcrow{6});
        end
        
    end
end