classdef markCounter
    %markCounter counts number of marked artifacts
    %   used in parallel processing
    
    properties
        count;
    end
    
    methods
        function self = markCounter ()
            self.count = 1;
        end
        
        function increment(self)
            %self.count = self.count + 1;
            self.count = self.count + 1;
        end
        
        function outcount = getCount(self)
            outcount = self.count;
        end
    end
    
end

