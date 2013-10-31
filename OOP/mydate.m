classdef mydate
    %MYDATE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        minute
        hour
        day
        month
        year
    end % properties
    properties(Constant = true)
        days_per_year = 365;
        month_per_year = 12;
        weeks_per_year = 52;
    end
    methods
        function obj = mydate(minute,hour,day,month,year) % constructor
            if (nargin > 0)
                obj.minute = minute;
                obj.hour = hour;
                obj.day = day;
                obj.month = month;
                obj.year = year;
            end
        end
        function obj2 = rollDay(obj,numdays)
           obj2 = obj.day + numdays;
        end
    end % methods
    methods(Static = true)
        function printCurrentDate()
            display(datestr(now));
        end
    end
end % Class