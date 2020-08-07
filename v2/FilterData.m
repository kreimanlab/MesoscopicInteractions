function z=FilterData(x,sr,varargin)
%Fs is sampling frequency
%d is raw signal
% notch filter at 60 Hz and harmonics
%            for f=60:60:Fs/2
%                d = FilterData(d,Fs,'notch',f);
%            end

notch_width_hz = 4;

high=0;
low=0;
notch=0;

for a=1:2:length(varargin)
   argument=lower(cell2mat(varargin(a)));
   value=cell2mat(varargin(a+1));
   eval([argument '=' num2str(value) ';']);
end

filterorder = 5;
if sr > 500
    filterorder = 4;
end
if sr > 5000
    filterorder = 3;
end

z=x;
if high
   [b,a]=butter(filterorder,2*high/sr,'high');
   z=filtfilt(b,a,z);
end
if low
   [b,a]=butter(filterorder,2*low/sr,'low');
   z=filtfilt(b,a,z);
end
if notch
   if (notch > (sr/2))
       notch = sr - notch;
      [b,a]=butter(filterorder,2*[notch-(notch_width_hz/2) notch+(notch_width_hz/2)]/sr,'stop');
   else
      %[b,a]=butter(filterorder,2*[notch-3 notch+3]/sr,'stop');
      range = 2*[notch-(notch_width_hz/2) notch+(notch_width_hz/2)]/sr;
      if (range(end) > 1)
          range(end) = 1-1e-6;
      end
      if (range(1) < 0)
          range(1) = 1e-6;
      end
      [b,a]=butter(filterorder,range,'stop');
   end
   z=filtfilt(b,a,z);
end

return
