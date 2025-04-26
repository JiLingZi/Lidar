function saber_yymmdd=Date_format_shfit(x)
Days_Month=[31 28 31 30 31 30 31 31 30 31 30 31 ];
Date=cast(x,'double');
y=floor(Date/1000);
t1=Date-y*1000;
saber_yymmdd=strings(size(y));
for i=1:length(t1)
    for j=1:12
        m=1;
        d=0;
        if t1(i)-sum(Days_Month(1:j))<=0
            m=m+j-1;
            d=t1(i)-sum(Days_Month(1:j-1));
            if d==0
                d=31;
            end
            
            str=strcat(num2str(y(i),'%02d'),num2str(m,'%02d'),num2str(d,'%02d'));
            saber_yymmdd(i)=str;
            break;
        end
    end
end
end

