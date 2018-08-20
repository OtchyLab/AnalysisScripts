function SAPAnalyze(data, text, DOBstr)
%Purpose of this code is to take the .csv files produced by SAP (via
%Excel...) and process them for Variability over time.  The files need to
%first be imported into MATLAB (as separate text and data files) so that 
%they are accessible to the function call.  The input arguments 'data' and
%'text' should point to these data sources; 'DOBstr' is the bird's DOB in
%'mm-dd-yyyy' format.

%Convert string DOB to univsersal time code
DOB=datenum(DOBstr, 'mm-dd-yyyy');                 %For Rd153, DOB=09-20-2008

%format inbound data for only 'useful' columns
text = [text(:,2) text(:,3)];       %
data=[data(:,2)];
 
% %Remove X-X similarity scores
% k=(1:length(data))';
% index(k,1)=strcmp(text(k,1),text(k,2));
% one=ones(length(data),1);
% index=(index-one).*-k;
% m=1;
% for k=1:length(data)
%      if index(k,1)~=0;
%          date(m,1)=text(index(k,1),1);
%          sim(m,1)=data(index(k,1),1);
%          m=m+1;
%      end
% end
dateraw=text(:,2);
sim=data;

% %Strip out unuseful info and collect file month
n=1:length(sim);
% date(n,1) = regexprep(date(n,1), '\s', '');
% date(n,1) = regexprep(date(n,1), '1Rd153_', '');
% date(n,1) = regexprep(date(n,1), '.wav', '');
% month(n,1) = regexprep(date(n,1), '-(\w*)_(\w*)', '');
% 
% %Strip out Syllable Type 
% type(n,1) = regexprep(date(n,1), '(\w*)-(\w*)_', '');
% type(n,1) = regexprep(type(n,1), '(\w*)_', '');
% %type=str2double(type);
% 
% %Recover MM-DD from filename
% date(n,1) = (regexprep(date(n,1), '_(\w*)', ''));
% 
% %Assign year based on month file recorded
% for k=1:length(sim)
%      if str2double(month(k,1))>5
%          year(k,1)=2008;
%      else
%          year(k,1)=2009;
%      end
% end
% year = num2str(year);
date = [];
type=cell(length(sim),1);
for n=1:length(sim)
    %Strip out Date
    q=regexp(char(dateraw(n)), '_', 'split');
    daten = datenum([char(q(3)) char('-') char(q(4)) char('-') char(q(5))],'yyyy-mm-dd');
    date(n)= daten-DOB;
    %Strip out Syl Type
    type(n) = regexprep(q(11), '.wav', '');
    
end

%Reassemble string date for calculating DPH
% dash(n,1)='-';
% date=char(date);
% date = [date, dash, year];
%daten = datenum(date, 'yyyy-mm-dd');
%DOBn(n,1) = DOB;
%nine(n,1) = 9;  At one time, the date seemed to be off for some reason.
%date = (daten-DOBn)-nine;  
%date = (daten-DOBn');


%Sort data by syllable type in to separate vectors (date,sim)
syltypes = sort(unique(type));
%syltypes = ['A'; 'B'; 'C'; 'D'; 'E', 'F'];
n=1:length(sim);
for k = 1:length(syltypes);
    curtype = char(genvarname(syltypes(k)));
    m=1;
    for n = 1:length(sim);
        if strcmp(type(n,1),syltypes(k)) == 1;
            if k==4 & date(n)<74
            else   
                eval([curtype '(m,1)=date(n);']);              %1st array is the date
                eval([curtype '(m,2)=sim(n);']);               %2nd array is the similarity
                m=m+1;
            end
        end
    end
    %And plot those vectors on individual axes (similarity vs dph)
    subplot(length(syltypes),2,k*2-1);
    plot(eval([curtype '(:,1)']),eval([curtype '(:,2)']),'.');
    %xlim([35 max(eval([curtype '(:,1)']))]);
    xlim([35 max(date)+15]);
    ylim([0 100]);
    ylabel('Similarity (ADU)');
    xlabel('Days Post Hatch');
    hold on;
    
    %Calculate and plot variability
    curvar = char(genvarname(['var' char(syltypes(k))]));
    eval([curvar '=sort(unique(' curtype '(:,1)));']);
    s=1;
    for p=1:length(eval([curvar '(:,1)']));
        filecnt=0;
        eval([curvar '(p,2)=0;']);
        for q=1:length(eval([curtype '(:,2);']));
            if eval([curvar '(p,1)']) == eval([curtype '(q,1)']);
                eval([curvar '(p,2)=' curvar '(p,2) + ' curtype '(q,2);']);
                filecnt = filecnt + 1;
            end
        end
        eval([curvar '(p,2) = ' curvar '(p,2) / filecnt;']);       %Calc mean similarity
        eval([curvar '(p,2)=((100-' curvar '(p,2))/50)*100;']);       %Calc variability
    
        %Load variability data into an array that can be used to calc average
        %varability for all syllables
        AllVar(eval([curvar '(p,1)']),1)=eval([curvar '(p,1)']);
        AllVar(eval([curvar '(p,1)']),k+2)=eval([curvar '(p,2)']);
    end

    %Plot variability on left-hand axes
    subplot(length(syltypes),2,k*2-1);
    plot(eval([curvar '(:,1)']),eval([curvar '(:,2)']));
    hold on
    
    %Plot variability (lightly)on right hand axis
    axes_sel=[1:length(syltypes)]*2;
    subplot(length(syltypes),2,axes_sel);
    plot(eval([curvar '(:,1)']),eval([curvar '(:,2)']), 'g');
    xlim([35 max(date)+15]);
    ylim([0 100]);
    ylabel('Similarity (ADU)');
    xlabel('Days Post Hatch');
    hold on;

end

%Calculate average variability (anf fitline) across syllables accros time.
%z=1;
%for x=1:length(AllVar);
%    if AllVar(x,1) ~= 0;
%        AllVar(z,:) = AllVar(x,:);
%        AllVar(z,(length(syltypes)+3)) = nnz(AllVar(z,:))-1;
%        AllVar(z,2) = (sum(AllVar(z,3:(length(syltypes)+2))))/AllVar(z,(length(syltypes)+3));
%        z=z+1;
%    end
%end

R = find(AllVar(:,1)>1)
AllVarStrp = AllVar(R,:);
for x=1:length(AllVarStrp)
    AllVarStrp(x,(length(syltypes)+3)) = nnz(AllVarStrp(x,:))-1;
    AllVarStrp(x,2) = (sum(AllVarStrp(x,3:(length(syltypes)+2))))/AllVarStrp(x,(length(syltypes)+3));
end

fitline = polyfit(AllVarStrp(:,1), AllVarStrp(:,2),1); % least squares fitting to a line
a1 = fitline(2) % y-intercept of the fitted line
a2 = fitline(1) % slope of fitted lines

%Plot it!
subplot(length(syltypes),2,axes_sel);
plot(AllVarStrp(:,1), AllVarStrp(:,2),'dr');
plot(AllVarStrp(:,1),fit,'k')
hold off
a=1;


    
   