function displayCells(allCells,cells, adult_active, ifr_bins, adult_IFR)

%get the adult from ISI_to_active routine
young=0;
old=0;

for i=1:length(cells)
    cellNo=cells(i);
    if (allCells{1,i}.dph<63)
        
      
        if(~isnan(allCells{1,cellNo}.activeNoDrug{1,1}(1)))
            young=young+1;
            activeNoDrug_young(young,1:length(allCells{1,cellNo}.activeNoDrug{1,1}))=allCells{1,cellNo}.activeNoDrug{1,1};
            ISINoDrug_young(young,1:length(allCells{1,cellNo}.ISINoDrug{1,1}))=allCells{1,cellNo}.ISINoDrug{1,1};
        end
            %semilogy(allCells{1,i}.activeNoDrug{1,2},allCells{1,i}.activeNoDrug{1,1}, 'color',color);
        %hold on;
        
    elseif (allCells{1,cellNo}.dph>77)
        
        
        if(~isnan(allCells{1,cellNo}.activeNoDrug{1,1}(1)))
            old=old+1;
            activeNoDrug_old(old,1:length(allCells{1,cellNo}.activeNoDrug{1,1}))=allCells{1,cellNo}.activeNoDrug{1,1};
            ISINoDrug_old(old,1:length(allCells{1,cellNo}.ISINoDrug{1,1}))=allCells{1,cellNo}.ISINoDrug{1,1};
        end
            %semilogy(allCells{1,i}.activeNoDrug{1,2},allCells{1,i}.activeNoDrug{1,1}, 'color',color);
        %hold on;
        
    end
end

activeNoDrug_old=sum(activeNoDrug_old,1)/old;
activeNoDrug_young=sum(activeNoDrug_young,1)/young;
ISINoDrug_old=sum(ISINoDrug_old,1)/old;
ISINoDrug_young=sum(ISINoDrug_young,1)/young;
figure;
semilogy(allCells{1,1}.activeNoDrug{1,2},activeNoDrug_young, 'color','r');
hold on;
semilogy(allCells{1,1}.activeNoDrug{1,2},activeNoDrug_old, 'color','k');
semilogy(ifr_bins,adult_active,'color','g');
hold off;
figure;
plot(allCells{1,1}.ISINoDrug{1,2},ISINoDrug_young, 'color','r');
hold on;
plot(allCells{1,1}.ISINoDrug{1,2},ISINoDrug_old, 'color','k');
plot(ifr_bins,adult_IFR,'color','g');
hold off;
young
old
