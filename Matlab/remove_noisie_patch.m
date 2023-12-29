function final_matrix = remove_noisie_patch(final_matrix,n, index)

Middle_label = final_matrix(n-1).label;
Middle_Area = final_matrix(n-1).Area;
Middle_xc = final_matrix(n-1).xc;
Middle_yc = final_matrix(n-1).yc;
edge_=final_matrix(n-1).edge;

% Removing information noise patches 
for i=1:length(index)
 
if length((Middle_Area(index(i)))==1)
% round(Middle_yc(index(i)))
% round(Middle_xc(index(i)))
edge_(round(Middle_yc(index(i))),round(Middle_xc(index(i))))= 0;
else 

edge_(round( Middle_yc(index(i))),(round(Middle_xc(index(i))-1)))=0;
edge_(round( Middle_yc(index(i))),(round(Middle_xc(index(i))+1)))=0;
edge_(round(Middle_yc(index(i))-1),round(Middle_xc(index(i))))=0;
edge_(round( Middle_yc(index(i))+1),round(Middle_xc(index(i))))=0;
edge_(round(Middle_yc(index(i))-1),round(Middle_xc(index(i))-1))=0;
edge_(round(Middle_yc(index(i))+1),round(Middle_xc(index(i))+1))=0;
edge_(round(Middle_yc(index(i))-1),round(Middle_xc(index(i))+1))=0;
edge_(round(Middle_yc(index(i))+1),round(Middle_xc(index(i))-1))=0;

end    
    
    
Middle_label{index(i)}=[];    
Middle_Area(index(i))=0;
Middle_xc(index(i))=0;
Middle_yc(index(i))=0;


end

w1= find(Middle_Area==0);
Middle_Area(w1)=[]; 

w2= find(Middle_xc==0);
Middle_xc(w2)=[];

w3= find(Middle_yc==0);
Middle_yc(w3)=[];

Middle_label = Middle_label(~cellfun('isempty',Middle_label));

final_matrix(n-1).label = Middle_label;
final_matrix(n-1).Area = Middle_Area;
final_matrix(n-1).xc = Middle_xc;
final_matrix(n-1).yc = Middle_yc;
final_matrix(n-1).edge =edge_;
end


