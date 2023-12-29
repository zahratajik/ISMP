function arr=pad_array(arr,all)


if all==1 %|| %(dim==1 && before==1)   
    arr1=zeros(size(arr,1)+1,size(arr,2));
    arr1(2:end,:) = arr ; %adds a row to top
    clear arr
    arr=arr1;
    clear arr1
end   
    
if all==1 %|| (dim==1 && after==1)
    arr1=zeros(size(arr,1)+1,size(arr,2));
    arr1(1:end-1,:) = arr ; %adds a row to bottom
    clear arr
    arr=arr1;
    clear arr1  ; 
end

if all==1 %|| dim==0 && before==1
     arr1=zeros(size(arr,1),size(arr,2)+1);
     arr1(:,2:end)= arr; %adds col before
     clear arr
     arr=arr1;
     clear arr1  ; 
end
                                               
if all==1 %|| dim==1 && after==1
    arr1=zeros(size(arr,1),size(arr,2)+1);
    arr1(:,1:end-1) = arr ; %adds col after
    clear arr
    arr=arr1;
    clear arr1  ; 
end

end