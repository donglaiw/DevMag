function d = U_lenGeo(geo_type,geo_param)
d = zeros(numel(geo_type),1);
for j=0:2
    ind = reshape(find(geo_type==j),1,[]);
    if ~isempty(ind)
        switch j
            case 0
                % line
                d(ind) = sqrt(sum((1+abs(geo_param(ind,[1 2])-geo_param(ind,[3 4]))).^2,2));
            case 1
                % circle
                d(ind) = geo_param(ind,3).*abs(geo_param(ind,7)-geo_param(ind,6));
            case 2
                % ellipse
                % approx elliptic intregal
                %rr2  = sort(geo_param(:,3:4),2,'ascend').^2;
                %m = (rr2(:,2)-rr2(:,1))./rr2(:,2);
                for i = ind
                    %d(i) = rr2(i,2)*abs(U_elliptic12(geo_param(i,7),m(i))-U_elliptic12(geo_param(i,6),m(i)));
                    d(i) = U_elliptic(geo_param(i,3),geo_param(i,4),geo_param(i,6),geo_param(i,7));
                end
                % check:
                %{
            [10*pi U_lenGeo(2,[0 0 10 10 0 0 pi])]
             
            ellp1 = @(ab) 2*pi*sqrt(sum(ab.^2)/2);
            ellp2 = @(ab) pi*(3*sum(ab)-sqrt((ab(1)+3*ab(2))*(3*ab(1)+ab(2))));
            [ellp1([5 10]) ellp2([5 10]) U_lenGeo(2,[0 0 5 10 0 0 2*pi])]
            2*pi*sqrt(sum([5 10].^2)/2)
                %}
        end
    end
end
