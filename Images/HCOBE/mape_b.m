function B = mape_b(B,site_reference)
% Mapping the coefficients(B) to have the same distribution
site_num = length(B);
feature_num = size(B{1,1},1);
b_reference = B{1,site_reference};
for i = 1:feature_num
    for j = 1:site_num
        if j~=site_reference
        mean_reference = mean(b_reference(i,:));
        std_reference = std(b_reference(i,:));
        b_site = B{1,j};
        mean_site = mean(b_site(i,:));
        std_site = std(b_site(i,:));
        b_site(i,:) = std_reference*(b_site(i,:)-mean_site)/std_site+mean_reference;
        B{1,j} = b_site;
        end
    end
end

end

