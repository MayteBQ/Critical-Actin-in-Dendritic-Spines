function A = area_S(S)
    S_p = [S(2:end,:);S(1,:)];
    A = sum(S(:,1).*S_p(:,2) - S(:,2).*S_p(:,1))/2;
end