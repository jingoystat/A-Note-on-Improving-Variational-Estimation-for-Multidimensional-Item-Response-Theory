function lb = lower_bound(w)

    lb = mean(w, 3);
    lb = log(lb);
    lb = mean(lb, 2);
    lb = sum(lb);
end