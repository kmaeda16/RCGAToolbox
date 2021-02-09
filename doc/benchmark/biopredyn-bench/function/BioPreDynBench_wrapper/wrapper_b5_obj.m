function f = wrapper_b5_obj(par)

[f,~,~] = b5_obj(par);

if isnan(f) || isinf(f)
    f = realmax;
end
