function f = wrapper_b6_obj(par)

[objective,~,~] = b6_obj(par);
f = objective;

if isnan(f) || isinf(f)
    f = realmax;
end
