function f = wrapper_b2_obj(par)

[objective,~,~] = b2_obj(par);
f = objective;

if isnan(f) || isinf(f)
    f = realmax;
end
