function f = wrapper_b4_obj(par)

[objective,~,~] = b4_obj(par);
f = objective;

if isnan(f) || isinf(f)
    f = realmax;
end
