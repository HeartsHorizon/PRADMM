function q=SO3_to_Q(R)
q=axisangle_to_Q(SO3_to_axisangle(R));
end