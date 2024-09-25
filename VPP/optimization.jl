using JuMP, GLPK

m=Model(GLPK.Optimizer)

@variable(m,x>=0)
@variable(m,y>=0)
@constraint(m,6x+8y>=100)
@constraint(m,7x+12y>=120)
@objective(m,Min,12x+20y)
optimize!(m)
@show value(x)
@show value(y)
@show objective_value(m)