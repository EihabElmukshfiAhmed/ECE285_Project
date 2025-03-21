using PowerModels,Ipopt, DataFrames, XLSX
using Gurobi, JuMP, DataFrames

 
df_price_rtm = DataFrame(XLSX.readtable("ECE285/prices_RTM.xlsx", "Sheet1")[1],XLSX.readtable("ECE285/prices_RTM.xlsx", "Sheet1")[2]) 

df_price = DataFrame(XLSX.readtable("ECE285/prices.xlsx", "Sheet1")[1],XLSX.readtable("ECE285/prices.xlsx", "Sheet1")[2]) 

price_matrix = Matrix(df_price[:,1:9])'

# Reshape the matrix to size (24, 6, 6)
price_3d = reshape(Matrix(df_price[:,1:9])', 3, 3, 24)

Δ_d = DataFrame(XLSX.readtable("ECE285/load_change.xlsx", "load")[1],XLSX.readtable("ECE285/load_change.xlsx", "load")[2]) 





case_file = "case6ww.m"
network_data = PowerModels.parse_file(case_file) 

# define bidding price

γ_G= reshape(Matrix(df_price[:,1:9])', 3, 3, 24) /network_data["baseMVA"]
γ_D= reshape(Matrix(df_price[:,10:18])', 3, 3, 24) /network_data["baseMVA"]
γ_CB= reshape(Matrix(df_price[:,37:42])', 6, 24) /network_data["baseMVA"]


γ_G=γ_G[1,:,:]
γ_D=γ_D[1,:,:]

num_buses = length(network_data["bus"]) 
function create_loads_dictionary(
    load_type::String,
    num_loads::Int,
    buses::Vector{Int},
    min_values::Vector{Float64},
    max_values::Vector{Float64},
    segments::Any
)
    # Initialize an empty dictionary to store the loads/EVs/etc
    loads_dict = Dict{String, Dict{String, Any}}()

    # Check if the input vectors have the correct size
    if length(buses) != num_loads || length(min_values) != num_loads || length(max_values) != num_loads
        throw(ArgumentError("The length of buses, min_values, and max_values must match num_loads."))
    end

    # Populate the dictionary with load information
    for i in 1:num_loads
        bus_id = buses[i]  # Current bus ID

        # Add the load information
        load_info = Dict(
            "load_type" => load_type,
            load_type * "_" * "bus" => bus_id,
            "min_value" => min_values[i]/network_data["baseMVA"],
            "max_value" => max_values[i]/network_data["baseMVA"]
        )

        loads_dict[string(i)] = load_info
    end

    return loads_dict
end












#----------------------------------------------ADD CB-------------------------
max_CB=5
load_type = "CB"
num_loads = length(network_data["bus"])
buses = [1, 2, 3,4,5,6]
min_values = ones(num_buses)*-max_CB
max_values = ones(num_buses)*max_CB

CB_dictionary = create_loads_dictionary(load_type, num_loads, buses, min_values, max_values,"None")

network_data["CB"]=CB_dictionary 



num_gens = length(network_data["gen"]) 
num_loads=length(network_data["load"])
num_CBs=length(network_data["CB"])
num_Brs=length(network_data["branch"])
T=1:24
D=1:num_loads
I=1:num_buses
G=1:num_gens

CB=1:num_CBs
Br=1:num_Brs



function get_indicident_matrix(component::String)
    if haskey(network_data,component)==false
        throw(ArgumentError("component Does not Exsit!"))
    end
    if component=="branch"
        incidence_matrix_line = zeros(Int, num_buses, num_buses)
        capacity_matrix = zeros(Float64, num_buses, num_buses)
        # Populate the incidence matrix
        branch_data = sort(network_data["branch"])
        for branch in values(branch_data)
            from_bus = branch["f_bus"]
            to_bus = branch["t_bus"]

            # Set the appropriate entries
            incidence_matrix_line[from_bus, to_bus] = 1   # From bus
            incidence_matrix_line[to_bus, from_bus] = -1   # To bus
            thermal_limit = get(branch, "rate_a", 0.0)  # Default to 0 if not available

            # Set the capacity in the matrix
            capacity_matrix[from_bus, to_bus] = thermal_limit
            capacity_matrix[to_bus, from_bus] = thermal_limit  # If considering both directions
        end
    return capacity_matrix
    
    #How to end this function here?
    else
        comp_data = sort(network_data[component])
        component_length=length(comp_data)
        incidence_matrix = zeros(Int,num_buses,component_length)
        

        for (i, x) in enumerate(values(comp_data))
    
            bus=x[component*"_"*"bus"]
            
        
            incidence_matrix[bus,i]=1
        
        end 
        
    end


    return Array(incidence_matrix)
end


#--------------------------------Optimization Code-----------------------------------------------

A_i_g = get_indicident_matrix("gen")
A_i_d = get_indicident_matrix("load")
A_i_cb = get_indicident_matrix("CB")




P_max_line  = get_indicident_matrix("branch")  
B=PowerModels.calc_basic_susceptance_matrix(network_data)



γ_G_up= Array(df_price_rtm[!,1])' /network_data["baseMVA"]
γ_G_down= Array(df_price_rtm[!,2])' /network_data["baseMVA"]





get_max_values = load -> [sum(load["segment"][segment]["max_value"] for segment in keys(load["segment"]))][1]  # get maximum gen/load from the segments

discretization_range = collect(-0.1:.005:0.1)

# Define the discretization range and the number of discrete values
# discretization_range = LinRange(-max_CB,max_CB,20)/network_data["baseMVA"]  # Discretization range
Q = length(discretization_range)          # Number of discrete values
β_values = ones(Q,num_buses,24)

# Fill the array with the discretization range
for i in 1:num_buses
    for t in T
        β_values[:,i, t] .= discretization_range
    end
end



#----------------------------------KKT-------------------------------------------------------------
M_g=1e1
M_d=1e1
M_θ=1e1
M_𝜱=1e1
C_B= Model(Gurobi.Optimizer)  # optimal bidding maximization problem
# set_optimizer_attribute(C_B,"MIPGap",0.000001)

@variable(C_B,P_g_kkt[ g in G,t in T]>=0)     #generation units
@variable(C_B,P_d_kkt[ d in D,t in T]>=0)     #demand units
@variable(C_B,P_cb_kkt[cb in CB,t in T])     #CB units
@variable(C_B,θ_kkt[i in I,t in T])     #Voltage Angle

## dual variables
@variable(C_B,μ_ub_g_kkt[ g in G,t in T]>=0)
@variable(C_B,μ_lb_g_kkt[ g in G,t in T]>=0)
@variable(C_B,μ_ub_d_kkt[ d in D,t in T]>=0)
@variable(C_B,μ_lb_d_kkt[ d in D,t in T]>=0)

@variable(C_B,μ_lb_cb_kkt[cb in CB,t in T]>=0)
@variable(C_B,μ_ub_cb_kkt[cb in CB,t in T]>=0)

@variable(C_B,λ_i_kkt[i in I, t in T])
@variable(C_B,π_i_t_ub_kkt[I,I,T]>=0)
@variable(C_B,π_i_t_lb_kkt[I,I,T]>=0)
@variable(C_B,z_g_o_kkt[G,T],Bin)
@variable(C_B,z_g_u_kkt[G,T],Bin)
@variable(C_B,z_d_u_kkt[D,T],Bin)
@variable(C_B,z_d_o_kkt[D,T],Bin)

@variable(C_B,z_cb_u_kkt[CB,T],Bin)
@variable(C_B,z_cb_o_kkt[CB,T],Bin)


@variable(C_B,z_line_o_kkt[I,I,T],Bin)
@variable(C_B,z_line_u_kkt[I,I,T],Bin)

 

@variable(C_B,P_g_up_kkt[G,T]>=0)
@variable(C_B,P_g_down_kkt[G,T]>=0)



@variable(C_B,θ_kkt_rtm[I,T])

@variable(C_B, π_rtm_ub_kkt[ I, I, T]>=0)
@variable(C_B, π_rtm_lb_kkt[ I, I, T]>=0)


@variable(C_B, λ_rtm_kkt[ I, T])

@variable(C_B, μ_G_ub_up_kkt[ G, T]>=0)
@variable(C_B, μ_G_ub_down_kkt[ G, T]>=0)
@variable(C_B, μ_G_lb_up_kkt[ G, T]>=0)
@variable(C_B, μ_G_lb_down_kkt[ G, T]>=0)




@variable(C_B, z_G_ub_up_kkt[G,T],Bin)
@variable(C_B, z_G_ub_down_kkt[G,T],Bin)
@variable(C_B, z_G_lb_up_kkt[G,T],Bin)
@variable(C_B, z_G_lb_down_kkt[G,T],Bin)


@variable(C_B, z_rtm_ub_kkt[I,I,T],Bin)
@variable(C_B, z_rtm_lb_kkt[I,I,T],Bin)


@variable(C_B, z_kkt[1:Q,I,T], Bin)  # Binary variables z_q

@variable(C_B, 𝜱_rtm[1:Q,I,T] >=0)  # Auxiliary variables φ_q



@objective(C_B, Max,
    sum(

        + sum(sum(γ_G[g, t] * P_g_kkt[g, t]  for g in G))
        - sum(sum(γ_D[d, t] * P_d_kkt[d, t]  for d in D))
        +sum(network_data["gen"][string(g)]["pg"]*μ_ub_g_kkt[g,t] for  g in G)
        +sum(network_data["load"][string(d)]["pd"]*μ_ub_d_kkt[d,t] for d in D)
        +sum(sum(P_max_line[i,j]*(π_i_t_ub_kkt[i,j,t]+π_i_t_lb_kkt[i,j,t]) for j in I) for i in I)
        # -sum( sum(β_values[q, cb, t] * 𝜱_rtm[q, cb, t] for q in 1:Q)  for cb in CB )   # Sum over q

        -sum(  λ_rtm_kkt[network_data["CB"][string(cb)]["CB_bus"],t] * P_cb_kkt[cb, t] for cb in CB )   # Sum over q

       for t in T
             )             # Sum over i and t
)

# # Linearization of λ_rtm_kkt*P_cb_kkt using binary expansion

# @constraint(C_B,[ cb in CB, t in T],  sum(β_values[q, cb, t]*z_kkt[q, cb, t] for q in 1:Q)-0.025 <=P_cb_kkt[cb,t] )   

# @constraint(C_B,[ cb in CB, t in T],  P_cb_kkt[cb,t] <=  sum(β_values[q, cb, t]*z_kkt[q, cb, t] for q in 1:Q)+0.025)   


# @constraint(C_B,[ cb in CB, t in T],sum(z_kkt[q, cb, t] for q in 1:Q)==1) 

# @constraint(C_B,[q in 1:Q, cb in CB, t in T],𝜱_rtm[q, cb, t]<=M_𝜱*z_kkt[q, cb, t])   

# @constraint(C_B,[q in 1:Q, cb in CB, t in T],λ_rtm_kkt[network_data["CB"][string(cb)]["CB_bus"],t]-𝜱_rtm[q, cb, t]<=M_𝜱*(1-z_kkt[q, cb, t]))

# @constraint(C_B,[q in 1:Q, cb in CB, t in T],  𝜱_rtm[q, cb, t] <= λ_rtm_kkt[network_data["CB"][string(cb)]["CB_bus"],t]  ) 

#--------End of Linearization Approach - ---------------------------------


@constraint(C_B, [ i in I , t in T], 
sum(P_cb_kkt[cb, t] for cb in CB if network_data["CB"][string(cb)]["CB_bus"] == i) 
+sum(P_g_kkt[g,t]  for g in G if network_data["gen"][string(g)]["gen_bus"]==i)

-sum(P_d_kkt[d,t]  for d in D  if network_data["load"][string(d)]["load_bus"]==i )

- sum(B[i, j] * (θ_kkt[i, t] - θ_kkt[j, t]) for j in I )

    ==0)

@constraint(C_B, [g in G, t in T], P_g_kkt[g, t]  <= network_data["gen"][string(g)]["pg"])


@constraint(C_B, [d in D, t in T], P_d_kkt[d, t]  <= network_data["load"][string(d)]["pd"])


@constraint(C_B,[i in I,j in I,t in T], -P_max_line[i,j] <= B[i,j]*(θ_kkt[i,t]-θ_kkt[j,t]) <= P_max_line[i,j] )   

@constraint(C_B,[cb in CB,t in T],network_data["CB"][string(cb)]["min_value"] <= P_cb_kkt[cb,t]<=network_data["CB"][string(cb)]["max_value"] )

@constraint(C_B,[ g in G, t in T],γ_G[g,t] + λ_i_kkt[network_data["gen"][string(g)]["gen_bus"],t] + μ_ub_g_kkt[g,t]-μ_lb_g_kkt[g,t] ==0)

@constraint(C_B,[ d in D, t in T],-γ_D[d,t] - λ_i_kkt[network_data["load"][string(d)]["load_bus"],t] + μ_ub_d_kkt[d,t]-μ_lb_d_kkt[d,t] ==0)


@constraint(C_B,[cb in CB, t in T],γ_CB[cb,t] +sum(λ_i_kkt[i,t] for i in I if network_data["CB"][string(cb)]["CB_bus"]==i)+ μ_ub_cb_kkt[cb,t]-μ_lb_cb_kkt[cb,t]==0)

@constraint(C_B, [ i in I, t in T ], sum( B[i,j] * (λ_i_kkt[j,t]-λ_i_kkt[i,t]) + B[i,j] * (π_i_t_ub_kkt[i,j,t]-π_i_t_ub_kkt[j,i,t]) - B[i,j] * (π_i_t_lb_kkt[i,j,t]-π_i_t_lb_kkt[j,i,t]) for j in I)==0) 



@constraint(C_B, [g in G,t in T],(network_data["gen"][string(g)]["pg"]-P_g_kkt[g,t])<=z_g_o_kkt[g,t]*M_g)
@constraint(C_B, [g in G,t in T],μ_ub_g_kkt[g,t]<=(1-z_g_o_kkt[g,t])*M_g)
@constraint(C_B, [g in G,t in T],P_g_kkt[g,t]<=z_g_u_kkt[g,t]*M_g)
@constraint(C_B, [g in G,t in T],μ_lb_g_kkt[g,t]<=(1-z_g_u_kkt[g,t])*M_g)

@constraint(C_B, [d in D, t in T],(network_data["load"][string(d)]["pd"]-P_d_kkt[d,t])<=z_d_o_kkt[d,t]*M_d)
@constraint(C_B, [d in D, t in T],μ_ub_d_kkt[d,t]<=(1-z_d_o_kkt[d,t])*M_d)
@constraint(C_B, [d in D, t in T],P_d_kkt[d,t]<=z_d_u_kkt[d,t]*M_d)
@constraint(C_B, [d in D, t in T],μ_lb_d_kkt[d,t]<=(1-z_d_u_kkt[d,t])*M_d)




@constraint(C_B, [cb in CB, t in T],(network_data["CB"][string(cb)]["max_value"]-P_cb_kkt[cb,t])<=z_cb_o_kkt[cb,t]*M_d)
@constraint(C_B, [cb in CB, t in T],μ_ub_cb_kkt[cb,t]<=(1-z_cb_o_kkt[cb,t])*M_d)
@constraint(C_B, [cb in CB, t in T],P_cb_kkt[cb,t] - network_data["CB"][string(cb)]["min_value"]<=z_cb_u_kkt[cb,t]*M_d)
@constraint(C_B, [cb in CB, t in T],μ_lb_cb_kkt[cb,t]<=(1-z_cb_u_kkt[cb,t])*M_d)

@constraint(C_B, [ i in I, j in I, t in T ], P_max_line[i,j]-B[i,j]*(θ_kkt[i,t]-θ_kkt[j,t])<=z_line_o_kkt[i,j,t]*M_θ)
@constraint(C_B, [ i in I, j in I, t in T ],π_i_t_ub_kkt[i,j,t]<=(1-z_line_o_kkt[i,j,t])*M_θ)
@constraint(C_B, [ i in I, j in I, t in T ],P_max_line[i,j]+B[i,j]*(θ_kkt[i,t]-θ_kkt[j,t])<=z_line_u_kkt[i,j,t]*M_θ)
@constraint(C_B, [ i in I, j in I, t in T ],π_i_t_lb_kkt[i,j,t]<=(1-z_line_u_kkt[i,j,t])*M_θ)


@constraint(C_B, [ i in I , t in T], 

 +sum( P_g_kkt[g,t] + P_g_up_kkt[g,t]- P_g_down_kkt[g,t] for g in G if network_data["gen"][string(g)]["gen_bus"]==i)

 -sum( P_d_kkt[d,t]  - Δ_d[d,t] for d in D if network_data["load"][string(d)]["load_bus"]==i ) 
  

 -sum(B[i, j] * (θ_kkt_rtm[i, t] - θ_kkt_rtm[j, t]) for j in I if j != i) 
  ==0) 

@constraint(C_B,[g in G,t in T],P_g_up_kkt[g,t]<=network_data["gen"][string(g)]["pg"] - P_g_kkt[g,t] )    
@constraint(C_B,[g in G,t in T],P_g_down_kkt[g,t]<=P_g_kkt[g,t] )  


@constraint(C_B,[i in I,j in I,t in T],B[i,j]*(θ_kkt_rtm[i,t]-θ_kkt_rtm[j,t])<=P_max_line[i,j] )   
@constraint(C_B,[i in I,j in I,t in T],-P_max_line[i,j]-B[i,j]*(θ_kkt_rtm[i,t]-θ_kkt_rtm[j,t])<=0 ) 

#Dual Feasibility

@constraint(C_B,[g in G,t in T],γ_G_up[t] +λ_rtm_kkt[network_data["gen"][string(g)]["gen_bus"],t] + μ_G_ub_up_kkt[g,t]-μ_G_lb_up_kkt[g,t]==0)
@constraint(C_B,[g in G,t in T],-γ_G_down[t] -λ_rtm_kkt[network_data["gen"][string(g)]["gen_bus"],t] + μ_G_ub_down_kkt[g,t]-μ_G_lb_down_kkt[g,t]==0)



@constraint(C_B, [ i in I, t in T ], sum(B[i,j]*(λ_rtm_kkt[j,t]-λ_rtm_kkt[i,t])+B[i,j]*(π_rtm_ub_kkt[i,j,t]-π_rtm_ub_kkt[j,i,t])-B[i,j]*(π_rtm_lb_kkt[i,j,t]-π_rtm_lb_kkt[j,i,t]) for j in I)==0)

#Complementary Slackness


@constraint(C_B, [g in G,t in T],network_data["gen"][string(g)]["pg"] - P_g_kkt[g,t]  -P_g_up_kkt[g,t] <=z_G_ub_up_kkt[g,t]*M_g)
@constraint(C_B, [g in G,t in T],μ_G_ub_up_kkt[g,t]<=(1-z_G_ub_up_kkt[g,t])*M_g)
@constraint(C_B, [g in G,t in T],P_g_up_kkt[g,t]<=z_G_lb_up_kkt[g,t]*M_g)
@constraint(C_B, [g in G,t in T],μ_G_lb_up_kkt[g,t]<=(1-z_G_lb_up_kkt[g,t])*M_g)


@constraint(C_B, [g in G,t in T],P_g_kkt[g,t] -P_g_down_kkt[g,t]<=z_G_ub_down_kkt[g,t]*M_g)
@constraint(C_B, [g in G,t in T],μ_G_ub_down_kkt[g,t]<=(1-z_G_ub_down_kkt[g,t])*M_g)
@constraint(C_B, [g in G,t in T],P_g_down_kkt[g,t]<=z_G_lb_down_kkt[g,t]*M_g)
@constraint(C_B, [g in G,t in T],μ_G_lb_down_kkt[g,t]<=(1-z_G_lb_down_kkt[g,t])*M_g)



@constraint(C_B, [ i in I, j in I, t in T ], P_max_line[i,j]-B[i,j]*(θ_kkt_rtm[i,t]-θ_kkt_rtm[j,t])<=z_rtm_ub_kkt[i,j,t]*M_θ)
@constraint(C_B, [ i in I, j in I, t in T ],π_rtm_ub_kkt[i,j,t]<=(1-z_rtm_ub_kkt[i,j,t])*M_θ)
@constraint(C_B, [ i in I, j in I, t in T ],P_max_line[i,j]+B[i,j]*(θ_kkt_rtm[i,t]-θ_kkt_rtm[j,t])<=z_rtm_lb_kkt[i,j,t]*M_θ)
@constraint(C_B, [ i in I, j in I, t in T ],π_rtm_lb_kkt[i,j,t]<=(1-z_rtm_lb_kkt[i,j,t])*M_θ)


JuMP.optimize!(C_B)
status = termination_status(C_B)
 

obj_value_KKT_primal= sum(
    + sum(γ_CB[ cb, t]* value(P_cb_kkt[ cb, t])  for cb in CB)
    + sum(γ_G[ g, t] * value(P_g_kkt[ g, t])  for g in G)
    - sum(γ_D[ d, t] * value(P_d_kkt[ d, t])   for d in D)

    for t in T
)  


obj_value_KKT_dual = sum(
    -sum(network_data["CB"][string(cb)]["max_value"]* value(μ_ub_cb_kkt[cb,t]) for cb in CB)
    +sum(network_data["CB"][string(cb)]["min_value"]* value(μ_lb_cb_kkt[cb,t]) for cb in CB)
    -sum(network_data["gen"][string(g)]["pg"]* value(μ_ub_g_kkt[g,t])  for g in G)
    -sum(network_data["load"][string(d)]["pd"]* value(μ_ub_d_kkt[d,t])  for d in D)

    -sum(sum(P_max_line[i,j]*( value(π_i_t_ub_kkt[i,j,t]) + value(π_i_t_lb_kkt[i,j,t])) for j in I) for i in I)

    for t in T) 


Objective_value_kkt_primal  = sum(  


+sum(γ_G_up[t]* value(P_g_up_kkt[g,t])-γ_G_down[t]*value(P_g_down_kkt[g,t] ) for g in G)

    for t in  T
) 
        
Objective_value_kkt_dual = sum(  

sum(value(λ_rtm_kkt[network_data["gen"][string(g)]["gen_bus"],t])*value(P_g_kkt[g,t]) for g in G)
-sum(value(λ_rtm_kkt[network_data["load"][string(d)]["load_bus"],t])*value(P_d_kkt[d,t]) for d in D)
+sum(value(λ_rtm_kkt[network_data["load"][string(d)]["load_bus"],t])*value(Δ_d[d,t]) for d in D)

-sum( (network_data["gen"][string(g)]["pg"]- value(P_g_kkt[g,t]))*value(μ_G_ub_up_kkt[g,t]) for g in G)
-sum( value(P_g_kkt[g,t])*value(μ_G_ub_down_kkt[g,t]) for g in G)


-sum(sum(P_max_line[i,j]*(value(π_rtm_ub_kkt[i,j,t])+value(π_rtm_lb_kkt[i,j,t])) for j in I) for i in I)

    for t in  T
)



sum_biinear = sum(value(P_cb_kkt[cb,t])*value(λ_rtm_kkt[cb,t]) for cb in CB for t in T)
println("Bilinear term =$sum_biinear ")



rtm_price=Array(value.(λ_rtm_kkt)')
dam_price= Array(value.(λ_i_kkt)')


dam_price = DataFrame(Array(value.(-λ_i_kkt)')*100, :auto) 

rtm_price = DataFrame(Array(value.(-λ_rtm_kkt)')*100, :auto) 

cb = DataFrame(Array(value.(-P_cb_kkt)')*100, :auto) 

# XLSX.openxlsx("ECE285/CB.xlsx", mode="w") do xf
#     # Create new sheets and assign them
#     # sheet1 = XLSX.addsheet!(xf, "dam_price")
#     # sheet2 = XLSX.addsheet!(xf, "rtm_price")
#     sheet1 = XLSX.addsheet!(xf, "CB")
    

#     XLSX.writetable!(sheet1, round.(cb, digits=3) ) 
#     # XLSX.writetable!(sheet2,round.(rtm_price, digits=3))

# end  
network_data["load"]