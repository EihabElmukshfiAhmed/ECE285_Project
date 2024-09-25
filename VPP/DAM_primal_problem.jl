using PowerModels,Ipopt, DataFrames
using Gurobi, JuMP, DataFrames, ExcelFiles 


case_file = "C:/Users/eihab/OneDrive/Desktop/SDSU/Julia/VPP/case6ww.m"
network_data = PowerModels.parse_file(case_file)

branch_data = network_data["branch"]

# Create a B matrix based on the branch data
n_buses = length(network_data["bus"])  # Number of buses
B_matrix = zeros(n_buses, n_buses)      # Initialize B matrix
 
B = PowerModels.calc_susceptance_matrix(network_data).matrix

gen_data = network_data["gen"] 




# B=DataFrame(load("susceptance_matrix.xlsx")) # susceptance_matrix
B=ones(length(bus),length(bus))
P_max_line=ones(length(bus),length(bus))   # max power flow in the lines
P_max_CB=ones(length(bus))  # max CB at each bus
P_min_CB=ones(length(bus))  # min CB at each bus

#define segment limits

P_max_Genco_Seg=ones(length(Gen_Seg),length(Gen))  #max Genco_Seg
P_max_Demand_Seg=ones(length(Demand_Seg),length(Demand))  #max Demand_Seg
P_max_EV_Seg=ones(length(EV_Seg),length(EV))  #max EV_Seg
P_max_D_C_Seg=ones(length( D_C_Seg),length( D_C))  #max  D_C_Seg

# define bidding price

γ_CB=ones(length(bus),length(T))
γ_Gen=ones(length(Gen_Seg),length(Gen),length(T))
γ_Demand=ones(length(Demand_Seg),length(Demand),length(T))
γ_EV=ones(length(EV_Seg),length(EV),length(T))
γ_D_C=ones(length( D_C_Seg),length( D_C),length(T))

#Define Incidence  matrix:
A_g = ones(length(bus),length(Gen)) #Gencos Incidence matrix
A_cb = ones(length(bus),length(bus)) #CB Incidence matrix
A_d = ones(length(bus),length(Demand)) #Demand Incidence matrix
A_c = ones(length(bus),length( D_C)) #Data_center Incidence matrix
A_e = ones(length(bus),length(EV)) #EV Incidence matrix


m=Model(Gurobi.Optimizer)

#Power variables

@variable(m,P_d[Demand,T]>=0) #demand
@variable(m,P_e[EV,T]>=0)     #EV
@variable(m,P_c[ D_C,T]>=0)  #Data_center 
@variable(m,P_g[Gen,T]>=0)  #Gencos
@variable(m,P_CB[bus,T]>=0) #Convergence bidding

@variable(m,θ_DA[bus,T])
#segment variables

@variable(m,P_g_w[Gen_Seg,Gen,T] >=0)
@variable(m,P_d_y[Demand_Seg,Demand,T]   >=0)
@variable(m,P_e_z[EV_Seg,EV,T]  >=0)
@variable(m,P_c_n[ D_C_Seg, D_C,T] >=0)




#----------------------------------Constraints-------------------------------------------------------------

@constraint(m,sum( 

    sum((γ_CB[i,t]*P_CB[i,t] for i in bus)  

    +sum(  sum( γ_Gen[w,g,t]*P_g_w[w,g,t] for w in Gen_Seg) for g in Gen )

    +sum( sum( γ_Demand[y,d,t]*P_d_y[y,d,t] for y in Demand_Seg) for d in Demand)

    +sum( γ_EV[z,e,t]*P_e_z[z,e,t] for z in EV_Seg) for e in EV)
    
    +sum( sum( γ_D_C[n,c,t]*P_c_n[n,c,t] for z in  D_C_Seg) for c in  D_C)
    

    for t in  T
)==0)

@constraint(m, sum(A_g[i,g]*P_g[g,t] for g in G) + P_CB[i,t] .- sum(A_d[i,d]*P_d[d,t] for d in Demand) .- sum(A_c[i,c]*P_c[c,t] for c in  D_C) .- sum(A_e[i,e]*P_e[e,t] for e in EV)
  + sum(B[i,j]*(θ_DA[i,t]-θ_DA[j,t])) == 0 for i in bus for t in T)





@constraint(m,[g in Gen, t in T],P_g[g,t].-sum(P_g_w[w,g,t] for w in Gen_Seg)==0 ) 

@constraint(m,[d in Demand, t in T],P_d[d,t].-sum(P_d_y[y,d,t] for y in Demand_Seg)==0)

@constraint(m,[e in EV , t in T],P_e[e,t].-sum(P_e_z[z,e,t] for z in EV_Seg)==0 )

@constraint(m, [c in D_C, t in T], P_c[c,t]-sum(P_c_n[n,c,t] for n in  D_C_Seg)==0 )  

@constraint(m,[i in bus,t in T],P_min_CB[i]<=P_CB[i,t]<=P_max_CB[i])       
    
@constraint(m,[w in Gen_Seg,g in Gen,t in T],P_g_w[w,g,t]<=P_max_Genco_Seg[w,g])  

@constraint(m,[y in Demand_Seg,d in Demand,t in T],P_d_y[y,d,t]<=P_max_Demand_Seg[y,d])  

@constraint(m,[z in EV_Seg,e in EV,t in T],P_e_z[z,e,t]<=P_max_EV_Seg[z,e]) 

@constraint(m,[n in D_c_Seg,c in D_C,t in T],P_c_n[n,c,t]<=P_max_D_C_Seg[z,e]) 

@constraint(m,[i in bus,j in bus,t in T],-P_max_line[i,j]<=B[i,j]*(θ_DA[i,t]-θ_DA[j,t])<=P_max_line[i,j] )   
