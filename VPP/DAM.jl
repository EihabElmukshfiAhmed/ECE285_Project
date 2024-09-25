using Gurobi, JuMP, DataFrames, ExcelFiles
#test

#--------------------------------------------------------------define sets--------------------------------------------------------------#

bus=1:10
Gen=1:10
Demand=1:10
EV=1:10
T=1:24
 D_C=1:10
Demand_Seg=1:3
Gen_Seg=1:3
EV_Seg=1:3
 D_C_Seg=1:3


#--------------------------------------------------------------define parameters--------------------------------------------------------------#

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

#Dual variables------------------

#Inequality dual variables

@variable(m,μ_CB_Upper[ bus, T])
@variable(m,μ_CB_Lower[ bus, T])

@variable(m,μ_Gen_Seg_Upper[ Gen_Seg, Gen, T])
@variable(m,μ_Gen_Seg_Lower[ Gen_Seg, Gen, T])

@variable(m,μ_Demand_Seg_Upper[ Demand_Seg, Demand, T])
@variable(m,μ_Demand_Seg_Lower[ Demand_Seg, Demand, T])

@variable(m,μ_EV_Seg_Upper[ EV_Seg, EV, T])
@variable(m,μ_EV_Seg_Lower[ EV_Seg, EV, T])

@variable(m,μ_D_C_Seg_Upper[  D_C_Seg,  D_C, T])
@variable(m,μ_D_C_Seg_Lower[  D_C_Seg,  D_C, T])

@variable(m,π_Upper_DA[ bus, bus, T])
@variable(m,π_Lower_DA[ bus, bus, T])

#Equality dual variables

@variable(m,λ_i_DA[ bus, T])
@variable(m,λ_g_DA[ Gen, T])
@variable(m,λ_d_DA[ Demand, T])
@variable(m,λ_e_DA[ EV, T])
@variable(m,λ_c_DA[  D_C, T])

#----------------------------------Constraints-------------------------------------------------------------

@constraint(m,sum( 

    sum((μ_CB_Upper[i,t]*P_max_CB[i].- μ_CB_Lower*P_min_CB[i]).-(γ_CB[i,t]*P_CB[i,t]) .+ sum(P_max_line[i,j]*(π_Upper_DA[i,j,t] .+π_Lower_DA[i,j,t]) for j in bus) for i in bus)  

    +sum(  sum( (μ_Gen_Seg_Upper[w,g,t]*P_max_Genco_Seg[w,g,t]).-(γ_Gen[w,g,t]*P_g_w[w,g,t]) for w in Gen_Seg) for g in Gen )

    +sum( sum( (μ_Demand_Seg_Upper[y,d,t]*P_max_Demand_Seg[y,d,t]).+(γ_Demand[y,d,t]*P_d_y[y,d,t]) for y in Demand_Seg) for d in Demand)

    +sum( sum( (μ_EV_Seg_Upper[z,e,t]*P_max_EV_Seg[z,e,t]) .+ (γ_EV[z,e,t]*P_e_z[z,e,t]) for z in EV_Seg) for e in EV)
    
    +sum( sum( (μ_D_C_Seg_Upper[n,c,t]*P_max_D_C_Seg[n,c,t]) .+ (γ_D_C[n,c,t]*P_c_n[n,c,t]) for z in  D_C_Seg) for c in  D_C)
    

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

@constraint(m,[w in Gen_Seg,g in Gen, t in T], γ_Gen[w,g,t]-λ_g_DA[g,t]+μ_Gen_Seg_Upper[w,g,t]-μ_Gen_Seg_Lower[w,g,t]==0)

@constraint(m,[y in Demand_Seg,d in Demand, t in T], -γ_Demand[y,d,t]-λ_d_DA[d,t]+μ_Demand_Seg_Upper[y,d,t]-μ_Demand_Seg_Lower[y,d,t]==0)

@constraint(m,[z in EV_Seg,e in EV, t in T], -γ_EV[z,e,t]-λ_d_DA[e,t]+μ_EV_Seg_Upper[z,e,t]-μ_EV_Seg_Lower[z,e,t]==0)

@constraint(m,[n in D_c_Seg,c in EV, t in T], -γ_D_c[n,c,t]-λ_c_DA[c,t]+μ_D_c_Seg_Upper[n,c,t]-μ_D_C_Seg_Lower[n,c,t]==0)

@constraint(m,[i in bus, t in T], γ_CB[i,t]+λ_i_DA[i,t]+μ_CB_Upper[i,t]-μ_CB_Lower[i,t]==0)

@constraint(m,[ g in Gen , i in bus, t in T ], A_g[g,i]*λ_i_DA[i,t]+λ_g_DA[g,t]==0) 

@constraint(m,[ d in Demand , i in bus, t in T ], -A_d[d,i]*λ_i_DA[i,t]+λ_d_DA[d,t]==0) 

@constraint(m,[ c in D_C , i in bus, t in T ], -A_c[c,i]*λ_i_DA[i,t]+λ_c_DA[c,t]==0) 

@constraint(m,[ e in EV , i in bus, t in T ], -A_e[e,i]*λ_i_DA[i,t]+λ_e_DA[e,t]==0)

@constraint(m,[ i in bus,  t in T ], sum(B[i,j]*(λ_i_DA[j,t]-λ_i_DA[i,t])+B[i,j]*(π_Upper_DA[i,j,t]-π_Upper_DA[j,i,t])+B[i,j]*(π_Lower_DA[j,i,t]-π_Lower_DA[i,j,t]) for j in bus)==0)
#
