using PowerModels,Ipopt, DataFrames, XLSX
using Gurobi, JuMP, DataFrames

 

df_price = DataFrame(XLSX.readtable("prices.xlsx", "Sheet1")[1],XLSX.readtable("prices.xlsx", "Sheet1")[2]) 

price_matrix = Matrix(df_price[:,1:9])'

# Reshape the matrix to size (24, 6, 6)
price_3d = reshape(Matrix(df_price[:,1:9])', 3, 3, 24)






case_file = "case6ww.m"
network_data = PowerModels.parse_file(case_file) 


# define bidding price

γ_G= reshape(Matrix(df_price[:,1:9])', 3, 3, 24) /network_data["baseMVA"]
γ_D= reshape(Matrix(df_price[:,10:18])', 3, 3, 24) /network_data["baseMVA"]
γ_EV= reshape(Matrix(df_price[:,19:27])', 3, 3, 24)/network_data["baseMVA"]
γ_DC= reshape(Matrix(df_price[:,28:36])', 3, 3, 24)/network_data["baseMVA"]
γ_CB= reshape(Matrix(df_price[:,37:42])', 6, 24) /network_data["baseMVA"]


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

        # Add segment information if segments are provided
        if segments !== "None"
            if haskey(segments, bus_id)
                segment_info = segments[bus_id]  # Retrieve the segment info for this bus
                load_info["segment"] = segment_info
            else
                throw(ArgumentError("Bus $bus_id does not exist in the segments dictionary."))
            end
        end

        loads_dict[string(i)] = load_info
    end

    return loads_dict
end



  function create_segment(component::String, buses::Vector{Int}, max_values::Vector{T}) where T
    if length(buses) != length(max_values)
        throw(ArgumentError("The length of buses and max_values must match."))
    end

    bus_dict = Dict{Int, Dict{String, Dict{String, Any}}}()  # Outer dictionary for buses with integer keys

    for (bus, max_value) in zip(buses, max_values)
        seg_dict = Dict{String, Dict{String, Any}}()
        segments = [max_value / 3,  max_value / 3, max_value/3]  # Divide into 3 segments

        for (i, seg_value) in enumerate(segments)
            seg_info = Dict(
                "component" => component,
                "bus" => bus,
                "max_value" => seg_value/network_data["baseMVA"]
            )
            seg_dict[string(i)] = seg_info
        end
        bus_dict[bus] = seg_dict  # Assign segments to the current bus
    end

    return bus_dict
end










# ---------------------------------Add EV load-----------------------------------
load_type = "EV"
num_loads = 3
buses = [1, 2, 3]
min_values = [0.0, 0.0, 0.0]
max_values = [10.0, 5.0, 4.0]
EV_segment=create_segment("EV",buses,max_values) 
EV_dictionary = create_loads_dictionary(load_type, num_loads, buses, min_values, max_values,EV_segment) 

network_data["EV"]=EV_dictionary

# -----------------------------------Add DC load---------------------------------
load_type = "DC"
num_loads = 3
buses = [1, 2, 6]
min_values = [0.0,0.0, 0.0]
max_values = [4.0, 5.0, 2.0]
DC_segment=create_segment("DC",buses,max_values)
DC_dictionary = create_loads_dictionary(load_type, num_loads, buses, min_values, max_values,DC_segment)

network_data["DC"]=DC_dictionary  




#----------------------------------------------ADD CB-------------------------
max_CB=15
load_type = "CB"
num_loads = length(network_data["bus"])
buses = [1, 2, 3,4,5,6]
min_values = ones(num_buses)*-max_CB
max_values = ones(num_buses)*max_CB

CB_dictionary = create_loads_dictionary(load_type, num_loads, buses, min_values, max_values,"None")

network_data["CB"]=CB_dictionary 







num_gens = length(network_data["gen"]) 
num_loads=length(network_data["load"])
num_EVs=length(network_data["EV"])
num_DCs=length(network_data["DC"])
num_CBs=length(network_data["CB"])
num_Brs=length(network_data["branch"])
T=1:24
D=1:num_loads
I=1:num_buses
G=1:num_gens
EV=1:num_EVs
DC=1:num_DCs
CB=1:num_CBs
Br=1:num_Brs






G_Seg = 1:3
D_Seg = 1:3
DC_Seg = 1:3
EV_Seg = 1:3
CB_Seg = 1:3


for i in G
    G_segment=create_segment("gen",[i],[network_data["baseMVA"]*network_data["gen"]["$i"]["pg"]]) 
    network_data["gen"][string(i)]["segment"]= G_segment[i]
end

for i in D
    D_segment=create_segment("load",[i],[network_data["baseMVA"]*network_data["load"]["$i"]["pd"]]) 
    network_data["load"][string(i)]["segment"]= D_segment[i]
  end





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
A_i_e = get_indicident_matrix("EV")
A_i_c = get_indicident_matrix("DC")
A_i_cb = get_indicident_matrix("CB")




P_max_line  = get_indicident_matrix("branch")  
B=PowerModels.calc_basic_susceptance_matrix(network_data)

