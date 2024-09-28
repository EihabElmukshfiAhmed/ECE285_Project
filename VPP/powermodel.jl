using PowerModels,Ipopt, DataFrames, XLSX
using Gurobi, JuMP, DataFrames, ExcelFiles 

# Load the network data (replace the path with your case file)
case_file = "C:/Users/eihab/OneDrive/Desktop/SDSU/Julia/VPP/case6ww.m"
network_data = PowerModels.parse_file(case_file)


function create_loads_dictionary(load_type::String, num_loads::Int, buses::Vector{Int}, min_values::Vector{Float64}, max_values::Vector{Float64})
    # Initialize an empty dictionary to store the loads/EVs/etc
    loads_dict = Dict{String, Dict{String, Any}}()

    # Check if the input vectors have the correct size
    if length(buses) != num_loads || length(min_values) != num_loads || length(max_values) != num_loads
        throw(ArgumentError("The length of buses, min_values, and max_values must match num_loads."))
    end

    # Populate the dictionary with load information
    for i in 1:num_loads
        load_info = Dict(
            "load_type" => load_type,
            load_type*"_"*"bus" => buses[i],
            "min_value" => min_values[i],
            "max_value" => max_values[i]
        )
        loads_dict[string(i)] = load_info
    end

    return loads_dict
end




# ---------------------------------Add EV load-----------------------------------
load_type = "EV"
num_loads = 3
buses = [1, 2, 3]
min_values = [10.0, 15.0, 5.0]
max_values = [50.0, 70.0, 20.0]

EV_dictionary = create_loads_dictionary(load_type, num_loads, buses, min_values, max_values)

network_data["EV"]=EV_dictionary

# -----------------------------------Add DV load---------------------------------
load_type = "DC"
num_loads = 3
buses = [1, 2, 6]
min_values = [0.0,0.0, 0.0]
max_values = [50.0, 70.0, 20.0]

DC_dictionary = create_loads_dictionary(load_type, num_loads, buses, min_values, max_values)

network_data["DC"]=DC_dictionary 

#----------------------------------------------

num_buses = length(network_data["bus"]) 
num_gens = length(network_data["gen"]) 
num_loads=length(network_data["load"])
num_EVs=length(network_data["EV"])
num_DCs=length(network_data["DC"])









gen_data = network_data["gen"]

load_data = network_data["load"]

DC_data = network_data["DC"]

EV_data = network_data["EV"]

# Initialize a 6x6 incidence matrix with zeros









function get_indicident_matrix(component::String)
    if haskey(network_data,component)==false
        throw(ArgumentError("component Does not Exsit!"))
    end
    if component=="branch"
        incidence_matrix_line = zeros(Int, num_buses, num_buses)
        capacity_matrix = zeros(Int, num_buses, num_buses)
        # Populate the incidence matrix
        branch_data = network_data["branch"]
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
        comp_data = network_data[component]
        component_length=length(comp_data)
        incidence_matrix = zeros(Int,num_buses,component_length)
        

        for (i, x) in enumerate(values(comp_data))
    
            bus=x[component*"_"*"bus"]
            println("------------------BUS-----------------------------")
            println(bus)
        
            incidence_matrix[bus,i]=1
        
        end 
        
    end


    return incidence_matrix
end


c=get_indicident_matrix("DC")

incidence_matrix_lune_capacity = get_indicident_matrix("brach")

incidence_matrix_gen = get_indicident_matrix("gen")
incidence_matrix_load = get_indicident_matrix("load")
incidence_matrix_EV = get_indicident_matrix("EV")
incidence_matrix_DC = get_indicident_matrix("DC")




# Print the incidence matrix
# println("Incidence Matrix:")

# incidence_matrix_transpose = incidence_matrix'

# file_path = "incidence_matrix_transpose.xlsx"

# df = DataFrame(incidence_matrix_transpose) 
# # Save the transposed matrix to Excel
# XLSX.writetable(file_path,  incidence_matrix_transpose) 


# data1 = [[1, 2, 3], [4, 5, 6]];
# data2 = [[1, 1, 1], [2, 2, 2]];    
# data3 = [[0, 1,-1], [-2,0, 2]];

# # write data1 to SHEET_A and data2 to SHEET_B
# XLSX.writetable("TEST12.xlsx", SHEET_A=(incidence_matrix_transpose),
#         overwrite=true)

# for i in keys(branch_data["1"])
#     println(i)
# end
