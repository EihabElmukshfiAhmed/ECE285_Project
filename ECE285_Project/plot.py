import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data from Excel files
file_without_vb = "With_out_CB.xlsx"
file_with_vb = "With_CB.xlsx"
CB = "CB.xlsx"

# Read data from Excel sheets
prices_da_without_vb = pd.read_excel(file_without_vb, sheet_name="dam_price")
prices_rt_without_vb = pd.read_excel(file_without_vb, sheet_name="rtm_price")

prices_da_with_vb = pd.read_excel(file_with_vb, sheet_name="dam_price")
prices_rt_with_vb = pd.read_excel(file_with_vb, sheet_name="rtm_price")

P_CB = pd.read_excel(CB, sheet_name="CB")
# Extract time (assuming first column is time)
CB_Profit =(( prices_da_with_vb - prices_rt_with_vb)*P_CB).sum()



# Define bus numbers (1 to 14)
# bus_numbers = np.arange(1, 7)

# # Example profit data for each bus (some positive, some negative)
# profits_black = np.random.randint(-500, 500, size=14)
# profits_white = profits_black + np.random.randint(-50, 50, size=14)  # Slight variation

# # Create bar chart
# fig, ax = plt.subplots(figsize=(6, 2))

# ax.bar(bus_numbers - 0.2, CB_Profit, width=0.4, color='black', label="Black Bars")
# # ax.bar(bus_numbers + 0.2, profits_white, width=0.4, edgecolor='black', facecolor='white', label="White Bars")

# # Labels and titles
# ax.set_xlabel("Bus Number")
# ax.set_ylabel("Profit from CB [$]")
# ax.axhline(0, color='black', linewidth=1)  # Zero reference line

# # Adjust tick labels
# ax.set_xticks(bus_numbers)

# # Remove top and right borders
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.set_xlabel("Bus Number", fontsize=20)  # Increase label font size
# ax.set_ylabel("Profit from CB [$]", fontsize=20)  

# # Increase tick label font size
# ax.tick_params(axis='both', labelsize=20)  

# # Adjust x-axis ticks
# ax.set_xticks(bus_numbers)

# plt.tight_layout()
# plt.show()

time = range(1, 25)  # 1 to 24

# Plot settings
fig, axes = plt.subplots(2, 1, figsize=(10, 20))  # Slightly larger figure size

# Plot Without Virtual Bidding
axes[0].set_title("Without Convergence Bidding", fontsize=20)
axes[0].set_xlabel("Time", fontsize=20)
axes[0].set_ylabel("Price ($/MWh)", fontsize=20)
axes[1].set_ylabel("Price ($/MWh)", fontsize=20)
# Plot With Virtual Bidding
axes[1].set_title("With Convergence Bidding", fontsize=20)
axes[1].set_xlabel("Time", fontsize=20)

# Assuming we only plot Bus 2 data
bus_index = 1  # Adjust based on column indexing

axes[0].plot(time, prices_da_without_vb.iloc[:, bus_index], label="DA LMP at Node 2", linestyle='solid')
axes[0].plot(time, prices_rt_without_vb.iloc[:, bus_index], label="RT LMP at Node 2", linestyle='solid')

axes[1].plot(time, prices_da_with_vb.iloc[:, bus_index],  linestyle='solid')
axes[1].plot(time, prices_rt_with_vb.iloc[:, bus_index],  linestyle='solid')

# Increase legend font size
axes[0].legend(loc="upper left", fontsize=20)
# axes[1].legend(loc="upper right", fontsize=20)

# Increase tick label font size
for ax in axes:
    ax.tick_params(axis='both', labelsize=20)

# Save and show plot
plt.tight_layout()
plt.savefig("ECE285/CB.png", dpi=300)  # Increase resolution
plt.show()
