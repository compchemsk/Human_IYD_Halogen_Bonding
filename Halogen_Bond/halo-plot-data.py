import re

# Define file paths
input_file_path = "halo.out"  # Replace with the actual path if needed
output_file_path = "d-theta-pathway.csv"

# Initialize variables
trajectory_number = None
extracted_data = []

# Regex patterns
trajectory_pattern = re.compile(r"For trajectory (\d+\.\d+)")
distance_angle_pattern = re.compile(r"Distance:\s+([\d.]+)\s+nm.*Angle:\s+([\d.]+)°")

# Read and process file
with open(input_file_path, "r") as file:
    for line in file:
        traj_match = trajectory_pattern.search(line)
        data_match = distance_angle_pattern.search(line)

        if traj_match:
            # Update trajectory number
            trajectory_number = float(traj_match.group(1))

        if data_match and trajectory_number is not None:
            # Extract distance and angle values
            distance = float(data_match.group(1)) * 10  # Convert nm to Å (multiplying by 10)
            angle = float(data_match.group(2))

            # Store the extracted data
            extracted_data.append((trajectory_number, distance, angle))

# Write data to CSV file
with open(output_file_path, "w") as file:
    file.write("Trajectory_Number,Distance_Angstrom,Angle_Degrees\n")  # Header
    for row in extracted_data:
        file.write(f"{row[0]},{row[1]},{row[2]}\n")

print(f"Data saved to {output_file_path}")
