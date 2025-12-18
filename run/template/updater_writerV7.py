import os
import numpy as np
import re
import reader_functionsV4 as read
#from stl import mesh
##08 FEB 2025
###V4 updating vertices as well ++ infowriter
##09 FEB 2025
###written probe location
##11 FEB 2025
###V5 writeFV2 writeMV2  is added

###V6 moodify_velocity_BC & writeV3 added. (son terimi almiyordu oncekiler)
## modify_box_coordinates & modify_wave_dict & body_locater

###V7 2112 versions are added

def body_locater(input_file,output_file,x,y,z):
    # Load the mesh
    mesh_data = mesh.Mesh.from_file(input_file)

    # Define the translation vector
    translation_vector = np.array([x, y, z])  # Replace dx, dy, dz with your desired translation values

    # Apply the translation
    mesh_data.vectors += translation_vector

    # Save the updated STL file
    mesh_data.save(output_file)

    print(f"Translated STL saved to {output_file}")    

def modify_wave_dict(file_path, new_values, output_path=None):
    """
    Modifies the values in an OpenFOAM dictionary file.

    Parameters:
        file_path (str): Path to the original dictionary file.
        new_values (dict): Dictionary with new values to replace.
        output_path (str): Optional path to save the modified file. 
                           If None, it overwrites the original.
    """
    # Read the original file
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Prepare modified lines
    modified_lines = []
    pattern = re.compile(r'^(\s*)(\w+)(\s+)([-+]?[0-9]*\.?[0-9]+)(\s*;.*)$')

    for line in lines:
        match = pattern.match(line)
        if match:
            indent, key, spacing, old_value, suffix = match.groups()
            if key in new_values:
                new_value = new_values[key]
                line = f"{indent}{key}{spacing}{new_value}{suffix}\n"
        modified_lines.append(line)

    # Write to file
    target_path = output_path if output_path else file_path
    with open(target_path, 'w') as f:
        f.writelines(modified_lines)



def modify_box_coordinates(input_file: str, output_file: str, box_name: str, new_min: str, new_max: str):
    """
    Modify the min and max coordinates for a specified box in a snappyHexMeshDict file.

    Args:
        input_file (str): Path to the original snappyHexMeshDict file.
        output_file (str): Path to save the modified file.
        box_name (str): The name of the box (e.g., "box1x1x1").
        new_min (str): New minimum coordinates as a string, e.g., "(1.0 2.0 3.0)".
        new_max (str): New maximum coordinates as a string, e.g., "(4.0 5.0 6.0)".
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()

    modified_lines = []
    inside_box = False

    for line in lines:
        if re.search(rf'\b{re.escape(box_name)}\b', line):
            inside_box = True

        if inside_box and "min" in line:
            line = re.sub(r'min\s+\(.*?\);', f'min    {new_min};', line)
        if inside_box and "max" in line:
            line = re.sub(r'max\s+\(.*?\);', f'max    {new_max};', line)
            inside_box = False  # done with this block

        modified_lines.append(line)

    with open(output_file, 'w') as f:
        f.writelines(modified_lines)

    print(f"Updated '{box_name}' coordinates in '{output_file}'")


def modify_velocity_bc(file_path, internal_velocity, inlet_velocity):
    """
    Modify internalField and inlet boundary condition velocity in OpenFOAM 0/U file.

    Args:
        file_path (str): Path to the OpenFOAM velocity dictionary file (typically "0/U").
        internal_velocity (tuple): Desired internal velocity, e.g., (1, 0, 0).
        inlet_velocity (tuple): Desired inlet velocity, e.g., (1, 0, 0).
    """
    def format_vector(vec):
        return f"({vec[0]} {vec[1]} {vec[2]})"

    with open(file_path, 'r') as f:
        content = f.read()

    # Update internalField
    content = re.sub(
        r"(internalField\s+uniform\s+)\([^\)]+\);",
        rf"\1{format_vector(internal_velocity)};",
        content
    )

    # Update inlet value
    content = re.sub(
        r"(inlet\s*\{[^}]*?value\s+uniform\s+)\([^\)]+\);",
        lambda m: re.sub(
            r"\([^\)]+\)",
            format_vector(inlet_velocity),
            m.group(0)
        ),
        content,
        flags=re.DOTALL
    )

    # Save back
    with open(file_path, 'w') as f:
        f.write(content)

    print(f"Modified '{file_path}' with internalField={internal_velocity} and inlet={inlet_velocity}")


def read_probe_locations(tr,dx_value,dt_value):
    probe_locations = {}
    
    fid1 = "/postProcessing/patchProbes1/0/U"
    file_path = tr + fid1
    with open(file_path, 'r') as file:
        for line in file:
            match = re.match(r"# Probe (\d+) \(([-\d\.e ]+)\)", line)
            if match:
                probe_number = int(match.group(1))
                coordinates = tuple(map(float, match.group(2).split()))
                probe_locations[probe_number] = coordinates
    for probe, coordinates in probe_locations.items():
        print(f"Probe {probe}: X={coordinates[0]}, Y={coordinates[1]}, Z={coordinates[2]}")

    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/Probloc_Nx_{dx_str}_nt_{dt_str}.dat')

    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for probe, coordinates in probe_locations.items():
            # Write the time value in scientific notation
            file.write(f"{coordinates[0]:.6e} {coordinates[1]:.6e} {coordinates[2]:.6e}")
            file.write('\n')
        
        print(f"Probloc zeta successfully saved to {output_file_path}")


    return probe_locations




def simwriter(output_file_path, parameters_list,geolist):
    """
    Writes a dictionary of simulation parameters to a .dat file.

    Args:
        output_file_path (str): Path to the output file.
        parameters_list (dict): Dictionary containing simulation parameters.
    """
    with open(output_file_path, 'w') as file:
        for key, value in parameters_list.items():
            file.write(f"{key} = {value}\n")  # Writing in key = value format
        for key, value in geolist.items():
            file.write(f"{key} = {value}\n")  # Writing in key = value format
    print(f"Simulation parameters written to {output_file_path}")




def update_deltaT(file_path, new_deltaT):
    """
    Update the deltaT value in the controlDict file.

    Args:
    file_path (str): Path to the controlDict file.
    new_deltaT (float): New value for deltaT.

    """
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: File '{file_path}' not found.")
        return
    
    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Update the deltaT value
    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('deltaT'):
                file.write(f"deltaT          {new_deltaT};\n")
            else:
                file.write(line)
    
    print(f"Updated deltaT in '{file_path}' to {new_deltaT}.")


def update_control(file_path, new_deltaT, endtime,writecont,probeloc):
    """
    Update the deltaT value in the controlDict file.

    Args:
    file_path (str): Path to the controlDict file.
    new_deltaT (float): New value for deltaT.

    """
    pro="("
    ncase=len(probeloc)
    for j in range(0,ncase):
        pro+=f"({probeloc[j]} 0 0) "
    pro+=")"


    print("probe locations written:",pro)
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: File '{file_path}' not found.")
        return
    
    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Update the deltaT value
    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('deltaT'):
                file.write(f"deltaT          {new_deltaT};\n")
            elif line.strip().startswith('endTime'):
                file.write(f"endTime          {endtime};\n")
            elif line.strip().startswith('writeInterval'):
                file.write(f"writeInterval       {writecont};\n")
            elif line.strip().startswith('probeLocations'):
                file.write(f"        probeLocations  {pro};\n")
            else:
                file.write(line)
    
    print(f"Control Dict. Updated in '{file_path}'.")    


def update_control2(file_path, new_deltaT, endtime,writecont,probelocx,probelocy):
    """
    Update the deltaT value in the controlDict file.

    Args:
    file_path (str): Path to the controlDict file.
    new_deltaT (float): New value for deltaT.

    """
    pro="("
    ncase=len(probelocx)
    for j in range(0,ncase):
        pro+=f"({probelocx[j]} {probelocy[j]} 0) "
    pro+=")"


    #print("probe locations written:",pro)
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: File '{file_path}' not found.")
        return
    
    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Update the deltaT value
    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('deltaT'):
                file.write(f"deltaT          {new_deltaT};\n")
            elif line.strip().startswith('endTime'):
                file.write(f"endTime          {endtime};\n")
            elif line.strip().startswith('writeInterval'):
                file.write(f"writeInterval       {writecont};\n")
            elif line.strip().startswith('probeLocations'):
                file.write(f"        probeLocations  {pro};\n")
            else:
                file.write(line)
    
    print(f"Control Dict. Updated in '{file_path}'.")    

def update_hex_block(file_path, new_hex_values):
    """
    Update the hex block values in the blockMeshDict file.
    
    new_hex_values = f"(0 1 2 3 4 5 6 7) (1000 1 {number[j]}) simpleGrading (1 1 1)"
    update_hex_block(file_path, new_hex_values)
    Args:
    
    file_path (str): Path to the blockMeshDict file.
    new_hex_values (str): New values for the hex block in the format:
                          "(0 1 2 3 4 5 6 7) (4000 1 150) simpleGrading (1 1 1)"
    """
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: File '{file_path}' not found.")
        return
    
    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Update the hex block values
    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('hex'):
                # Replace the entire hex line with the new values
                file.write(f"    hex {new_hex_values}\n")
            else:
                file.write(line)
    
    print(f"Updated hex block in '{file_path}' to {new_hex_values}.")

def update_vertices(file_path, new_vertices):
    """
    Replace the vertices section in the blockMeshDict file.

    Args:
        file_path (str): Path to the blockMeshDict file.
        new_vertices (list of tuples): New vertex coordinates in the format:
                                       [(x1, y1, z1), (x2, y2, z2), ...]
    """
    if not os.path.isfile(file_path):
        print(f"Error: File '{file_path}' not found.")
        return
    
    # Convert new vertices to the required format
    new_vertices_str = "vertices (" + " ".join(f"({x} {y} {z})" for x, y, z in new_vertices) + ");"
    
    # Read the file contents
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Replace the vertices section completely
    updated_content = re.sub(r"vertices\s*\([^;]*\);", new_vertices_str, content, count=1)

    # Write back to file
    with open(file_path, 'w') as file:
        file.write(updated_content)
    
    print(f"Successfully replaced vertices in '{file_path}'.")



def update_block_Mesh(file_path, new_hex_values,domainlen):
    """
    Update the hex block values in the blockMeshDict file.
    
    new_hex_values = f"(0 1 2 3 4 5 6 7) (1000 1 {number[j]}) simpleGrading (1 1 1)"
    update_hex_block(file_path, new_hex_values)
    Args:
    
    file_path (str): Path to the blockMeshDict file.
    new_hex_values (str): New values for the hex block in the format:
                          "(0 1 2 3 4 5 6 7) (4000 1 150) simpleGrading (1 1 1)"
    """
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: File '{file_path}' not found.")
        return
    
    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Update the hex block values
    with open(file_path, 'w') as file:
        for line in lines:
            if line.strip().startswith('hex'):
                # Replace the entire hex line with the new values
                file.write(f"    hex {new_hex_values}\n")
            else:
                file.write(line)
    
    print(f"Updated hex block in '{file_path}' to {new_hex_values}.")

def writeP(tr,dt_value):
    fid1 = "/postProcessing/patchProbes1/0/p"
    input_file_path = tr + fid1
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        numberofprobes = sum(1 for line in data_lines if line.startswith('#') and 'Probe' in line)
        data_lines = [line for line in data_lines if not line.startswith('#')]

    print(f"Number of probes: {numberofprobes}")

    datainprobes = np.zeros((numberofprobes, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            coordinates = [tuple(map(float, element.strip('()').split())) for element in elements[1:]]
            flattened_list = [value for coord in coordinates for value in coord]
            datainprobes[:, i] = flattened_list 

    dt_str = f"{dt_value:d}"
    output_file_path = os.path.join(tr, f'yedek/probes_P_nz_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Data P successfully saved to {output_file_path}")

def writeP_gh(tr,dt_value):
    fid1 = "/postProcessing/patchProbes1/0/p_gh"
    input_file_path = tr + fid1
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        numberofprobes = sum(1 for line in data_lines if line.startswith('#') and 'Probe' in line)
        data_lines = [line for line in data_lines if not line.startswith('#')]

    print(f"Number of probes: {numberofprobes}")

    datainprobes = np.zeros((numberofprobes, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            coordinates = [tuple(map(float, element.strip('()').split())) for element in elements[1:]]
            flattened_list = [value for coord in coordinates for value in coord]
            datainprobes[:, i] = flattened_list 

    dt_str = f"{dt_value:d}"
    output_file_path = os.path.join(tr, f'yedek/probes_P_gh_nz_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Data P_gh successfully saved to {output_file_path}")    

def writeU(tr,dx_value,dt_value):
    fid1 = "/postProcessing/patchProbes1/0/U"
    input_file_path = tr + fid1

    probe_locs = []
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        numberofprobes = sum(1 for line in data_lines if line.startswith('#') and 'Probe' in line)
        data_lines = [line for line in data_lines if not line.startswith('#')]

    

    datainprobes = np.zeros((numberofprobes * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            coordinates = [tuple(map(float, element.strip('()').split())) for element in elements[1:]]
            flattened_list = [value for coord in coordinates for value in coord]
            datainprobes[:, i] = flattened_list 
    
    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/U_nx_{dx_str}_nt_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Data U successfully saved to {output_file_path}")


def writeU_2112(tr,dx_value,dt_value):
    fid1 = "/postProcessing/patchProbes1/0/U"
    input_file_path = tr + fid1

    probe_locs = []
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        numberofprobes = sum(1 for line in data_lines if line.startswith('#') and 'Probe' in line)
        data_lines = [line for line in data_lines if not line.startswith('#')]

    

    datainprobes = np.zeros(((numberofprobes-1) * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            coordinates = [tuple(map(float, element.strip('()').split())) for element in elements[1:]]
            flattened_list = [value for coord in coordinates for value in coord]
            datainprobes[:, i] = flattened_list 
    
    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/U_nx_{dx_str}_nt_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Data U successfully saved to {output_file_path}")


def writeF(tr,dx_value,dt_value):
    fid1 = "/postProcessing/forces1/0/force.dat"
    input_file_path = tr + fid1

   
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        data_lines = [line for line in data_lines if not line.startswith('#')]

    
    
    datainprobes = np.zeros((3 * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            for j in range(1, 9):
                datainprobes[j-1,i]=float(elements[j])
                

    
    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/F_nx_{dx_str}_nt_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Force data successfully saved to {output_file_path}")    


def writeFV2(tr,st2,dx_value,dt_value):
    fid1 = "/postProcessing/"+ st2 +"/0/force.dat"
    input_file_path = tr + fid1

   
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        data_lines = [line for line in data_lines if not line.startswith('#')]

    
    
    datainprobes = np.zeros((3 * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            for j in range(1, 9):
                datainprobes[j-1,i]=float(elements[j])
                

    
    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/{st2}_nx_{dx_str}_nt_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Force data successfully saved to {output_file_path}")        

def writeFV3(tr, name):
    fid1 = "/postProcessing/forces1/0/force.dat"
    input_file_path = tr + fid1

    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        data_lines = [line for line in data_lines if not line.startswith('#')]

    if not data_lines:
        raise ValueError("No valid data lines found in file.")

    num_data_columns = len(data_lines[0].split()) - 1
    datainprobes = np.zeros((num_data_columns, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        time[i] = float(elements[0])
        for j in range(num_data_columns):
            datainprobes[j, i] = float(elements[j + 1])

    output_file_path = os.path.join(tr, name)

    with open(output_file_path, 'w') as file:
        for i in range(len(time)):
            file.write(f"{time[i]:.6e} ")
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            file.write('\n')

    print(f"Force data successfully saved to {output_file_path}")

   
def writeFV3_2112(tr, name):
    fid1 = "/postProcessing/forces1/0/force.dat"
    input_file_path = tr + fid1

    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        data_lines = [line for line in data_lines if not line.startswith('#')]

    if not data_lines:
        raise ValueError("No valid data lines found in file.")

    # Detect if vectors are inside parentheses
    is_vector_format = '(' in data_lines[0]

    time = np.zeros(len(data_lines))
    datainprobes = []

    for i, line in enumerate(data_lines):
        elements = line.split()
        time[i] = float(elements[0])
        if is_vector_format:
            values = []
            for group in elements[1:]:
                if group.startswith('('):
                    current = group.strip('(')
                    if current:
                        values.append(float(current))
                elif group.endswith(')'):
                    current = group.strip(')')
                    if current:
                        values.append(float(current))
                else:
                    values.append(float(group))
            datainprobes.append(values)
        else:
            datainprobes.append([float(x) for x in elements[1:]])

    datainprobes = np.array(datainprobes).T  # Shape: (numVars, numTimesteps)

    output_file_path = os.path.join(tr, name)

    with open(output_file_path, 'w') as file:
        for i in range(len(time)):
            file.write(f"{time[i]:.6e} ")
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            file.write('\n')

    print(f"Force data successfully saved to {output_file_path}")



    

def writeM(tr,dx_value,dt_value):
    fid1 = "/postProcessing/forces1/0/moment.dat"
    input_file_path = tr + fid1

   
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        data_lines = [line for line in data_lines if not line.startswith('#')]

    
    
    datainprobes = np.zeros((3 * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            for j in range(1, 9):
                datainprobes[j-1,i]=float(elements[j])
                

    
    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/M_nx_{dx_str}_nt_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Moment Data successfully saved to {output_file_path}")    


def writeMV2(tr,st2,dx_value,dt_value):
    fid1 = "/postProcessing/"+st2+"/0/moment.dat"
    input_file_path = tr + fid1

   
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        data_lines = [line for line in data_lines if not line.startswith('#')]

    
    
    datainprobes = np.zeros((3 * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            for j in range(1, 9):
                datainprobes[j-1,i]=float(elements[j])
                

    
    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/M{st2}_nx_{dx_str}_nt_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Moment Data successfully saved to {output_file_path}")    


def writePhi(tr,dx_value,dt_value):
    fid1 = "/postProcessing/patchProbes1/0/Phi"
    input_file_path = tr + fid1
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        numberofprobes = sum(1 for line in data_lines if line.startswith('#') and 'Probe' in line)
        data_lines = [line for line in data_lines if not line.startswith('#')]

    print(f"Number of probes: {numberofprobes}")

    datainprobes = np.zeros((numberofprobes, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            coordinates = [tuple(map(float, element.strip('()').split())) for element in elements[1:]]
            flattened_list = [value for coord in coordinates for value in coord]
            datainprobes[:, i] = flattened_list 

    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/Phi_nx_{dx_str}_nt_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Data Phi successfully saved to {output_file_path}")

def writeZetas(tr,dx_value,dt_value):
    fid1 = "/postProcessing/patchProbes1/0/zeta"
    input_file_path = tr + fid1

    probe_locs = []
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        numberofprobes = sum(1 for line in data_lines if line.startswith('#') and 'Probe' in line)
        data_lines = [line for line in data_lines if not line.startswith('#')]

    

    datainprobes = np.zeros((numberofprobes * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            coordinates = [tuple(map(float, element.strip('()').split())) for element in elements[1:]]
            flattened_list = [value for coord in coordinates for value in coord]
            datainprobes[:, i] = flattened_list 
    
    dx_str = f"{dx_value:.3f}"
    dt_str = f"{dt_value:.3f}"
    output_file_path = os.path.join(tr, f'yedek/Zeta_L_{dx_str}_cour_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Data Zeta successfully saved to {output_file_path}")

def writeZetas_2112(tr,dx_value,dt_value):
    fid1 = "/postProcessing/patchProbes1/0/zeta"
    input_file_path = tr + fid1

    probe_locs = []
   
    # Read the file and process lines
    with open(input_file_path, 'r') as file:
        data_lines = [line.strip() for line in file if line.strip()]
        numberofprobes = sum(1 for line in data_lines if line.startswith('#') and 'Probe' in line)
        data_lines = [line for line in data_lines if not line.startswith('#')]

    

    datainprobes = np.zeros(((numberofprobes-1) * 3, len(data_lines)))
    time = np.zeros(len(data_lines))

    for i, line in enumerate(data_lines):
        elements = line.split()
        if elements:
            time[i] = float(elements[0])
            coordinates = [tuple(map(float, element.strip('()').split())) for element in elements[1:]]
            flattened_list = [value for coord in coordinates for value in coord]
            datainprobes[:, i] = flattened_list 
    
    dx_str = f"{dx_value:.3f}"
    dt_str = f"{dt_value:.3f}"
    output_file_path = os.path.join(tr, f'yedek/Zeta_L_{dx_str}_cour_{dt_str}.dat')
    
    # Save the data to a .dat file
    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(time)):
            # Write the time value in scientific notation
            file.write(f"{time[i]:.6e} ")
            
            # Write the probe data for the corresponding time step in scientific notation
            file.write(" ".join(f"{value:.6e}" for value in datainprobes[:, i]))
            
            # Newline for the next time step
            file.write('\n')

    print(f"Data Zeta successfully saved to {output_file_path}")

def writeZeta(time,tr,dx_value,dt_value):
    
    freepath = f"{time:f}/zeta"
    data = read.extract_upper_vector_field(freepath) 
    zeta = data['third_components']
    
    file_path3= f"{time:f}/C"
    c_upper_data=read.parse_upper_section(file_path3)

   
   
    
    
    dx_str = f"{dx_value:d}"
    dt_str = f"{dt_value:.5f}"
    output_file_path = os.path.join(tr, f'yedek/zeta{time:f}_nx_{dx_str}_nt_{dt_str}.dat')

    with open(output_file_path, 'w') as file:
        # Write the header
        #file.write('# Time and Probe Data\n')
        
        # Write the time and probe data
        for i in range(len(c_upper_data[:,0])):
            # Write the time value in scientific notation
            file.write(f"{c_upper_data[i,0]:.6e} {c_upper_data[i,1]:.6e} {c_upper_data[i,2]:.6e} {zeta[i]}")
            file.write('\n')
        
        print(f"Data zeta successfully saved to {output_file_path}")

def writeZeta2(tr,dx_value,dt_value,timeshot):
    data = {}
    fid1 = f"{timeshot}/zeta"
    input_file_path = tr + fid1
    with open(fid1, 'r') as file:
        content = file.read()

    # Extract `upper` boundary nonuniform list
    upper_pattern = r"upper\s*\{\s*type\s+calculated;\s*value\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;"
    upper_match = re.search(upper_pattern, content, re.DOTALL)

    if upper_match:
        num_values = int(upper_match.group(1))
        vector_strings = upper_match.group(2).splitlines()
        third_components = [
            float(re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?", vector)[2])
            for vector in vector_strings if vector.strip()
        ]
        if len(third_components) != num_values:
            raise ValueError("Number of extracted values does not match the count specified.")
        
        node = num_values
        zeta = np.array(third_components)
    else:
        raise ValueError("Upper boundary block not found in the file content.")
    

def modify_openfoam_PHIdict(file_path, new_T, new_Ls, new_xl):

    """
    Update the T, Ls, and xl values in the Phi file.

    Args:
    file_path (str): Path to the Phi file.
    new_T (float): New value for T in the outlet.
    new_Ls (float): New value for Ls in the upper boundary.
    new_xl (float): New value for xl in the upper boundary.
    """
    # Check if the file exists
    if not os.path.isfile(file_path):
        print(f"Error: File '{file_path}' not found.")
        return
    
    # Read the file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Update the values
    with open(file_path, 'w') as file:
        for line in lines:
            if re.match(r'\s*T\s+', line):
                file.write(f"        T               {new_T};\n")
            elif re.match(r'\s*Ls\s+', line):
                file.write(f"            Ls              {new_Ls};\n")
            elif re.match(r'\s*xl\s+', line):
                file.write(f"            xl              {new_xl};\n")
            else:
                file.write(line)
    
    print(f"Phi file updated in '{file_path}'.")   
    