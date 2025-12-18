import re
import numpy as np
import os

##08 FEB 2025
##V3 Update: change in oscparameter to read amp ## added dispersion osman 08 FEB 2025
##08 FEB 2025

##27 March 2025 V4 
# #extract_waveCur_parameters_from_file
##28 March
##read_wavecon_dict & update_radius_in_topoSetDict


def extract_waveCur_parameters_from_file(file_path):
    """
    Extracts steepness and wavelength from a given C++ file.

    Parameters:
    - file_path (str): Path to the C++ file.

    Returns:
    - dict: Extracted wave parameters.
    """
    with open(file_path, 'r') as file:
        code_string = file.read()

    # Regular expressions to extract values
    steepness_match = re.search(r'const\s+scalar\s+steepness\(([\d.]+)\);', code_string)
    wavelength_match = re.search(r'const\s+scalar\s+wavelength\s+\(([\d.]+)\);', code_string)
    current_match = re.search(r'const\s+scalar\s+currentspeed\s+\(([\d.]+)\);', code_string)
    hdepth_match = re.search(r'const\s+scalar\s+hdepth\s+\(([\d.]+)\);', code_string)
    
    # Parse values
    steepness = float(steepness_match.group(1)) if steepness_match else None
    wavelength = float(wavelength_match.group(1)) if wavelength_match else None
    currentspeede = float(current_match.group(1)) if current_match else None
    depth = float(hdepth_match.group(1)) if hdepth_match else None
    return {
        "steepness": steepness,
        "wavelength": wavelength,
        "current" : currentspeede,
        "depth" : depth
    }

def update_radius_in_topoSetDict(tr, new_radius):
    """
    Reads an OpenFOAM topoSetDict file, modifies the radius value,
    and writes the updated file back.

    Parameters:
    - filename (str): Path to the topoSetDict file.
    - new_radius (float): The new radius to set.
    """
    filename = os.path.join(tr, "system", "topoSetDict")
    
    with open(filename, 'r') as f:
        lines = f.readlines()

    updated_lines = []
    for line in lines:
        if 'radius' in line.strip().split():  # safer check for the word "radius"
            parts = line.strip().split()
            indent = line[:len(line) - len(line.lstrip())]
            updated_line = f"{indent}radius      {new_radius};\n"
            updated_lines.append(updated_line)
        else:
            updated_lines.append(line)

    with open(filename, 'w') as f:
        f.writelines(updated_lines)

    print(f"Updated radius in {filename} to {new_radius}")


def read_wavecon_dict(file_path):
    data = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        line = line.strip()

        # Skip comments and empty lines
        if not line or line.startswith('//') or line.startswith('/*') or line.startswith('FoamFile'):
            continue
        
        # Stop parsing when the header block ends
        if line == '}':
            continue
        
        # Parse key-value pairs
        if ';' in line:
            parts = line.split(';')[0].split()
            if len(parts) == 2:
                key, value = parts
                try:
                    data[key] = float(value)
                except ValueError:
                    data[key] = value
                    
    return data

def parse_internal_field_from_file(file_path):
    """
    Parse the internal field and dimensions from an OpenFOAM file.

    Args:
        file_path (str): Path to the OpenFOAM file.

    Returns:
        dict: A dictionary containing dimensions and internal field values.
    """
    data = {}

    with open(file_path, 'r') as file:
        content = file.read()

    # Extract dimensions
    dimensions_pattern = r"dimensions\s+\[([^\]]+)\];"
    dimensions_match = re.search(dimensions_pattern, content)
    if dimensions_match:
        data['dimensions'] = [float(x) for x in dimensions_match.group(1).split()]
    
    # Extract internal field
    internal_field_pattern = r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s+\((.*?)\)\s*;"
    internal_field_match = re.search(internal_field_pattern, content, re.DOTALL)
    if internal_field_match:
        num_values = int(internal_field_match.group(1))
        values = [float(x) for x in internal_field_match.group(2).split()]
        if len(values) != num_values:
            raise ValueError("Number of values does not match the count specified.")
        data['internalField'] = np.array(values)  # Convert to NumPy array
    
    return data



def parse_internal_vector_field(file_path):
    """
    Parse the internal vector field and dimensions from an OpenFOAM file.

    Args:
        file_path (str): Path to the OpenFOAM file.

    Returns:
        dict: A dictionary containing dimensions and internal vector field values.
    """
    data = {}

    with open(file_path, 'r') as file:
        content = file.read()

    # Extract dimensions
    dimensions_pattern = r"dimensions\s+\[([^\]]+)\];"
    dimensions_match = re.search(dimensions_pattern, content)
    if dimensions_match:
        data['dimensions'] = [float(x) for x in dimensions_match.group(1).split()]
    
    # Extract internal vector field
    internal_field_pattern = r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s+\((.*?)\)\s*;"
    internal_field_match = re.search(internal_field_pattern, content, re.DOTALL)
    if internal_field_match:
        num_values = int(internal_field_match.group(1))
        vector_strings = internal_field_match.group(2).splitlines()
        vectors = [
            [float(component) for component in re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?", vector)]
            for vector in vector_strings if vector.strip()
        ]
        if len(vectors) != num_values:
            raise ValueError("Number of vectors does not match the count specified.")
        data['internalField'] = np.array(vectors)  # Convert to NumPy array
    
    return data


def extract_upper_vector_field(file_path: str):
    """
    Extract the numerical values from the `upper` boundary field in an OpenFOAM file.

    Args:
        file_path (str): Path to the OpenFOAM file.

    Returns:
        dict: A dictionary containing the number of values and the third vector component values.
    """
    data = {}

    with open(file_path, 'r') as file:
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
        
        data['num_values'] = num_values
        data['third_components'] = np.array(third_components)
    else:
        raise ValueError("Upper boundary block not found in the file content.")
    
    return data

### 3 FUNCTIONS BELOWS READ WAVE PARAMETERS

def extract_osc_or_wave_parameters(file_pathosc,file_wave):

    if os.path.isfile(file_wave):
        print(f"Wave File exists: {file_wave}")
        if os.path.isfile(file_pathosc):
            raise ValueError("Both wave and forced oscillation exist")

        parameters_list=extract_wave_parameters_from_file(file_wave)
        wavelength = parameters_list["wavelength"]  # Wavelength (in meters)
        wavenumber = 2.0 * np.pi / wavelength  # Wavenumber (2*pi / wavelength)
        omega = np.sqrt(9.81 * wavenumber)  # Angular frequency from dispersion relation
        steepness=parameters_list["steepness"]
        amp=0.5*steepness*wavelength
        period=2.0 * np.pi / omega
        celerity=wavelength/period
        Cg=0.5*celerity

        simtype="Wave Scattering"
        print("Simulation type: Wave Scattering")
        print("w angular frequency:", omega,    "T period:", period,"       L wavelength:", wavelength)
        print("Phase velocity:",celerity,"     Group velocity:",Cg) 
        print("Steepnes:",steepness,"     ampilutude:",amp)   


    elif os.path.isfile(file_pathosc):
        print(f"Forced File exists: {file_pathosc}")
        parameters_list=extract_osc_parameters_from_file(file_pathosc)        
        period = parameters_list["T"]  # Wavelength (in meters)
        omega = 2*np.pi/period  # Angular frequency from dispersion relation
        wavenumber=omega*omega/9.81
        wavelength= 2*np.pi/wavenumber 
        amp=parameters_list["amp"]

        period=2.0 * np.pi / omega
        celerity=wavelength/period
        Cg=0.5*celerity
        simtype="Forced Oscillation"
        print("Simulation type: Forced Oscillation")
        print("w angular frequency:", omega,    "T period:", period,"       L wavelength:", wavelength)
        print("Phase velocity:",celerity,"     Group velocity:",Cg)
        print("Oscillation amp:", parameters_list["amp"])
    
    
    else:
        print(f"Force or wave File not found: {file_pathosc} or {file_wave}")

    return {
        "celerity": celerity,
        "wavelength": wavelength,
        "period": period,
        "amplitude":amp,
        "simtype":simtype
    }

def extract_wave_parameters_from_file(file_path):
    """
    Extracts steepness and wavelength from a given C++ file.

    Parameters:
    - file_path (str): Path to the C++ file.

    Returns:
    - dict: Extracted wave parameters.
    """
    with open(file_path, 'r') as file:
        code_string = file.read()

    # Regular expressions to extract values
    steepness_match = re.search(r'const\s+scalar\s+steepness\(([\d.]+)\);', code_string)
    wavelength_match = re.search(r'const\s+scalar\s+wavelength\s+\(([\d.]+)\);', code_string)

    # Parse values
    steepness = float(steepness_match.group(1)) if steepness_match else None
    wavelength = float(wavelength_match.group(1)) if wavelength_match else None

    return {
        "steepness": steepness,
        "wavelength": wavelength
    }

def dispersionosman(T, h):
    g = 9.81
    omega = 2 * np.pi / T
    y = (omega ** 2) * h / g

    # Step 1: Initial approximation of kh
    kh0 = (y + (y ** 1.986) * np.exp(-1.863 - 1.198 * (y ** 1.366))) / np.sqrt(np.tanh(y))

    # Step 2: Improved approximation
    kh = (kh0 ** 2 + y * (np.cosh(kh0) ** 2)) / (kh0 + 0.5 * np.sinh(2 * kh0))
    k=kh/h
    return k

def checks(period, tankdepth,tankwidth):
    wavenumber=dispersionosman(period,tankdepth)
    wavelengthfromdispersion = 2*np.pi/wavenumber
    if tankdepth < 0.5*wavelengthfromdispersion:
        raise ValueError("Deepwater limit exceeded.")
    else:
        print("Deep Water Wave")
    pi = np.pi
    
    # Check primary resonance condition (kw = 2 * n * pi)
    n_resonance = (wavenumber * tankwidth) / (2 * pi)
    resonance = np.isclose(n_resonance, round(n_resonance))  # Check if n is an integer
    print("RES CHECK",n_resonance,resonance)

    # Check near-resonance condition (gap = n * lambda or (n + 1/2) * lambda) BUNA BAK ELONGATED BODIES
    #n_gap = gap / wavelength
    #near_resonance = np.isclose(n_gap, round(n_gap)) or np.isclose(n_gap, round(n_gap - 0.5) + 0.5)

def extract_osc_parameters_from_file(file_path):
    """
    Extracts steepness and wavelength from a given C++ file.

    Parameters:
    - file_path (str): Path to the C++ file.

    Returns:
    - dict: Extracted wave parameters.
    """
    with open(file_path, 'r') as file:
        code_string = file.read()

    # Regular expressions to extract values
    period_match = re.search(r'const\s+scalar\s+T\(([\d.]+)\);', code_string)
    amp_match = re.search(r'const\s+scalar\s+amp\s*\(\s*([\d.]+)\s*\);', code_string)

    
    # Parse values
    period = float(period_match.group(1)) if period_match else None
    amp = float(amp_match.group(1)) if amp_match else None

    return {
        "T": period,
        "amp": amp
        #"wavelength": wavelength
    }

def parse_probes(file_path: str):
    """
    Parses probe coordinates and velocities from the file.

    Args:
        file_path (str): Path to the probe data file.

    Returns:
        dict: A dictionary with coordinates and velocities for each time step.
    """
    probe_data = {
        'coordinates': [],
        'time_series': []
    }

    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Extract coordinates
    coordinate_pattern = re.compile(r"# Probe \d+ \((.*?)\)")
    for line in lines:
        coord_match = coordinate_pattern.match(line)
        if coord_match:
            coords = tuple(map(float, coord_match.group(1).split()))
            probe_data['coordinates'].append(coords)

    # Extract time and velocities
    time_vel_pattern = re.compile(r"^(\d+\.\d+)\s+(.*)$")
    for line in lines:
        time_match = time_vel_pattern.match(line)
        if time_match:
            time = float(time_match.group(1))
            velocities = re.findall(
                r"\(([-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?)\s+"  # x-component
                r"([-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?)\s+"    # y-component
                r"([-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?)\)",   # z-component
                time_match.group(2)
            )
            velocities = [(float(v[0]), float(v[1]), float(v[2])) for v in velocities]  # Extract all components
            probe_data['time_series'].append((time, velocities))

    return probe_data


def parse_upper_section(file_path):
    """
    Parse the internal vector field values for the 'upper' section from an OpenFOAM file.

    Args:
        file_path (str): Path to the OpenFOAM file.

    Returns:
        np.ndarray: A NumPy array containing the internal vector field values for the 'upper' section.
    """
    data = []

    with open(file_path, 'r') as file:
        content = file.read()

    # Regex to locate the "upper" section with vectors
    upper_section_pattern = (
        r"upper\s*\{\s*type\s*calculated;\s*value\s*nonuniform\s+List<vector>\s+"
        r"(\d+)\s*\(\s*(.*?)\s*\)\s*;"
    )
    upper_match = re.search(upper_section_pattern, content, re.DOTALL)

    if upper_match:
        num_values = int(upper_match.group(1))  # Number of vectors
        vector_strings = upper_match.group(2).strip().splitlines()  # Split vectors into lines

        # Parse each vector string into a list of floats
        vectors = [
            [float(component) for component in re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:e[-+]?[0-9]+)?", vector)]
            for vector in vector_strings
        ]

        # Ensure the number of parsed vectors matches the expected count
        if len(vectors) != num_values:
            raise ValueError("Number of vectors does not match the count specified for the 'upper' section.")
        
        data = np.array(vectors)  # Convert to NumPy array

    return data