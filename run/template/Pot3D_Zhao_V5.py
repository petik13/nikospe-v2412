import os
import subprocess
import numpy as np
import reader_functionsV4 as read
import updater_writerV7 as up
import re
import shutil

import math
##import glob

def find_latest_time_folder(path, target_time, tol=0.4):
    all_dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    time_dirs = []
    for d in all_dirs:
        try:
            t = float(d)
            if abs(t - target_time) < tol:
                time_dirs.append((t, d))
        except ValueError:
            pass
    if not time_dirs:
        raise FileNotFoundError(f"No folder close to time {target_time}")
    # Return the closest match
    return os.path.join(path, sorted(time_dirs)[-1][1])

def real_roots_only(a, b, c):
    D = b**2 - 4*a*c
    if a==0:
        return -c/b
    
    if D < 0:
        return None
    elif D == 0:
        root = -b / (2*a)
        return root
    else:
        root1 = (-b + math.sqrt(D)) / (2*a)
        root2 = (-b - math.sqrt(D)) / (2*a)
        if root1 > 0:
            return root1
        elif root2 > 0:
            return root2




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



def modify_sphere(file_path, new_centre, new_radius):
    """
    Modify centre and radius values of a searchableSphere block in a 
    snappyHexMeshDict-style OpenFOAM dictionary.

    Args:
        file_path (str): Path to the dictionary file.
        new_centre (tuple or list of float): New coordinates for sphere centre (x, y, z).
        new_radius (float): New radius value.
    """
    with open(file_path, 'r') as f:
        content = f.read()

    # Modify the 'centre' line
    centre_pattern = r'centre\s*\(\s*[^ ]+\s+[^ ]+\s+[^ )]+\s*\);'
    new_centre_str = f'centre          ({new_centre[0]} {new_centre[1]} {new_centre[2]});'
    content = re.sub(centre_pattern, new_centre_str, content)

    # Modify the 'radius' line
    radius_pattern = r'radius\s+[^;]+;'
    new_radius_str = f'radius          {new_radius};'
    content = re.sub(radius_pattern, new_radius_str, content)

    with open(file_path, 'w') as f:
        f.write(content)


def modify_distance_levels(file_path, new_val1, new_level1, new_val2, new_level2):
    with open(file_path, 'r') as f:
        content = f.read()

    pattern = re.compile(
        r'(refinementRegions\s*\{\s*cylinder\s*\{\s*[^{}]*?levels\s*\()\s*\([^)]+\)\s*\([^)]+\)\s*\)(\s*;\s*)',
        re.DOTALL
    )

    replacement = (
        r'\1\n'
        f'        ({new_val1} {new_level1})\n'
        f'        ({new_val2} {new_level2})\n    )' + r'\2'
    )

    new_content = re.sub(pattern, replacement, content)

    with open(file_path, 'w') as f:
        f.write(new_content)


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

def extract_and_save_data(input_file, output_file):
    """
    Reads numerical time/value pairs from a text file (skipping headers)
    and writes them to a .dat file in tab-separated format.

    Args:
        input_file (str): Path to the input text file.
        output_file (str): Path for the output .dat file.
    """
    data = []

    with open(input_file, 'r') as f:
        for line in f:
            # Skip comments, blank lines, and header
            if line.strip().startswith('#') or line.strip() == '' or line.strip().startswith('Time'):
                continue
            parts = line.split()
            if len(parts) == 2:
                try:
                    time = float(parts[0])
                    value = float(parts[1])
                    data.append((time, value))
                except ValueError:
                    continue

    with open(output_file, 'w') as f:
        for time, value in data:
            f.write(f"{time:.6f}\t{value:.12e}\n")

    print(f"Data saved to {output_file}")

def readF(tr,str):
    fid1 = "/postProcessing/forces"+ str +"/0/force.dat"
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


    return datainprobes,time   



def writeF(tr,name,datainprobes,time):
    
    output_file_path = os.path.join(tr, name)
    
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

def ceil_to4(n):
    return ((n + 3) // 4) * 4

runfoam=1

R=0.5

#w2Rg=[1.58]
w2Rg=[1.58,1.50,1.525,1.48] #1.48

#w2Rg=[1.21,1.24,1.25,1.18] ##1.21 #5VB de kullandik


#w2Rg=[0.90,0.89,1.18] ##V5BMeshde oldu method          hierarchical;     n       (9 4 1);




#w2Rg=[1.3,1.325,1.46,1.47] #w2Rg=[1.48,1.50,1.525] hiearchi de oldu

#w2Rg=[1.17,1.18,1.35,1.40,1.45] #1.14

#w2Rg=[1.15,1.24,1.2,1.35,1.40,1.45] #1.14



##simpleda calisabilenler  ## 1.14 1.24,1.35,

#0.82,0.825 1.15 1.185 1.19 1.2 1.21 1.25  1.45hata cikti paralell no neighbours error


#w2Rg=[0.85,0.90,0.95,0.98,1.1,0.875,1.03,1.17,1.18] #worked simple method

#w2Rg=[1.15,1.17,1.18,1.19,1.2,1.35,1.40,1.45] #worked on heararchical

#w2Rg=[0.85,0.90,0.95,0.98,1.1]
#w2Rg=[0.9,0.95,0.98,1.1]
#w2Rg=[0.95,0.98,1.1,1.15,1.18,1.21,1.25]
#w2Rg=[1.21,1.25,1.3,1.44,1.47,1.5,1.55,1.58]
#w2Rg=np.linspace(0.1,1.6,16)  #[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.  1.1 1.2 1.3 1.4 1.5 1.6]


print(w2Rg)

##UsqrtgR=0.128
UsqrtgR=0.064
steepness=0.0167


Nx=40

cour=0.1




domainlen=20           
model_pos=10
print("Ucur :",UsqrtgR*np.sqrt(9.81*R))



U0 =  UsqrtgR*np.sqrt(9.81*R)


tr = os.getcwd()
blockMesh_path = tr + "/system/blockMeshDict"
controlDict_path = tr + "/system/controlDict"
file_path0 = os.path.join(tr, "constant", "waveCurConditions")
snappyHexDict = os.path.join(tr, "system", "snappyHexMeshDict")
##new_vals = {'steepness': steepness,'wavelength': wavelength,'currentspeed': U0,'waterdepth': wdepth,'v0': 3.0,'Ls': 5*wavelength,'xl': 15*wavelength}
##modify_wave_dict(file_path0, new_vals)


##modify_box_coordinates( tr + "/system/snappyHexMeshDict_orig", tr + "/system/snappyHexMeshDict_mod", "box1x1x1", "(1 -100 -0.01)", "(4 100 0.01)")

##alttakini de degistirmek gerekebilir
##new_hex_values = f"(0 1 2 3 4 5 6 7) ({number*domainlen} 1 {int(np.round(wdepth/dx))}) simpleGrading (1 1 0.5)"
    ##up.update_hex_block(blockMesh_path, new_hex_values)
    
ncase=len(w2Rg)
tankwidth=10.5

for c in range(0,ncase):
    
    

    
    
    omega0=real_roots_only(U0/9.81, 1, -np.sqrt(w2Rg[c]*9.81/(R)))
    wavenumber=omega0*omega0/9.81
    wavelength= 2.0 * np.pi / wavenumber 
     
    tank_length = round(domainlen*wavelength,2)
    wdepth=round(1*wavelength,2)   ##depth 1 oldu
    nwidth=int(tankwidth/wavelength)


    omega = U0*wavenumber + np.sqrt(9.81*wavenumber*np.tanh(wavenumber*wdepth)) 
    period=2.0 * np.pi / omega
    
    probelocx=[5*wavelength,10*wavelength-1.1*R,10*wavelength+1.1*R,15*wavelength,10*wavelength,10*wavelength]
    probelocy=[0,0,0,0,+1.1*R,-1.1*R]
    
        

    Number_x=int(Nx*domainlen)
    #Number_y=4*int(np.floor(0.25*tankwidth*Nx/wavelength))
    Number_y=ceil_to4(int(tankwidth*Nx/wavelength))
    
    Number_z=int(Nx*1)
    
    celerity=wavelength/period
    Cg=0.5*celerity
    dt=round(cour*(wavelength/Nx)/(Cg+U0),5)
    ##dt=round(cour*(wavelength/Nx)/(U0),5)

    print("w2Rg[c] :",w2Rg[c])
    print("omega :",omega)
    print("period :",period)
    print("wavelength :",wavelength)
    print("wavenumber :",wavenumber)
    
    
    print("time step: ",dt)
    print("tank length:",tank_length)
    print("water depth:",wdepth)
    print("CG",Cg)
    print("Tank width/wavelength =",tankwidth/wavelength )
    print("Approximate number of Cells =",Number_x*Number_y*Number_z)
    
    

    new_vertices = [( 0.0, -tankwidth/2, -wdepth), ( tank_length, -tankwidth/2, -wdepth),( tank_length, tankwidth/2, -wdepth),( 0.0, tankwidth/2, -wdepth),
                ( 0.0, -tankwidth/2, 0),( tank_length, -tankwidth/2, 0),( tank_length, tankwidth/2, 0),( 0.0, tankwidth/2, 0)]
    new_hex_values = f"(0 1 2 3 4 5 6 7) ({Number_x:d} {Number_y:d} {Number_z:d}) simpleGrading (1 1 0.5)"
    up.update_vertices(blockMesh_path,new_vertices )
    up.update_hex_block(blockMesh_path, new_hex_values)
    
    new_vals = {'steepness': steepness,'wavelength': wavelength,'currentspeed': U0,'waterdepth': wdepth,'v0': 3.0,'Ls': round(5*wavelength,2),'xl': round(15*wavelength,2),'R0': R,'xc' : model_pos} ###buraya xc ekeleme
    modify_wave_dict(file_path0, new_vals)
    
    
    modify_sphere(snappyHexDict, new_centre=(model_pos*wavelength, 0.0, 0.0), new_radius=R)

    
    
    
    #up.update_control2(controlDict_path, 0.001, 0.005,1,probelocx,probelocy) ##dt,int(np.ceil(45*period)),int(np.ceil(15*period))
    up.update_control2(controlDict_path, dt, int(np.ceil(60*period)),int(np.ceil(30*period)),probelocx,probelocy)
    #up.update_control(controlDict_path, dt,3*my_int,my_int,[0])
    
    

    if runfoam==1:
        cmd = 'foamCleanTutorials'
        completed_process = subprocess.run(cmd, shell=True)
        #cmd = './runMesh'
        cmd = './Allrun'
        completed_process = subprocess.run(cmd, shell=True)
        
        #for dir in range(0,2):
        #    src_folder = find_latest_time_folder(tr, (dir+1)*my_int)
        #    dst_folder = tr +'/yedek/' +str((dir+1)*my_int) + f'_wavelength_{wavelength:.3f}'
        #    shutil.copytree(src_folder, dst_folder)
        up.writeZetas(tr,wavelength,cour)
        output_file_path = os.path.join(tr, f'yedek/Awavelength_{wavelength:.3f}.dat')
        data,time=readF(tr,"1")
        writeF(tr,output_file_path,data,time) ## up.writeF(tr,radius,dt)
        output_file_path = os.path.join(tr, f'yedek/Bwavelength_{wavelength:.3f}.dat')
        data,time=readF(tr,"2")
        writeF(tr,output_file_path,data,time) 
        output_file_path = os.path.join(tr, f'yedek/Cwavelength_{wavelength:.3f}.dat')
        data,time=readF(tr,"3")
        writeF(tr,output_file_path,data,time) 
		

##The End
#vertices (( 0.0 -0.6 -1) ( 20.0 -0.6 -1)( 20.0 0.6 -1)( 0.0 0.6 -1)( 0.0 -0.6 0)( 20.0 -0.6 0)( 20.0 0.6 0)( 0.0 0.6 0));