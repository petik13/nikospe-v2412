import subprocess

def mesh(lam):
    subprocess.run(['surfaceFeatureExtract'])
    subprocess.run(['blockMesh'])

    subprocess.run(['decomposePar'])

    # Build topoSets and refineMesh in parallel (creates processor*/constant/polyMesh/sets/*)
    subprocess.run(['foamJob', '-parallel', '-screen', 'topoSet', '-dict', 'system/topoSetDict.1'])
    subprocess.run(['foamJob', '-parallel', '-screen', 'refineMesh', '-dict', 'system/refineMeshDict.1', '-overwrite'])
    subprocess.run(['foamJob', '-parallel', '-screen', 'topoSet', '-dict', 'system/topoSetDict.2'])
    subprocess.run(['foamJob', '-parallel', '-screen', 'refineMesh', '-dict', 'system/refineMeshDict.2', '-overwrite'])
    if lam >= 3.0:
        subprocess.run(['foamJob', '-parallel', '-screen', 'topoSet', '-dict', 'system/topoSetDict.3'])
        subprocess.run(['foamJob', '-parallel', '-screen', 'refineMesh', '-dict', 'system/refineMeshDict.3', '-overwrite'])
        subprocess.run(['foamJob', '-parallel', '-screen', 'topoSet', '-dict', 'system/topoSetDict.4'])
        subprocess.run(['foamJob', '-parallel', '-screen', 'refineMesh', '-dict', 'system/refineMeshDict.4', '-overwrite'])
    if lam >= 6.0:
        subprocess.run(['foamJob', '-parallel', '-screen', 'topoSet', '-dict', 'system/topoSetDict.5'])
        subprocess.run(['foamJob', '-parallel', '-screen', 'refineMesh', '-dict', 'system/refineMeshDict.5', '-overwrite'])

    subprocess.run(['./snappy.sh'])
    

def run_case(Nproc):
    subprocess.run(['rm', '-r', '0'])
    subprocess.run(['cp', '-r', '0.orig', '0'])
    subprocess.run(['topoSet', '-dict', 'system/topoSetDict'])
    subprocess.run(['topoSet', '-dict', 'system/topoSetDict_2'])
    subprocess.run(['renumberMesh', '-overwrite'])
    subprocess.run(['decomposePar'])
    subprocess.run(['foamJob', '-s', '-p', 'renumberMesh', '-overwrite'])

    with open("log", "w") as log:
        proc = subprocess.Popen(
            ['mpirun', '-np', str(Nproc),
             'shipFlow', '-parallel', '-withFunctionObjects'],
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        rc = proc.wait()   # <-- THIS is “detect finish”
    return rc

def post_process():
    subprocess.run(['./prepPost.sh'])