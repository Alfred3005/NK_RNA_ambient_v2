import h5py
try:
    f = h5py.File('referencias/scanvi_sin_adultos.h5ad', 'r')
    print("Keys in H5AD:", list(f.keys()))
    if 'obs' in f:
        print("\nObs datasets:", list(f['obs'].keys()))
        # Check for age_group
        if 'age_group' in f['obs']:
            # For categorical data in H5AD, it might be in __categories
            pass
    if 'var' in f:
        print("\nVar datasets:", list(f['var'].keys()))
except Exception as e:
    print(f"Error: {e}")
