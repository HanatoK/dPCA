# dPCA
Principle component analysis of dihedral angles.

# Compilation
```
g++ dPCA.cpp -lfmt -I/usr/include/eigen3 -O2 -o dPCA
```

# Input file format
- Columns: dihedral angles.
- Lines: time series or simulation frames.
Columns are separated by spaces.

# Running
```
./dPCA data.dat dihedral
```

# Output
- dihedral.transformed: the sine and cosine values of dihedral angles.
- dihedral.projected: projected values.
