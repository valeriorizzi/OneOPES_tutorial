#!/bin/bash

filename="fes_phi.dat"
folders=(0 1 2 3)
n=${#folders[@]}

# Construct paths
files=()
for f in "${folders[@]}"; do
    files+=("$f/$filename")
done

# Paste files and process with awk
paste "${files[@]}" | awk -v n="$n" '
# Skip lines starting with #
/^#!/ { next }

# Process data lines
{
    x = $1
    sum = 0
    sq_sum = 0
    
    # Free Energy values are in columns 2, 4, 6, 8, 10, 12, 14, 16
    # Index formula: i * 2
    for (i=1; i<=n; i++) {
        val = $(i * 2)
        sum += val
        sq_sum += (val * val)
    }

    avg = sum / n
    # Standard deviation: sqrt( <x^2> - <x>^2 )
    var = (sq_sum / n) - (avg * avg)
    # Handle potential tiny negative numbers from floating point precision
    std = (var > 0) ? sqrt(var) : 0
    
    printf "%15.6f %15.6f %15.6f\n", x, avg, std
}' > fes_phi_average_std.dat

echo "Done! Results saved in fes_average_std.dat"
