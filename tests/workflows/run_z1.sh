#!/bin/bash
# Run Z1+ on every .Z1 file in a directory, one subdirectory each.
# Z1+ writes its output into the working directory, so each input needs its
# own or the second run overwrites the first.
SRC="$1"
rm -rf ~/z1run && mkdir -p ~/z1run && cd ~/z1run || exit 1
for f in "$SRC"/*.Z1; do
    name=$(basename "$f" .Z1)
    mkdir -p "$name" && cp "$f" "$name/" && cd "$name" || continue
    ~/z1/Z1+ "$(basename "$f")" > z1.log 2>&1
    echo "=== $name ==="
    if [ -f Z1+summary.dat ]; then
        # column 2 = chains, 6 = mean entanglements per chain, 5 = <Lpp>
        awk '!/^#/ && NF>6 {printf "  chains %s   mean Z per chain %s   <Lpp> %s\n", $2, $6, $5}' Z1+summary.dat
    else
        echo "  no summary; tail of log:"; tail -5 z1.log | sed 's/^/    /'
    fi
    [ -f Z_values.dat ] && echo "  per-chain Z: $(tr '\n' ' ' < Z_values.dat)"
    cd ..
done
