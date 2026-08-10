#!/bin/bash
# Z1+ with partner attribution on every .Z1 in a directory.
#
# Z1+ writes its output into the working directory, so each configuration
# gets its own. The shortest-path file is copied back beside the input as
# SP_<name>.dat, since that is the one carrying which chain is responsible
# for each entanglement -- the summary only gives counts.
SRC="$1"
WORK="${2:-$HOME/z1batch}"

rm -rf "$WORK" && mkdir -p "$WORK" && cd "$WORK" || exit 1
for f in "$SRC"/*.Z1; do
    [ -e "$f" ] || continue
    name=$(basename "$f" .Z1)
    mkdir -p "$name" && cp "$f" "$name/" && cd "$name" || continue
    ~/z1/Z1+ + "$(basename "$f")" > z1.log 2>&1
    [ -f "Z1+SP.dat" ] && cp "Z1+SP.dat" "$SRC/SP_$name.dat"
    cd ..
done
