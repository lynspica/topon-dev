#!/bin/bash
# Z1+ with partner attribution on every .Z1 in a directory.
#
# Z1+ writes its output into the working directory, so each configuration
# gets its own. The shortest-path file is copied back beside the input as
# SP_<name>.dat, since that is the one carrying which chain is responsible
# for each entanglement -- the summary only gives counts.
SRC="$1"

# A private working directory per invocation, not a fixed one.
#
# This used to be $HOME/z1batch, wiped with rm -rf on entry. Two batches
# running at once -- a search in the background and anything else in the
# foreground -- both claimed it, and the second one's rm -rf destroyed the
# first one's working directory mid-run. Configurations still being
# processed lost their output, so those candidates silently vanished from
# the batch and the search scored a different set than it had written.
#
# That was the whole of the "same case, same seed, different answer" the
# gallery reported. Everything else is deterministic: the network, the
# conformation, the proposal RNG, Z1+ itself on identical input, and the
# search loop end to end all reproduce exactly when checked in isolation.
WORK="${2:-$(mktemp -d "${TMPDIR:-/tmp}/z1batch.XXXXXX")}"
mkdir -p "$WORK" && cd "$WORK" || exit 1
trap 'rm -rf "$WORK"' EXIT
for f in "$SRC"/*.Z1; do
    [ -e "$f" ] || continue
    name=$(basename "$f" .Z1)
    mkdir -p "$name" && cp "$f" "$name/" && cd "$name" || continue
    ~/z1/Z1+ + "$(basename "$f")" > z1.log 2>&1
    [ -f "Z1+SP.dat" ] && cp "Z1+SP.dat" "$SRC/SP_$name.dat"
    cd ..
done
