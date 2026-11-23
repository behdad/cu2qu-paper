#!/bin/bash
# Instructions per approx_spline call on the font cubics, for a kurbo tree.
# Usage: bench.sh <label> <path to kurbo crate>
S=/tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad
W=/home/behdad/hb/cu2qu-work
D="$W/kdiff-out"
label=$1; crate=$2
B="$W/kbench-$label"
mkdir -p "$B/src"
cp "$S/kdiff/bench.rs" "$B/src/main.rs"
cat > "$B/Cargo.toml" <<EOF
[package]
name = "kbench"
version = "0.0.0"
edition = "2021"
publish = false

[dependencies]
kurbo = { path = "$crate" }

[profile.release]
debug = false

[workspace]
EOF
[ -s "$D/font_inputs.txt" ] || sed -n '380001,535066p' "$D/inputs.txt" > "$D/font_inputs.txt"
(cd "$B" && CARGO_TARGET_DIR="$B/target" cargo +1.99.0 build --release -q -j 8) || exit 1
count () { perf stat -x, -e instructions:u "$B/target/release/kbench" "$D/font_inputs.txt" "$1" 2>&1 >/dev/null | awk -F, '/instructions/ {print $1}'; }
base=$(count 0)
n=$(count 20)
calls=$(( $(wc -l < "$D/font_inputs.txt") * 20 ))
echo "$label: $(( (n - base) / calls )) instructions per call ($("$B/target/release/kbench" "$D/font_inputs.txt" 1))"
