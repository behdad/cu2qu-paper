#!/bin/bash
# The checks of kurbo's CI that can run here without cargo-hack / nextest /
# taplo / typos: fmt, copyright headers, clippy (std, all features, tests;
# no_std + libm on the host target), tests, doc tests, MSRV check, docs with
# private items.
export CARGO_TARGET_DIR=/home/behdad/hb/cu2qu-work/kurbo-target
cd /home/behdad/kurbo || exit 1
J="-j 8"
step () { echo "=== $1"; shift; "$@" > /tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/kurbo_step.log 2>&1; rc=$?; tail -${TAIL:-3} /tmp/claude-1000/-home-behdad-harfbuzz/7a6a948d-c79c-4ca1-8b42-b271f66712b8/scratchpad/kurbo_step.log; echo "    exit $rc"; }
step "cargo fmt --all --check" cargo +1.99.0 fmt --all --check
step "copyright headers" bash .github/copyright.sh
step "clippy (std, all features)" cargo +1.99.0 clippy $J --workspace --locked --all-features -- -D warnings
step "clippy (std, default features, tests benches examples)" cargo +1.99.0 clippy $J --workspace --locked --tests --benches --examples -- -D warnings
step "clippy (all features, tests benches examples)" cargo +1.99.0 clippy $J --workspace --locked --all-features --tests --benches --examples -- -D warnings
step "clippy (no_std + libm, host target)" cargo +1.99.0 clippy $J -p kurbo --locked --no-default-features --features libm -- -D warnings
TAIL=4 step "tests (workspace, all features)" cargo +1.99.0 test $J --workspace --locked --all-features --no-fail-fast --lib --bins --tests
step "doc tests" cargo +1.99.0 test $J --doc --workspace --locked --all-features --no-fail-fast
step "msrv 1.85 check (std)" cargo +1.85 check $J -p kurbo -p polycool --locked --features std
step "msrv 1.85 check (no_std + libm)" cargo +1.85 check $J -p kurbo --locked --no-default-features --features libm
step "cargo doc (private items)" env RUSTDOCFLAGS="-D warnings" cargo +nightly doc $J --workspace --locked --all-features --no-deps --document-private-items
