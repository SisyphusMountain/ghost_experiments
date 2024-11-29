#!/bin/bash

# Exit script if any command fails
set -e

# Step 1: Initialize and update submodules
echo "Cloning submodules..."
git submodule update --init --recursive

# Step 2: Compile the Rust script in gene_transfer_script
echo "Compiling Rust script in gene_transfer_script..."
cd gene_transfer_script

# Check if Rust is installed
if ! command -v cargo &> /dev/null
then
    echo "Rust is not installed. Please install Rust before proceeding."
    exit 1
fi

# Build the Rust project
cargo build --release

# Step 3: Copy and rename the binary to ./gene_transfer_script_
echo "Copying and renaming the binary to ./gene_transfer_script_..."
cd ..
cp gene_transfer_script/target/release/gene_transfer_script ./gene_transfer_script_

echo "Installation complete."
