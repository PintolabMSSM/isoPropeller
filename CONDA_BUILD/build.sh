#!/bin/bash
set -e  # Exit on error

# Ensure target directory exists
mkdir -p "$PREFIX/bin"

# Copy scripts
cp -ra bin/* "$PREFIX/bin/"

# Make sure all scripts are executable
chmod +x "$PREFIX/bin/"*
