# Custom installation script to install mxfold2 package from a wheel file

# Define the path to the wheel file
#change to mxfold2-0.1.2-cp310-cp310-manylinux_2_17_x86_64
wheel_path="setting/mxfold2-0.1.2-cp310-cp310-macosx_13_0_arm64.whl"

# Print the current directory for debugging
echo "Current directory: $(pwd)"

# Print the full path to the wheel file for debugging
echo "Wheel file path: $(pwd)/$wheel_path"

# Check if the wheel file exists
if [[ -f "$wheel_path" ]]; then
    # Install mxfold2 package
    pip3 install "$wheel_path"
else
    echo "Wheel file not found: $wheel_path"
    exit 1
fi