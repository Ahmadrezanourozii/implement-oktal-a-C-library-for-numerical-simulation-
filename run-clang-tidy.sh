#!/bin/bash
# Script to run clang-tidy on the entire project

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${PROJECT_DIR}/build"

# Check if compile_commands.json exists
if [ ! -f "${BUILD_DIR}/compile_commands.json" ]; then
    echo "Error: compile_commands.json not found. Please build the project first."
    echo "Run: cmake -S . -B build && cmake --build build"
    exit 1
fi

echo "Running clang-tidy-19 on entire project..."
echo ""

cd "${PROJECT_DIR}"

# Find all .hpp and .cpp files in include and src directories
# Exclude test files and build directories
# Find all .hpp and .cpp files in include, src, and apps directories
# Exclude test files and build directories
HEADER_FILES=($(find include -name "*.hpp" -type f | grep -v "/test" | sort))
if [ -d "apps" ]; then
    SOURCE_FILES=($(find src apps -name "*.cpp" -type f | grep -v "/test" | sort))
else
    SOURCE_FILES=($(find src -name "*.cpp" -type f | grep -v "/test" | sort))
fi

# Combine into single array
ALL_FILES=("${HEADER_FILES[@]}" "${SOURCE_FILES[@]}")

TOTAL_FILES=${#ALL_FILES[@]}

if [ "${TOTAL_FILES}" -eq 0 ]; then
    echo "No files found to check."
    exit 1
fi

echo "Found ${TOTAL_FILES} files to check"
echo ""

# Determine clang-tidy command
if command -v clang-tidy-19 &> /dev/null; then
    CLANG_TIDY=clang-tidy-19
elif command -v clang-tidy &> /dev/null; then
    CLANG_TIDY=clang-tidy
elif [ -x "/opt/homebrew/opt/llvm/bin/clang-tidy" ]; then
    CLANG_TIDY="/opt/homebrew/opt/llvm/bin/clang-tidy"
else
    echo "Error: clang-tidy not found."
    exit 1
fi

# Run clang-tidy
CURRENT=0

# Detect SDK path for macOS to fix missing headers
SDK_PATH=$(xcrun --show-sdk-path 2>/dev/null)
EXTRA_ARGS=()
if [ -n "$SDK_PATH" ]; then
    EXTRA_ARGS+=("-extra-arg=-isysroot" "-extra-arg=$SDK_PATH")
fi

for file in "${ALL_FILES[@]}"; do
    if [ -f "${file}" ]; then
        CURRENT=$((CURRENT + 1))
        echo "[${CURRENT}/${TOTAL_FILES}] Checking ${file}..."
        
        # Run clang-tidy
        OUTPUT=$($CLANG_TIDY \
            --config-file=.clang-tidy \
            -p "${BUILD_DIR}" \
            "${EXTRA_ARGS[@]}" \
            "${file}" \
            2>&1 | grep -E "(error|warning):" || true)
        
        if [ -n "${OUTPUT}" ]; then
            echo "  Issues found:"
            echo "${OUTPUT}" | sed 's/^/    /'
        else
            echo "  ✓ No issues found"
        fi
        echo ""
    fi
done

echo "Done checking all files!"

