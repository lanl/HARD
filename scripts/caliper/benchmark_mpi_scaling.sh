#!/bin/bash

# List of process counts
PROCESS_COUNTS=(1 2 4 8 16)

# Caliper environment configuration
export CALI_CONFIG_FILE=../configs/caliper/hard.conf
export CALI_LOG_VERBOSITY=0

# Command template (use \$X to delay expansion)
COMMAND_TEMPLATE="mpirun -np \$X ./apps/hydro/hydro ./configs/sod.py -d 2"

echo "=== Starting executions ==="
echo "Process counts: ${PROCESS_COUNTS[@]}"
echo ""

# Loop over each process count
for X in "${PROCESS_COUNTS[@]}"; do
    echo ">>> Running with $X processes... [$COMMAND_TEMPLATE]"
    eval "$COMMAND_TEMPLATE"

    # Rename result files
    if [ -f hard-report.cali ]; then
        mv hard-report.cali hard-report-$X.cali
        echo "    File hard-report-$X.cali created"
    fi
    
    if [ -f hard-raw.cali ]; then
        mv hard-raw.cali hard-raw-$X.cali
        echo "    File hard-raw-$X.cali created"
    fi
    echo ""
done

echo "=== Analyzing scaling results ==="

# Build command line for analyze_scaling.py
ANALYZE_ARGS=""
for X in "${PROCESS_COUNTS[@]}"; do
    if [ -f hard-report-$X.cali ]; then
        ANALYZE_ARGS="$ANALYZE_ARGS hard-report-$X.cali=$X"
    else
        echo "WARNING: hard-report-$X.cali does not exist"
    fi
done

# Run analysis
if [ -n "$ANALYZE_ARGS" ]; then
    echo "Command: python ../scripts/caliper/analyze_scaling.py$ANALYZE_ARGS"
    python ../scripts/caliper/analyze_scaling.py$ANALYZE_ARGS
else
    echo "ERROR: No result files found"
fi

echo ""
echo "=== Completed ==="