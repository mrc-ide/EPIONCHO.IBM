#!/bin/bash
set -euo pipefail
IFS=$'\n\t'


USER_TO_CHECK="${1:-$USER}"

qstat -t -u "$USER_TO_CHECK" | awk '
NR>5 {
    # Extract base job ID (strip any suffix after first dot)
    job=$1
    sub(/\[.*\]/,"",job)
    sub(/\..*/,"",job)

    # Sum counts by job ID
    if ($10=="R") running[job]++
    else if ($10=="Q") waiting[job]++
    else if ($10=="H") held[job]++
    else finished[job]++

    # Track total tasks
    total[job]++
}
END {
    # Print header
    printf "%-12s %-8s %-8s %-8s %-8s %-10s\n", "Job ID", "Running", "Waiting", "Held", "Finished", "Total"
    
    # Print counts per job ID
    for (j in total) {
        r=running[j]+0
        w=waiting[j]+0
        h=held[j]+0
        f=finished[j]+0
        t=total[j]+0
        printf "%-12s %-8d %-8d %-8d %-8d %-10d\n", j, r, w, h, f, t
    }
}'

exit 0