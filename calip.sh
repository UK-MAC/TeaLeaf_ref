num_ranks=4

# e.g. to display percentage total time spent in each function for each rank
for rank in $(seq 0 $((num_ranks-1))); do
    echo "rank ${rank}"
    cali-query -t \
        -q "select annotation, loop, function, mpi.function, percent_total(sum#time.duration)" \
        caliper/caliper-${rank}.cali
    echo
done
