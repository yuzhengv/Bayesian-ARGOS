shell_start=4.7 # 2 Control the number of observations
shell_end=5 # 4.4
shell_by=0.3
shell_seq_end=$(awk "BEGIN {print $shell_end - $shell_by}")
shell_loop_seq=$(awk "BEGIN {for (i=$shell_start; i<=$shell_seq_end; i+=$shell_by) printf i\" \"}")
echo $shell_loop_seq