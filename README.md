# LinCapR

## Compile
```
$ make
```

## Run
```
$ ./LinCapR <input_file> <output_file> <beam_size> [-e]
```

input_file: File name of the input RNA sequence (in FASTA format).
output_file: File name of the output structural profile.
beam_size: Parameter for beam_pruning. Set 0 for infinite beam size.
options:
	-e: Output the value of g_ensemble to stdout.

The format of output structural profile is consistent with CapR.
The first line is the name of the sequence stated in the input file.
Succeeding six lines contain the probability that each base position belongs to each structural context, delimited by spaces.