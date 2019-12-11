
Very simple simulation of sequences with random mutations, and simulating linked-read (https://www.10xgenomics.com/linked-reads/) sequencing FASTQ data using `LRTK-SIM`.

To simulate sequences with random mutations, first a FASTA sequence is mutated randomly with `bin/mutate_fasta.py`, and then `bin/make_haps.py` shuffles these mutations randomly, creating FASTA seqences with random mutations.

The `LRTK-SIM/LRTK-SIM.py` script from `zhanglu295/LRTK-SIM` is used to simulate linked reads from the FASTA haplotypes.

> zhanglu295/LRTK-SIM
>
> https://github.com/zhanglu295/LRTK-SIM/

`LRTK-SIM` is included as a submodule, so use `--recursive` when cloning this repo:

```

git clone --recursive https://github.com/olavurmortensen/linkedsim.git
```

