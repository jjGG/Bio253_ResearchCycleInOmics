new_sequence = None
old_sequence = None
with open('GCF_000462955.1_ASM46295v1_genomic.fna', 'r') as f:
    new_sequence = "".join([line.strip() for line in f.readlines()[1:]])
    print(f"length of sequence 1: {len(new_sequence)}")

with open('SA_6850_ncbi_sequence.fasta', 'r') as f:
    old_sequence = "".join([line.strip() for line in f.readlines()[1:]])
    print(f"length of sequence 2: {len(old_sequence)}")

if len(new_sequence) != len(old_sequence):
    print("sequence lengths differ")
    exit()

differences = []
for i, n in enumerate(old_sequence):
    if new_sequence[i] != n:
        differences.append(i)
print(differences)
print(f"{len(differences)} differences found")
