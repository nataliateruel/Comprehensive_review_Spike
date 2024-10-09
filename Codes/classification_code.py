import os
import pandas as pd

"""
This code will describe a tutorial for classifying antibodies among the 14 predefined epitopes.
This pipeline requires the use of Surfaces (https://doi.org/10.1093/bioinformatics/btad608).
We will use as an example the structure 7KXK, to evaluate Fab 15033-7.
"""

# Step 1: Install Surfaces, as described in https://surfaces-tutorial.readthedocs.io/en/latest/index.html

# Step 2: Run Surfaces, making sure that these files are in the Surfaces directory

pdb_file = '7KXK'
chain_Spike = 'B'
chains_AB = 'MI'

try:
    print ('Cleaning structure...')
    os.system('python clean_structure.py -f '+pdb_file+'.pdb -def AMINO_FlexAID.def')
    print ('Running Surfaces...')
    os.system('python surface_cont.py -f clean_'+pdb_file+'.pdb -c1 '+chain_Spike+' -c2 '+chains_AB+' -o '+pdb_file+'_output.csv -def AMINO_FlexAID.def -dat FlexAID.dat')
except:
    print ('ERROR with Surfaces')

# Step 3: Verify alignment

print ('Verifying alignment...')

# from 14 to 1165
seq_ref = '-------------QCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVD'
aa = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

Surfaces_results = pd.read_csv(pdb_file + "_output.csv", index_col=0)
index = Surfaces_results.index.values.tolist()

def build_sequence(index, start=1, end=1165):
    numbers = []
    aminoacids = []
    for res in index:
        numbers.append(int(res[3:-1]))
        aminoacids.append(aa[res[:3]])
    seq = ''
    for i in range(start, end + 1):
        if i in numbers:
            idx = numbers.index(i)
            seq += aminoacids[idx]
        else:
            seq += '-'
    return seq

def compare_sequence(reference_sequence, generated_sequence):
    max_length = max(len(generated_sequence), len(reference_sequence))
    generated_sequence = generated_sequence.ljust(max_length, '-')
    reference_sequence = reference_sequence.ljust(max_length, '-')
    matches = 0
    valid_positions = 0
    for gen_char, ref_char in zip(generated_sequence, reference_sequence):
        if gen_char != '-' and ref_char != '-':
            valid_positions += 1
            if gen_char == ref_char:
                matches += 1
    if valid_positions == 0:
        return 0
    match_percentage = (matches / valid_positions) * 100
    return match_percentage

seq = build_sequence(index)
k = compare_sequence(seq_ref, seq)

if k < 70:
    print ('WARNING: VERIFY RESIDUE NUMBERING')

# Step 4: Build a vector of interactions

print ('Building vector of interactions...')

def build_vector(df_cut):
    AB_residues = df_cut.columns.values.tolist()
    Spike_residues = df_cut.index.values.tolist()
    Sum = df_cut.sum(axis=1).values.tolist()
    Full_interaction = sum(Sum)
    return (Spike_residues, Sum, Full_interaction)

def order_vector(residues, interactions, file, chain):
    name = file[:4] + "_" + chain
    num_residues = []
    net_interaction = []
    for i in range(14, 1166):
        num_residues.append(i)
        net_interaction.append(0)
    for j in range(len(residues)):
        res_num = int(residues[j][3:-1])
        if res_num in num_residues:
            id = num_residues.index(res_num)
            net_interaction[id] = interactions[j]
    return (name, net_interaction)

selected_residues, selected_interactions, full_interaction = build_vector(Surfaces_results)

if sum(selected_interactions) != 0:
    name, vector = order_vector(selected_residues, selected_interactions, pdb_file, chains_AB)
    print ('Vector of interactions DONE')
    print ('Total interaction:', full_interaction)
else:
    print ('No interaction')

# Step 5: Calculate reference point of interaction

with open('reference_residues.txt', 'r') as file:
    contents = file.read()
    reference_residue_positions = eval(contents)

def weighted_average(positions, values):
    X_total = 0
    Y_total = 0
    Z_total = 0
    for i in range(len(positions)):
        X_total = X_total + positions[i][0]*values[i]
        Y_total = Y_total + positions[i][1]*values[i]
        Z_total = Z_total + positions[i][2]*values[i]
    average_posit = [X_total/sum(values), Y_total/sum(values), Z_total/sum(values)]
    return (average_posit)

print ('Calculating reference point...')
point1 = weighted_average(reference_residue_positions, vector)
print ('Reference point', point1)

# Step 6: Calculate distances to the reference points of each epitope

with open('reference_epitopes.txt', 'r') as file:
    contents = file.read()
    reference_epitope_positions = eval(contents)

all_distances = []

for point2 in reference_epitope_positions:
    distance = ((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2 + (point2[2] - point1[2])**2)**0.5
    all_distances.append(distance)

print ('Assigning epitope...')
id = all_distances.index(min(all_distances))
print ('Closest Epitope:', id + 1)
print ('Distance to reference:', min(all_distances))