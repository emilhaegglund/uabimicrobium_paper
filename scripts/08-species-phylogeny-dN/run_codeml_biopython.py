## Run the codeml using bio python (this alow us to set a working directory useful for parallelisation)

from Bio.Phylo.PAML import codeml
import os
### create the codeml object with input and output

def create_codeml_object(workdir, seqfile, outfile):
    cml = codeml.Codeml()
    cml.alignment = str(seqfile)
    cml.tree = 'empty.tree' # Pairewise dNdS computation do not require a tree file in CODEML however the biopython object must have one
    cml.out_file = str(outfile)
    cml.working_dir = str(workdir) + '/' + str(seqfile).split('aa_msa_converted_to_nt/')[-1].split('.')[0]

    return cml

### set the options for the codeml run

def add_codeml_options(cml):
    cml.set_options(noisy=0)        # Minimum text on the screen
    cml.set_options(verbose=False)  # Set to true to view what is happening on the screen
    cml.set_options(runmode=-2)     # pairwise comparison
    cml.set_options(seqtype=1)      # codons
    cml.set_options(CodonFreq=2)    # calculated from the average nucleotide frequencies at the three codon positions
    cml.set_options(model=1)        # branch model = 1, means one omega ratio per branch
    #cml.set_options(NSsites=[0, 1, 2])      # (model M0 in Yang et al. 2000b)
    cml.set_options(icode=0)        # universial code table in genebank translation table
    cml.set_options(fix_kappa=1)    # fixed kappa
    cml.set_options(kappa=1)        # initial or fixed kappa
    cml.set_options(fix_omega=0)    # estimate omega
    cml.set_options(omega=0.5)      # initial omega value
    return cml

### create a codeml object and add settings:
if not os.path.exists(snakemake.params.outdir):
    os.mkdir(snakemake.params.outdir)
cml = create_codeml_object(snakemake.params.outdir, snakemake.input[0], snakemake.output[0])
cml = add_codeml_options(cml)

### run CODEML
cml.run(verbose=True, parse=False)
