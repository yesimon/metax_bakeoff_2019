MEGAN = config.get('MEGAN', 'MEGAN')
MEGAN_LICENSE = config.get('MEGAN_LICENSE')
MEGAN_SPECIES_TEMPLATE = config.get('MEGAN_SPECIES_TEMPLATE', 'megan_species_template.txt')
MEGAN_SPECIES_FILTERED_TEMPLATE = config.get(
    'MEGAN_SPECIES_FILTERED_TEMPLATE', 'megan_species_filtered_template.txt')

# rule megan_all:
#     input: expand('results/{sample}.diamond.megan', sample=samples)

# rule megan:
#   input: 'data/{seq}.diamond.{db}.daa'
#   output: 'reports/{seq}.diamond_megan.{db}.tsv'
#   resources: megan=1
#   params: commands_file='{seq}.megan_commands.txt'
#   shadow: 'shallow'
#   run:
#       with open(MEGAN_SPECIES_TEMPLATE) as f:
#           commands = f.read().format(blast_file=input, megan_file='megan_file.rma', species_file=output)
#       with open(params.commands_file, 'w') as f:
#           f.write(commands)
#       shell('env INSTALL4J_ADD_VM_PARAMS=-Xmx30000M xvfb-run -a -n 100 {MEGAN} -g -L {MEGAN_LICENSE} -c {commands_file} | tee megan_run.log',
#             commands_file=params.commands_file)

# rule megan_filtered_all:
#   input: expand('reports/{sample}.diamond.filtered.megan', sample=samples)

# rule megan_filtered:
#   input: 'data/{seq}.m8.gz'
#   output: 'reports/{seq}.diamond.filtered.megan'
#   params: commands_file='{seq}.megan_commands.txt'
#   shadow: 'shallow'
#   run:
#       shell('echo {input}')
#       with open(MEGAN_SPECIES_FILTERED_TEMPLATE) as f:
#           commands = f.read().format(blast_file=input, megan_file='megan_file.rma', species_file=output)
#       with open(params.commands_file, 'w') as f:
#           f.write(commands)
#       shell('env INSTALL4J_ADD_VM_PARAMS=-Xmx30000M xvfb-run -a -n 100 {MEGAN} -g -L {MEGAN_LICENSE} -c {commands_file} | tee megan_run.log',
#             commands_file=params.commands_file)
