import os 

# sample prefixes:
SAMPLES_CONTROL = ['C2', 'C3', ]
SAMPLES_TREATED = ['T1', 'T2', 'T3']

# sample suffix for paired reads
SUFFIX1 = "_1.fq"
SUFFIX2 = "_2.fq"

INDIR= os.path.dirname(os.getcwd())+"/test_data"
OUTDIR= os.path.dirname(os.getcwd())+"/test_output"


SGADIR=""

# SGA settings:
QUALFILTER=5
CORRECTK=21
OVERLAP=31

# other settings:
threads=4
READLEN=100

# simplification and filtering
MIN_LENGTH = 200
DEADENDS_REMOVE_ROUNDS = 10
DEADENDS_MIN_LENGTH = 200

# heuristics to run:
HEURISTICS = {"longest", "longestfc", "bestfc"}
MIN_FC = 2
MIN_CONTIG_LENGTH=200
