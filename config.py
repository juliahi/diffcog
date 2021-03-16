
# sample prefixes:
#SAMPLES_CONTROL = ['6685_04-06-2015', '6690_04-06-2015', '6695_04-06-2015', '6704_04-06-2015']
#SAMPLES_TREATED = ['6683_16-06-2015', '6685_16-06-2015', '6690_16-06-2015', '6695_16-06-2015', '6704_16-06-2015']
SAMPLES_CONTROL = ['6685_04-06-2015', '6690_04-06-2015', ]
SAMPLES_TREATED = ['6683_16-06-2015', '6685_16-06-2015', '6690_16-06-2015']

# sample suffix for paired reads
#SUFFIX1 = "_depl_1.fq.gz"
#SUFFIX2 = "_depl_2.fq.gz"
SUFFIX1 = "_1.fq"
SUFFIX2 = "_2.fq"

#OUTDIR="/mnt/chr7/data/julia"
#INDIR="/mnt/chr4/mikrobiomy-2/Wyniki_sekwencjonowania/demultiplexed"
#INDIR="/home/julia/chr4/mikrobiomy-2/Wyniki_sekwencjonowania/demultiplexed"
INDIR="/home/julia/studia/test_data"
OUTDIR="/home/julia/studia/heurystyki_test"

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
