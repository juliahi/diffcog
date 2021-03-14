
import glob
#SAMPLES_CONTROL = ['6685_04-06-2015', '6690_04-06-2015', '6695_04-06-2015', '6704_04-06-2015']
#SAMPLES_TREATED = ['6683_16-06-2015', '6685_16-06-2015', '6690_16-06-2015',
#             '6695_16-06-2015', '6704_16-06-2015']
SAMPLES_CONTROL = ['6685_04-06-2015', '6690_04-06-2015', ]
SAMPLES_TREATED = ['6683_16-06-2015', '6685_16-06-2015', '6690_16-06-2015']


#SAMPLES=glob.glob(INDIR+"/*_depl_1.fq.gz")+glob.glob(INDIR+"/*_depl_2.fq.gz")
#SUFFIX1 = "_depl_1.fq.gz"
#SUFFIX2 = "_depl_2.fq.gz"
SUFFIX1 = "_1.fq"
SUFFIX2 = "_2.fq"

#OUTDIR="/mnt/chr7/data/julia"
#INDIR="/mnt/chr4/mikrobiomy-2/Wyniki_sekwencjonowania/demultiplexed"
#INDIR="/home/julia/chr4/mikrobiomy-2/Wyniki_sekwencjonowania/demultiplexed"
INDIR="/home/julia/studia/test_data"
OUTDIR="/home/julia/studia/heurystyki_test"


#SGADIR="/usr/bin/"
SGADIR=""



# SGA settings:
QUALFILTER=5
CORRECTK=21
OVERLAP=31

# other settings:
MIN_FC = 2
threads=2
READLEN=100

# simplification and filtering
#MINCONTIGLENGTH=200
MIN_LENGTH = 200
DEADENDS_REMOVE_ROUNDS = 1
DEADENDS_MIN_LENGTH = 200

#
#
#
# ### parameters of experiment and file paths
# conditions = {'6683_16-06-2015': 1, '6685_04-06-2015': 0, '6685_16-06-2015': 1,
#              '6690_04-06-2015': 0, '6690_16-06-2015': 1, '6695_04-06-2015': 0,
#              '6695_16-06-2015': 1, '6704_04-06-2015': 0, '6704_16-06-2015': 1}
# test_name = 'sga_test_full_notrim_paired_reversed'
# folder = "/mnt/chr7/data/julia/" + test_name
# suf = ".preprocessed_qf5.ec.filter.pass.rmdup"
#
# K=31
# graph_filename = folder + "/merged" + suf + "_" + str(K) + ".asqg"
#
# # simplification and filtering
# MIN_LENGTH = 200
# DEADENDS_REMOVE_ROUNDS = 1
# DEADENDS_MIN_LENGTH = 200
#
# OUTPUTDIR="/home/julia/mikrobiomy_results/"
# outname = "mysimplified%d_deadends%d_minlen%d"%(DEADENDS_REMOVE_ROUNDS, DEADENDS_MIN_LENGTH, MIN_LENGTH)
# my_simplified_graph_filename = "%s%s.asqg" % (OUTPUTDIR, outname)
# my_simplified_graph_pickle = "%s%s.pickle" % (OUTPUTDIR, outname)
# renamed_nodenames = "%s%s_nodes_renamed.tsv" % (OUTPUTDIR, outname)
# my_renamed_graph_pickle = "%s%s_renamed.pickle" % (OUTPUTDIR, outname)
# my_renamed_graph_filename = "%s%s_renamed.asqg" % (OUTPUTDIR, outname)
# simplify_stats_filename = "%s%s_stats.pickle" % (OUTPUTDIR, outname)
#
# #import seaborn as sns
# #color_palette = sns.color_palette('muted6')
# # labels = ["longest", "longestfc", "bestfc", "sga", "megahit"]
# # colors = {"longest": color_palette[0], "longestfc": color_palette[5], "bestfc": color_palette[2],
# #           "sga": color_palette[4], "megahit": color_palette[1]}
#
# # for heuristics
# MIN_FC = 2
# longest_fasta = "%sheuristics_%s/longest_filtered_fc%d.fa" % (OUTPUTDIR, outname, MIN_FC)
# longestfc_fasta = "%sheuristics_%s/longestfc_filtered_fc%d.fa" % (OUTPUTDIR,  outname, MIN_FC)
# bestfc_fasta = "%sheuristics_%s/bestfc_filtered_fc%d.fa" % (OUTPUTDIR,  outname, MIN_FC)
# sgafile = "%ssga_scaffold_%s/sga-scaffolds.fa" % (OUTPUTDIR, outname)
# sga_fasta = "%sheuristics_%s/sga_filtered_fc%d.fa" % (OUTPUTDIR, outname, MIN_FC)
#
# NPAIRS = 90582474
