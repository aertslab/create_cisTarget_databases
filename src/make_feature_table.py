from commands import getoutput
import numpy as np
from collections import defaultdict
import pandas as pd
from collections import OrderedDict
import optparse
import sys

fmt = optparse.IndentedHelpFormatter(indent_increment=2, max_help_position=9, width=79, short_first=1)
usage = "Usage: %prog [options]"
parser = optparse.OptionParser(usage = usage, version = "%prog v1.0", formatter = fmt)

parser.add_option("-o", "--feature_order", action = "store", type = "string", dest = "feature_order_path", help = 'Name of the TF to process')
parser.add_option("-f", "--fasta_file", action = "store", type = "string", dest = "fasta_file", help = 'Path to fasta file with DNA sequences of the samples')
parser.add_option("-s", "--save_path", action = "store", type = "string", dest = "path_to_save_results", help = 'Path to save results')
parser.add_option("-m", "--motifs_path", action = "store", type = "string", dest = "motifs_path", help = 'Path to folder with singletons')
parser.add_option("-c", "--cbust_path", action = "store", type = "string", dest = "cbust_path", help = 'Path to cluster-buster tool')
(options, args) = parser.parse_args()


# Check if we have an expression matrix filea FASTA or twobit file is given as input.
if ( (options.feature_order_path is None) or (options.fasta_file is None) or (options.path_to_save_results is None)  \
             or (options.motifs_path is None) or (options.cbust_path is None)):
    parser.print_help()
    print >> sys.stderr, '\nERROR: minimum required options not satisfied:\n'
    sys.exit(1)

PATH_TO_FEATURE_ORDER = options.feature_order_path
PATH_TO_FASTA = options.fasta_file
PATH_TO_SAVE = options.path_to_save_results
PATH_TO_FOLDER_WITH_SINGLETONS = options.motifs_path
PATH_TO_CBUST = options.cbust_path


def readFetureOrder(path_to_file):
    FeatureOrder_list=[]
    with open(path_to_file) as my_file:
        for line in my_file:
            motif_name = line.split()[0]
            FeatureOrder_list.append(motif_name)
    return FeatureOrder_list


def read_fasta_IDs(path_to_file):
    fasta_ID_dict=OrderedDict()
    with open(path_to_file, 'r') as fh:
        for line in fh:
            if line[0] == '>':
                regionID=line[1:].rstrip('\n')
                fasta_ID_dict[regionID] = 1

    ids = fasta_ID_dict.keys()
    return ids


def merge2oneFileSingletonScanWG(Text, oldText):

    """
    Input:
    # Score	Start	End	Sequence	transfac_pro-M01238
    23.2	436	1285	chrX-reg7487	47.1
    """


    """
    Output is in format:
    chr1    247279147    247282336    transfac_pro-M00801    136.0
    """
    Merged_file_RAM = ''
    ####Parsing of stdout####
    Text = Text.split("\n")
    motif_name = ''
    for line in Text:
        chr = ''
        start = ''
        end = ''
        crm_score = ''
        if line.startswith("#"):
            tmp_line = line.split()
            motif_name = tmp_line[5]
        line = line.split()
        # if len(line) == 5 and str(line[3]).startswith('chr'):
        if len(line) == 5:
            try:
                crm_score = str(float(line[0]))
            except:
                continue
            chr = str(line[3])   ###chr:start-end
            start = str(line[1])
            end = str(line[2])
            regionID = line[3]
        if chr != '' and str(start) != '' and str(end) != '' and motif_name != '' and str(crm_score) != '':
            Merged_file_RAM += regionID + "\t" + start  + "\t" + end + "\t" + motif_name + "\t" + str(crm_score) + "\n"
    ###---Merge perocessed results with existing---###
    oldText = oldText + Merged_file_RAM
    return oldText


def run_cbust_scan(singletone_list, path_to_fasta):
    ###---Variable to save results of scanning---###
    initScanningData=''
    for motif in singletone_list:
        Cbustscann_singl = getoutput(PATH_TO_CBUST + " " + '-c 0 -m 0 -f 3' + " " + PATH_TO_FOLDER_WITH_SINGLETONS + motif + "  " + path_to_fasta)
        initScanningData = merge2oneFileSingletonScanWG(Cbustscann_singl,initScanningData)
    return initScanningData



def CRMmaxScore_Homo_Cluster(Merged_file_RAM, motif_order, regionOrder):


    """
    input file should have the format:
    chr16:262550-262750    262678    262697    transfac_pro-M01652    0.747
    chr16:262550-262750    262678    262697    transfac_pro-M01652    0.747
    chr16:262550-262750    262678    262697    transfac_pro-M01652    0.747
    """

    MaxScoredict = defaultdict(dict)

    ###---put regions to the dictionary:
    reg_idx=dict()
    i=0
    for region_id in regionOrder:
        reg_idx[region_id] = i
        i=i+1
    motif_idx=dict()
    i=0
    for motif in motif_order:
        ###----remove cb extension----###
        motif_idx[motif[:-3]]=i
        i=i+1
    ###-------------------------------###
    ####-----Read merged results for cbust scoring------####
    Merged_file_RAM = Merged_file_RAM.split("\n")

    for line in Merged_file_RAM:
        line=line.split()
        if len(line) ==5:
            regionID=line[0]
            motifID=line[3]
            score=float(line[4])
            if regionID not in MaxScoredict:
                MaxScoredict[regionID][motifID]=score
            else:
                if motifID not in MaxScoredict[regionID]:
                    MaxScoredict[regionID][motifID]=score
                else:
                    if MaxScoredict[regionID][motifID]<score:
                        MaxScoredict[regionID][motifID]=score
    ###---initialize the matrix with values
    featureMatrix_POS = np.zeros((len(regionOrder), len(motif_order)))
    row_idx=0
    col_idx=0
    for region in reg_idx:
        row_idx=reg_idx[region]
        for motif in motif_idx:
            col_idx=motif_idx[motif]
            if region in MaxScoredict and motif in MaxScoredict[region]:
                score=MaxScoredict[region][motif]
                featureMatrix_POS[row_idx,col_idx] = score
    return featureMatrix_POS


if __name__ == '__main__':
    ###---Run cbust scoring---###
    FeatureOrder_list = readFetureOrder( PATH_TO_FEATURE_ORDER )
    region_order_ref = read_fasta_IDs(PATH_TO_FASTA)
    cbust_allmotifs_merged_results_ref = run_cbust_scan(FeatureOrder_list, PATH_TO_FASTA)
    featureMatrix_region_ref = CRMmaxScore_Homo_Cluster(cbust_allmotifs_merged_results_ref, FeatureOrder_list, region_order_ref)

    ###---Save FT in pandas format----###
    FM_data_frame = pd.DataFrame(featureMatrix_region_ref, index=region_order_ref, columns=FeatureOrder_list)
    print "Saving feature-table"
    FM_data_frame.to_csv(PATH_TO_SAVE , sep='\t')
    print "Done"