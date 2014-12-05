from pymol import cmd
from Bio import SeqIO, pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from collections import OrderedDict
import sys


matrix = matlist.blosum62
#gp120 ends at 511
_hxbc_seq_gp120 = "MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPNPQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFNISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTWSTEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWRSELYKYKVVKIEPLGVAPTKAKRRVVQREKR"

#gp41 starts at 512
_hxbc_seq_gp41 = "AVGIGALFLGFLGAAGSTMGAASMTLTVQARQLLSGIVQQQNNLLRAIEAQQHLLQLTVWGIKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNASWSNKSLEQIWNHTTWMEWDREINNYTSLIHSLIEESQNQQEKNEQELLELDKWASLWNWFNITNWLWYIKLFIMIVGGLVGLRIVFAVLSIVNRVRQGYSPLSFQTHLPTPRGPDRPEGIEEEGGERDRDRSIRLVNGSLALIWDDLRSLCLFSYHRLRDLLLIVTRIVELLGRRGWEALKYWWNLLQYWSQELKNSAVSLLNATAIAVAEGTDRVIEVVQGACRAIRHIPRRIRQGLERILL"

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
          'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
          'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
          'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
          'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
          'TYS': 'Y'}


class HXBCAlignment():
    def __init__(self,alignement_object,interest_numbering):
        self.alignement_object = alignement_object
        self.best_alignment = alignement_object[0]
        self.hxbc_part = self.best_alignment[1]
        self.interest_part = self.best_alignment[0]
        self.score = self.best_alignment[2]

        self.gp120_epitopes = [{'V1':range(131,157)},
                                  {'V2':range(157,197)},
                                  {'V3':range(296,331)},
                                  {'V4':range(385,419)},
                                  {'V5':range(461,472)},
                                  {'CD4-binding_loop':range(364,375)},
                                  {'b12_epitope':[257,280,281,282,364,365,366,
                                                  367,368,369,370,371,372,373,
                                                  384,385,386,417,418,519,530,
                                                  431,432,455,456,457,469,472,
                                                  473,474]}]

        self.gp41_epitopes = [{'Fusion_peptide':range(512,533)},
                              {'MPER':range(662,684)}]

        self.interest_numbering = interest_numbering

    def get_hxbc2(self):
        return self.hxbc_part

    def get_interest_part(self):
        return self.interest_part

    def get_score(self):
        return self.score

    def get_alignment(self,epitope,epi="gp120"):
        name = epitope.keys()[0]
        location = epitope.values()[0]
        self.interest_epitope[name] = set()
        interest_counter = 0
        if epi == "gp41":
            hx_counter = 511
        else:
            hx_counter = 0

        for align_hx_index, align_interest_index in zip(
            range(0,len(self.hxbc_part)),
            range(0,len(self.interest_part))):
            
            if self.hxbc_part[align_hx_index] != "-":
                hx_counter += 1

            if self.interest_part[align_interest_index] != "-":
                interest_counter += 1

            if hx_counter in location:
                self.interest_epitope[name].add(self.interest_numbering[interest_counter])


    def match_epitopes(self, epi="gp120"):
        self.interest_epitope = {}
        if epi == "gp120":
            for epitope in self.gp120_epitopes:
                self.get_alignment(epitope)
        if epi == "gp41":
            for epitope in self.gp41_epitopes:
                self.get_alignment(epitope,epi="gp41")

    def make_objects(self,chain):
        for epitope in self.interest_epitope:
             if self.interest_epitope[epitope]:
                name = epitope
                selection = ",".join(str(x) for x in self.interest_epitope[epitope])
                cmd.select("{}_{}".format(chain.lower(),name),"chain {} and resi {}".format(chain,selection))
                cmd.group("chain_{}".format(chain),"{}_*".format(chain.lower()))


def align(key_value_pairing,chain,ascutoff):
    check_string = "".join(key_value_pairing.values())
    pw_gp120 = pairwise2.align.globaldx(check_string,_hxbc_seq_gp120.strip(),matrix)
    pw_gp41 = pairwise2.align.globaldx(check_string,_hxbc_seq_gp41.strip(),matrix)
    
    

    if pw_gp120:
        alignment = HXBCAlignment(pw_gp120,key_value_pairing.keys())
        if alignment.get_score() > ascutoff:
            alignment.match_epitopes(epi="gp120")
            alignment.make_objects(chain)
            return 
    
    if pw_gp41:
        alignment = HXBCAlignment(pw_gp41,key_value_pairing.keys())
     
        if alignment.get_score() > 550:
            alignment.match_epitopes(epi="gp41")
            alignment.make_objects(chain)

def get_sane_pairing(pairing):
    l = {}
    for position in pairing:
        try:
            position_n = int(position[0])
        except ValueError:
            position_n = position[0]
        try:
            position_l = longer_names[position[1]]
        except KeyError:
            print 'Warning: {} not found'.format(position[1])
            continue
        l[position_n] = position_l
    return l

def select_gp160(alignment_score_cutoff = 800):
    '''
        DESCRIPTION

            makes alignement groups for gp160 groups
            
        ARGUMENTS
            
            alignment_score_cutoff - default 1000, the alignment score that all the chains must pass
            in order to be considered as a gp160 molecule. If you get no alignments, try to lower
            this number. 
        USAGE

            select_gp160 , [alignnment=int]
    '''
    all_chains = cmd.get_chains()
    for chain in all_chains:
    #chain = "B"
        myspace = {'myset':set()}
        cmd.iterate('chain {}'.format(chain), 'myset.add((resi,resn))',space=myspace)
        chain_value_pair = myspace['myset']
        chain_value_pair = get_sane_pairing(chain_value_pair)
        align(chain_value_pair,chain,ascutoff=alignment_score_cutoff)
cmd.extend('select_gp160',select_gp160)