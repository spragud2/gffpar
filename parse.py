import numpy as np 
import collections

from collections import defaultdict


class GFF(object):

    '''
    A GFF3 PARSER OBJECT DESIGNED TO CREATE
    A HIERARCHAL STRUCTURE WITH PARENTAL GENES
    CONTAINING A DICTIONARY OF ALL 'CHILD' FEATURES
    INCLUDING TRANSCRIPTS,EXONS,CDS,ETC

    ATTRIBUTES:
    gff_file: a file object of a specified gff3 file path
    species: the species belonging to the gff3 file
    version: the annotation version of the file (for example: gencode_m18)


    '''

    def __init__(self,gff_file,species,version):    
        self.gff_file = open(gff_file)
        self.species = species
        self.version = version
        self.par_map = defaultdict(list)
        self.chrom_map = defaultdict(list)

    def format_gff(self):
        # Convert the gff3 file to a 
        # list of lists with elements seperated by \t and
        # split attributes element by ; character

        # I: self
        # O: list of lists containing seperated elements of gff3 file
        lines = [line.strip() for line in self.gff_file]
        split_tabs = list(map(lambda i: i.split('\t'),lines ))
        split_attributes = list(map(lambda i: i[:7]+i[8].split(';'),split_tabs))

        return split_attributes

    def parents(self,formatted_gff):
        # return the 'parent' or 'canonical' gene 

        # I: list of lists (formatted_gff file)
        # O: list of lists containing only 'parent' genes
        return [line[:11] for line in formatted_gff if not any('Parent=' in i for i in line)]

    def children(self,formatted_gff):
        # return the 'children' of 'parent' genes above
        # including transcripts,exons,cds 

        # I: list of lists (formatted_gff file)
        # O: list of lists containing only transcripts,exons,etc
        return [line[:17] for line in formatted_gff if any('Parent=' in i for i in line)] 

    def idxer(self,opt = None):
        # convert the ordereddict to a regular dictionary (so list indexing is not necessary)
        # and add key:value pairs for each child of a parent gene

        # I: self,optional ordereddict
        # O: dict
        if opt is None:
            opt = self.par_map
        d = dict(opt)
        for k in d:
            x = list(range(len(d[k])))
            d[k] = dict(zip(x,d[k]))
        return d

    def parser(self,formatted_gff):
        parents = self.parents(formatted_gff)
        children = self.children(formatted_gff)
        for i,line in enumerate(parents):
            line = [i.split('=')[1] if '=' in i else i for i in line]
            chro,src,typ,start,stop,_,strand,pid,gene_id,gene_type,gene_name = line
            self.par_map[pid].append({'PARENT':True,'chr':chro,'source':src,'type':typ,
                'start':start,'stop':stop,'strand':strand,
                'gene_id':gene_id,'gene_type':gene_type,
                'gene_name':gene_name,'seq_range':range(int(start),int(stop)),'peaks':0})

        for i,line in enumerate(children):
            line = [i.split('=')[1] if '=' in i else i for i in line]
            chro,src,typ,start,stop,_,strand,pid,child_of,gene_id,transcript_id,gene_type,gene_name,transcript_type, transcript_name,*_ = line
            self.par_map[child_of].append({'PARENT':False,'chr':chro,'source':src,'type':typ,
                'start':start,'stop':stop,'strand':strand,'ID':pid[3:],
                'child_of':child_of,'gene_id':gene_id,'transcript_id':transcript_id,
                'gene_type':gene_type,'gene_name':gene_name,'transcript_type':transcript_type,
                'transcript_name':transcript_name,'seq_range':range(int(start),int(stop)),'peaks':0})

        return self.idxer()

    def get_genes_transcript001(self,par_dict):
        of_interest = defaultdict(list)
        for k in par_dict:
            for i in par_dict[k]:
                if par_dict[k][i]['PARENT'] or '001' in par_dict[k][i]['transcript_name']:
                    of_interest[k].append(par_dict[k][i])

        return self.idxer(of_interest)

    def rna_ip_overlap(self):
        pass 

    def chrom_mapper(self,of_interest):
        for k in of_interest:
            loc = of_interest[k][0]['chr']
            self.chrom_map[loc].append({k:of_interest[k]})
            
        self.chrom_map = dict(self.chrom_map)
        gene_map = collections.defaultdict(dict)
        for k in self.chrom_map:
            for i,val in enumerate(self.chrom_map[k]):
                gene_map[list(val.keys())[0]] = self.chrom_map[k][i][list(val.keys())[0]]
            self.chrom_map[k] = dict(gene_map) 

        return self.chrom_map



