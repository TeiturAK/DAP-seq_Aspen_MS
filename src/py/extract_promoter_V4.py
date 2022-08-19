#!/usr/bin/env python

import sys
import re

'''
Extract promoter regions for genes in maker annotation gff3 file.
Promoter region defines as: num_bases_upstream_promoterEnd upstream from first 5'UTR to first CDS in each gene / num_bases_upstream_promoterEnd upstream from first CDS in each gene.
Complete genes are defined as the whole transcribed region, 5'UTRs, exons, CDSs, 3'UTRs. If a gene is nested, the gene flanking the nested gene is ignored.

IMPORTANT: It is assumed that genes are sorted within contigs in the maker_gff file.
'''

maker_gff = sys.argv[1]
contiglengths = sys.argv[2]
num_bases_upstream_promoterEnd = int(sys.argv[3])
output_filename = sys.argv[4]

#maker_gff = "Potra02_genes_short.gff"
#contiglengths = "contiglengths"
#num_bases_upstream_promoterEnd = 3000
#output_filename = "promoters.gff"

with open(maker_gff) as f:
    data = f.readlines()

#contiglengths is a two columnn tab separated file, first column with
#chromosome/scaffold id, second with the length of that chromosome/scaffold
with open(contiglengths) as f:
    chromosome_sizes = f.readlines()

#Separate lines by 't'
data = [line.strip('\n') for line in data]
data = [line.split('\t') for line in data]
chromosome_sizes = [line.split('\t') for line in chromosome_sizes]

#Create chromosome to chromosome size dictionary
chromosomes = [line[0] for line in chromosome_sizes]
sizes = [line[1] for line in chromosome_sizes]
chromosomesize_dictionary = dict(zip(chromosomes, sizes))

#Remove any rows not containing annotation data (9 columns)
data = [line for line in data if len(line) == 9]

#Save indexes of lines that have line[3] == gene ("type" in gff3).
#Also add index of the last line in the file (needed for complete gene parsing later on).
completegene_line_indexes = [i for i in range(0, len(data)) if data[i][2] == "gene"]
completegene_line_indexes.append((len(data) + 1))

#Evaluate information in complete genes. Get line where the whole gene is defined and all the rows
#up until the next gene using indexes for lines. Store genes with all information in lists.
completegenes_indexes = [[completegene_line_indexes[i], completegene_line_indexes[i + 1]] for i in range(0, (len(completegene_line_indexes)) - 1)]
#completegenes = [data[line[0]:line[1]] for line in completegenes_indexes]


#for gene in completegenes:
#    print (gene)



promoters = []
for i in range(0, len(completegenes_indexes)):
    
    gene_start_index = completegenes_indexes[i][0]
    gene_end_index = completegenes_indexes[i][1]
    
    completegene = data[gene_start_index:gene_end_index]

    #Create copy of the gene row 
    gene = completegene[0][:]
    
    gene_start = int(gene[3])
    gene_end = int(gene[4])
    
    #Create copy of the gene row to be modified for output
    promoter_output_row = gene
    
    gene_id = gene[8]
    gene_name = re.search('ID=(.+?)(;|$)', gene_id).group(1)
    promoter_output_row[8] = gene_name

    #testing to remove this as function is unclear
    #gene[8]

    gene_chromosome = gene[0]
    chromosome_size = int(chromosomesize_dictionary[gene_chromosome])

    strand = gene[6]

    #Extract only CDS rows
    codingfeatures = [[feat[0], feat[2], int(feat[3]), int(feat[4])] for feat in completegene if feat[2] == 'CDS']
    #Extract only 5'UTR rows
    untranslatedfeatures = [[feat[0], feat[2], int(feat[3]), int(feat[4])] for feat in completegene if feat[2] == 'five_prime_UTR']

    #Calculate promoter regions for genes
    #+/- blocks for readability, should probably be a function though
    if strand == '+':
        codingfeatures.sort(key=lambda x: x[2])

        #Promoter end should be the first coordinate of the first CDS in the gene
        promoterEnd_coordinate = codingfeatures[0][2]

        #On the + strand the promoter start coordinate is < end coordinate
        promoter_output_row[4] = str(promoterEnd_coordinate)

        if untranslatedfeatures:
            #Promoter start should be num_bases_upstream_promoterEnd upstream first 5'UTR coordinate
            untranslatedfeatures.sort(key=lambda x: int(x[2]))

            fivePrUTR_coordinate = untranslatedfeatures[0][2]

            #if first 5'UTR is placed after first CDS, then dont use fivePrUTR to calculate promoterStart_coordinate
            if fivePrUTR_coordinate > promoterEnd_coordinate:
                promoterStart_coordinate = promoterEnd_coordinate - num_bases_upstream_promoterEnd
            #Should be if first 5'UTR > first CDS, then dont use that as the start
            else:
                promoterStart_coordinate = fivePrUTR_coordinate - num_bases_upstream_promoterEnd

        else:
            #Promoter start should be num_bases_upstream_promoterEnd upstream first CDS coordinate
            promoterStart_coordinate = promoterEnd_coordinate - num_bases_upstream_promoterEnd

        #Account for contig/chromosome length, promoter region may not exceed the boundaries of the
        #chromosome/scaffold
        if promoterStart_coordinate < 1:
            promoterStart_coordinate = 1
        
        
        #Check end coordinate of upstream gene
        if i == 0:
            pass
        else:
            
            #check if upstream gene overlaps the entire promoter/or the entire gene
            #j = 1
            #ignore = False
            #while not ignore:
                
                
            upstream_gene_start_index = completegenes_indexes[i - 1][0]
            upstream_gene_end_index = completegenes_indexes[i - 1][1]            
                
            upstream_gene = data[upstream_gene_start_index:upstream_gene_end_index]
            
            upstream_gene_chromosome = upstream_gene[0][0]
            
            if upstream_gene_chromosome != gene_chromosome:
                pass
            
            else:
                upstream_gene_start = int(upstream_gene[0][3])
                upstream_gene_end = int(upstream_gene[0][4])                
                
                if upstream_gene_start <= promoterEnd_coordinate and upstream_gene_end >= promoterEnd_coordinate:
                    
                    if i - 2 >= 0:
                    
                        upstream_gene_start_index = completegenes_indexes[i - 2][0]
                        upstream_gene_end_index = completegenes_indexes[i - 2][1]            
                        
                        upstream_gene = data[upstream_gene_start_index:upstream_gene_end_index]
                    
                        upstream_gene_chromosome = upstream_gene[0][0]
                    
                        if upstream_gene_chromosome != gene_chromosome:
                            pass                    
                        else:
                            upstream_gene_start = int(upstream_gene[0][3])
                            upstream_gene_end = int(upstream_gene[0][4])                              
                            if upstream_gene_end > promoterStart_coordinate and upstream_gene_end < promoterEnd_coordinate:
                                promoterStart_coordinate = upstream_gene_end
                            else:
                                pass
                    else:
                        pass
                else:
                    if upstream_gene_end > promoterStart_coordinate and upstream_gene_end < promoterEnd_coordinate:
                        promoterStart_coordinate = upstream_gene_end
                    else:
                        pass                    

        #On the + strand the promoter start coordinate is < end coordinate
        promoter_output_row[3] = str(promoterStart_coordinate)

    elif strand == '-':
        codingfeatures.sort(key=lambda x: x[3], reverse = True)
        

        #Promoter end should be the second coordinate of the first CDS in the gene
        promoterEnd_coordinate = codingfeatures[0][3]

        #On the - strand the promoter start coordinate is > end coordinate
        promoter_output_row[3] = str(promoterEnd_coordinate)

        if untranslatedfeatures:
            #Promoter start should be num_bases_upstream_promoterEnd upstream second 5'UTR coordinate
            untranslatedfeatures.sort(key=lambda x: int(x[3]), reverse = True)

            fivePrUTR_coordinate = untranslatedfeatures[0][3]

            #if first 5'UTR is placed after first CDS (- strand perspective), then dont use fivePrUTR to calculate promoterStart_coordinate
            if fivePrUTR_coordinate < promoterEnd_coordinate:
                promoterStart_coordinate = promoterEnd_coordinate + num_bases_upstream_promoterEnd
            else:
                promoterStart_coordinate = fivePrUTR_coordinate + num_bases_upstream_promoterEnd

        else:
            #Promoter start should be num_bases_upstream_promoterEnd upstream second CDS coordinate
            promoterStart_coordinate = promoterEnd_coordinate + num_bases_upstream_promoterEnd

        #Account for contig/chromosome length, promoter region may not exceed the boundaries of the
        #chromosome/scaffold
        if promoterStart_coordinate > chromosome_size:
            promoterStart_coordinate = chromosome_size
            
        #Check start coordinate of downstream gene
        if i == len(completegenes_indexes) - 1:
            pass
        else:
            
            downstream_gene_start_index = completegenes_indexes[i + 1][0]
            downstream_gene_end_index = completegenes_indexes[i + 1][1]            
            
            downstream_gene = data[downstream_gene_start_index:downstream_gene_end_index]
            
            downstream_gene_chromosome = downstream_gene[0][0]
            
            if downstream_gene_chromosome != gene_chromosome:
                pass
            
            else:
                downstream_gene_start = int(downstream_gene[0][3])
                downstream_gene_end = int(downstream_gene[0][4])
                #if downstream_gene_end > gene_end and downstream_gene_start < gene_start:
                    #continue
                
                if downstream_gene_end <= promoterEnd_coordinate and downstream_gene_start >= promoterEnd_coordinate:
                    
                    if i + 2 <= len(completegenes_indexes): 
                    
                        downstream_gene_start_index = completegenes_indexes[i + 2][0]
                        downstream_gene_end_index = completegenes_indexes[i + 2][1]            
                        
                        downstream_gene = data[downstream_gene_start_index:downstream_gene_end_index]
                    
                        downstream_gene_chromosome = downstream_gene[0][0]
                    
                        if downstream_gene_chromosome != gene_chromosome:
                            pass                    
                        else:
                            downstream_gene_start = int(downstream_gene[0][3])
                            downstream_gene_end = int(downstream_gene[0][4])                            
                            if downstream_gene_start < promoterStart_coordinate and downstream_gene_start > promoterEnd_coordinate:
                                promoterStart_coordinate = downstream_gene_start
                            else:
                                pass   
                    else:
                        pass
                    
                    
                else:
                
                    if downstream_gene_start < promoterStart_coordinate and downstream_gene_start > promoterEnd_coordinate:
                        promoterStart_coordinate = downstream_gene_start
                    else:
                        pass        
        
                
        #On the - strand the promoter start coordinate is > end coordinate
        promoter_output_row[4] = str(promoterStart_coordinate)

    else:
        #Maybe add an exit code because that could be part of a larger problem/wrong input.
        print ("No strand information, check data")
        continue
    promoter_output_row[2] = 'promoter'+str(num_bases_upstream_promoterEnd/1000)+'kb'

    #Control, print out all of the gene information and chromosome length and the computed promoter
    #to check that it is correct
    #for line in gene:
        #print (line)
    #print (chromosome_size)
    #print (promoter)
    #print ('\n')


    promoter_row = "\t".join(promoter_output_row)
    #print (promoter_row)
    promoters.append(promoter_row)


#Write promoters to file in gff3 format
with open(output_filename, 'w') as output:
    output.write('##gff-version\t3\n')
    for promoter in promoters:
        output.write(promoter + '\n')
        #for feature in promoterfeatures:
            ##Do not add tab to last
            #output.write(str(feature) + '\t')
        #output.write('\n')

        #output.write('\n')
        #'\n' included in ID

