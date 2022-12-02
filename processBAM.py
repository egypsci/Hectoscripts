# SYNOPSIS
# When mapping RNAseq reads to a *transcriptome* using STAR, how does one identify unique mapping reads? 
# A read mapping to two different transcripts of the same gene would be considered a multimapper, which is not accurate.
# For the salamander genome project, using the TRINITY assembled transcriptome ID and counting based on the ‘gene’ level ID instead of the ‘isoform’ level ID is not helpful, 
# since transcripts with different ‘gene’ IDs can have the same exact annotation. 
# In the salamander genome project, protein coding transcripts from the assembled transcriptome are organized into “orthology groups” (OGs) 
# based on orthology with protein coding genes from several model organisms. Each OG can consist of several transcripts belonging to different TRINITY ´genes´. 
# This python script processes BAM files and asks if a read is mapping to various transcripts within the same OG or not. 
# If it is, then the mapping is considered as unique even if it mapped to different TRINITY genes.
# Note that BBL or bbl is used in the script to refer to OGs since I named the annotation framework BABEL (since it allows for easy comparative genomics, i.e. multiple organisms or languages)

import os
import pysam
from collections import defaultdict

# annotate ids and make dictionaries
f1 = open('/Users/ahmele/Documents/Postdoc_KI/1-PwGenome/datasets/Clean/Pw_annotationX','r')  # A file with OGs and their annotation (species: Pleurodeles waltl)
f2 = open('/Users/ahmele/Documents/Postdoc_KI/1-PwGenome/datasets/Pw_RNAME','r') # A file corresponding the START index ID (RNAME, numeric) with its TRINITY transcript (characters).

dict_RNAME = defaultdict(str) # RNAME as key, Annotation as value
dict_BBLAnn = defaultdict(str) # Orthology Group as key and Annotation as value
dict_AnnBBL = defaultdict(list)# Annotation as key and Orthology Groups as values

for a in f1:
  a		= a.strip()                # For each entry,
  bbl = a.split(' ')[0]          # Get orthology group,
	ann	= a.split(' ')[1]+'_'+bbl  # and annotation. Tag annotation with orthology group ID to make things easier downstream,
	dict_BBLAnn[bbl]=ann


# Reverse the dictionary so that Annotation is key and Orthology Groups are values
for k,v in dict_BBLAnn.items():
	dict_AnnBBL[v].append(k)


i=0
for a in f2:
	a		= a.strip()
	if a.startswith('ER'):
		dict_RNAME[str(i)]='ERCC::'+a
	elif a.startswith('T'):
                dict_RNAME[str(i)]=a+'::'+dict_BBLAnn[a.split(' ')[0]]
	elif a.split('_')[0] in dict_BBLAnn:
		dict_RNAME[str(i)]=a.split('_')[0]+'::'+dict_BBLAnn[a.split('_')[0]]
	else:
		dict_RNAME[str(i)]=a.split('_')[0]+'::UNKNOWN'
	i+=1

# Process BAM files in folder
for filename in os.listdir('/bamfiles/SS2_19_152'):
	if filename.endswith('bam'):
		well = filename.split('.')[0] # Get the name of the well/cell
        
		# open bam file
		samfile = pysam.AlignmentFile('/bamfiles/SS2_19_152/'+well+'.Aligned.sortedByCoord.out.bam', 'rb')

		dict_uniqueD = defaultdict(int) # Create a dictionary for unique reads
		dict_temp1 = defaultdict(list)  # Create a temporary dictionary for putative multimappers (more than 1 and less than 41 targets)


		for read in samfile:
			# first process the cigar tuples
			M = M1 = M2 = S = S1 = S2 = s = m =0
			cigar = read.cigar
			for o in cigar:
				if o[0] == 0:  # some CIGARs have two Ms with an insertion btw. Add Ms
					if m == 0:
						M1 = o[1]
					if m == 1:
						M2 = o[1]
				m+=1
			if o[0] == 4:  # some CIGARs have S before and after M, so choose longer one
				if s == 0:
					S1 = o[1]
				if s == 1:
					S2 = o[1]
				s+=1
			M = M1+M2
			S = max([S1,S2])
    
			read = str(read)
			qname= read.split('\t')[0]
			rname= read.split('\t')[2]
			cigar= read.split('\t')[5] # write over cigar the string version
			star = read.split('\t')[-1][1:-1]
			star = star.split('), ')
			NH = star[0].split(', ')[1]

      # Get the OG and annotation for the target that the read mapped to
			bbl    = dict_RNAME[rname].split('::')[0]      # get orthology group
			annotation = dict_RNAME[rname].split('::')[1]  # get Annotation
	
			if int(NH) == 1 and int(M) >= 40:  # unique read (class A)
				dict_uniqueD[annotation]+=1

			if int(NH) >  1 and int(M) >= 40:  # criteria for putative multimapper
				dict_temp1[qname].append(annotation)
		
    # Check and see if the targets of the putative multimappers all have the same annotation. 
    # If the answer is yes, then upgrade the read to "Unique"
    # Check and see if the targets of the putative multimappers all have the same annotation.
    
		for read,annotation in dict_temp1.items():
			if len(set(annotation)) == 1:  # then read maps to the same annotation and it's unique (class B)
				dict_uniqueD[annotation[0]]+=1

    # Output the tally of unique reads (Class A and Class B) for each cell.
		oD = open('/Volumes/Corazon 1/bamfilecounts/'+well+'_uniqueD','w')

		for k,v in dict_uniqueD.items():
			oD.write(k)
			oD.write('\t')
			oD.write(str(v))
			oD.write('\n')
		oD.close()


