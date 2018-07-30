#!/usr/bin/python2.7
# -*-coding:Utf-8 -*
import re
import sys
annotation_file=sys.argv[1]			# first arg: path or name of blast2go file
colnum_of_gene=int(sys.argv[2])		# second arg: gene column number
colnum_of_GOID=int(sys.argv[3])		# third arg: GO ID column number
back="\n"
tab="\t"
nb_gene=0
annotated_genes_MF=0
annotated_genes_BP=0
annotated_genes_CC=0
annotated_genes_MF_BP=0
annotated_genes_MF_CC=0
annotated_genes_BP_CC=0
totally_annotated=0
unannotated_genes=0
with open("annot_BP.tsv", "w") as GOBP:					# file creation with "Biological Process" annotation
	with open("annot_MF.tsv", "w") as GOMF:				# file creation with "Molecular Function" annotation
		with open("annot_CC.tsv", "w") as GOCC:			# file creation with "Cellular Component" annotation
			with open(annotation_file, "r") as data:	# opening the blast2go file
				for line in data:
					nb_gene+=1														# gene count
					MF=False
					BP=False
					CC=False
					result=""
					splitted_line=line.split("\t")									# line splitting with tabulation 
					gene_name = splitted_line[colnum_of_gene-1]						# gene name recovery
					go_id_list = splitted_line[colnum_of_GOID-1].split(";")			# GOterms recovery
					##print(go_id_list)
					if len(go_id_list) > 1:											# if there are several GOterms 
					 	for n in range(0,len(go_id_list)):
					 		go=str(go_id_list[n])
					 		go=go.replace("_", "")									# removal "_"
					 		go=go.replace(" ", "")									# removal " "
					 		go=go.replace('"','')									# removal '"'
					 		
					 		if go[0:2]=="F:":										# if GOterm takes part of "Molecular Function" domain
					 			go_id=go.replace("F:", "")							# code domain removal
					 			result = "%s%s%s%s" % (gene_name, tab, go_id, back)	# string values concatenation
					 			GOMF.write(result)									# writting of result line into "Molecular Function" file
					 			MF=True 
					 			#print(result)
					 			
					 		elif go[0:2]=="P:":										# if GOterm takes part of "Biological Process" domain
					 			go_id=go.replace("P:", "")							# code domain removal
					 			result = "%s%s%s%s" % (gene_name, tab, go_id, back)	# string values concatenation
					 			GOBP.write(result)									# writting of result line into "Biological Process" file
					 			BP=True
					 			#print(result)
					 			
					 		elif go[0:2]=="C:":										# if GOterm takes part of "Cellular Component" domain
					 			go_id=go.replace("C:", "")							# code domain removal
					 			result = "%s%s%s%s" % (gene_name, tab, go_id, back)	# string values concatenation
					 			GOCC.write(result)									# writting of result line into "Cellular Component" file
					 			CC=True
					 			#print(result)
					 				
					else :															# if there is an unique GOterm
					 	go = splitted_line[colnum_of_GOID-1]
					 	if go == "NA":
					 		unannotated_genes+=1
					 		#result = "%s%s%s" % (gene_name, tab, "No GOterm")
					 		#print(result)
					 	else:
					 		go=go.replace("_", "")									# removal "_"
					 		go=go.replace(" ", "")									# removal ""
					 		go=go.replace('"','')									# removal '"'
					 		
					 		if go[0:2]=="F:":										# if GOterm takes part of "Molecular Function" domain
					 			go_id=go.replace("F:", "")							# code domain removal
					 			result = "%s%s%s%s" % (gene_name, tab, go_id, back)	# string values concatenation
					 			GOMF.write(result)									# writting of result line into "Molecular Function" file
					 			MF=True
					 			#print(result)
					 			
					 		elif go[0:2]=="P:":										# if GOterm takes part of "Biological Process" domain
					 			go_id=go.replace("P:", "")							# code domain removal
					 			result = "%s%s%s%s" % (gene_name, tab, go_id, back)	# string values concatenation
					 			GOBP.write(result)									# writting of result line into "Biological Process" file
					 			BP=True
					 			#print(result)
					 			
					 		elif go[0:2]=="C:":										# if GOterm takes part of "Cellular Component" domain
					 			go_id=go.replace("C:", "")							# code domain removal
					 			result = "%s%s%s%s" % (gene_name, tab, go_id, back)	# string values concatenation
					 			GOCC.write(result)									# writting of result line into "Cellular Component" file
					 			CC=True
					 			#print(result)
					if MF==True:												#
						annotated_genes_MF+=1									#
					if BP==True:												#
						annotated_genes_BP+=1									#
					if CC==True:												#
						annotated_genes_CC+=1									#
					if MF==True and BP==True:									# gene counts for each GO domains
						annotated_genes_MF_BP+=1								#
					if MF==True and CC==True:									#
						annotated_genes_MF_CC+=1								#
					if BP==True and CC==True:									#
						annotated_genes_BP_CC+=1								#
					if MF==True and BP==True and CC==True:						#
						totally_annotated+=1									#
s1 = "%s%s" % ("number of annotated genes MF: ", annotated_genes_MF)
s2 = "%s%s" % ("number of annotated genes BP: ", annotated_genes_BP)
s3 = "%s%s" % ("number of annotated genes CC: ", annotated_genes_CC)
s4 = "%s%s" % ("number of annotated genes MF + BP: ", annotated_genes_MF_BP)
s5 = "%s%s" % ("number of annotated genes MF + CC: ", annotated_genes_MF_CC)
s6 = "%s%s" % ("number of annotated genes BP + CC: ", annotated_genes_BP_CC)
s7 = "%s%s" % ("number of totally annotated genes: ", totally_annotated)
s8 = "%s%s" % ("number of unannotated genes: ", unannotated_genes)
s9 = "%s%s" % ("total number of genes: ", nb_gene)
print(s1)
print(s2)
print(s3)
print(s4)
print(s5)
print(s6)
print(s7)
print(s8)
print(s9)