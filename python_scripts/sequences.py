# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 11:13:22 2016

@author: Joana
"""

from Bio import Entrez
Entrez.email = "joana_smar@hotmail.com"


class sequence():

#Objetivo do Script:
#Obter a sequência da proteína com recurso ao NCBI
#Guardar a sequência em formato FASTA
#Criar péptidos para serem analisados pelo netMHCpan
    
    def __init__(self, id_prot = 0, id_gene = 0):
        self.id_prot = id_prot #id da proteína
        self.id_gene = id_gene #id do gene
        self.sequence = "" #sequência da proteína
        self.mut_sequence = "" #sequência da proteína mutada
        self.definition = "" #informação sobre o nome da proteína e respetivo organismo associado
    
    def edit(self, feature):
        #Função que permite editar/processar texto proveniente de ficheiros.
        
        f = str(feature)
        f = f.replace("[", "")
        f = f.replace("]", "")
        f = f.replace("'", "")
        f = f.replace(", ", "\n")
        return f


#Procurar a sequência por ID da proteína:

    def protein_ncbi(self):
        #Função que acede ao NCBI a partir do id de uma proteína, retirando a respetiva sequência.
    
        handle = Entrez.efetch(db = "protein", rettype = "gb", retmode = "xml", id = self.id_prot)
        results = Entrez.read(handle)
        self.sequence = results[0]["GBSeq_sequence"].upper()      
        self.definition = results[0]["GBSeq_definition"]
        #print("Protein: " + str(self.definition[0:4]) + "\nOrganism: " + str(self.edit(self.definition[6::])) + "\nGI: " + self.id_prot + "\nSequence: " + self.sequence)
        return self.sequence, self.definition


#Procurar a sequência por ID do gene e conversão na respetiva proteína:

    def gene_ncbi(self):
        #Função que acede ao NCBI a partir do id de um gene, retirando a respetiva sequência.
        
        handle = Entrez.efetch(db = "nucleotide", rettype="gb", retmode="xml", id = self.id_gene)      
        results = Entrez.read(handle)
        self.sequence = results[0]["GBSeq_feature-table"][2]["GBFeature_quals"][11]["GBQualifier_value"] #permite aceder diretamente à sequência proteica  
        #print("GI: " + self.id_gene + "\nSequence: " + self.sequence)
        return self.sequence
        

#Induzir as mutações:

    def induce_mutation(self, pos, aa_wt, aa_mut):
        #Funcão que insere a mutação na sequência wildtype, dando a posição, o aa da sequência wildtype e o aa da sequência mutada.
        
        mutated_sequence = ""
        if self.sequence[int(pos)-1] == aa_wt: #primeiro elemento do tuplo (aa. wt, aa. mut)
            mutated_sequence += self.sequence[:int(pos)-1]
            mutated_sequence += str(aa_mut) #aminoacido de substituição
            mutated_sequence += self.sequence[int(pos):]
            self.mut_sequence = mutated_sequence
        return self.mut_sequence
        

#Guardar a sequência num ficheiro:
    
    def protein_fasta(self):
        #Função que cria um ficheiro em formato fasta com a sequência da proteína obtida anteriormente.       
        
        if self.id_gene == 0:
            if self.mut_sequence != "":
                protein = open("mut_" + str(self.id_prot) + ".fasta", "w")
                protein.write(">mutation|" + self.id_prot + "| " + self.definition +  "\n" + self.mut_sequence)
                protein.close() 
            else:
                protein = open("wt_" + str(self.id_prot) + ".fasta", "w")
                protein.write(">gi|" + self.id_prot + "| " + self.definition +  "\n" + self.sequence)
                protein.close() 
        else:
            if self.mut_sequence != "":
                protein = open("mut_" + str(self.id_gene) + ".fasta", "w")
                protein.write(">mutation|" + self.id_gene + "\n" + self.mut_sequence)
                protein.close() 
            else:
                protein = open("wt_" + str(self.id_gene) + ".fasta", "w")
                protein.write(">gi|" + self.id_gene + "\n" + self.sequence)
                protein.close()   


#Criar fragmentos:

    def create_fragments(self, length, pos_initial, pos_final):
        #Função que cria fragmentos de tamanho x, entre uma posição inicial e final, a partir de uma sequência guardada em formato fasta.
        if self.sequence != "":
            wt_frags = []
            for pos in range(pos_initial - 1, pos_final - length + 1):
                wt_frags.append(self.sequence[pos:pos+length])
            if self.mut_sequence == "":
                save_wt_frags = open("wt_frags_" + str(self.id_prot) + ".fasta", "w")
                save_wt_frags.write(">WT \n" + str(self.edit(wt_frags)))
                save_wt_frags.close() 
        if self.mut_sequence != "":
            mut_frags = []
            for pos in range(pos_initial - 1, pos_final - length + 1):
                mut_frags.append(self.mut_sequence[pos:pos+length])
            if self.sequence == "":
                save_mut_frags = open("mut_frags_" + str(self.id_prot) + ".fasta", "w")
                save_mut_frags.write(">MUT \n" + str(self.edit(mut_frags)))
                save_mut_frags.close() 
        if self.sequence != "" and self.mut_sequence != "":
            save_frags = open("fragments_" + str(self.id_prot) + ".fasta", "w")
            save_frags.write(">WT \n" + str(self.edit(wt_frags)) + "\n\n>MUT \n" + str(self.edit(mut_frags)))
            save_frags.close() 
    

#Criar péptidos:
    
    def create_peptide(self, len_frag, aa):
        #Função que cria a sequência peptídica para ser fragmentada e analisada pelo NetMHCpan, apartir do tamanho de péptido e do aa que sofre mutação.
                
        for n in len_frag:
            if self.sequence != "":
                wt_peptide = []
                wt_peptide.append(self.sequence[aa - n : aa + n - 1])
                if self.mut_sequence == "":
                    save_wt_peptide = open(str(n) + "wt_peptide_" + str(self.id_prot) + ".fasta", "w")
                    save_wt_peptide.write(">WT \n" + str(self.edit(wt_peptide)))
                    save_wt_peptide.close() 
            if self.mut_sequence != "":
                mut_peptide = []
                mut_peptide.append(self.mut_sequence[aa - n : aa + n - 1])
                if self.sequence == "":
                    save_mut_peptide = open(str(n) + "mut_peptide_" + str(self.id_prot) + ".fasta", "w")
                    save_mut_peptide.write(">MUT \n" + str(self.edit(mut_peptide)))
                    save_mut_peptide.close() 
            if self.sequence != "" and self.mut_sequence != "":
                save_peptide = open(str(n) + "peptide_IDH1.fasta", "w")
                save_peptide.write(">WT \n" + str(self.edit(wt_peptide)) + "\n\n>MUT \n" + str(self.edit(mut_peptide)))
                save_peptide.close()

        
if __name__=="__main__":
    
    def test1(): #ID de uma proteína
        s = sequence("49168486", 0) #IDH1
        s.protein_ncbi()
        s.protein_fasta()
        s.create_fragments(9,126,140)
    
    def test2(): #ID de um gene
        s = sequence(0, "49456350") #IDH1
        s.gene_ncbi()
        s.protein_fasta()

    def test3(): 
        #s = sequence(0, 0) #IDH1
        #s.getSequencefromFasta("protein_49168486.fasta")
        s = sequence("49168486", 0) #IDH1
        s.protein_ncbi()  
        s.induce_mutation(132,"R","H") #mutation
        #s.create_fragments(9,124,140)
        s.create_peptide([8,9,10,11,12,13,14,15,16,17,18,19,20],132)


    test3()
