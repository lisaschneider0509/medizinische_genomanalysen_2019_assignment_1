import mysql.connector
import pysam

__author__ = 'Lisa Schneider'

##
## Concept:
## TODO
##


class Assignment1:
    
    def __init__(self):
        ## Your gene of interest
        self.gene = 'CRYZL1'
        self.alignment_file = pysam.AlignmentFile("chr21.bam", "rb")
        self.download_gene_coordinates("hg38", "gene_coordinates.txt")
        self.gene_list = self.read_gene_coordinates()
        self.gene_coordinates = self.get_coordinates_of_gene()
        self.reads = self.get_all_reads()



    
    def download_gene_coordinates(self, genome_reference, file_name):
        ## TODO concept
        
        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query
        query = "SELECT DISTINCT {} from refGene WHERE name2 = '{}'".format(",".join(query_fields), self.gene)
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        ## TODO this may need some work 
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
    
            
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print("Done fetching data\n")

    def read_gene_coordinates(self):
        with open("gene_coordinates.txt", "r") as fh:
            for row in fh:
                if self.gene in row:
                    gene_list = row.replace(")", "").replace("(", "").replace("'", "")
                    gene_list = gene_list.split(", ")
                    break
        return gene_list

    def get_coordinates_of_gene(self):
        gene = self.gene_list[0]
        chromosome = self.gene_list[2]
        start = self.gene_list[3]
        stop = self.gene_list[4]
        print(f"Coordinates of gene {gene}: \n"
              f"\tchromosome - {chromosome} \n"
              f"\tstart - {start} \n"
              f"\tstop - {stop}")
        return [chromosome, int(start), int(stop)]
        
    def get_gene_symbol(self):
        print(f"Gene symbol: {self.gene}")
                        
    def get_sam_header(self):
        print("todo")

    def get_all_reads(self):
        reads = list(self.alignment_file.fetch(self.gene_coordinates[0],
                                               self.gene_coordinates[1],
                                               self.gene_coordinates[2]))
        return reads

    def get_properly_paired_reads_of_gene(self):
        proper_reads = len([i for i in self.reads if i.is_proper_pair])
        print(f"Number of properly paired reads: {proper_reads}")
        
    def get_gene_reads_with_indels(self):
        print("todo")
        
    def calculate_total_average_coverage(self):
        print("todo")
        
    def calculate_gene_average_coverage(self):
        print("todo")
        
    def get_number_mapped_reads(self):
        print("todo")

    def get_region_of_gene(self):
        print("todo")
        
    def get_number_of_exons(self):
        exons = self.gene_list[6]
        print(f"exons : {exons}")
    
    
    def print_summary(self):
        self.read_gene_coordinates()
        self.get_gene_symbol()
        self.get_number_of_exons()
        self.get_properly_paired_reads_of_gene()

        print("Print all results here")
    
    
def main():
    print("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()
    
    
    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()
