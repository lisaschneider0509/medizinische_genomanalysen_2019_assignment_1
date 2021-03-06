import mysql.connector
import pysam
import numpy

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
        self.chromosome = self.gene_coordinates[0]
        self.start = self.gene_coordinates[1]
        self.stop = self.gene_coordinates[2]
        self.reads = self.get_all_reads()

    def download_gene_coordinates(self, genome_reference, file_name):
        '''
        :param genome_reference:
        :param file_name:
        :return:
        '''
        ## TODO concept
        
        print("Connecting to UCSC to fetch data")
        
        # Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        # Get cursor
        cursor = cnx.cursor()
        
        # Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        # Build query
        query = "SELECT DISTINCT {} from refGene WHERE name2 = '{}'".format(",".join(query_fields), self.gene)
        
        # Execute query
        cursor.execute(query)
        
        ## Write to file
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
        chromosome = self.gene_list[2]
        start = self.gene_list[3]
        stop = self.gene_list[4]
        return [chromosome, int(start), int(stop)]
        
    def get_gene_symbol(self):
        # todo check if gene symbol is right
        print(f"Gene symbol: {self.gene}")
                        
    def get_sam_header(self):
        # todo change format
        print("\nSam Header:")
        for k, v in self.alignment_file.header["HD"].items():
            if k == "SO":
                print(f"\tSO - Sorting order of Alignments: {v}")
            if k == "VN":
                print(f"\tVN - Format version: {v}")
            if k == "GO":
                print(f"\tGO - Grouping of alignments: {v} \n")

    def get_all_reads(self):
        reads = list(self.alignment_file.fetch(self.chromosome,
                                               self.start,
                                               self.stop))
        return reads

    def get_properly_paired_reads_of_gene(self):
        proper_reads = len([i for i in self.reads if i.is_proper_pair])
        print(f"Properly paired reads: {proper_reads}")
        
    def get_gene_reads_with_indels(self):
        indel_list = []
        for i in self.reads:
            if not i.is_unmapped:
                cig = i.cigartuples
                for (operation, length) in cig:
                    if (operation == 1) or (operation == 2): # cigar operation 1 = insertion, 2 = deletion
                        indel_list.append(i)
        indel_number = len(indel_list)
        print(f"Reads with indels: {indel_number}")
        
    def calculate_total_average_coverage(self):
        total_coverage = []
        for pileup in self.alignment_file.pileup(self.chromosome):
            total_coverage.append(pileup.n)
        average_total_coverage = sum(total_coverage) / float(len(total_coverage))
        print(f"Average total coverage: {average_total_coverage} \n")

    def calculate_gene_average_coverage(self):
        coverage = self.alignment_file.count_coverage(self.chromosome,
                                                      start=self.start,
                                                      stop=self.stop)
        average_gene_coverage = round(numpy.mean(coverage), 2)
        print(f"Average gene coverage: {average_gene_coverage}")
        
    def get_number_mapped_reads(self):
        mapped_reads = 0
        for i in self.reads:
            if not i.is_unmapped:
                mapped_reads += 1
        print(f"Mapped reads: {mapped_reads}")

    def get_region_of_gene(self):
        # todo actually gene locus i think
        print(f"Genomic Region:\n"
              f"\tChromosome: {self.chromosome}\n"
              f"\tStart: {self.start}\n"
              f"\tStop: {self.stop} \n")
        
    def get_number_of_exons(self):
        exons = self.gene_list[6]
        print(f"Exons : {exons}")

    def print_summary(self):
        self.read_gene_coordinates()
        print(f"Coordinates of gene {self.gene}: \n"
              f"\tChromosome - {self.chromosome} \n"
              f"\tStart - {self.start} \n"
              f"\tStop - {self.stop} \n")
        self.get_gene_symbol()
        self.get_number_of_exons()
        self.get_properly_paired_reads_of_gene()
        self.get_number_mapped_reads()
        self.get_sam_header()
        self.calculate_total_average_coverage()
        self.calculate_gene_average_coverage()
        self.get_region_of_gene()
        self.get_gene_reads_with_indels()


def main():
    assignment1 = Assignment1()
    assignment1.print_summary()


if __name__ == '__main__':
    main()
