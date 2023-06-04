#!/usr/bin/env python3
"""
Author : Conrad Izydorczyk <conrad.izydorczyk@ucalgary.ca>
Date   : 2023-06-01
Purpose: Filter VCF file from snp-sites to encode genotypes properly
"""

import argparse


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Filter VCF file from snp-sites to encode genotypes properly',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--input_vcf',
                        help='VCF file from snp-sites to recode',
                        metavar='str',
                        type=str,
                        required=True)
    
    parser.add_argument('-o',
                        '--output_vcf',
                        help='VCF file to write',
                        metavar='str',
                        type=str,
                        required=True)


    return parser.parse_args()


# --------------------------------------------------
def main():
    """Filter VCF"""

    args = get_args()

    output_vcf_lines = []
    
    with open(args.input_vcf, "r") as infile1:
        for line in infile1:
            if line.startswith("#"):
                output_vcf_lines.append(line)
            
            elif not line.startswith("#"):
                line_elements = line.strip().split("\t") 
                alt_alleles = line_elements[4].split(",")
                genotypes = line_elements[9:]

                if "*" in alt_alleles:
                        
                    if len(alt_alleles) == 1:
                        star_pos = alt_alleles.index("*") # asterisk position in alt alleles list
                        genotypes = ["." if i == str(star_pos+1) else i for i in genotypes] # +1 since 0 index = 1st alt allele = coded as 1 in GT
                        alt_allele = "."
                        
                    elif len(alt_alleles) == 2:
                        star_pos = alt_alleles.index("*")
                        genotypes = ["." if i == str(star_pos+1) else i for i in genotypes]
                            
                        if star_pos == 0:
                            alt_pos = 1
                            alt_allele = alt_alleles[alt_pos]
                            genotypes = [str(alt_pos) if i == str(alt_pos + 1) else i for i in genotypes]
                            
                        else:
                            alt_allele = alt_alleles[0]
                        
                    elif len(alt_alleles) == 3:
                        star_pos = alt_alleles.index("*")
                        genotypes = ["." if i == str(star_pos+1) else i for i in genotypes]
                            
                        if star_pos == 0:
                            alt_pos_1 = 1
                            alt_pos_2 = 2
                                
                            alt_allele = ",".join(alt_alleles[1:])

                            genotypes = [str(alt_pos_1) if i == str(alt_pos_1 + 1) else i for i in genotypes]
                            genotypes = [str(alt_pos_2) if i == str(alt_pos_2 + 1) else i for i in genotypes]
                            
                        elif star_pos == 1:
                            alt_pos = 2

                            alt_allele = alt_alleles[0] + "," + alt_alleles[2]

                            genotypes = [str(alt_pos) if i == str(alt_pos + 1) else i for i in genotypes]
                            
                        else:
                            alt_allele = ",".join(alt_alleles[0:2])
                        
                    elif len(alt_alleles) == 4:
                        star_pos = alt_alleles.index("*")
                        genotypes = ["." if i == str(star_pos+1) else i for i in genotypes]

                        if star_pos == 0:
                            alt_pos_1 = 1
                            alt_pos_2 = 2
                            alt_pos_3 = 3

                            alt_allele = ",".join(alt_alleles[1:])

                            genotypes = [str(alt_pos_1) if i == str(alt_pos_1 + 1) else i for i in genotypes]
                            genotypes = [str(alt_pos_2) if i == str(alt_pos_2 + 1) else i for i in genotypes]
                            genotypes = [str(alt_pos_3) if i == str(alt_pos_3 + 1) else i for i in genotypes]

                        elif star_pos == 1:
                            alt_pos_1 = 2
                            alt_pos_2 = 3

                            alt_allele = alt_alleles[0] + "," + alt_alleles[2] + "," + alt_alleles[3]

                            genotypes = [str(alt_pos_1) if i == str(alt_pos_1 + 1) else i for i in genotypes]
                            genotypes = [str(alt_pos_2) if i == str(alt_pos_2 + 1) else i for i in genotypes]
                            
                        elif star_pos == 2:
                            alt_pos = 3

                            alt_allele = alt_alleles[0] + "," + alt_alleles[1] + "," + alt_alleles[3]

                            genotypes = [str(alt_pos) if i == str(alt_pos + 1) else i for i in genotypes]
                            
                        else:
                            alt_allele = ",".join(alt_alleles[0:3])
                
                line_elements[4] = alt_allele

                genotypes = ["0" if i == "00" else i for i in genotypes]
                line_elements = line_elements[0:9] + genotypes

                out_line = "\t".join(line_elements) + "\n"
                output_vcf_lines.append(out_line)
    
    with open(args.output_vcf, "w") as outfile1:
        outfile1.write("".join(output_vcf_lines))
                        






                            



# --------------------------------------------------
if __name__ == '__main__':
    main()
